#!/usr/bin/env python
"""
7a — Run reninness-score bedspace training in batch over multiple param configs.

Self-contained integration of three originally-separate files from
``geniml_dev/geniml/bedspace/tests/``:

    1. wandb_integration.py           (wandb wrapper + streaming StarSpace trainer)
    2. test_training_evaluation.py    (embedding stats, t-SNE/UMAP, eval report)
    3. experiment_3cluster_sweep.py   (SWEEP_CONFIGS catalog + orchestration)

Pipeline per config:
    preprocess (once, shared)  →  StarSpace train  →  evaluate (stats + plots)
    (optional) multi-k clustering diagnostic  →  wandb logging if enabled

For the LOOCV-flavored pipeline (preprocess → train → distances → score, per
held-out sample), see ``7b_run_reninness_score_loo_in_batch.py``.

------------------------------------------------------------------
USAGE
------------------------------------------------------------------

    source /home/bingjiexue/Documents/code/reninness_workspace/.venv/bin/activate
    pip install wandb        # optional, only if you want wandb tracking
    wandb login              # optional

    # Run all SWEEP_CONFIGS:
    python 7a_run_reninness_score_in_batch.py \\
        --data-path  /path/to/bed/files/ \\
        --metadata   /path/to/metadata.csv \\
        --universe   /path/to/universe.bed \\
        --labels     "label" \\
        --starspace  /path/to/Starspace/ \\
        --output-dir /path/to/output/ \\
        --wandb-project bedspace-3cluster        # optional

    # Run a subset of configs:
    python 7a_run_reninness_score_in_batch.py --configs d3_e10_m0.1 d10_e10_m1.0 ...

    # List available configs and exit:
    python 7a_run_reninness_score_in_batch.py --list-configs

Output (under <output-dir>/<experiment_name>/):
    preprocessed/                          (shared across all configs)
    <config_name>/
        model_<config>.tsv                 trained StarSpace model
        <config>_stats.json                statistical metrics
        <config>_report.txt                summary report
        <config>_tsne.png, _umap.png       embedding visualizations
        <config>_silhouette_vs_k.png       multi-k clustering diagnostic
        <config>_elbow.png
    sweep_comparison.txt                   side-by-side comparison table

This file is bundled into renin_atac/src/ so it remains usable even if the
upstream tests/ folder is moved or removed. It only depends on the installed
``geniml`` package (``geniml.bedspace.preprocess`` and ``geniml.bedspace.const``).
"""

from __future__ import annotations

import argparse
import datetime
import json
import logging
import os
import re
import subprocess
import sys
import threading
import time
from datetime import datetime as _datetime
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import gaussian_kde
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

try:
    import wandb
    HAS_WANDB = True
except ImportError:
    wandb = None
    HAS_WANDB = False

from geniml.bedspace.preprocess import main as preprocess_main
from geniml.bedspace.const import MODE_BULK

_LOGGER = logging.getLogger("bedspace.batch")


# =============================================================================
# SECTION 1 — wandb integration  (from wandb_integration.py)
# =============================================================================

WANDB_PROGRESS_THROTTLE_SECONDS = 5

# stdout: "---+++ Epoch 1 Train error : 0.123"
_EPOCH_PATTERN = re.compile(
    r"---\+\+\+\s+Epoch\s+(\d+)\s+Train error\s*:\s*([\d.]+(?:e[+-]?\d+)?)"
)
# stderr progress: "Epoch:  12.3% lr: 0.0100 loss: 0.1234"
_PROGRESS_PATTERN = re.compile(
    r"Epoch:\s*([\d.]+)%\s+lr:\s*([\d.]+(?:e[+-]?\d+)?)\s+loss:\s*([\d.]+(?:e[+-]?\d+)?)"
)


def parse_starspace_line(line: str) -> Optional[Dict]:
    """Parse a single line of StarSpace output into structured metrics."""
    m = _EPOCH_PATTERN.search(line)
    if m:
        return {"type": "epoch", "epoch": int(m.group(1)), "train_error": float(m.group(2))}
    m = _PROGRESS_PATTERN.search(line)
    if m:
        return {"type": "progress", "pct": float(m.group(1)),
                "lr": float(m.group(2)), "loss": float(m.group(3))}
    return None


def _flatten_dict(d: Dict, prefix: str = "", sep: str = "/") -> Dict:
    flat = {}
    for key, value in d.items():
        full_key = f"{prefix}{sep}{key}" if prefix else key
        if isinstance(value, dict):
            flat.update(_flatten_dict(value, full_key, sep))
        else:
            flat[full_key] = value
    return flat


class WandbTracker:
    """No-op safe wrapper around the wandb API."""

    def __init__(self, project: Optional[str] = None, config: Optional[Dict] = None,
                 run_name: Optional[str] = None, tags: Optional[List[str]] = None,
                 group: Optional[str] = None, enabled: bool = True):
        self.active = False
        self._last_progress_time = 0.0
        self._run = None
        if not enabled or project is None:
            return
        if not HAS_WANDB:
            _LOGGER.warning("wandb not installed. Install with: pip install wandb.")
            return
        try:
            self._run = wandb.init(
                project=project, config=config or {}, name=run_name,
                tags=tags, group=group, reinit=True,
            )
            self.active = True
            _LOGGER.info(f"wandb run initialized: {self._run.url}")
        except Exception as e:
            _LOGGER.warning(f"Failed to initialize wandb: {e}")

    def log_epoch(self, epoch: int, train_error: float, lr: Optional[float] = None):
        if not self.active:
            return
        metrics = {"train/loss": train_error, "train/epoch": epoch}
        if lr is not None:
            metrics["train/lr"] = lr
        wandb.log(metrics, step=epoch)

    def log_progress(self, pct: float, lr: float, loss: float):
        if not self.active:
            return
        now = time.time()
        if now - self._last_progress_time < WANDB_PROGRESS_THROTTLE_SECONDS:
            return
        self._last_progress_time = now
        wandb.log({"train/progress_pct": pct, "train/lr": lr, "train/loss": loss})

    def log_evaluation(self, eval_results: Dict):
        if not self.active:
            return
        flat = _flatten_dict(eval_results, prefix="eval")
        for key, value in flat.items():
            if isinstance(value, (int, float)):
                wandb.summary[key] = value

    def log_images(self, tsne_path: Optional[str] = None, umap_path: Optional[str] = None):
        if not self.active:
            return
        images = {}
        if tsne_path:
            images["viz/tsne"] = wandb.Image(tsne_path)
        if umap_path:
            images["viz/umap"] = wandb.Image(umap_path)
        if images:
            wandb.log(images)

    def finish(self):
        if not self.active:
            return
        try:
            wandb.finish()
        except Exception as e:
            _LOGGER.warning(f"Error finishing wandb run: {e}")
        self.active = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.finish()
        return False


def _stream_pipe(pipe, tracker: WandbTracker, stream_out, is_stderr: bool = False):
    """Read lines from a subprocess pipe, print to terminal, and log metrics."""
    last_epoch_lr = None
    for raw_line in iter(pipe.readline, ""):
        if is_stderr:
            sub_lines = raw_line.replace("\r", "\n").split("\n")
        else:
            sub_lines = [raw_line]
        for line in sub_lines:
            stripped = line.rstrip("\n\r")
            if not stripped:
                continue
            if is_stderr:
                print(f"\r{stripped}", end="", file=stream_out, flush=True)
            else:
                print(stripped, file=stream_out, flush=True)
            parsed = parse_starspace_line(stripped)
            if parsed is None:
                continue
            if parsed["type"] == "epoch":
                tracker.log_epoch(epoch=parsed["epoch"], train_error=parsed["train_error"],
                                  lr=last_epoch_lr)
            elif parsed["type"] == "progress":
                last_epoch_lr = parsed["lr"]
                tracker.log_progress(pct=parsed["pct"], lr=parsed["lr"], loss=parsed["loss"])
    pipe.close()


def run_train_with_wandb(starspace_path: str, train_input: str, output_dir: str,
                         experiment_name: str, params: Dict[str, str],
                         tracker: WandbTracker) -> Optional[str]:
    """Run StarSpace training with streaming output and wandb metric logging."""
    print(f"\n  Training with config: {experiment_name}")
    print(f"    Parameters: {params}")

    model_output = os.path.join(output_dir, f"model_{experiment_name}")
    starspace_bin = os.path.join(starspace_path, "starspace")

    cmd = [
        starspace_bin, "train",
        "-trainFile", train_input,
        "-model", model_output,
        "-trainMode", "0",
        "-dim", params["dim"],
        "-epoch", params["epoch"],
        "-negSearchLimit", params["negSearchLimit"],
        "-lr", params["lr"],
        "-minCount", params["minCount"],
        "-thread", "20",
    ]
    if params.get("adagrad"):
        cmd.extend(["-adagrad", "true"])
    if "margin" in params:
        cmd.extend(["-margin", params["margin"]])

    print(f"    Command: {' '.join(cmd)}")
    start_time = datetime.datetime.now()

    if tracker.active:
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            text=True, bufsize=1,
        )
        stdout_thread = threading.Thread(
            target=_stream_pipe, args=(process.stdout, tracker, sys.stdout, False),
            daemon=True,
        )
        stderr_thread = threading.Thread(
            target=_stream_pipe, args=(process.stderr, tracker, sys.stderr, True),
            daemon=True,
        )
        stdout_thread.start()
        stderr_thread.start()
        process.wait()
        stdout_thread.join()
        stderr_thread.join()
        print(file=sys.stderr)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    ERROR: Training failed!\n    stderr: {result.stderr}")
            return None

    duration = (datetime.datetime.now() - start_time).total_seconds()
    model_tsv = model_output + ".tsv"
    if os.path.exists(model_tsv):
        print(f"    Training complete in {duration:.1f}s\n    Model saved to: {model_tsv}")
        return model_tsv
    print(f"    ERROR: Model file not found: {model_tsv}")
    return None


def run_train_fallback(starspace_path, train_input, output_dir, experiment_name, params):
    """Train without wandb (original subprocess.run approach)."""
    model_output = os.path.join(output_dir, f"model_{experiment_name}")
    starspace_bin = os.path.join(starspace_path, "starspace")
    cmd = [
        starspace_bin, "train",
        "-trainFile", train_input,
        "-model", model_output,
        "-trainMode", "0",
        "-dim", params["dim"],
        "-epoch", params["epoch"],
        "-negSearchLimit", params["negSearchLimit"],
        "-lr", params["lr"],
        "-minCount", params["minCount"],
        "-thread", "20",
    ]
    if params.get("adagrad"):
        cmd.extend(["-adagrad", "true"])
    if "margin" in params:
        cmd.extend(["-margin", params["margin"]])

    print(f"    Command: {' '.join(cmd)}")
    start = datetime.datetime.now()
    result = subprocess.run(cmd, capture_output=True, text=True)
    duration = (datetime.datetime.now() - start).total_seconds()

    if result.returncode != 0:
        print(f"    ERROR: {result.stderr}")
        return None

    model_tsv = model_output + ".tsv"
    if os.path.exists(model_tsv):
        print(f"    Training complete in {duration:.1f}s → {model_tsv}")
        return model_tsv
    print(f"    ERROR: Model not found: {model_tsv}")
    return None


# =============================================================================
# SECTION 2 — Embedding evaluation  (from test_training_evaluation.py)
# =============================================================================

FOLDCHANGE_CMAP = mcolors.LinearSegmentedColormap.from_list(
    "BuYlRd", ["#2166ac", "#fee08b", "#b2182b"]
)

OUTPUT_DIR = Path(os.environ.get("BEDSPACE_OUTPUT_DIR", "test_outputs"))


def load_model_embeddings(model_tsv_path: str, label_prefix: str = "__label__"):
    """Load region and label embeddings from a StarSpace .tsv model."""
    df = pd.read_csv(model_tsv_path, sep="\t", header=None)
    label_mask = df[0].str.startswith(label_prefix)
    label_df = df[label_mask].reset_index(drop=True)
    label_names = [name.replace(label_prefix, "") for name in label_df[0].tolist()]
    label_embeddings = label_df.iloc[:, 1:].values.astype(float)
    region_df = df[~label_mask].reset_index(drop=True)
    region_names = region_df[0].tolist()
    region_embeddings = region_df.iloc[:, 1:].values.astype(float)
    return region_embeddings, region_names, label_embeddings, label_names


def compute_statistics(region_embeddings, label_embeddings, label_names):
    """Compute embedding stats (norms, region-to-label distances, label separation)."""
    stats = {
        "embedding_dim": region_embeddings.shape[1],
        "num_regions": region_embeddings.shape[0],
        "num_labels": label_embeddings.shape[0],
    }
    region_norms = np.linalg.norm(region_embeddings, axis=1)
    label_norms = np.linalg.norm(label_embeddings, axis=1)
    stats["avg_region_norm"] = float(np.mean(region_norms))
    stats["std_region_norm"] = float(np.std(region_norms))
    stats["avg_label_norm"] = float(np.mean(label_norms))

    if len(label_embeddings) > 0 and len(region_embeddings) > 0:
        distances = cdist(region_embeddings, label_embeddings, metric="cosine")
        stats["region_label_distances"] = {
            label_names[i]: {
                "mean": float(np.mean(distances[:, i])),
                "std": float(np.std(distances[:, i])),
                "min": float(np.min(distances[:, i])),
                "max": float(np.max(distances[:, i])),
            }
            for i in range(len(label_names))
        }

    if len(label_embeddings) > 1:
        label_distances = cdist(label_embeddings, label_embeddings, metric="cosine")
        upper_tri = label_distances[np.triu_indices(len(label_names), k=1)]
        stats["label_separation"] = {
            "mean_pairwise_distance": float(np.mean(upper_tri)),
            "min_pairwise_distance": float(np.min(upper_tri)),
            "max_pairwise_distance": float(np.max(upper_tri)),
        }
        stats["label_pairwise"] = {}
        for i, n1 in enumerate(label_names):
            for j, n2 in enumerate(label_names):
                if i < j:
                    stats["label_pairwise"][f"{n1}_vs_{n2}"] = float(label_distances[i, j])
    return stats


def _label_colors(label_names):
    """Color list: 'Renin' → orange, 'Non_Renin' → blue, others → tab10."""
    tab10 = plt.cm.tab10(np.linspace(0, 1, max(len(label_names), 1)))
    return [
        "darkorange" if n.lower() == "renin"
        else "royalblue" if "non" in n.lower() and "renin" in n.lower()
        else tab10[i]
        for i, n in enumerate(label_names)
    ]


def _plot_label_markers(ax, label_coords, label_names):
    colors = _label_colors(label_names)
    for i, (name, color) in enumerate(zip(label_names, colors)):
        ax.scatter(label_coords[i, 0], label_coords[i, 1], c=[color], s=200,
                   marker="^", edgecolors="black", linewidths=1.5,
                   label=name, zorder=10)
        ax.annotate(name, (label_coords[i, 0], label_coords[i, 1]),
                    xytext=(5, 5), textcoords="offset points",
                    fontsize=10, fontweight="bold")


def create_tsne_visualization(region_embeddings, label_embeddings, label_names,
                              output_path, perplexity=30, random_state=42,
                              max_regions=5000):
    """t-SNE: regions = grey circles, labels = colored triangles."""
    if len(region_embeddings) > max_regions:
        idx = np.random.choice(len(region_embeddings), max_regions, replace=False)
        region_embeddings_plot = region_embeddings[idx]
    else:
        region_embeddings_plot = region_embeddings

    all_embeddings = np.vstack([region_embeddings_plot, label_embeddings])
    print(f"Running t-SNE on {len(all_embeddings)} embeddings...")
    tsne = TSNE(n_components=2, perplexity=min(perplexity, len(all_embeddings) - 1),
                random_state=random_state, metric="cosine")
    coords = tsne.fit_transform(all_embeddings)
    region_coords = coords[: len(region_embeddings_plot)]
    label_coords = coords[len(region_embeddings_plot):]

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.scatter(region_coords[:, 0], region_coords[:, 1], c="lightgrey",
               s=10, alpha=0.5, marker="o", label="Regions")
    _plot_label_markers(ax, label_coords, label_names)
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_title("t-SNE: Regions (circles) vs Labels (triangles)")
    ax.legend(loc="best", fontsize=10)
    plt.tight_layout(); plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved t-SNE plot to {output_path}")


def create_umap_visualization(region_embeddings, label_embeddings, label_names,
                              output_path, n_neighbors=15, min_dist=0.1,
                              random_state=42, max_regions=5000):
    """UMAP: regions = grey circles, labels = colored triangles."""
    if not HAS_UMAP:
        print("UMAP not available, skipping UMAP visualization")
        return
    if len(region_embeddings) > max_regions:
        idx = np.random.choice(len(region_embeddings), max_regions, replace=False)
        region_embeddings_plot = region_embeddings[idx]
    else:
        region_embeddings_plot = region_embeddings

    all_embeddings = np.vstack([region_embeddings_plot, label_embeddings])
    print(f"Running UMAP on {len(all_embeddings)} embeddings...")
    reducer = umap.UMAP(n_neighbors=min(n_neighbors, len(all_embeddings) - 1),
                        min_dist=min_dist, metric="cosine", random_state=random_state)
    coords = reducer.fit_transform(all_embeddings)
    region_coords = coords[: len(region_embeddings_plot)]
    label_coords = coords[len(region_embeddings_plot):]

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.scatter(region_coords[:, 0], region_coords[:, 1], c="lightgrey",
               s=10, alpha=0.5, marker="o", label="Regions")
    _plot_label_markers(ax, label_coords, label_names)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title("UMAP: Regions (circles) vs Labels (triangles)")
    ax.legend(loc="best", fontsize=10)
    plt.tight_layout(); plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved UMAP plot to {output_path}")


def _plot_value_distribution(ax, values, norm):
    """Horizontal histogram + KDE aligned to a colorbar y-range."""
    normed = np.asarray(norm(np.asarray(values, dtype=float)))
    ax.hist(normed, bins=40, orientation="horizontal", color="gray",
            alpha=0.4, density=True, edgecolor="black", linewidth=0.3)
    if len(normed) > 1 and float(np.std(normed)) > 0:
        try:
            kde = gaussian_kde(normed)
            yy = np.linspace(0.0, 1.0, 200)
            ax.plot(kde(yy), yy, color="black", linewidth=1.5)
        except Exception:
            pass
    zero_norm = float(norm(0.0))
    if 0.0 <= zero_norm <= 1.0:
        ax.axhline(zero_norm, color="dimgray", linewidth=0.6, linestyle="--", alpha=0.7)
    ax.set_ylim(0.0, 1.0); ax.set_yticks([])
    ax.set_xlabel("Density", fontsize=9); ax.tick_params(axis="x", labelsize=8)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)


def _foldchange_plot(coords, fc_values, label_coords, label_names, output_path, title):
    fig = plt.figure(figsize=(16, 10), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[20, 0.5, 3], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])
    dax = fig.add_subplot(gs[0, 2])
    fc_min, fc_max = float(np.min(fc_values)), float(np.max(fc_values))
    if fc_min < 0 < fc_max:
        norm = mcolors.TwoSlopeNorm(vmin=fc_min, vcenter=0, vmax=fc_max)
    else:
        norm = mcolors.Normalize(vmin=fc_min, vmax=fc_max)
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=fc_values, cmap=FOLDCHANGE_CMAP,
                    norm=norm, s=10, alpha=0.7, marker="o")
    fig.colorbar(sc, cax=cax, label="log2FoldChange")
    _plot_value_distribution(dax, fc_values, norm)
    _plot_label_markers(ax, label_coords, label_names)
    ax.set_xlabel(title.split(":")[0] + " 1"); ax.set_ylabel(title.split(":")[0] + " 2")
    ax.set_title(title); ax.legend(loc="best", fontsize=10)
    plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved {title} plot to {output_path}")


def create_tsne_foldchange(region_embeddings, region_names, foldchange_map,
                           label_embeddings, label_names, output_path,
                           perplexity=30, random_state=42, max_regions=5000):
    """t-SNE with regions colored by log2FoldChange."""
    if len(region_embeddings) > max_regions:
        idx = np.random.choice(len(region_embeddings), max_regions, replace=False)
        region_embeddings_plot = region_embeddings[idx]
        region_names_plot = [region_names[i] for i in idx]
    else:
        region_embeddings_plot = region_embeddings
        region_names_plot = region_names

    fc_values = np.array([foldchange_map.get(name, 0.0) for name in region_names_plot])
    all_embeddings = np.vstack([region_embeddings_plot, label_embeddings])
    print(f"Running t-SNE (foldchange) on {len(all_embeddings)} embeddings...")
    tsne = TSNE(n_components=2, perplexity=min(perplexity, len(all_embeddings) - 1),
                random_state=random_state, metric="cosine")
    coords = tsne.fit_transform(all_embeddings)
    region_coords = coords[: len(region_embeddings_plot)]
    label_coords = coords[len(region_embeddings_plot):]
    _foldchange_plot(region_coords, fc_values, label_coords, label_names,
                     output_path, "t-SNE: Regions colored by log2FoldChange")


def create_umap_foldchange(region_embeddings, region_names, foldchange_map,
                           label_embeddings, label_names, output_path,
                           n_neighbors=15, min_dist=0.1, random_state=42,
                           max_regions=5000):
    """UMAP with regions colored by log2FoldChange."""
    if not HAS_UMAP:
        print("UMAP not available, skipping UMAP foldchange visualization")
        return
    if len(region_embeddings) > max_regions:
        idx = np.random.choice(len(region_embeddings), max_regions, replace=False)
        region_embeddings_plot = region_embeddings[idx]
        region_names_plot = [region_names[i] for i in idx]
    else:
        region_embeddings_plot = region_embeddings
        region_names_plot = region_names

    fc_values = np.array([foldchange_map.get(name, 0.0) for name in region_names_plot])
    all_embeddings = np.vstack([region_embeddings_plot, label_embeddings])
    print(f"Running UMAP (foldchange) on {len(all_embeddings)} embeddings...")
    reducer = umap.UMAP(n_neighbors=min(n_neighbors, len(all_embeddings) - 1),
                        min_dist=min_dist, metric="cosine", random_state=random_state)
    coords = reducer.fit_transform(all_embeddings)
    region_coords = coords[: len(region_embeddings_plot)]
    label_coords = coords[len(region_embeddings_plot):]
    _foldchange_plot(region_coords, fc_values, label_coords, label_names,
                     output_path, "UMAP: Regions colored by log2FoldChange")


def generate_summary_report(stats, model_path):
    """Generate a short embedding-evaluation summary."""
    lines = ["=" * 70, "BEDSPACE TRAINING EVALUATION REPORT", "=" * 70,
             f"Generated: {_datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
             f"Model: {model_path}", "",
             "=" * 70, "EMBEDDING STATISTICS", "=" * 70,
             f"  Embedding dimension: {stats['embedding_dim']}",
             f"  Number of regions: {stats['num_regions']}",
             f"  Number of labels: {stats['num_labels']}",
             f"  Avg region norm: {stats['avg_region_norm']:.4f} "
             f"(+/- {stats['std_region_norm']:.4f})",
             f"  Avg label norm: {stats['avg_label_norm']:.4f}"]
    if "label_separation" in stats:
        lines += ["", "  Label Separation (cosine distance):",
                  f"    Mean pairwise: {stats['label_separation']['mean_pairwise_distance']:.4f}",
                  f"    Min pairwise:  {stats['label_separation']['min_pairwise_distance']:.4f}",
                  f"    Max pairwise:  {stats['label_separation']['max_pairwise_distance']:.4f}"]
    if "region_label_distances" in stats:
        lines += ["", "  Region-to-Label Distances (cosine):"]
        for label, dist in stats["region_label_distances"].items():
            lines.append(f"    {label}: mean={dist['mean']:.4f}, std={dist['std']:.4f}")
    lines += ["", "=" * 70, "QUALITY ASSESSMENT", "=" * 70]
    if "label_separation" in stats:
        sep = stats["label_separation"]["mean_pairwise_distance"]
        if sep > 0.5:
            q = "GOOD - Labels are well separated"
        elif sep > 0.3:
            q = "MODERATE - Labels have some separation"
        else:
            q = "POOR - Labels are too close together"
        lines.append(f"  Label separation: {q} (distance={sep:.4f})")
    lines.append("=" * 70)
    return "\n".join(lines)


def evaluate_model(model_tsv_path: str, output_name: str = None,
                   output_dir=None, deseq_path=None):
    """Full evaluation on a trained StarSpace model."""
    if output_name is None:
        output_name = Path(model_tsv_path).stem
    out_dir = Path(output_dir) if output_dir else OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'='*60}\nEvaluating model: {model_tsv_path}\n{'='*60}\n")
    region_embeddings, region_names, label_embeddings, label_names = load_model_embeddings(model_tsv_path)
    print(f"  Loaded {len(region_names)} region embeddings")
    print(f"  Loaded {len(label_names)} label embeddings: {label_names}")
    print(f"  Embedding dimension: {region_embeddings.shape[1]}")

    print("\nComputing statistics...")
    stats = compute_statistics(region_embeddings, label_embeddings, label_names)

    stats_path = out_dir / f"{output_name}_stats.json"
    with open(stats_path, "w") as f:
        json.dump(stats, f, indent=2)
    print(f"  Saved statistics to {stats_path}")

    print("\nGenerating summary report...")
    report = generate_summary_report(stats, model_tsv_path)
    report_path = out_dir / f"{output_name}_report.txt"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"  Saved report to {report_path}")
    print("\n" + report)

    print("\nCreating visualizations...")
    tsne_path = out_dir / f"{output_name}_tsne.png"
    create_tsne_visualization(region_embeddings, label_embeddings, label_names, tsne_path)
    umap_path = out_dir / f"{output_name}_umap.png"
    create_umap_visualization(region_embeddings, label_embeddings, label_names, umap_path)

    if deseq_path is not None:
        try:
            print(f"\nLoading DESeq data from {deseq_path}")
            deseq_df = pd.read_csv(deseq_path)
            region_keys = (deseq_df["seqnames"].astype(str) + "_" +
                           deseq_df["start"].astype(str) + "_" +
                           deseq_df["end"].astype(str))
            foldchange_map = dict(zip(region_keys, deseq_df["log2FoldChange"]))
            create_tsne_foldchange(region_embeddings, region_names, foldchange_map,
                                   label_embeddings, label_names,
                                   out_dir / f"{output_name}_tsne_foldchange.png")
            create_umap_foldchange(region_embeddings, region_names, foldchange_map,
                                   label_embeddings, label_names,
                                   out_dir / f"{output_name}_umap_foldchange.png")
        except Exception as e:
            print(f"    WARNING: DESeq foldchange plots failed: {e}")

    print(f"\n{'='*60}\nEvaluation complete!\n{'='*60}")
    return stats


# =============================================================================
# SECTION 3 — Experiment orchestration  (from experiment_3cluster_sweep.py)
# =============================================================================

EXPERIMENT_DIR = (
    Path(os.environ.get("BEDSPACE_OUTPUT_DIR", "test_outputs")) / "experiments_3cluster"
)
DEFAULT_STARSPACE_PATH = os.environ.get("STARSPACE_PATH", "")

# Sweep configs: comprehensive parameter search for 3-cluster region embeddings.
# See the upstream experiment_3cluster_sweep.py for the full strategy comments.
SWEEP_CONFIGS = {
    # --- REFERENCE: previous 2-cluster winners ---
    "ref_d10_e10_m0.3": {"description": "REFERENCE: 2-cluster winner dim=10 epoch=10 margin=0.3",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "ref_d50_e10_m0.3": {"description": "REFERENCE: 2-cluster winner dim=50 epoch=10 margin=0.3",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "ref_d10_e50_m0.3": {"description": "REFERENCE: 2-cluster winner dim=10 epoch=50 margin=0.3",
        "epoch": "50", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    # --- LOWER MARGINS at winning dim/epoch combos ---
    "d10_e10_m0.05": {"description": "dim=10, epoch=10, margin=0.05 — very gentle",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.05", "adagrad": False},
    "d10_e10_m0.1": {"description": "dim=10, epoch=10, margin=0.1",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    "d10_e10_m0.15": {"description": "dim=10, epoch=10, margin=0.15",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    "d10_e10_m0.2": {"description": "dim=10, epoch=10, margin=0.2",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    "d50_e10_m0.05": {"description": "dim=50, epoch=10, margin=0.05 — very gentle",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.05", "adagrad": False},
    "d50_e10_m0.1": {"description": "dim=50, epoch=10, margin=0.1",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    "d50_e10_m0.15": {"description": "dim=50, epoch=10, margin=0.15",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    "d50_e10_m0.2": {"description": "dim=50, epoch=10, margin=0.2",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    # --- VERY EARLY STOPPING ---
    "d10_e3_m0.3": {"description": "dim=10, epoch=3, margin=0.3 — very early stop",
        "epoch": "3", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d10_e5_m0.3": {"description": "dim=10, epoch=5, margin=0.3 — early stop",
        "epoch": "5", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d50_e3_m0.3": {"description": "dim=50, epoch=3, margin=0.3 — very early stop",
        "epoch": "3", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d50_e5_m0.3": {"description": "dim=50, epoch=5, margin=0.3 — early stop",
        "epoch": "5", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d10_e5_m0.15": {"description": "dim=10, epoch=5, margin=0.15",
        "epoch": "5", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    "d50_e5_m0.15": {"description": "dim=50, epoch=5, margin=0.15",
        "epoch": "5", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    # --- DIM x MARGIN GRID at epoch=10, lr=0.0005 ---
    "d3_e10_m0.1": {"description": "dim=3, epoch=10, margin=0.1",
        "epoch": "10", "negSearchLimit": "50", "dim": "3", "lr": "0.0005",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    "d3_e10_m0.2": {"description": "dim=3, epoch=10, margin=0.2",
        "epoch": "10", "negSearchLimit": "50", "dim": "3", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    "d3_e10_m0.3": {"description": "dim=3, epoch=10, margin=0.3",
        "epoch": "10", "negSearchLimit": "50", "dim": "3", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d5_e10_m0.1": {"description": "dim=5, epoch=10, margin=0.1",
        "epoch": "10", "negSearchLimit": "50", "dim": "5", "lr": "0.0005",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    "d5_e10_m0.2": {"description": "dim=5, epoch=10, margin=0.2",
        "epoch": "10", "negSearchLimit": "50", "dim": "5", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    "d5_e10_m0.3": {"description": "dim=5, epoch=10, margin=0.3",
        "epoch": "10", "negSearchLimit": "50", "dim": "5", "lr": "0.0005",
        "minCount": "1", "margin": "0.3", "adagrad": False},
    "d10_e10_m0.5": {"description": "dim=10, epoch=10, margin=0.5 — strong separation",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.5", "adagrad": False},
    "d50_e10_m0.5": {"description": "dim=50, epoch=10, margin=0.5 — strong separation",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.5", "adagrad": False},
    "d3_e10_m1.0": {"description": "dim=3, epoch=10, margin=1.0 — high-margin, low-dim",
        "epoch": "10", "negSearchLimit": "50", "dim": "3", "lr": "0.0005",
        "minCount": "1", "margin": "1.0", "adagrad": False},
    "d5_e10_m1.0": {"description": "dim=5, epoch=10, margin=1.0 — high-margin, low-dim",
        "epoch": "10", "negSearchLimit": "50", "dim": "5", "lr": "0.0005",
        "minCount": "1", "margin": "1.0", "adagrad": False},
    "d50_e10_m1.0": {"description": "dim=50, epoch=10, margin=1.0 — high-margin, high-dim",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "1.0", "adagrad": False},
    # --- HARDER NEGATIVES ---
    "d10_e10_m0.2_neg100": {"description": "dim=10, epoch=10, margin=0.2, negSearch=100",
        "epoch": "10", "negSearchLimit": "100", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    "d50_e10_m0.2_neg100": {"description": "dim=50, epoch=10, margin=0.2, negSearch=100",
        "epoch": "10", "negSearchLimit": "100", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.2", "adagrad": False},
    "d10_e10_m0.15_neg100": {"description": "dim=10, epoch=10, margin=0.15, negSearch=100",
        "epoch": "10", "negSearchLimit": "100", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    # --- LONGER TRAINING + VERY LOW LR ---
    "d10_e30_m0.15_slowlr": {"description": "dim=10, epoch=30, margin=0.15, lr=0.0001",
        "epoch": "30", "negSearchLimit": "50", "dim": "10", "lr": "0.0001",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    "d50_e30_m0.15_slowlr": {"description": "dim=50, epoch=30, margin=0.15, lr=0.0001",
        "epoch": "30", "negSearchLimit": "50", "dim": "50", "lr": "0.0001",
        "minCount": "1", "margin": "0.15", "adagrad": False},
    "d10_e50_m0.1_slowlr": {"description": "dim=10, epoch=50, margin=0.1, lr=0.0001",
        "epoch": "50", "negSearchLimit": "50", "dim": "10", "lr": "0.0001",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    "d50_e50_m0.1_slowlr": {"description": "dim=50, epoch=50, margin=0.1, lr=0.0001",
        "epoch": "50", "negSearchLimit": "50", "dim": "50", "lr": "0.0001",
        "minCount": "1", "margin": "0.1", "adagrad": False},
    # --- MINCOUNT + ADAGRAD at promising combos ---
    "d10_e10_m0.15_mc5": {"description": "dim=10, epoch=10, margin=0.15, minCount=5",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "5", "margin": "0.15", "adagrad": False},
    "d10_e10_m0.2_mc5": {"description": "dim=10, epoch=10, margin=0.2, minCount=5",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "5", "margin": "0.2", "adagrad": False},
    "d10_e10_m0.15_ada": {"description": "dim=10, epoch=10, margin=0.15, adagrad",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": True},
    "d50_e10_m0.15_ada": {"description": "dim=50, epoch=10, margin=0.15, adagrad",
        "epoch": "10", "negSearchLimit": "50", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.15", "adagrad": True},
    # --- UC GROUP (universe_cc, bulk ATAC) ---
    "uc_d10_e10_m0.7": {"description": "UC: dim=10, epoch=10, margin=0.7 — margin push",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.7", "adagrad": False},
    "uc_d10_e10_m1.0": {"description": "UC: dim=10, epoch=10, margin=1.0 — harder margin",
        "epoch": "10", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "1.0", "adagrad": False},
    "uc_d10_e10_m0.7_neg100": {"description": "UC: dim=10, epoch=10, margin=0.7, neg=100",
        "epoch": "10", "negSearchLimit": "100", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.7", "adagrad": False},
    "uc_d10_e10_m1.0_neg200": {"description": "UC: dim=10, epoch=10, margin=1.0, neg=200",
        "epoch": "10", "negSearchLimit": "200", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "1.0", "adagrad": False},
    "uc_d10_e30_m0.5_ada": {"description": "UC: dim=10, epoch=30, margin=0.5, adagrad",
        "epoch": "30", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.5", "adagrad": True},
    "uc_d10_e50_m0.5_ada": {"description": "UC: dim=10, epoch=50, margin=0.5, adagrad",
        "epoch": "50", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "1", "margin": "0.5", "adagrad": True},
    "uc_d10_e20_m0.5_mc5": {"description": "UC: dim=10, epoch=20, margin=0.5, minCount=5",
        "epoch": "20", "negSearchLimit": "50", "dim": "10", "lr": "0.0005",
        "minCount": "5", "margin": "0.5", "adagrad": False},
    "uc_d10_e20_m0.7_mc5_neg100": {"description": "UC: dim=10, epoch=20, margin=0.7, mc=5, neg=100",
        "epoch": "20", "negSearchLimit": "100", "dim": "10", "lr": "0.0005",
        "minCount": "5", "margin": "0.7", "adagrad": False},
    "uc_d50_e20_m0.7_neg100": {"description": "UC: dim=50, epoch=20, margin=0.7, neg=100",
        "epoch": "20", "negSearchLimit": "100", "dim": "50", "lr": "0.0005",
        "minCount": "1", "margin": "0.7", "adagrad": False},
    "uc_d50_e30_m0.7_mc5_ada_neg100": {"description": "UC: dim=50, epoch=30, margin=0.7, mc=5, ada, neg=100",
        "epoch": "30", "negSearchLimit": "100", "dim": "50", "lr": "0.0005",
        "minCount": "5", "margin": "0.7", "adagrad": True},
}

K_RANGE = range(2, 7)


# ---- Clustering diagnostics --------------------------------------------------

def multi_k_clustering(region_embeddings: np.ndarray, k_range=K_RANGE,
                       max_samples: int = 5000) -> Dict:
    """Run KMeans for multiple k values and return silhouette + inertia per k."""
    if len(region_embeddings) > max_samples:
        print(f"    Subsampling {len(region_embeddings)} -> {max_samples} regions", flush=True)
        idx = np.random.RandomState(42).choice(len(region_embeddings), max_samples, replace=False)
        emb = region_embeddings[idx]
    else:
        emb = region_embeddings
    results = {}
    for k in k_range:
        if k >= len(emb):
            break
        print(f"    KMeans k={k}...", end=" ", flush=True)
        km = KMeans(n_clusters=k, random_state=42, n_init=10)
        labels = km.fit_predict(emb)
        sil = float(silhouette_score(emb, labels, metric="cosine"))
        print(f"silhouette={sil:.4f}", flush=True)
        results[k] = {"silhouette": sil, "inertia": float(km.inertia_)}
    return results


def create_silhouette_plot(multi_k_results: Dict, output_path: str, config_name: str):
    ks = sorted(multi_k_results.keys())
    sils = [multi_k_results[k]["silhouette"] for k in ks]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(ks, sils, color="steelblue", edgecolor="black")
    best_k = ks[int(np.argmax(sils))]
    bars[ks.index(best_k)].set_color("darkorange")
    ax.set_xlabel("k (number of clusters)"); ax.set_ylabel("Silhouette Score (cosine)")
    ax.set_title(f"{config_name}: Silhouette vs k  (best k={best_k})")
    ax.set_xticks(ks)
    plt.tight_layout(); plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()


def create_elbow_plot(multi_k_results: Dict, output_path: str, config_name: str):
    ks = sorted(multi_k_results.keys())
    inertias = [multi_k_results[k]["inertia"] for k in ks]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(ks, inertias, "o-", color="steelblue", linewidth=2, markersize=8)
    ax.set_xlabel("k (number of clusters)"); ax.set_ylabel("Inertia (within-cluster SS)")
    ax.set_title(f"{config_name}: Elbow Plot"); ax.set_xticks(ks)
    plt.tight_layout(); plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()


# ---- Preprocessing (shared across configs) ----------------------------------

def run_preprocess(data_path, metadata, universe, output_dir, labels):
    """Run preprocessing (shared across all configs)."""
    print("\n" + "=" * 60)
    print("STEP 1: PREPROCESSING")
    print("=" * 60)
    os.makedirs(output_dir, exist_ok=True)
    preprocess_main(
        data_path=data_path, metadata=metadata, universe=universe,
        output=output_dir, labels=labels, mode=MODE_BULK,
    )
    train_input_path = os.path.join(output_dir, "train_input.txt")
    if not os.path.exists(train_input_path):
        raise FileNotFoundError(f"Preprocessing failed: {train_input_path} not found")
    with open(train_input_path, "r") as f:
        n = sum(1 for _ in f)
    print(f"  Preprocessing complete: {n} samples → {train_input_path}")
    return train_input_path


# ---- Main experiment orchestration ------------------------------------------

def run_experiment(
    data_path: str, metadata: str, universe: str, labels: str,
    starspace_path: str, configs: Optional[List[str]] = None,
    experiment_name: Optional[str] = None,
    wandb_project: Optional[str] = None,
    wandb_group: Optional[str] = None,
    output_dir: Optional[str] = None,
    skip_clustering: bool = False,
    deseq_path: Optional[str] = None,
):
    universe_name = Path(universe).stem
    if experiment_name is None:
        experiment_name = f"{universe_name}_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"

    base_dir = Path(output_dir) if output_dir else EXPERIMENT_DIR
    exp_dir = base_dir / experiment_name
    exp_dir.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 70)
    print("BEDSPACE PARAMETER SWEEP")
    print("=" * 70)
    print(f"Experiment: {experiment_name}")
    print(f"Universe:   {universe_name}")
    print(f"Output:     {exp_dir}")

    configs_to_run = configs or list(SWEEP_CONFIGS.keys())
    print(f"Configs:    {configs_to_run}")

    preprocess_dir = exp_dir / "preprocessed"
    train_input = run_preprocess(data_path, metadata, universe, str(preprocess_dir), labels)

    wb_group = wandb_group or experiment_name if wandb_project else None
    results = {}

    for name in configs_to_run:
        if name not in SWEEP_CONFIGS:
            print(f"\n  WARNING: Unknown config '{name}', skipping")
            continue

        cfg = SWEEP_CONFIGS[name]
        print(f"\n{'='*60}\nConfig: {name}\n  {cfg['description']}\n{'='*60}")

        cfg_dir = exp_dir / name
        cfg_dir.mkdir(parents=True, exist_ok=True)

        wb_config = {
            "universe": universe_name,
            "dim": int(cfg["dim"]),
            "epochs": int(cfg["epoch"]),
            "lr": float(cfg["lr"]),
            "negSearchLimit": int(cfg["negSearchLimit"]),
            "minCount": int(cfg["minCount"]),
            "margin": float(cfg["margin"]),
            "adagrad": cfg["adagrad"],
        }
        star_params = {
            "dim": cfg["dim"], "epoch": cfg["epoch"],
            "negSearchLimit": cfg["negSearchLimit"], "lr": cfg["lr"],
            "minCount": cfg["minCount"], "margin": cfg["margin"],
        }
        if cfg.get("adagrad"):
            star_params["adagrad"] = "true"

        run_name = f"{universe_name}/{name}"
        tracker = WandbTracker(
            project=wandb_project, config=wb_config, run_name=run_name,
            tags=[universe_name, name, f"dim_{cfg['dim']}", f"margin_{cfg['margin']}"],
            group=wb_group, enabled=wandb_project is not None,
        )

        # ---- Train ----
        if tracker.active:
            model_path = run_train_with_wandb(
                starspace_path=starspace_path, train_input=train_input,
                output_dir=str(cfg_dir), experiment_name=name,
                params=star_params, tracker=tracker,
            )
        else:
            model_path = run_train_fallback(
                starspace_path=starspace_path, train_input=train_input,
                output_dir=str(cfg_dir), experiment_name=name, params=star_params,
            )

        if model_path is None:
            results[name] = {"status": "FAILED"}
            tracker.finish()
            continue

        # ---- Evaluate ----
        try:
            stats = evaluate_model(model_path, output_name=name,
                                   output_dir=str(cfg_dir), deseq_path=deseq_path)
        except Exception as e:
            print(f"    Evaluation failed: {e}")
            results[name] = {"status": "EVAL_FAILED"}
            tracker.finish()
            continue

        if tracker.active:
            try:
                tracker.log_evaluation(stats)
                for suffix in ["tsne", "umap", "tsne_foldchange", "umap_foldchange"]:
                    img_path = cfg_dir / f"{name}_{suffix}.png"
                    if img_path.exists():
                        img = wandb.Image(str(img_path))
                        wandb.summary[f"viz/{suffix}"] = img
                        wandb.log({f"viz/{suffix}": img})
                if "label_separation" in stats:
                    wandb.summary["label_sep_mean"] = (
                        stats["label_separation"].get("mean_pairwise_distance")
                    )
            except Exception as e:
                print(f"    wandb image logging failed: {e}", flush=True)

        # ---- Multi-k clustering diagnostic ----
        if not skip_clustering:
            try:
                print("    Running multi-k clustering...", flush=True)
                region_emb, _, _, _ = load_model_embeddings(model_path)
                print(f"    Loaded {len(region_emb)} region embeddings for clustering", flush=True)
                mk = multi_k_clustering(region_emb)
                best_k = max(mk, key=lambda k: mk[k]["silhouette"])
                stats["multi_k"] = mk
                stats["best_k"] = best_k
                stats["silhouette_at_k3"] = mk.get(3, {}).get("silhouette")
                sil_summary = ", ".join("k=%d: %.3f" % (k, v["silhouette"])
                                        for k, v in sorted(mk.items()))
                print(f"    Multi-k silhouette: {sil_summary}", flush=True)
                print(f"    Best k={best_k} (silhouette={mk[best_k]['silhouette']:.3f})", flush=True)

                sil_path = str(cfg_dir / f"{name}_silhouette_vs_k.png")
                elbow_path = str(cfg_dir / f"{name}_elbow.png")
                create_silhouette_plot(mk, sil_path, name)
                create_elbow_plot(mk, elbow_path, name)

                if tracker.active:
                    for k, v in mk.items():
                        tracker.log_evaluation({f"cluster_k{k}": v})
                    wandb.summary["best_k"] = best_k
                    wandb.summary["silhouette_at_k3"] = stats["silhouette_at_k3"]
                    wandb.summary["silhouette_best"] = mk[best_k]["silhouette"]
                    for img_name, img_file in [("viz/silhouette_vs_k", sil_path),
                                                ("viz/elbow", elbow_path)]:
                        if os.path.exists(img_file):
                            img = wandb.Image(img_file)
                            wandb.summary[img_name] = img
                            wandb.log({img_name: img})
            except Exception as e:
                print(f"    Multi-k diagnostic failed: {e}", flush=True)

        results[name] = {
            "status": "SUCCESS", "config": wb_config,
            "model_path": model_path, "stats": stats,
        }
        tracker.finish()

    print_comparison(results, configs_to_run, exp_dir)
    return results


def print_comparison(results, config_names, exp_dir):
    """Print and save a comparison table."""
    lines = ["\n" + "=" * 100, "PARAMETER SWEEP COMPARISON", "=" * 100]
    header = (f"{'Config':<25} {'ep':<5} {'dim':<5} {'lr':<9} {'margin':<8} {'mc':<4} {'ada':<5} "
              f"{'best_k':<7} {'sil@3':<8} {'sil_best':<10} {'label_sep':<10} {'status':<8}")
    lines.append(header); lines.append("-" * 120)
    for name in config_names:
        r = results.get(name)
        if r is None or r["status"] != "SUCCESS":
            lines.append(f"{name:<25} {'':88} {r['status'] if r else 'SKIP':<8}")
            continue
        c = r["config"]
        s = r.get("stats", {})
        best_k = s.get("best_k", "?")
        sil3 = s.get("silhouette_at_k3")
        sil3_s = f"{sil3:.4f}" if sil3 is not None else "N/A"
        mk = s.get("multi_k", {})
        sil_best = mk.get(best_k, {}).get("silhouette")
        sil_best_s = f"{sil_best:.4f}" if sil_best is not None else "N/A"
        lsep = s.get("label_separation", {}).get("mean_pairwise_distance")
        lsep_s = f"{lsep:.4f}" if lsep is not None else "N/A"
        lines.append(
            f"{name:<25} {c['epochs']:<5} {c['dim']:<5} {c['lr']:<9} {c['margin']:<8} "
            f"{c['minCount']:<4} {str(c['adagrad']):<5} {str(best_k):<7} {sil3_s:<8} "
            f"{sil_best_s:<10} {lsep_s:<10} SUCCESS"
        )
    lines.append("-" * 100)
    three_cluster = [n for n in config_names
                     if results.get(n, {}).get("status") == "SUCCESS"
                     and results[n].get("stats", {}).get("best_k") == 3]
    if three_cluster:
        ranked = sorted(three_cluster,
                        key=lambda n: results[n]["stats"].get("silhouette_at_k3") or 0,
                        reverse=True)
        lines.append(f"\nConfigs with best_k=3: {', '.join(three_cluster)}")
        lines.append(f"Best 3-cluster config: {ranked[0]}")
    else:
        lines.append("\nNo config achieved best_k=3. Inspect silhouette plots.")
    lines.append("=" * 100)
    report = "\n".join(lines)
    print(report)
    report_path = exp_dir / "sweep_comparison.txt"
    with open(report_path, "w") as f:
        f.write(report)
    print(f"\nSaved to {report_path}")


# =============================================================================
# CLI
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Sweep bedspace params for region embeddings (bulk ATAC-seq).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--data-path", "-d", required=False, help="BED files directory")
    parser.add_argument("--metadata", "-m", required=False, help="Metadata CSV")
    parser.add_argument("--universe", "-u", required=False, help="Universe BED file")
    parser.add_argument("--labels", "-l", required=False, help="Label column in metadata")
    parser.add_argument("--starspace", "-s", default=DEFAULT_STARSPACE_PATH,
                        help="StarSpace directory")
    parser.add_argument("--configs", "-c", nargs="+", default=None, help="Subset of configs")
    parser.add_argument("--name", "-n", default=None, help="Experiment name")
    parser.add_argument("--output-dir", "-o", default=None,
                        help="Base output directory (default: BEDSPACE_OUTPUT_DIR/experiments_3cluster)")
    parser.add_argument("--wandb-project", default=None, help="wandb project (None=disabled)")
    parser.add_argument("--wandb-group", default=None, help="wandb group")
    parser.add_argument("--skip-clustering", action="store_true",
                        help="Skip multi-k clustering diagnostic")
    parser.add_argument("--deseq", default=None,
                        help="Path to DESeq result file (CSV with seqnames, start, end, log2FoldChange).")
    parser.add_argument("--list-configs", action="store_true",
                        help="List available configs and exit")
    args = parser.parse_args()

    if args.list_configs:
        print("\nAvailable configs:")
        print("-" * 70)
        for name, cfg in SWEEP_CONFIGS.items():
            print(f"  {name:<30} lr={cfg['lr']:<6} margin={cfg['margin']:<5} "
                  f"mc={cfg['minCount']:<2} ada={cfg['adagrad']!s:<5}  {cfg['description']}")
        return 0

    # Required args when actually running:
    for path, label, kind in [
        (args.data_path, "Data path", "dir"),
        (args.starspace, "StarSpace", "dir"),
        (args.metadata, "Metadata", "file"),
        (args.universe, "Universe", "file"),
    ]:
        if path is None:
            print(f"Error: {label} (--{label.lower().replace(' ', '-')}) is required.")
            return 1
        if kind == "dir" and not os.path.isdir(path):
            print(f"Error: {label} not found (or not a directory): {path}")
            return 1
        if kind == "file" and not os.path.isfile(path):
            print(f"Error: {label} not found: {path}")
            return 1
    if args.labels is None:
        print("Error: --labels is required.")
        return 1

    run_experiment(
        data_path=args.data_path, metadata=args.metadata, universe=args.universe,
        labels=args.labels, starspace_path=args.starspace, configs=args.configs,
        experiment_name=args.name, wandb_project=args.wandb_project,
        wandb_group=args.wandb_group, output_dir=args.output_dir,
        skip_clustering=args.skip_clustering, deseq_path=args.deseq,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

