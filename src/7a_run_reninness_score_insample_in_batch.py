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
    python 7a_run_reninness_score_insample_in_batch.py \\
        --data-path  /path/to/bed/files/ \\
        --metadata   /path/to/metadata.csv \\
        --universe   /path/to/universe.bed \\
        --labels     "label" \\
        --starspace  /path/to/Starspace/ \\
        --output-dir /path/to/output/ \\
        --wandb-project bedspace-3cluster        # optional

    # Run a subset of configs:
    python 7a_run_reninness_score_insample_in_batch.py --configs d3_e10_m0.1 d10_e10_m1.0 ...

    # List available configs and exit:
    python 7a_run_reninness_score_insample_in_batch.py --list-configs

Output (under <output-dir>/<experiment_name>/):
    preprocessed/                          (shared across all configs)
    <config_name>/
        model_<config>.tsv                 trained StarSpace model
        <config>_stats.json                statistical metrics
        <config>_report.txt                summary report
        <config>_tsne.png, _umap.png       embedding visualizations
        <config>_silhouette_vs_k.png       multi-k clustering diagnostic
        <config>_elbow.png
    <prefix>_comparison.txt                side-by-side comparison table

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
import shutil
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


def _fit_projection(region_embeddings, label_embeddings, fit_fn, n_plot,
                    random_state=42, fit_max_regions=None):
    """Fit a 2D projection on ALL regions, return coordinates for a subsample.

    FITTING AND PLOTTING ARE SEPARATE. `fit_fn` sees every region (or
    `fit_max_regions` of them, for a cheaper fit); only `n_plot` regions come
    back to be drawn. Thinning the plot therefore cannot change the topology.

    Under the old code one `max_regions` did both jobs, so raising it re-fit
    the projection on a different point cloud -- that, not sampling noise, is
    why 7f/7g saw the apparent islands at 5,000 regions merge into one
    connected sheet at 50,000. With a full-vocabulary fit the layout is a
    property of the embedding and the plotted count is purely cosmetic.

    The fitted and plotted subsets are prefixes of the SAME seeded
    permutation, so they nest and the plotted regions are always inside the
    fitted ones; `n_plot > fit_max_regions` is rejected rather than silently
    drawing points that were never fitted. Setting `fit_max_regions == n_plot`
    reproduces the old behaviour exactly (same regions, same row order into
    the fit, so identical coordinates).

    COST. A full fit is ~400k points, not 5k, and `umap.UMAP(random_state=...)`
    forces `n_jobs=1`. Budget tens of minutes per UMAP and rather more per
    t-SNE; pass `fit_max_regions` for a cheaper fit that is still independent
    of the plotted count.

    MEASURED, 2026-08-28: at the full 422,488-region vocabulary one config's
    t-SNE had not finished after 10 minutes, against 3.2 min for a WHOLE config
    in the 2026-08-22 run -- which was fitting on 5,003 points, because the
    full-vocabulary fit landed after it. A per-config rate taken from that run
    does not transfer to this code path.

    fit_max_regions : int or float
        > 1  -- an absolute number of regions.
        0 < x <= 1 -- a FRACTION of the vocabulary, resolved per config. 0.1 is
            ~42k regions on the universe_cc_073126 universe. A fraction is the
            safer thing to pass across configs whose vocabularies differ
            (minCount > 1 drops rare regions), since it keeps the fit a
            constant share of each model rather than a constant count that
            might exceed a smaller vocabulary.
        None -- the whole vocabulary.

    NOTHING IN THE SELECTION DEPENDS ON THIS. t-SNE and UMAP here produce
    per-config PNGs only. `sil_k2` is computed by multi_k_clustering on the
    full embeddings, `ctt`/`gdst` by geniml.eval on the full embeddings, and
    `nl_auprc` from raw_cosdist_rl.csv -- none of them see a projection. So the
    cap changes the pictures and no number in the summary table or Fig1.

    fit_fn : callable(X) -> array of shape (len(X), 2).

    Returns (region_coords_plot, label_coords, plot_idx, n_fit).
    """
    n_total = len(region_embeddings)
    n_plot = min(int(n_plot), n_total)
    order = np.random.RandomState(random_state).permutation(n_total)

    # A value in (0, 1] is a FRACTION of this config's vocabulary. Resolved
    # here rather than by the caller so every projection in the run uses the
    # same share of whatever vocabulary its own model has.
    if fit_max_regions is not None and 0 < float(fit_max_regions) <= 1:
        frac = float(fit_max_regions)
        fit_max_regions = max(int(round(frac * n_total)), n_plot)
        print(f"    fit on {frac:.0%} of the vocabulary = "
              f"{fit_max_regions} of {n_total} regions", flush=True)

    if fit_max_regions is None or fit_max_regions >= n_total:
        fit_idx = np.arange(n_total)          # full vocabulary, file order
    else:
        if n_plot > fit_max_regions:
            raise ValueError(
                f"plotted regions ({n_plot}) > fit_max_regions "
                f"({fit_max_regions}): cannot draw regions that were never "
                "fitted. Raise fit_max_regions, or leave it None to fit the "
                "whole vocabulary.")
        fit_idx = order[:fit_max_regions]
    n_fit = len(fit_idx)

    pos = np.full(n_total, -1, dtype=np.int64)   # universe index -> fitted row
    pos[fit_idx] = np.arange(n_fit)
    plot_idx = order[:n_plot]

    coords = np.asarray(fit_fn(np.vstack([region_embeddings[fit_idx],
                                          label_embeddings])))
    return coords[:n_fit][pos[plot_idx]], coords[n_fit:], plot_idx, n_fit


def create_tsne_visualization(region_embeddings, label_embeddings, label_names,
                              output_path, perplexity=30, random_state=42,
                              max_regions=5000, fit_max_regions=None):
    """t-SNE: regions = grey circles, labels = colored triangles.

    Fit on the whole vocabulary, draw `max_regions` of it -- see
    `_fit_projection`. The plotted subsample is now seeded, so repeated calls
    draw the same regions (it used to be the unseeded global RNG).
    """
    def _fit(X):
        print(f"Running t-SNE on {len(X)} embeddings...")
        return TSNE(n_components=2, perplexity=min(perplexity, len(X) - 1),
                    random_state=random_state, metric="cosine").fit_transform(X)

    region_coords, label_coords, _, n_fit = _fit_projection(
        region_embeddings, label_embeddings, _fit, max_regions,
        random_state=random_state, fit_max_regions=fit_max_regions)

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.scatter(region_coords[:, 0], region_coords[:, 1], c="lightgrey",
               s=10, alpha=0.5, marker="o", label="Regions")
    _plot_label_markers(ax, label_coords, label_names)
    ax.set_xlabel("t-SNE 1"); ax.set_ylabel("t-SNE 2")
    ax.set_title("t-SNE: Regions (circles) vs Labels (triangles)\n"
                 f"drawn {len(region_coords):,} of {n_fit:,} fitted",
                 fontsize=11)
    ax.legend(loc="best", fontsize=10)
    plt.tight_layout(); plt.savefig(output_path, dpi=150, bbox_inches="tight"); plt.close()
    print(f"Saved t-SNE plot to {output_path}")


def create_umap_visualization(region_embeddings, label_embeddings, label_names,
                              output_path, n_neighbors=15, min_dist=0.1,
                              random_state=42, max_regions=None, frac=0.10,
                              alpha=0.10, point_size=10, fit_max_regions=None):
    """UMAP: regions = grey circles, labels = colored triangles.

    FIT vs PLOT. The projection is fit on the FULL vocabulary
    (`fit_max_regions=None`); `frac` / `max_regions` only choose how many
    regions are drawn. See `_fit_projection` -- this is what makes the layout
    independent of the plotted count.

    SUBSAMPLE SIZE. `frac` (default 10% of the vocabulary) replaces the old
    fixed max_regions=5000. At 5,000 of 350k-420k regions the plot showed 1-2%
    of the data, and on a sparse random sample of a continuous manifold UMAP
    fragments the projection into apparent islands -- gaps that denser sampling
    fills in. Established three ways: the max_regions grid (7f), the 10%
    re-plot (7i), and a neighbourhood-fraction-matched test (7j). Those three
    all varied the FIT; under a full-vocabulary fit the fragmentation question
    no longer depends on this knob at all, and `frac` is a readability choice.
    Pass max_regions explicitly to override `frac`.

    ALPHA 0.10, down from 0.5: at ~40k points the old value saturates every
    group into solid fill and density stops being readable.

    SEEDED SUBSAMPLE: RandomState(random_state) rather than the unseeded global
    RNG, so repeated calls plot the SAME regions and panels are comparable. The
    previous behaviour re-drew a different sample every call.

    NOTE: n_neighbors is a COUNT, not a fraction, and it now spans the full
    vocabulary rather than the subsample, so a fixed 15 is a much more local
    neighbourhood than it was at 5,000 regions. Left at 15 for continuity.
    """
    if not HAS_UMAP:
        print("UMAP not available, skipping UMAP visualization")
        return
    n_plot = max_regions if max_regions is not None else int(round(frac * len(region_embeddings)))

    def _fit(X):
        print(f"Running UMAP on {len(X)} embeddings...")
        return umap.UMAP(n_neighbors=min(n_neighbors, len(X) - 1),
                         min_dist=min_dist, metric="cosine",
                         random_state=random_state).fit_transform(X)

    region_coords, label_coords, _, n_fit = _fit_projection(
        region_embeddings, label_embeddings, _fit, n_plot,
        random_state=random_state, fit_max_regions=fit_max_regions)

    fig, ax = plt.subplots(figsize=(12, 10))
    ax.scatter(region_coords[:, 0], region_coords[:, 1], c="lightgrey",
               s=point_size, alpha=alpha, marker="o", label="Regions",
               linewidths=0, rasterized=True)
    _plot_label_markers(ax, label_coords, label_names)
    ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    ax.set_title(f"UMAP: Regions (circles) vs Labels (triangles)\n"
                 f"fit on {n_fit:,} of {len(region_embeddings):,} regions, "
                 f"drawn n={len(region_coords):,} "
                 f"({100*len(region_coords)/len(region_embeddings):.0f}%), "
                 f"n_neighbors={n_neighbors}, min_dist={min_dist}, alpha={alpha}",
                 fontsize=11)
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
                           perplexity=30, random_state=42, max_regions=5000,
                           fit_max_regions=None):
    """t-SNE with regions colored by log2FoldChange.

    Fit on the whole vocabulary, color and draw `max_regions` of it
    (`_fit_projection`). The plotted subsample is seeded now, so this figure
    and `create_tsne_visualization` show the same regions.
    """
    def _fit(X):
        print(f"Running t-SNE (foldchange) on {len(X)} embeddings...")
        return TSNE(n_components=2, perplexity=min(perplexity, len(X) - 1),
                    random_state=random_state, metric="cosine").fit_transform(X)

    region_coords, label_coords, plot_idx, n_fit = _fit_projection(
        region_embeddings, label_embeddings, _fit, max_regions,
        random_state=random_state, fit_max_regions=fit_max_regions)

    fc_values = np.array([foldchange_map.get(region_names[i], 0.0) for i in plot_idx])
    _foldchange_plot(region_coords, fc_values, label_coords, label_names,
                     output_path,
                     "t-SNE: Regions colored by log2FoldChange\n"
                     f"drawn {len(region_coords):,} of {n_fit:,} fitted")


def create_umap_foldchange(region_embeddings, region_names, foldchange_map,
                           label_embeddings, label_names, output_path,
                           n_neighbors=15, min_dist=0.1, random_state=42,
                           max_regions=5000, fit_max_regions=None):
    """UMAP with regions colored by log2FoldChange.

    Fit on the whole vocabulary, color and draw `max_regions` of it
    (`_fit_projection`). The plotted subsample is seeded now, so this figure
    and `create_umap_visualization` show the same regions.
    """
    if not HAS_UMAP:
        print("UMAP not available, skipping UMAP foldchange visualization")
        return

    def _fit(X):
        print(f"Running UMAP (foldchange) on {len(X)} embeddings...")
        return umap.UMAP(n_neighbors=min(n_neighbors, len(X) - 1),
                         min_dist=min_dist, metric="cosine",
                         random_state=random_state).fit_transform(X)

    region_coords, label_coords, plot_idx, n_fit = _fit_projection(
        region_embeddings, label_embeddings, _fit, max_regions,
        random_state=random_state, fit_max_regions=fit_max_regions)

    fc_values = np.array([foldchange_map.get(region_names[i], 0.0) for i in plot_idx])
    _foldchange_plot(region_coords, fc_values, label_coords, label_names,
                     output_path,
                     "UMAP: Regions colored by log2FoldChange\n"
                     f"drawn {len(region_coords):,} of {n_fit:,} fitted")


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
                   output_dir=None, deseq_path=None, fit_max_regions=None):
    """Full evaluation on a trained StarSpace model.

    `fit_max_regions` caps how many regions the t-SNE/UMAP figures FIT on;
    None (the default) fits the whole vocabulary. The plotted counts are
    separate knobs -- see `_fit_projection`. A full fit is ~400k points and
    single threaded, and the t-SNE one is the expensive half; set this if the
    figures are holding up an unattended sweep.
    """
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
    create_tsne_visualization(region_embeddings, label_embeddings, label_names, tsne_path,
                              fit_max_regions=fit_max_regions)
    umap_path = out_dir / f"{output_name}_umap.png"
    create_umap_visualization(region_embeddings, label_embeddings, label_names, umap_path,
                              fit_max_regions=fit_max_regions)

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
                                   out_dir / f"{output_name}_tsne_foldchange.png",
                                   fit_max_regions=fit_max_regions)
            create_umap_foldchange(region_embeddings, region_names, foldchange_map,
                                   label_embeddings, label_names,
                                   out_dir / f"{output_name}_umap_foldchange.png",
                                   fit_max_regions=fit_max_regions)
        except Exception as e:
            print(f"    WARNING: DESeq foldchange plots failed: {e}")

    print(f"\n{'='*60}\nEvaluation complete!\n{'='*60}")
    return stats


# =============================================================================
# SECTION 2b — Shared evaluation metrics
#
# Defined here because 7b (LOOCV) and 7c (holdout) already load this file via
# importlib to read SWEEP_CONFIGS; keeping the metric definitions in one place
# means the three drivers cannot drift apart on what "accuracy" means.
#
# TWO accuracy families and TWO clusterability families, matching
# geniml_dev/geniml/bedspace/tests/{analyze_sweep,compute_clusterability}.py:
#
#   accuracy        acc_sign / bal_acc  -- true label = file_label, predicted =
#                                          sign of score (notebook convention)
#                   prec_rank_*         -- threshold-free precision at rank
#   clusterability  silhouette          -- multi-k KMeans (see multi_k_clustering)
#                   ctt / gdst / npt / rct -- geniml.eval statistical tests
# =============================================================================

LABEL_PREFIX = "__label__"


def precision_at_rank(y, s):
    """Precision-at-rank accuracy: grade the ranking against the true class sizes.

    Sign-based accuracy asks "is the score on the right side of 0?", which
    depends on where the score's zero point lands. This drops the threshold and
    uses the class sizes already known from `file_label`:

        x = number of true Renin, y = number of true Non-Renin

    Rank by score descending: the top x are what the score claims are Renin, the
    bottom y what it claims are Non-Renin. Then count how many actually are.
    Because the two blocks exactly partition the samples, prec_rank_renin is
    simultaneously the recall and the F1 of the Renin class at that cut.

    Ties are broken by descending score then original order, so the number is
    reproducible; a tie straddling the x boundary is resolved consistently.
    """
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    order = np.lexsort((np.arange(len(s)), -s))
    y_sorted = y[order]
    x = int(y.sum())
    n = len(y)
    n_neg = n - x
    if x == 0 or n_neg == 0:
        return dict(prec_rank_renin=np.nan, prec_rank_nonrenin=np.nan,
                    prec_rank_overall=np.nan, n_rank_renin=x, n_rank_nonrenin=n_neg)
    hit_pos = int(y_sorted[:x].sum())
    hit_neg = int((y_sorted[x:] == 0).sum())
    return dict(prec_rank_renin=hit_pos / x,
                prec_rank_nonrenin=hit_neg / n_neg,
                prec_rank_overall=(hit_pos + hit_neg) / n,
                n_rank_renin=x, n_rank_nonrenin=n_neg)


def score_table_metrics(path, target="Renin"):
    """Both accuracy families for one reninness score CSV.

    Anchor rows (filename == file_label, the label-to-label rows that
    `--include-label-distances` adds) are dropped first, per the notebook
    convention. Returns None when the table has only one class, which is the
    normal case for a single LOOCV fold.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score

    df = pd.read_csv(path)
    df = df[df["filename"].astype(str) != df["file_label"].astype(str)].copy()
    if not len(df):
        return None
    y = (df["file_label"].astype(str) == target).astype(int).to_numpy()
    s = df["score"].to_numpy(dtype=float)
    if y.sum() == 0 or y.sum() == len(y):
        return None

    tp = int(((y == 1) & (s > 0)).sum()); fn = int(((y == 1) & (s <= 0)).sum())
    tn = int(((y == 0) & (s < 0)).sum()); fp = int(((y == 0) & (s >= 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    return dict(
        n_samples=len(df), n_pos=int(y.sum()), n_neg=int((y == 0).sum()),
        acc_sign=float((((y == 1) & (s > 0)) | ((y == 0) & (s < 0))).mean()),
        bal_acc=float((sens + spec) / 2), sens=sens, spec=spec,
        tp=tp, fn=fn, tn=tn, fp=fp,
        auc=float(roc_auc_score(y, s)),
        auprc=float(average_precision_score(y, s)),
        **precision_at_rank(y, s),
    )


def region_to_geniml(token: str) -> str:
    """bedspace `chr1_100_200` -> geniml.eval `chr1:100-200`.

    geniml.eval parses regions as chr:start-end; bedspace tokenizes to
    underscore-joined triples. rsplit keeps chromosome names that themselves
    contain underscores (chrUn_GL456239) intact.
    """
    chrom, start, end = token.rsplit("_", 2)
    return f"{chrom}:{start}-{end}"


def bedspace_model_to_base_embeddings(tsv: str, out_path: str) -> str:
    """Cache a bedspace model.tsv as a geniml.eval "base" BaseEmbeddings pickle.

    Drops the trailing `__label__*` rows: those are label embeddings, not
    regions, and would corrupt every region-level test.
    """
    import pickle

    if os.path.exists(out_path):
        return out_path
    from geniml.eval.utils import BaseEmbeddings

    df = pd.read_csv(tsv, sep="\t", header=None)
    names = df[0].astype(str)
    is_region = ~names.str.startswith(LABEL_PREFIX)
    vocab = [region_to_geniml(r) for r in names[is_region]]
    emb = df.loc[is_region.values, 1:].to_numpy(dtype=float)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        pickle.dump(BaseEmbeddings(emb, vocab), f)
    return out_path


def geniml_eval_metrics(model_tsv, cache_path, seed=42, num_data=10000,
                        num_workers=10, npt_k=50, npt_samples=100,
                        npt_resolution=10, bin_embed=None, rct_out_dim=64,
                        rct_cv=5, keep_base=False):
    """CTT / GDST / NPT (+ optional RCT) for one bedspace model.

    CTT is the clusterability headline (cluster tendency in [0, 1]; 0.5 =
    uniformly distributed). GDST grades genome-distance agreement and NPT
    neighborhood preservation. RCT needs a binary embedding of the same universe
    and trains an MLP per CV fold, so it only runs when `bin_embed` is given.

    Each test is guarded independently: a model that is degenerate for one test
    should not cost the others.
    """
    from geniml.eval.ctt import get_ctt_score
    from geniml.eval.gdst import get_gdst_score
    from geniml.eval.npt import get_npt_score

    base = bedspace_model_to_base_embeddings(model_tsv, cache_path)
    out = {}
    try:
        out["ctt"] = float(get_ctt_score(base, "base", seed=seed,
                                         num_data=num_data, num_workers=num_workers))
    except Exception as e:
        print(f"    CTT failed: {e}", flush=True)
        out["ctt"] = np.nan
    try:
        out["gdst"] = float(get_gdst_score(base, "base", num_samples=num_data,
                                           seed=seed))
    except Exception as e:
        print(f"    GDST failed: {e}", flush=True)
        out["gdst"] = np.nan
    try:
        r = get_npt_score(base, "base", K=npt_k, num_samples=npt_samples,
                          seed=seed, resolution=npt_resolution,
                          num_workers=num_workers)
        snpr = np.asarray(r["SNPR"], dtype=float).ravel()
        out["npt_snpr_mean"] = float(np.nanmean(snpr))
        out["npt_snpr_at_k"] = float(snpr[-1])
        out["npt_k"] = int(r["K"])
    except Exception as e:
        print(f"    NPT failed: {e}", flush=True)
        out["npt_snpr_mean"] = out["npt_snpr_at_k"] = np.nan
        out["npt_k"] = npt_k
    if bin_embed:
        try:
            from geniml.eval.rct import get_rct_score

            out["rct"] = float(get_rct_score(base, "base", bin_embed,
                                             out_dim=rct_out_dim, cv_num=rct_cv,
                                             seed=seed, num_workers=num_workers))
        except Exception as e:
            print(f"    RCT failed: {e}", flush=True)
            out["rct"] = np.nan
    # The converted pickle is the full region matrix again (~675 MB at dim=200),
    # so it is deleted once the tests have read it.
    if not keep_base:
        try:
            os.remove(base)
        except OSError:
            pass
    return out


# --- Scoring + figures, shared by 7a / 7b / 7c -------------------------------
# The four scoring methods, named exactly as the sweep's analyze_sweep.py and
# plot_param_*.py expect. `weighted_multi_learned_post` is the CURRENT learned
# variant: the two weights are fit AFTER `distances`, by `score --learn-weights`.
SCORING_METHODS = ["pairwise", "weighted_multi_uniform",
                   "weighted_multi_learned_post", "projection", "line"]

# NOTE ordering: `projection` is inserted BEFORE `line` deliberately, so the
# resume sentinel (SCORING_METHODS[-1]) stays "line" and existing finished runs
# are not re-triggered on a re-submit.
# The five above are what a THREE-label model supports. A two-label model
# (e.g. Group = Renin / Non_Renin) supports only `pairwise`: `weighted_multi`
# raises "needs >=2 negatives", and `line` needs two negatives to define the
# line. Restrict the list with `--methods pairwise` in that case -- not merely
# to silence two guaranteed failures per config, but because the resume
# sentinel below is SCORING_METHODS[-1], and a sentinel naming a method that
# is never written would make every config look unfinished on a re-submit.
ALL_SCORING_METHODS = tuple(SCORING_METHODS)


def _run(cmd):
    """Run a subprocess, stream output, raise on failure."""
    print(f"    $ {' '.join(str(c) for c in cmd)}", flush=True)
    res = subprocess.run([str(c) for c in cmd], check=False)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {' '.join(map(str, cmd))}")


def run_distances_and_score(
    model_prefix, cfg_dir, starspace, universe, data_path, label_col,
    metadata_train, metadata_test, config_name, target="Renin",
    negatives="Control,Tumoral", line_labels="Control,Tumoral",
    weights_from=None, weight_constraint="nonneg-scaled", class_weight="balanced",
    analytic=True,
):
    """`distances` once, then `score` for all four methods, into `cfg_dir`.

    `weights_from` selects where the two weighted_multi weights come from:
      None      -- learn them from THIS table (`score --learn-weights`). Correct
                   when train and test are the same samples, as in 7a.
      <path>    -- a weights.json fit elsewhere (7c fits on train, applies here).

    Truth for the weight fit is the binary Renin / Non-Renin label, derived from
    the Label_1 / Label_2 columns the distance table carries -- the same rule
    that sets `file_label` in the score output (Renin if EITHER label is Renin).
    """
    cfg_dir = Path(cfg_dir)
    cfg_dir.mkdir(parents=True, exist_ok=True)
    dist_cmd = [
        "geniml", "bedspace", "distances", "-i", model_prefix, "-s", starspace,
        "--metadata-train", metadata_train, "--metadata-test", metadata_test,
        "-u", universe, "-p", config_name, "-f", data_path, "-l", label_col,
        "-o", str(cfg_dir), "--include-label-distances",
    ]
    if analytic:
        dist_cmd.append("--analytic")
    if line_labels:
        dist_cmd += ["--line-labels", line_labels]
    _run(dist_cmd)

    dist = cfg_dir / "raw_cosdist_rl.csv"
    wm = ["--target", target, "--method", "weighted_multi"]
    if negatives:
        wm += ["--negatives", negatives]

    out = {}
    failed = []
    for method in SCORING_METHODS:
        dest = cfg_dir / f"reninness_score_{method}.csv"
        cmd = ["geniml", "bedspace", "score", "-i", str(dist), "-o", str(dest)]
        if method == "pairwise":
            cmd += ["--target", target, "--method", "pairwise"]
            if negatives:
                cmd += ["--negatives", negatives.split(",")[0]]
        elif method == "line":
            if not line_labels:
                continue
            cmd += ["--target", target, "--method", "line"]
        elif method == "weighted_multi_uniform":
            cmd += wm + ["--weight-mode", "uniform"]
        elif method == "projection":
            # Signed projection onto the (negative-line -> target) direction.
            # Needs both negatives (they define the line) but no __line__
            # column, so it works wherever weighted_multi does.
            if not negatives:
                continue
            cmd += ["--target", target, "--method", "projection",
                    "--negatives", negatives]
        else:
            # Sentinel from a caller whose out-of-sample weight fit failed:
            # skip this method rather than silently learning in-sample.
            if str(weights_from) == "SKIP":
                print("    weighted_multi_learned_post skipped "
                      "(no out-of-sample weights available)", flush=True)
                continue
            cmd += wm
            if weights_from:
                cmd += ["--weights-file", str(weights_from)]
            else:
                cmd += ["--learn-weights",
                        "--weight-constraint", weight_constraint,
                        "--class-weight", class_weight,
                        "--weights-out", str(cfg_dir / "weights_learned_post.json"),
                        "--subscores-out", str(cfg_dir / "subscores.csv")]
        # A method that cannot be scored must not cost the others. The usual
        # case is a degenerate model (e.g. dim=3, epoch=1, lr=1e-5, which barely
        # trains): `--learn-weights` then fails loudly because NNLS returns
        # all-zero weights, since no negative's sub-score is positively
        # associated with the target. That is a real signal about THAT model,
        # not a reason to lose its pairwise / uniform / line scores too. The
        # missing method simply has no row in the metric CSV.
        try:
            _run(cmd)
            out[method] = dest
        except Exception as e:
            print(f"    scoring method '{method}' failed: {e}", flush=True)
            failed.append(method)
    if failed:
        print(f"    NOTE: {len(failed)}/{len(SCORING_METHODS)} method(s) failed "
              f"for this config: {failed}; kept {sorted(out)}", flush=True)
    if not out:
        raise RuntimeError(
            f"every scoring method failed for this config: {failed}"
        )
    return out


def params_row(cfg):
    """Param columns in the schema the figure scripts colour by."""
    v = SWEEP_CONFIGS.get(cfg, {})
    if not v:
        return dict(dim=np.nan, epoch=np.nan, neg=np.nan, lr=np.nan, mc=np.nan,
                    margin=np.nan, ada=False, in_sweep=False)
    return dict(dim=int(v["dim"]), epoch=int(v["epoch"]),
                neg=int(v["negSearchLimit"]), lr=float(v["lr"]),
                mc=int(v["minCount"]), margin=float(v["margin"]),
                ada=bool(v["adagrad"]), in_sweep=True)


def collect_config_metrics(cfg, cfg_dir, target="Renin", silhouette=True,
                           geniml_eval=True, geniml_eval_workers=10,
                           bin_embed=None):
    """(accuracy rows per scoring method, clusterability row) for one config."""
    cfg_dir = Path(cfg_dir)
    base = dict(config=cfg, **params_row(cfg))
    scored = []
    for method in SCORING_METHODS:
        p = cfg_dir / f"reninness_score_{method}.csv"
        if not p.exists():
            continue
        try:
            m = score_table_metrics(str(p), target=target)
        except Exception as e:
            print(f"    metrics failed for {cfg}/{method}: {e}", flush=True)
            continue
        if m is None:
            continue
        scored.append(dict(base, method=method, **m))

    crow = {"config": cfg}
    model_tsv = next(iter(sorted(cfg_dir.glob("*.tsv"))), None)
    if model_tsv is not None:
        if silhouette:
            try:
                emb, _, _, _ = load_model_embeddings(str(model_tsv))
                mk = multi_k_clustering(emb)
                for k, v in mk.items():
                    crow[f"sil_k{k}"] = v["silhouette"]
                    crow[f"inertia_k{k}"] = v["inertia"]
                crow["best_k"] = max(mk, key=lambda k: mk[k]["silhouette"])
                crow["sil_best"] = max(v["silhouette"] for v in mk.values())
                crow["sil_at3"] = mk.get(3, {}).get("silhouette")
            except Exception as e:
                print(f"    silhouette failed for {cfg}: {e}", flush=True)
        if geniml_eval:
            try:
                crow.update(geniml_eval_metrics(
                    str(model_tsv), cache_path=str(cfg_dir / f"{cfg}_base_embed.pt"),
                    num_workers=geniml_eval_workers, bin_embed=bin_embed))
            except Exception as e:
                print(f"    geniml.eval failed for {cfg}: {e}", flush=True)
    return scored, crow


def write_run_summary(out_dir, scored_rows, cluster_rows, prefix="run"):
    """The two tidy CSVs every figure script consumes."""
    out_dir = Path(out_dir)
    s_path = out_dir / f"{prefix}_scored.csv"
    c_path = out_dir / f"{prefix}_clusterability.csv"
    pd.DataFrame(scored_rows).to_csv(s_path, index=False)
    pd.DataFrame(cluster_rows).to_csv(c_path, index=False)
    print(f"\nwrote {s_path}  ({len(scored_rows)} config x method rows)")
    print(f"wrote {c_path}  ({len(cluster_rows)} configs)")
    return s_path, c_path


def make_figures(scored_csv, cluster_csv, prefix, tests_dir=None, design="ofat10"):
    """The 12 line figures + 8 scatter figures, from the sweep's own scripts.

    Invoked as subprocesses rather than reimplemented, so a run driven from here
    and one driven from geniml_dev/.../tests produce byte-identical figures.
    """
    tests_dir = tests_dir or BEDSPACE_TESTS_DIR
    line = os.path.join(tests_dir, "plot_param_figure.py")
    scatter = os.path.join(tests_dir, "plot_param_scatter.py")
    # The joint-only scatter is emitted alongside the combined one: mixing the
    # OAT configs back in re-dilutes the colour gradient the joint design exists
    # to provide (most OAT points sit at the baseline on any axis but the swept
    # one), so both views are worth having.
    #
    # The line figures come in TWO flavours, for the same reason:
    #   default (ofat10)  -- addresses each panel point by config NAME off the
    #                        OAT table. Clean local response curve, but it can
    #                        only ever draw the 56 OAT configs; on a joint run,
    #                        or on a hand-picked shortlist, it draws nothing.
    #   marginal          -- addresses points by parameter VALUE, so it plots
    #                        whatever configs it is given (7, 56, 120). Each
    #                        point is a mean over the configs at that level,
    #                        i.e. a main effect rather than an OAT curve.
    # `marginal` is run twice, all-configs and joint-only: the joint design is
    # the one where every axis varies at once, so it is the only one whose
    # marginal is a genuine main effect rather than "the baseline plus noise".
    #
    # `design` selects WHICH by-name table the OAT line figures are laid out
    # from: 'ofat10' is the original baseline, 'ofat10b' the one re-centred on
    # the best-scoring config. A run of one design has no configs from the
    # other, so drawing both would emit a figure of empty panels.
    jobs = [(line, ["--design", design], str(prefix))]
    # A run whose primary layout IS marginal (the upsampled full factorial) would
    # otherwise draw the identical figure twice under two names.
    if design != "marginal":
        jobs.append((line, ["--design", "marginal"], f"{prefix}_marginal"))
    # The joint-only views only exist for a run that HAS joint configs. An OAT-B
    # run is 56 one-at-a-time configs, so both would be empty.
    if design == "ofat10":
        jobs += [(line, ["--design", "marginal", "--config-filter", "joint"],
                  f"{prefix}_marginal_jointonly"),
                 (scatter, [], str(prefix)),
                 (scatter, ["--design", "joint"], f"{prefix}_jointonly")]
    else:
        jobs += [(scatter, [], str(prefix))]
    for script, extra, out_prefix in jobs:
        if not os.path.exists(script):
            print(f"  figure script missing: {script}; skipping", flush=True)
            continue
        try:
            _run(["python3", script, "--scored", str(scored_csv),
                  "--clusterability", str(cluster_csv),
                  "--prefix", out_prefix] + extra)
        except Exception as e:
            print(f"  {os.path.basename(script)} failed: {e}", flush=True)


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

# ---------------------------------------------------------------------------
# JOINT (space-filling) design
#
# The OAT configs above move ONE parameter at a time off a fixed baseline. That
# is the right shape for the line figures, but it makes the SCATTER figures
# nearly unreadable: colouring the points by `dim` leaves most of them at the
# baseline dim=50, so there is no colour gradient to read a trend from.
#
# This design varies EVERY parameter at once, as a Latin hypercube over the same
# level grids -- each level of each axis used an equal number of times, permuted
# independently per axis -- so every scatter panel spans the full colour range.
#
# MUST STAY IN SYNC with geniml_dev/geniml/bedspace/tests/experiment_3cluster_sweep.py,
# which carries an identical seeded copy. `--check-joint-sync` compares them.
# ---------------------------------------------------------------------------
OFAT10_BASELINE = "d50_e10_m0.2_neg100"
OFAT10_LEVELS = {
    "dim": ["2", "3", "5", "10", "15", "25", "50", "75", "100", "200"],
    "epoch": ["1", "2", "3", "5", "10", "15", "20", "30", "50", "100"],
    "margin": ["0.05", "0.1", "0.15", "0.2", "0.3", "0.4", "0.5", "0.7", "1.0", "1.5"],
    "negSearchLimit": ["2", "5", "10", "20", "35", "50", "75", "100", "150", "200"],
    "minCount": ["1", "2", "3", "4", "5", "7", "10", "15", "20", "30"],
    "lr": ["0.00001", "0.00005", "0.0001", "0.0002", "0.0005",
           "0.001", "0.002", "0.005", "0.01", "0.02"],
    "adagrad": [False, True],  # binary -- 2 levels is the maximum that exists
}
JOINT_N = 64
JOINT_SEED = 20260811
JOINT_PREFIX = "jt"


def joint_design_configs(n: int = JOINT_N, seed: int = JOINT_SEED) -> dict:
    """Deterministic Latin-hypercube configs in which every parameter varies."""
    rng = np.random.default_rng(seed)
    axes = ["dim", "epoch", "margin", "negSearchLimit", "minCount", "lr"]
    picks = {}
    for ax in axes:
        lv = OFAT10_LEVELS[ax]
        pool = np.tile(np.arange(len(lv)), int(np.ceil(n / len(lv))))[:n]
        rng.shuffle(pool)
        picks[ax] = [lv[int(i)] for i in pool]
    ada = np.tile([False, True], int(np.ceil(n / 2)))[:n]
    rng.shuffle(ada)

    out = {}
    for i in range(n):
        c = {ax: picks[ax][i] for ax in axes}
        c["adagrad"] = bool(ada[i])
        name = (f"{JOINT_PREFIX}{i:03d}_d{c['dim']}_e{c['epoch']}_m{c['margin']}"
                f"_neg{c['negSearchLimit']}_mc{c['minCount']}"
                f"_lr{c['lr']}_ada{int(c['adagrad'])}")
        out[name] = {"description": f"JOINT[{i:03d}] space-filling (all params vary)",
                     **c}
    return out


JOINT_CONFIGS = joint_design_configs()
SWEEP_CONFIGS.update(JOINT_CONFIGS)

# ---------------------------------------------------------------------------
# Canonical config table.
#
# `experiment_3cluster_sweep.py` in the geniml tests dir is the source of truth
# for config NAMES: plot_param_figure.py lays its line-figure panels out by
# calling that module's ofat10_config_for(), so a run whose configs are named
# differently produces empty panels. The local table above is a stale partial
# copy (it predates the 10-level OAT group), so when the canonical module is
# importable we take its definitions wholesale. The local copy remains the
# fallback for a machine without the geniml checkout.
# ---------------------------------------------------------------------------
BEDSPACE_TESTS_DIR = os.environ.get(
    "BEDSPACE_TESTS_DIR", "/home/bx2ur/code/geniml_dev/geniml/bedspace/tests"
)
CANONICAL_CONFIGS = False
try:
    if BEDSPACE_TESTS_DIR not in sys.path:
        sys.path.insert(0, BEDSPACE_TESTS_DIR)
    import experiment_3cluster_sweep as _canon

    SWEEP_CONFIGS.update(_canon.SWEEP_CONFIGS)
    JOINT_CONFIGS = dict(_canon.JOINT_CONFIGS)
    OFAT10_LEVELS = dict(_canon.OFAT10_LEVELS)
    OFAT10_BASELINE = _canon.OFAT10_BASELINE
    ofat10_config_for = _canon.ofat10_config_for
    ofat10_configs = _canon.ofat10_configs
    # The re-centred OAT design exists only in the canonical module -- the local
    # fallback table below predates it, so --ofat10b-configs is gated on this.
    OFAT10B_BASELINE = _canon.OFAT10B_BASELINE
    ofat10b_configs = _canon.ofat10b_configs
    # Third OAT design, re-centred on jt042 -- also canonical-only.
    OFAT10C_BASELINE = _canon.OFAT10C_BASELINE
    ofat10c_configs = _canon.ofat10c_configs
    # Fourth OAT design, re-centred on the joint-criterion winner of all 221.
    OFAT10D_BASELINE = _canon.OFAT10D_BASELINE
    ofat10d_configs = _canon.ofat10d_configs
    # Fifth OAT design, centred on oc_d75_e20 -- the epoch=20 rung of OAT-C.
    OFAT10E_BASELINE = _canon.OFAT10E_BASELINE
    ofat10e_configs = _canon.ofat10e_configs
    # The 36-point FULL FACTORIAL for the class-balanced run. Canonical-only:
    # unlike the OAT families it has no local fallback, because it exists to be
    # read against a specific training file and a stale copy of it would be
    # worse than an error.
    UPSAMPLED_CONFIGS = dict(_canon.UPSAMPLED_CONFIGS)
    UPSAMPLED_LEVELS = dict(_canon.UPSAMPLED_LEVELS)
    upsampled_configs = _canon.upsampled_configs
    upsampled_configs_for_dim = _canon.upsampled_configs_for_dim
    CANONICAL_CONFIGS = True
except Exception as _e:  # noqa: BLE001 - any import problem falls back cleanly
    print(f"NOTE: canonical config table unavailable ({_e}); using the local "
          f"copy. Line-figure panels may be incomplete.", file=sys.stderr)


def ofat10_config_for(axis: str, level):
    """Name of the OAT config at `level` on `axis`, everything else at baseline."""
    base = SWEEP_CONFIGS[OFAT10_BASELINE]
    want = {**base, axis: level}
    for name, cfg in SWEEP_CONFIGS.items():
        if name.startswith(JOINT_PREFIX):
            continue
        if all(str(cfg[k]) == str(want[k]) for k in
               ("dim", "epoch", "margin", "negSearchLimit", "minCount", "lr", "adagrad")):
            return name
    return None


def ofat10_configs() -> List[str]:
    """The 56 config names of the 10-level OAT design, baseline first."""
    out, missing = [OFAT10_BASELINE], []
    for axis, levels in OFAT10_LEVELS.items():
        for lv in levels:
            name = ofat10_config_for(axis, lv)
            if name is None:
                missing.append(f"{axis}={lv}")
            elif name not in out:
                out.append(name)
    if missing:
        raise KeyError(f"10-level OAT design has no config for: {', '.join(missing)}")
    return out


def joint_configs() -> List[str]:
    """The joint-design config names, in index order."""
    return list(JOINT_CONFIGS.keys())


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

def _parse_balance_factors(spec):
    """'GROUP=N,GROUP=N' -> {group: n}. Used by --balance-strategy factors."""
    if not spec:
        return None
    out = {}
    for part in str(spec).split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise ValueError(f"--balance-factors entry '{part}' is not GROUP=N")
        key, _, val = part.rpartition("=")
        out[key.strip()] = int(val)
    return out


def balance_train_input(train_input, strategy=None, factors=None, min_ratio=None,
                        max_rounds=None, seed=42):
    """Oversample the minority label groups of the shared sweep corpus, in place.

    WHY IT LIVES HERE AND NOT IN THE PER-CONFIG TRAIN CALL. The holdout and
    LOOCV drivers preprocess per config, so they balance inside
    `geniml bedspace train --balance`. The sweep preprocesses ONCE and every
    config trains from that one file, so there is nothing per-config to balance
    -- it has to happen here, after preprocess and before the config loop.

    THE RULE IS geniml's, not a copy. `geniml.bedspace.balance` is the single
    implementation both paths call, so the two arms provably apply one
    transform. Default strategy 'labels': count how often each LABEL appears,
    counting a multi-label document once per label; take the commonest label's
    count as the target; duplicate every document carrying the SCARCEST label
    until it reaches it.

        sweep, 252 documents
          occurrences  Control 178 | Tumoral  63 | Renin  23    ratio 7.74
          x8 on the 23 Renin-carrying documents
          occurrences  Control 178 | Tumoral 147 | Renin 184    ratio 1.25
          252 -> 413 rows

    Tumoral rises too, because half the duplicated documents carry it. That is
    correct, and it is why ONE pass brings every label close instead of fixing
    one and stranding another.

    The unbalanced corpus is kept beside it as train_input_original.txt: RCS
    needs it (duplicate documents become duplicate columns and inflate the
    score), and it is the provenance record for the tokenization.

    Returns the path to use for training -- balanced when geniml decided to
    balance, the original when it decided the file was already within
    --balance-min-ratio.
    """
    from geniml.bedspace.balance import (DEFAULT_MAX_ROUNDS, DEFAULT_MIN_RATIO,
                                         DEFAULT_STRATEGY, balance_training_data)

    strategy = strategy or DEFAULT_STRATEGY
    min_ratio = DEFAULT_MIN_RATIO if min_ratio is None else min_ratio
    max_rounds = DEFAULT_MAX_ROUNDS if max_rounds is None else max_rounds

    out_dir = os.path.dirname(os.path.abspath(train_input))
    original = os.path.join(out_dir, "train_input_original.txt")
    receipt_path = os.path.join(out_dir, "BALANCE_RECEIPT.json")

    # Already balanced by an earlier job: train_input_original.txt exists, so
    # train_input.txt is the balanced copy. Re-balancing it would oversample
    # the oversampled file.
    if os.path.exists(original):
        with open(train_input) as f:
            n = sum(1 for _ in f)
        print(f"  already balanced ({n} rows); original kept at "
              f"{os.path.basename(original)}")
        return train_input

    if not os.path.exists(original):
        shutil.copy2(train_input, original)

    # geniml writes output_path only when it decides to balance, so build into
    # a temp name and move: train_input.txt is then either fully balanced or
    # untouched, never half-written.
    tmp = train_input + ".partial"
    receipt = balance_training_data(
        input_path=original, output_path=tmp, strategy=strategy,
        factors=_parse_balance_factors(factors), min_ratio=min_ratio,
        max_rounds=max_rounds, random_seed=seed,
    )
    if not receipt.get("balanced"):
        if os.path.exists(tmp):
            os.remove(tmp)
        os.remove(original)          # nothing was changed; no provenance to keep
        print(f"  NOT balanced: {receipt.get('reason')} "
              f"(within --balance-min-ratio {min_ratio})")
        return train_input

    shutil.move(tmp, train_input)
    receipt["output_path"] = os.path.abspath(train_input)
    with open(receipt_path, "w") as f:
        json.dump(receipt, f, indent=2)
    print(f"  BALANCED {receipt['total_original']} -> "
          f"{receipt['total_balanced']} rows (strategy '{strategy}')")
    print(f"  unbalanced corpus kept at {os.path.basename(original)}")
    print(f"  receipt -> {os.path.basename(receipt_path)}")
    return train_input


def run_preprocess(data_path, metadata, universe, output_dir, labels, reuse=False,
                   balance=False, balance_strategy=None, balance_factors=None,
                   balance_min_ratio=None, balance_max_rounds=None,
                   balance_seed=42):
    """Run preprocessing (shared across all configs).

    `reuse` returns an existing train_input.txt untouched instead of rebuilding
    it. That matters for more than the ~20 minutes it saves: when a run resumes
    models trained by an EARLIER job -- which is the whole point of training a
    second design into the same directory -- rebuilding the tokenization risks
    putting two different tokenizations inside one figure, the same way a
    metadata change would. Reusing the file makes "the resumed models and the
    new ones saw identical input" a fact rather than an assumption. It is opt-in
    because it is only safe when data/metadata/universe are unchanged, which the
    caller knows and this function cannot check.
    """
    print("\n" + "=" * 60)
    print("STEP 1: PREPROCESSING")
    print("=" * 60)
    os.makedirs(output_dir, exist_ok=True)
    existing = os.path.join(output_dir, "train_input.txt")
    if reuse and os.path.exists(existing) and os.path.getsize(existing) > 0:
        with open(existing) as f:
            n = sum(1 for _ in f)
        print(f"  REUSING existing tokenization: {n} samples -> {existing}")
        print("  (--reuse-preprocessed; pass nothing to rebuild it)")
        if balance:
            return balance_train_input(existing, balance_strategy, balance_factors,
                                       balance_min_ratio, balance_max_rounds,
                                       balance_seed)
        return existing
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
    if balance:
        return balance_train_input(train_input_path, balance_strategy,
                                   balance_factors, balance_min_ratio,
                                   balance_max_rounds, balance_seed)
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
    skip_geniml_eval: bool = False,
    geniml_eval_workers: int = 10,
    bin_embed: Optional[str] = None,
    score_after_train: bool = True,
    skip_eval: bool = False,
    force: bool = False,
    metadata_test: Optional[str] = None,
    target: str = "Renin",
    negatives: str = "Control,Tumoral",
    line_labels: str = "Control,Tumoral",
    weight_constraint: str = "nonneg-scaled",
    class_weight: str = "balanced",
    summary_prefix: Optional[str] = None,
    line_design: str = "ofat10",
    reuse_preprocessed: bool = False,
    balance: bool = False,
    balance_strategy: Optional[str] = None,
    balance_factors: Optional[str] = None,
    balance_min_ratio: Optional[float] = None,
    balance_max_rounds: Optional[int] = None,
    balance_seed: int = 42,
    rescore: bool = False,
    clusterability_from: Optional[str] = None,
    fit_max_regions: Optional[int] = None,
    arm_figures: bool = True,
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
    train_input = run_preprocess(data_path, metadata, universe, str(preprocess_dir),
                                 labels, reuse=reuse_preprocessed,
                                 balance=balance,
                                 balance_strategy=balance_strategy,
                                 balance_factors=balance_factors,
                                 balance_min_ratio=balance_min_ratio,
                                 balance_max_rounds=balance_max_rounds,
                                 balance_seed=balance_seed)

    wb_group = wandb_group or experiment_name if wandb_project else None
    results = {}
    resumed = 0

    for name in configs_to_run:
        if name not in SWEEP_CONFIGS:
            print(f"\n  WARNING: Unknown config '{name}', skipping")
            continue

        cfg = SWEEP_CONFIGS[name]
        print(f"\n{'='*60}\nConfig: {name}\n  {cfg['description']}\n{'='*60}")

        cfg_dir = exp_dir / name
        cfg_dir.mkdir(parents=True, exist_ok=True)

        # ---- Resume ----
        # A 120-config run is many hours; without this, a job that hits its wall
        # clock has to retrain everything. A config counts as done when its
        # model and the last scoring method's output are both present.
        sentinel = cfg_dir / f"reninness_score_{SCORING_METHODS[-1]}.csv"
        existing_model = next(iter(sorted(cfg_dir.glob("*.tsv"))), None)
        if not force and existing_model is not None and (
            sentinel.exists() or not score_after_train or rescore
        ):
            # --rescore: the model on disk is fine, the SCORING changed. Re-run
            # distances + score against it and keep going, rather than skipping
            # the config (plain resume) or retraining it (--force, which would
            # replace the models every published number was computed from).
            if rescore and score_after_train:
                print("  model on disk -> re-running distances + score "
                      "(no retraining)", flush=True)
                try:
                    run_distances_and_score(
                        model_prefix=str(existing_model)[:-4]
                        if str(existing_model).endswith(".tsv") else str(existing_model),
                        cfg_dir=cfg_dir, starspace=starspace_path, universe=universe,
                        data_path=data_path, label_col=labels,
                        metadata_train=metadata,
                        metadata_test=metadata_test or metadata,
                        config_name=name, target=target, negatives=negatives,
                        line_labels=line_labels,
                        weight_constraint=weight_constraint,
                        class_weight=class_weight,
                    )
                except Exception as e:
                    print(f"    distances/score failed: {e}", flush=True)
            else:
                print(f"  already done, skipping (use --force to redo)")
            # Same shape as the training path's wb_config below. It used to be
            # None, which crashed print_comparison -- and print_comparison runs
            # AFTER the whole training loop and BEFORE the metrics and figures,
            # so on a resumed run that one None threw away hours of completed
            # work at the last step. The parameters are known from the config
            # table, so there is no reason to leave them out.
            results[name] = {"status": "SUCCESS", "config": comparison_config(cfg),
                             "model_path": str(existing_model), "stats": {},
                             "resumed": True}
            resumed += 1
            continue

        wb_config = {"universe": universe_name, **comparison_config(cfg)}
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
                                   output_dir=str(cfg_dir), deseq_path=deseq_path,
                                   fit_max_regions=fit_max_regions)
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

        # ---- geniml.eval clusterability (second clusterability family) ----
        # Independent of the silhouette above: these grade the embedding as a
        # genomic-region embedding rather than against an imposed k.
        if not skip_geniml_eval:
            try:
                print("    Running geniml.eval region-embedding tests...", flush=True)
                ge = geniml_eval_metrics(
                    model_path,
                    cache_path=str(cfg_dir / f"{name}_base_embed.pt"),
                    num_workers=geniml_eval_workers,
                    bin_embed=bin_embed,
                )
                stats["geniml_eval"] = ge
                print("    " + "  ".join(f"{k}={v:.4f}" for k, v in ge.items()
                                         if isinstance(v, float)), flush=True)
                if tracker.active:
                    for k, v in ge.items():
                        if isinstance(v, (int, float)):
                            wandb.summary[f"geniml_eval/{k}"] = v
            except Exception as e:
                print(f"    geniml.eval diagnostic failed: {e}", flush=True)

        # ---- distances + score (all four methods) ----
        # Runs here rather than in a separate driver so one 7a invocation
        # produces the same CSVs and figures as the geniml tests pipeline.
        if score_after_train:
            try:
                print("    Running distances + score (4 methods)...", flush=True)
                run_distances_and_score(
                    model_prefix=str(model_path)[:-4] if str(model_path).endswith(".tsv")
                    else model_path,
                    cfg_dir=cfg_dir, starspace=starspace_path, universe=universe,
                    data_path=data_path, label_col=labels,
                    metadata_train=metadata, metadata_test=metadata_test or metadata,
                    config_name=name, target=target, negatives=negatives,
                    line_labels=line_labels, weight_constraint=weight_constraint,
                    class_weight=class_weight,
                )
            except Exception as e:
                print(f"    distances/score failed: {e}", flush=True)

        results[name] = {
            "status": "SUCCESS", "config": wb_config,
            "model_path": model_path, "stats": stats,
        }
        tracker.finish()

    if resumed:
        print(f"\nResumed (skipped) {resumed}/{len(configs_to_run)} configs "
              f"already on disk; pass --force to retrain them.")
    print_comparison(results, configs_to_run, exp_dir,
                     prefix=summary_prefix or experiment_name)

    # ---- run-level metric CSVs + figures ----
    if score_after_train and not skip_eval:
        print(f"\n{'='*70}\nEVALUATION METRICS + FIGURES\n{'='*70}")
        scored_rows, cluster_rows = [], []
        # The clusterability families (silhouette, geniml.eval CTT/GDST/NPT)
        # grade the MODEL, not the score, so a --rescore run must not recompute
        # them: the models are untouched, the tests are expensive, and a fresh
        # run of them would move published CTT numbers for reasons unrelated to
        # the scoring change. --clusterability-from reuses the table as written.
        reused_cluster = {}
        if clusterability_from:
            _cl = pd.read_csv(clusterability_from)
            reused_cluster = {str(r["config"]): dict(r)
                              for _, r in _cl.iterrows()}
            print(f"reusing clusterability rows for {len(reused_cluster)} configs "
                  f"from {clusterability_from}")
        ok = [n for n in configs_to_run if results.get(n, {}).get("status") == "SUCCESS"]
        for i, name in enumerate(ok, 1):
            print(f"[{i}/{len(ok)}] metrics for {name}", flush=True)
            have = name in reused_cluster
            s, c = collect_config_metrics(
                name, exp_dir / name, target=target,
                silhouette=(not skip_clustering) and not have,
                geniml_eval=(not skip_geniml_eval) and not have,
                geniml_eval_workers=geniml_eval_workers, bin_embed=bin_embed,
            )
            scored_rows.extend(s)
            cluster_rows.append(reused_cluster[name] if have else c)
        # The summary CSVs and the figures are named by summary_prefix, NOT by
        # experiment_name. They cover only `configs_to_run`, so a second design
        # trained into the same models directory (which is what makes the ten
        # shared configs reusable) would otherwise overwrite the first design's
        # tables and figures with a table that does not contain its configs.
        prefix = summary_prefix or experiment_name
        s_csv, c_csv = write_run_summary(exp_dir, scored_rows, cluster_rows,
                                         prefix=prefix)
        if arm_figures:
            make_figures(s_csv, c_csv, exp_dir / prefix, design=line_design)
        else:
            print("  --no-arm-figures: summary CSVs written, per-arm parameter "
                  "figures skipped", flush=True)

    return results


def comparison_config(cfg: dict) -> dict:
    """Parameters of one config, in the shape print_comparison and wandb expect.

    NB the key is `epochs`, not the config table's `epoch` -- that spelling is
    what wandb logged from the first sweep, and print_comparison reads it.
    """
    return {
        "dim": int(cfg["dim"]),
        "epochs": int(cfg["epoch"]),
        "lr": float(cfg["lr"]),
        "negSearchLimit": int(cfg["negSearchLimit"]),
        "minCount": int(cfg["minCount"]),
        "margin": float(cfg["margin"]),
        "adagrad": cfg["adagrad"],
    }


def print_comparison(results, config_names, exp_dir, prefix="sweep"):
    """Print and save a comparison table.

    `prefix` names the file for the same reason the summary CSVs and figures
    are named: it describes only `config_names`, so a second design trained
    into the same models directory would otherwise overwrite the first
    design's table with one that does not contain its configs.
    """
    lines = ["\n" + "=" * 100, "PARAMETER SWEEP COMPARISON", "=" * 100]
    header = (f"{'Config':<25} {'ep':<5} {'dim':<5} {'lr':<9} {'margin':<8} {'mc':<4} {'ada':<5} "
              f"{'best_k':<7} {'sil@3':<8} {'sil_best':<10} {'label_sep':<10} {'status':<8}")
    lines.append(header); lines.append("-" * 120)
    for name in config_names:
        r = results.get(name)
        if r is None or r["status"] != "SUCCESS":
            lines.append(f"{name:<25} {'':88} {r['status'] if r else 'SKIP':<8}")
            continue
        # Never subscript this blind. A run dir written by an older version of
        # this script carries config=None, and the whole point of the table is
        # to report on work that is already done -- it must not be able to
        # destroy it.
        c = r.get("config")
        if not c and name in SWEEP_CONFIGS:
            c = comparison_config(SWEEP_CONFIGS[name])
        s = r.get("stats", {})
        if not c:
            # Unknown config with no recorded parameters: report the status and
            # move on. This table runs after the training loop and before the
            # metrics and figures, so it must not be able to abort a run whose
            # work is already done.
            lines.append(f"{name:<25} {'':88} {'SUCCESS':<8}")
            continue
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
    report_path = exp_dir / f"{prefix}_comparison.txt"
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
                        help="Skip multi-k clustering diagnostic (silhouette family)")
    parser.add_argument("--skip-geniml-eval", action="store_true",
                        help="Skip the geniml.eval region-embedding tests "
                             "(CTT/GDST/NPT), the second clusterability family")
    parser.add_argument("--geniml-eval-workers", type=int, default=10,
                        help="Parallel workers for the geniml.eval tests (default: 10)")
    parser.add_argument("--bin-embed", default=None,
                        help="Path to a binary-embedding pickle over the same "
                             "universe; enables the geniml.eval RCT test. "
                             "Omitted -> RCT is skipped (it trains an MLP per "
                             "CV fold, so it belongs in its own job).")
    parser.add_argument("--no-score", dest="score_after_train",
                        action="store_false", default=True,
                        help="Train only: skip distances + score and the "
                             "metric CSVs / figures.")
    parser.add_argument("--no-arm-figures", dest="arm_figures",
                        action="store_false", default=True,
                        help="Write the run-level summary CSVs but SKIP the old per-arm parameter figures. Those are superseded by the five deliverable figures (7z_upsampled_final_figures.py); drawing them here just refills the arm directory with plots nobody reads.")
    parser.add_argument("--skip-eval", action="store_true",
                        help="Score, but skip the run-level metric CSVs and figures.")
    parser.add_argument("--metadata-test", default=None,
                        help="Test metadata for the distances/score step "
                             "(default: --metadata, i.e. score the training "
                             "samples themselves).")
    parser.add_argument("--fit-frac", type=float, default=None,
                        help="Fit the t-SNE / UMAP projections on this FRACTION "
                             "of each config's vocabulary (e.g. 0.1 for 10%%, "
                             "~42k of 422k regions). Mutually exclusive with "
                             "--fit-max-regions. The projections are per-config "
                             "PNGs only -- silhouette, CTT/GDST and nl_auprc are "
                             "all computed on the FULL embeddings, so this "
                             "changes the pictures and no number in the summary "
                             "table or Fig1.")
    parser.add_argument("--fit-max-regions", type=int, default=None,
                        help="cap how many regions the t-SNE/UMAP figures are FIT on "
                             "(default: the whole vocabulary). The plotted counts are "
                             "separate; a full fit is ~400k points and single threaded, "
                             "so cap this if the figures hold up a sweep.")
    parser.add_argument("--target", default="Renin", help="Positive class (default: Renin)")
    parser.add_argument("--negatives", default="Control,Tumoral",
                        help="Negatives for weighted_multi (default: Control,Tumoral). "
                             "The first is also the single negative for --method pairwise.")
    parser.add_argument("--line-labels", default="Control,Tumoral",
                        help="Two negatives forming the line for --method line; "
                             "empty string disables that method.")
    parser.add_argument("--weight-constraint", default="nonneg-scaled",
                        choices=["nonneg-scaled", "simplex"])
    parser.add_argument("--class-weight", default="balanced",
                        choices=["none", "balanced"])
    parser.add_argument("--force", action="store_true",
                        help="Retrain every config, even those already on disk. "
                             "Default: resume — a config with a model and score "
                             "outputs is skipped.")
    parser.add_argument("--joint-configs", action="store_true",
                        help="Run the 64-config joint (space-filling) design.")
    parser.add_argument("--ofat10-configs", action="store_true",
                        help="Run the 56-config 10-level OAT design.")
    parser.add_argument("--ofat10b-configs", action="store_true",
                        help="Run the 56-config 10-level OAT design RE-CENTRED "
                             "on the best-scoring config of the 120-config run "
                             "(d50_e10_m0.2_neg100_lr02). 10 of the 56 are the "
                             "ofat10 lr axis and resume from its models; 46 are "
                             "new. Requires the canonical config table.")
    parser.add_argument("--ofat10c-configs", action="store_true",
                        help="The 56-config 10-level OAT design re-centred on "
                             "jt042_d75_e2_m0.3_neg10_mc1_lr0.0005_ada1, the "
                             "config nothing beats on both CTT and "
                             "weighted_multi_learned_post accuracy. 55 of the "
                             "56 are new; the baseline itself is reused.")
    parser.add_argument("--ofat10d-configs", action="store_true",
                        help="The 56-config 10-level OAT design re-centred on "
                             "ob_d50_e10_m0.7_neg100_mc1_lr0.02_ada0, the best "
                             "config of all 221 on CTT + learned_post accuracy. "
                             "45 of 56 are new; the margin axis is shared with "
                             "ofat10b and is reused.")
    parser.add_argument("--ofat10e-configs", action="store_true",
                        help="The 56-config 10-level OAT design centred on "
                             "oc_d75_e20_m0.3_neg10_mc1_lr0.0005_ada1. 46 of 56 "
                             "are new. The EPOCH axis is entirely reused: that "
                             "centre differs from the ofat10c baseline only in "
                             "epoch, so sweeping epoch from it retraces the line "
                             "ofat10c already traced and adds nothing. What is "
                             "new is whether the other five responses look "
                             "different at epoch=20 than at epoch=2.")
    parser.add_argument("--upsampled-configs", action="store_true",
                        help="The 36-config FULL FACTORIAL of the "
                             "class-balanced run: dim x epoch x margin x "
                             "negSearchLimit = 3x3x2x2, with lr=0.001, "
                             "minCount=1 and adagrad=true held constant. "
                             "Unlike every OAT design this one crosses its "
                             "axes completely, so it can show interactions; "
                             "it is meant to be read against the SAME 36 "
                             "points trained on the unbalanced input, not "
                             "against the 312-config sweep (different "
                             "train_input.txt). Requires the canonical config "
                             "table.")
    parser.add_argument("--upsampled-dim", default=None,
                        help="Restrict --upsampled-configs to ONE dim level "
                             "(25/50/100/150), i.e. 36 of the 144. The sweep is "
                             "split across four concurrent jobs this way "
                             "because per-config cost is driven almost entirely "
                             "by dim -- measured 239s/354s/457s at dim "
                             "10/50/100 -- so one job per dim has a predictable "
                             "runtime, and the serial 144 would exceed the 12 h "
                             "wall.")
    parser.add_argument(
        "--balance", action="store_true",
        help="Oversample the minority label groups of the shared sweep corpus "
        "after preprocessing, before any config trains. The rule is geniml's "
        "(geniml.bedspace.balance), the same one `geniml bedspace train "
        "--balance` applies per config in the holdout / LOOCV arms, so the arms "
        "provably share one transform. The unbalanced corpus is kept beside it "
        "as train_input_original.txt -- RCS needs it.")
    parser.add_argument(
        "--balance-strategy", default=None, choices=["labels", "to-majority", "factors"],
        help="'labels' (geniml default) equalises LABEL OCCURRENCES, counting a "
        "multi-label document once per label.")
    parser.add_argument("--balance-factors", default=None,
                        help="--balance-strategy factors: 'GROUP=N,...'")
    parser.add_argument("--balance-min-ratio", type=float, default=None,
                        help="detection threshold; below it nothing is written")
    parser.add_argument("--balance-max-rounds", type=int, default=None)
    parser.add_argument("--balance-seed", type=int, default=42)
    parser.add_argument("--reuse-preprocessed", action="store_true",
                        help="Reuse an existing preprocessed/train_input.txt "
                             "instead of rebuilding it. Use when resuming into "
                             "a directory whose models were trained by an "
                             "earlier job, so old and new models provably share "
                             "one tokenization. Only valid when --data-path, "
                             "--metadata and --universe are unchanged.")
    parser.add_argument("--rescore", action="store_true",
                        help="Re-run distances + score (and the metric CSVs / "
                             "figures) against the models already on disk, "
                             "without retraining any of them. Use after a "
                             "scoring-code change; --force would retrain "
                             "instead and replace the models.")
    parser.add_argument("--clusterability-from", default=None,
                        help="Reuse the silhouette / geniml.eval rows from an "
                             "existing <prefix>_clusterability.csv instead of "
                             "recomputing them. They grade the model, so a "
                             "--rescore run should reuse them verbatim.")
    parser.add_argument("--summary-prefix", default=None,
                        help="Basename for the run summary CSVs and figures "
                             "(default: --name). Set it when training a second "
                             "design into an existing models directory, so the "
                             "first design's tables and figures survive.")
    parser.add_argument("--deseq", default=None,
                        help="Path to DESeq result file (CSV with seqnames, start, end, log2FoldChange).")
    parser.add_argument("--methods", default=None,
                        help="Comma-separated subset of scoring methods to run "
                             f"(default: all of {','.join(ALL_SCORING_METHODS)}). "
                             "A two-label model supports only 'pairwise'.")
    parser.add_argument("--list-configs", action="store_true",
                        help="List available configs and exit")
    args = parser.parse_args()

    if args.fit_frac is not None:
        if args.fit_max_regions is not None:
            print("Error: pass --fit-frac or --fit-max-regions, not both.",
                  file=sys.stderr)
            return 1
        if not 0 < args.fit_frac <= 1:
            print(f"Error: --fit-frac must be in (0, 1], got {args.fit_frac}",
                  file=sys.stderr)
            return 1
        # _fit_projection reads a value in (0, 1] as a fraction.
        args.fit_max_regions = args.fit_frac

    # Design shortcuts, so a full run does not need 120 names on the command line.
    if args.methods:
        chosen = [m.strip() for m in args.methods.split(",") if m.strip()]
        unknown = [m for m in chosen if m not in ALL_SCORING_METHODS]
        if unknown:
            print(f"Error: unknown scoring method(s) {unknown}; "
                  f"choose from {list(ALL_SCORING_METHODS)}", file=sys.stderr)
            return 1
        # Slice-assign so every module-level reader (the score dispatch, the
        # metric loop, and the resume sentinel) sees the restriction.
        SCORING_METHODS[:] = chosen
        print(f"Scoring methods restricted to: {SCORING_METHODS}")

    design = []
    if args.ofat10_configs:
        design += ofat10_configs()
    if args.ofat10b_configs:
        if not CANONICAL_CONFIGS:
            print("Error: --ofat10b-configs needs the canonical config table "
                  f"(experiment_3cluster_sweep.py under {BEDSPACE_TESTS_DIR}); "
                  "the local fallback table predates that design.",
                  file=sys.stderr)
            return 1
        design += [c for c in ofat10b_configs() if c not in design]
    if args.ofat10c_configs:
        if not CANONICAL_CONFIGS:
            print("Error: --ofat10c-configs needs the canonical config table "
                  f"(experiment_3cluster_sweep.py under {BEDSPACE_TESTS_DIR}); "
                  "the local fallback table predates that design.",
                  file=sys.stderr)
            return 1
        design += [c for c in ofat10c_configs() if c not in design]
    if args.ofat10d_configs:
        if not CANONICAL_CONFIGS:
            print("Error: --ofat10d-configs needs the canonical config table "
                  f"(experiment_3cluster_sweep.py under {BEDSPACE_TESTS_DIR}); "
                  "the local fallback table predates that design.",
                  file=sys.stderr)
            return 1
        design += [c for c in ofat10d_configs() if c not in design]
    if args.ofat10e_configs:
        if not CANONICAL_CONFIGS:
            print("Error: --ofat10e-configs needs the canonical config table "
                  f"(experiment_3cluster_sweep.py under {BEDSPACE_TESTS_DIR}); "
                  "the local fallback table predates that design.",
                  file=sys.stderr)
            return 1
        design += [c for c in ofat10e_configs() if c not in design]
    if args.upsampled_configs:
        if not CANONICAL_CONFIGS:
            print("Error: --upsampled-configs needs the canonical config table "
                  f"(experiment_3cluster_sweep.py under {BEDSPACE_TESTS_DIR}); "
                  "there is no local fallback for the full factorial.",
                  file=sys.stderr)
            return 1
        _ups = (upsampled_configs_for_dim(args.upsampled_dim)
                if args.upsampled_dim else upsampled_configs())
        if args.upsampled_dim and not _ups:
            print(f"Error: --upsampled-dim {args.upsampled_dim} matches no "
                  f"config; levels are {UPSAMPLED_LEVELS['dim']}", file=sys.stderr)
            return 1
        design += [c for c in _ups if c not in design]
    if args.joint_configs:
        design += [c for c in JOINT_CONFIGS if c not in design]
    if design:
        args.configs = design + [c for c in (args.configs or []) if c not in design]
        print(f"Design selection: {len(args.configs)} configs "
              f"(ofat10={args.ofat10_configs}, ofat10b={args.ofat10b_configs}, "
              f"ofat10c={args.ofat10c_configs}, ofat10d={args.ofat10d_configs}, "
              f"ofat10e={args.ofat10e_configs}, upsampled={args.upsampled_configs}, "
              f"joint={args.joint_configs})")

    # The by-name line figures are laid out from ONE table. A run that is purely
    # the re-centred design gets that design's panels; anything else keeps the
    # original, whose baseline the mixed run still contains.
    _pure = lambda flag, *others: flag and not any(others)
    # The upsampled grid has no by-name panel table and cannot have one: 'marginal'
    # addresses panel points by parameter VALUE, which is the only layout that
    # works for a fully crossed design (each point is then a main effect, the
    # mean over the 12 configs sharing that level).
    if _pure(args.upsampled_configs, args.ofat10_configs, args.ofat10b_configs,
             args.ofat10c_configs, args.ofat10d_configs, args.ofat10e_configs,
             args.joint_configs):
        line_design = "marginal"
    elif _pure(args.ofat10e_configs, args.ofat10_configs, args.ofat10b_configs,
             args.ofat10c_configs, args.ofat10d_configs, args.joint_configs):
        line_design = "ofat10e"
    elif _pure(args.ofat10d_configs, args.ofat10_configs, args.ofat10b_configs,
             args.ofat10c_configs, args.ofat10e_configs, args.joint_configs):
        line_design = "ofat10d"
    elif _pure(args.ofat10c_configs, args.ofat10_configs, args.ofat10b_configs,
               args.ofat10d_configs, args.ofat10e_configs, args.joint_configs):
        line_design = "ofat10c"
    elif _pure(args.ofat10b_configs, args.ofat10_configs, args.ofat10c_configs,
               args.ofat10d_configs, args.ofat10e_configs, args.joint_configs):
        line_design = "ofat10b"
    else:
        line_design = "ofat10"


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
        skip_geniml_eval=args.skip_geniml_eval,
        geniml_eval_workers=args.geniml_eval_workers,
        bin_embed=args.bin_embed,
        score_after_train=args.score_after_train,
        skip_eval=args.skip_eval, force=args.force,
        metadata_test=args.metadata_test,
        target=args.target, negatives=args.negatives,
        line_labels=args.line_labels or None,
        weight_constraint=args.weight_constraint, class_weight=args.class_weight,
        summary_prefix=args.summary_prefix, line_design=line_design,
        reuse_preprocessed=args.reuse_preprocessed,
        balance=args.balance,
        balance_strategy=args.balance_strategy,
        balance_factors=args.balance_factors,
        balance_min_ratio=args.balance_min_ratio,
        balance_max_rounds=args.balance_max_rounds,
        balance_seed=args.balance_seed,
        rescore=args.rescore, clusterability_from=args.clusterability_from,
        fit_max_regions=args.fit_max_regions,
        arm_figures=args.arm_figures,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

