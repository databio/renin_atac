#!/usr/bin/env python
"""
7c  Leave-one-out cross-validation (LOOCV) driver for the bedspace pipeline.

Per fold:
  - hold out one sample; train on the remaining N-1
  - preprocess -> train -> distances(train) -> fit weights -> distances(test) -> score
  - the held-out sample never enters the fit that scores it

Weights
  - `weighted_multi` weights are fit on each fold's TRAIN split (train-fit).
  - `--learn-weights` adds a pooled OUT-OF-FOLD table as a comparison: each
    fold's weights fit on the other folds' held-out sub-scores.
  - the two appear side by side in the metrics as `weighted_multi_learned_post`
    (train-fit) and `weighted_multi_learned_out_of_fold` (pooled).

Balancing
  - `--balance` oversamples the minority label inside `geniml bedspace train`.
  - each fold is preprocessed separately, so it cannot reach the held-out sample.
  - required whenever the sweep/holdout arms being compared against were balanced.

Usage
    python 7c_run_reninness_score_loo_in_batch.py \\
        --config     dim100_ep50_mar0.4_neg50 \\
        --metadata   metadata/renin_annotation_sheet_test_l3.csv \\
        --data-path  /path/to/bed/files/ \\
        --universe   /path/to/universe.bed \\
        --starspace  /path/to/Starspace/ \\
        --label-col  Label_1,Label_2 \\
        --target     Renin --negatives Control,Tumoral \\
        --line-labels Control,Tumoral \\
        --balance --learn-weights --class-weight balanced \\
        --output-dir <run>/loocv/<config>/

Output
    <output-dir>/
        fold_<sample>/
            reninness_score_<method>.csv   one per scoring method
            reninness_score.csv            copy of the --method table
            subscores.csv                  held-out per-negative sub-scores
            weights_trainfit.json          weights fit on this fold's train split
        reninness_score_loocv_<method>.csv concatenated over folds
        reninness_score_loocv_learned.csv  out-of-fold comparison
        loocv_trainfit_weights.csv         per-fold weights + fit diagnostics
        loocv_eval_metrics.csv             accuracy + clusterability
        model_<config>.tsv                 first-fold regions, per-label means
        <config>_starspace_embed.txt       held-out doc embeddings, all folds
        _loocv_artifacts/                  per-fold cached model + doc embed

Resume
  - a fold is done if its reninness_score.csv is readable; re-run to retry the rest
  - `--force` wipes every fold and redoes it
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Per-fold dump of the held-out sample's per-negative sub-scores, written by
# `bedspace score --subscores-out`. Survives the fold cleanup because
# `learn_weights_loocv` pools these across folds after every fold has run.
SUBSCORES_NAME = "subscores.csv"

# Weights fit on a fold's TRAIN split; applied to that fold's held-out sample.
TRAINFIT_WEIGHTS = "weights_trainfit.json"

# StarSpace marks label embeddings with this prefix in model.tsv. Label rows are
# found by the prefix everywhere, never by position/count -- the label count is
# 2 for a single `label` column and 3 for `--label-col Label_1,Label_2`.
LABEL_PREFIX = "__label__"

# Balancing defaults come from geniml, so this parser cannot offer a strategy
# geniml does not implement.
try:
    from geniml.bedspace.balance import (
        DEFAULT_MAX_ROUNDS as _BAL_ROUNDS, DEFAULT_MIN_RATIO as _BAL_RATIO,
        DEFAULT_STRATEGY as _BAL_STRAT, STRATEGIES as _BAL_STRATS)
except Exception:
    _BAL_STRATS = ("labels", "to-majority", "factors")
    _BAL_STRAT, _BAL_RATIO, _BAL_ROUNDS = "labels", 2.0, 1


# Load SWEEP_CONFIGS from the co-located 7a file. We can't `import` it
# directly because the module name starts with a digit, so use importlib.
_THIS_DIR = Path(__file__).resolve().parent
_SWEEP_FILE = _THIS_DIR / "7a_run_reninness_score_insample_in_batch.py"

if not _SWEEP_FILE.exists():
    print(f"ERROR: Expected sibling config file not found: {_SWEEP_FILE}", file=sys.stderr)
    print("This script reads SWEEP_CONFIGS from 7a_run_reninness_score_insample_in_batch.py.", file=sys.stderr)
    sys.exit(1)

_spec = importlib.util.spec_from_file_location("_sweep_configs_module", _SWEEP_FILE)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
SWEEP_CONFIGS = _mod.SWEEP_CONFIGS
# 7a also owns the shared evaluation metrics (both accuracy families and both
# clusterability families), so the three drivers cannot drift on definitions.
_sweep = _mod

# Same four scoring methods, under the same names, as 7a and 7c -- so the three
# drivers' outputs are directly comparable and feed the same figure scripts.
METHODS = _sweep.SCORING_METHODS


def method_alias(method: str) -> str:
    """Map the legacy --method value onto the per-method output name.

    `--method` used to select the single method a run scored with; it now only
    chooses which of the four tables is copied to the historical
    reninness_score.csv filename.
    """
    return {"pairwise": "pairwise",
            "weighted_multi": "weighted_multi_uniform"}.get(method, method)



def run(cmd: list[str]) -> None:
    """Run a subprocess, stream output, raise on failure."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {' '.join(cmd)}")


def write_fold_metadata(
    meta_df: pd.DataFrame, held_out_idx: int, train_path: Path, test_path: Path
) -> None:
    """Write per-fold train (N-1) and test (1) metadata CSVs."""
    train_df = meta_df.drop(index=held_out_idx)
    test_df = meta_df.loc[[held_out_idx]]
    train_df.to_csv(train_path, index=False)
    test_df.to_csv(test_path, index=False)


def run_fold(
    config_name: str,
    params: dict,
    held_out_name: str,
    train_meta_path: Path,
    test_meta_path: Path,
    data_path: str,
    universe: str,
    starspace: str,
    label_col: str,
    fold_dir: Path,
    target: str,
    method: str = "weighted_multi",
    negatives: str | None = None,
    line_labels: str | None = None,
    weight_constraint: str = "nonneg-scaled",
    class_weight: str = "balanced",
    balance: bool = False,
    balance_strategy: str = _BAL_STRAT,
    balance_factors: str | None = None,
    balance_min_ratio: float = _BAL_RATIO,
    balance_max_rounds: int = _BAL_ROUNDS,
    balance_seed: int = 42,
) -> tuple[Path | None, dict | None]:
    """One fold, weights fit on its TRAIN split.

    Returns (reninness_score.csv path, train-fit weights dict or None).
    """
    fold_dir.mkdir(parents=True, exist_ok=True)
    preprocess_dir = fold_dir / "preprocessed"

    # Train split only; the held-out BED is tokenized inside `distances`.
    run([
        "geniml", "bedspace", "preprocess",
        "-i", data_path, "-m", str(train_meta_path), "-u", universe,
        "-o", str(preprocess_dir), "-l", label_col, "--mode", "bulk",
    ])
    train_input = preprocess_dir / "train_input.txt"
    if not train_input.exists():
        raise FileNotFoundError(f"preprocess did not produce {train_input}")

    train_cmd = [
        "geniml", "bedspace", "train",
        "-s", starspace, "-i", str(train_input), "-o", str(fold_dir),
        "-n", str(params["epoch"]), "-d", str(params["dim"]), "-l", str(params["lr"]),
        "--neg-search-limit", str(params["negSearchLimit"]),
        "--margin", str(params["margin"]), "--min-count", str(params["minCount"]),
    ]
    if balance:
        # This fold's train_input holds the N-1 training samples only, so the
        # oversampling cannot reach the held-out sample.
        train_cmd += ["--balance",
                      "--balance-strategy", balance_strategy,
                      "--balance-min-ratio", str(balance_min_ratio),
                      "--balance-max-rounds", str(balance_max_rounds),
                      "--balance-seed", str(balance_seed)]
        if balance_factors:
            train_cmd += ["--balance-factors", balance_factors]
    train_cmd.append("--adagrad" if params.get("adagrad") else "--no-adagrad")
    run(train_cmd)

    candidates = sorted(fold_dir.glob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No model .tsv found in {fold_dir}")
    model_prefix = str(candidates[0])[:-4]

    # Distances on TRAIN -> `score --learn-weights` fits on training samples only.
    train_dir = fold_dir / "_train_distances"
    weights_info = None
    try:
        _sweep.run_distances_and_score(
            model_prefix=model_prefix, cfg_dir=train_dir, starspace=starspace,
            universe=universe, data_path=data_path, label_col=label_col,
            metadata_train=str(train_meta_path), metadata_test=str(train_meta_path),
            config_name=f"{config_name}_train", target=target,
            negatives=negatives, line_labels=line_labels,
            weight_constraint=weight_constraint, class_weight=class_weight,
        )
    except Exception as e:
        # A degenerate model can fail the fit (NNLS returns all-zero weights).
        # That must not cost this fold its other scorings.
        print(f"  WARN: train-side weight fit failed: {e}")

    weights_json = train_dir / "weights_learned_post.json"
    if weights_json.exists():
        shutil.copy2(weights_json, fold_dir / TRAINFIT_WEIGHTS)
        with open(weights_json) as fh:
            weights_info = json.load(fh)
    else:
        # No fallback to fitting on the held-out table: one sample carries one
        # class, so nothing can be fit, and uniform weights reported under the
        # `learned` name is what the train-fit path exists to remove.
        print("  WARN: no weights.json from the train pass; SKIPPING "
              "weighted_multi_learned_post for this fold")
        weights_json = "SKIP"

    # Distances on the HELD-OUT sample, scored with the train-fit weights.
    outputs = _sweep.run_distances_and_score(
        model_prefix=model_prefix, cfg_dir=fold_dir, starspace=starspace,
        universe=universe, data_path=data_path, label_col=label_col,
        metadata_train=str(train_meta_path), metadata_test=str(test_meta_path),
        config_name=config_name, target=target, negatives=negatives,
        line_labels=line_labels, weights_from=weights_json,
        weight_constraint=weight_constraint, class_weight=class_weight,
    )

    # Held-out sub-scores, for the out-of-fold comparison. Re-runs the uniform
    # scoring (identical output, ~1s) purely to get --subscores-out:
    # run_distances_and_score dumps sub-scores only on the learn-weights branch,
    # which this pass no longer takes.
    dest = fold_dir / "reninness_score_weighted_multi_uniform.csv"
    try:
        cmd = ["geniml", "bedspace", "score",
               "-i", str(fold_dir / "raw_cosdist_rl.csv"), "-o", str(dest),
               "--target", target, "--method", "weighted_multi",
               "--weight-mode", "uniform",
               "--subscores-out", str(fold_dir / SUBSCORES_NAME)]
        if negatives:
            cmd += ["--negatives", negatives]
        run(cmd)
    except Exception as e:
        print(f"  WARN: sub-score dump failed (out-of-fold comparison will "
              f"skip this fold): {e}")

    # Canonical per-fold name, kept so the resume check and the notebook work.
    score_out = fold_dir / "reninness_score.csv"
    primary = outputs.get(method_alias(method)) or next(iter(outputs.values()), None)
    if primary is not None and Path(primary).exists():
        shutil.copy2(primary, score_out)
    if not score_out.exists():
        print(f"  WARN: {score_out} was not produced for fold {held_out_name}")
        return None, weights_info
    return score_out, weights_info


def _flatten_weights(info: dict) -> dict:
    """weights.json -> one flat row (w_<negative> + fit diagnostics)."""
    row = {f"w_{k}": v for k, v in info.get("weights", {}).items()}
    for k in ("n_positive", "n_negative", "mse", "auc"):
        if k in info:
            row[k] = info[k]
    return row


def cache_fold_artifacts(fold_dir: Path, config_name: str, artifact_dir: Path) -> None:
    """Stash the per-fold model.tsv and held-out doc-embed outside the fold dir
    so they survive the post-fold cleanup. These are the inputs the notebook's
    ``get_embedding_matrix`` consumes to reconstruct sample-vs-label embeddings
    for LOOCV-wide plots (heatmap / scatter / FC).

    `geniml bedspace train` writes the model as ``starspace_trained_model.tsv``;
    `geniml bedspace distances -p <config>` writes the held-out doc embed as
    ``<config>_starspace_embed.txt``.
    """
    artifact_dir.mkdir(parents=True, exist_ok=True)
    model_tsv = next(fold_dir.glob("*.tsv"), None)
    if model_tsv is not None:
        shutil.copy2(model_tsv, artifact_dir / "model.tsv")
    doc_embed = fold_dir / f"{config_name}_starspace_embed.txt"
    if doc_embed.exists():
        shutil.copy2(doc_embed, artifact_dir / "held_doc_embed.txt")


def aggregate_loocv_artifacts(
    artifact_root: Path,
    config_name: str,
    output_dir: Path,
) -> None:
    """Build LOOCV-wide model.tsv + concatenated held-out doc-embed for plotting.

    Writes
    ------
    <output_dir>/<config_name>_starspace_embed.txt
        StarSpace ``embed_doc`` 4-row header (copied from the first fold) +
        ``(input_row, embedding_row)`` pairs concatenated from every fold's
        single held-out sample. Parsed by ``get_sample_embedding`` in the
        notebook (`skiprows=4` then `iloc[1::2]`).
    <output_dir>/model_<config_name>.tsv
        First-fold's header + region rows, with the trailing 2 ``__label__*``
        rows replaced by per-label means across all folds. Compatible with
        ``get_label_embedding(path, no_labels=2)`` (reads only the last 2 rows).

    Region rows are taken from the first fold's model and are only meaningful
    for plots that consume per-region embeddings (FC scatter, region+sample
    background). The heatmap path uses only the sample + label rows, so it is
    unaffected by which fold the regions come from.
    """
    fold_dirs = sorted(d for d in artifact_root.iterdir() if d.is_dir())
    if not fold_dirs:
        print(f"  No LOOCV artifact folds in {artifact_root}; nothing to aggregate.")
        return

    # ---- concatenate held-out doc embeddings ----
    docs_out = output_dir / f"{config_name}_starspace_embed.txt"
    header_lines: list[str] | None = None
    data_lines: list[str] = []
    n_doc_folds = 0
    for d in fold_dirs:
        f = d / "held_doc_embed.txt"
        if not f.exists():
            continue
        with open(f) as fh:
            lines = fh.readlines()
        if len(lines) < 6:
            print(f"    skip {d.name}: doc embed has only {len(lines)} lines")
            continue
        if header_lines is None:
            header_lines = lines[:4]
        data_lines.extend(lines[4:])
        n_doc_folds += 1

    if header_lines is None:
        print(f"  No held-out doc embeddings found in {artifact_root}.")
    else:
        with open(docs_out, "w") as fh:
            fh.writelines(header_lines + data_lines)
        print(
            f"  Wrote LOOCV doc embed: {docs_out}  "
            f"(4-row header + {n_doc_folds} held-out samples × 2 rows)"
        )

    # ---- average label rows across folds + reuse fold-0 region rows ----
    base_model_path = fold_dirs[0] / "model.tsv"
    if not base_model_path.exists():
        print(f"  Base model.tsv missing in {fold_dirs[0]}; skipping model aggregation.")
        return
    base_df = pd.read_csv(base_model_path, sep="\t", header=None)

    label_acc: dict[str, list[np.ndarray]] = {}
    label_order: list[str] = []
    for d in fold_dirs:
        m = d / "model.tsv"
        if not m.exists():
            continue
        df = pd.read_csv(m, sep="\t", header=None)
        # Select the label rows by their `__label__` prefix, NOT by position.
        # The model has one label row per distinct label, which is 2 for a single
        # `label` column but 3 for `--label-col Label_1,Label_2`
        # (Renin/Control/Tumoral). A hardcoded tail(2) averaged only two of the
        # three and left the third sitting in the region block of the aggregated
        # model. This is the same prefix mask 7a's `load_model_embeddings` uses.
        labels = df[df[0].astype(str).str.startswith(LABEL_PREFIX)]
        for _, row in labels.iterrows():
            name = str(row.iloc[0])
            vec = row.iloc[1:].astype(float).to_numpy()
            label_acc.setdefault(name, []).append(vec)
            if name not in label_order:
                label_order.append(name)

    if not label_acc:
        print(f"  No label rows in fold models; skipping model aggregation.")
        return

    n_value_cols = base_df.shape[1] - 1
    avg_rows = []
    for name in label_order:
        vecs = np.stack(label_acc[name])
        avg = vecs.mean(axis=0)
        if len(avg) != n_value_cols:
            raise ValueError(
                f"label '{name}' avg vector dim {len(avg)} != "
                f"model.tsv value cols {n_value_cols}"
            )
        avg_rows.append([name] + list(avg))

    # Region rows are everything that is NOT a label row -- again by prefix, so
    # the region block stays clean whatever the label count is.
    preserved = base_df[~base_df[0].astype(str).str.startswith(LABEL_PREFIX)]
    avg_df = pd.DataFrame(avg_rows, columns=base_df.columns)
    combined = pd.concat([preserved, avg_df], ignore_index=True)
    model_out = output_dir / f"model_{config_name}.tsv"
    combined.to_csv(model_out, sep="\t", index=False, header=False)
    print(
        f"  Wrote LOOCV model: {model_out}  "
        f"({len(preserved)} region rows + {len(avg_df)} averaged label rows "
        f"{[n.replace(LABEL_PREFIX, '') for n in label_order]}; "
        f"averaged across {len(fold_dirs)} folds)"
    )


def learn_weights_loocv(
    fold_subscore_paths: list[tuple[str, Path]],
    negatives: list[str],
    out_dir: Path,
    constraint: str = "nonneg-scaled",
    class_weight: str = "none",
    min_train: int = 3,
) -> Path | None:
    """Out-of-fold analogue of `bedspace score --learn-weights`, for LOOCV.

    `bedspace score --learn-weights` fits the two `weighted_multi` weights on the
    same table it scores, which is in-sample. Under LOOCV we can do better for
    free: each fold already wrote its held-out sample's per-negative sub-scores
    (comparable across folds, since each is divided by its own fold's reference
    distance `dl`), so for fold *i* we fit the weights on the other N-1 folds'
    held-out sub-scores and apply them only to fold *i*. No sample ever
    contributes to the weights used to score it.

    Writes
    ------
    <out_dir>/reninness_score_loocv_learned.csv
        filename, file_label, truth, <one column per negative>, w_<negative>...,
        score  -- one row per fold, scored with that fold's out-of-fold weights.
    <out_dir>/loocv_learned_weights.json
        Per-fold weights plus the all-folds fit (the weights you would get from
        `score --learn-weights` on the pooled table) and the LOOCV AUC.
    """
    from geniml.bedspace.weights import learn_weights_from_subscores

    frames = []
    for fold_name, path in fold_subscore_paths:
        df = pd.read_csv(path)
        df["fold"] = fold_name
        frames.append(df)
    if not frames:
        print("  No per-fold sub-scores found; skipping out-of-fold weight learning.")
        return None

    pooled = pd.concat(frames, ignore_index=True)
    missing = [n for n in negatives if n not in pooled.columns]
    if missing:
        print(f"  Sub-score files missing negative column(s) {missing}; skipping.")
        return None

    # Keep only labeled, real samples (anchor rows never make it into the
    # sub-score dump, but a fold whose truth could not be resolved will be NaN).
    pooled = pooled[pooled["truth"].notna()].reset_index(drop=True)
    if pooled["filename"].duplicated().any():
        dups = pooled.loc[pooled["filename"].duplicated(), "filename"].tolist()
        print(f"  WARN: duplicate held-out filenames across folds: {dups[:5]}")
    n_pos = int((pooled["truth"] > 0).sum())
    n_neg = int((pooled["truth"] < 0).sum())
    print(f"  Pooled held-out sub-scores: {len(pooled)} folds "
          f"({n_pos} positive, {n_neg} negative)")
    if n_pos < min_train + 1 or n_neg < min_train + 1:
        print(f"  Not enough labeled folds to leave one out and still fit "
              f"(need > {min_train} per class); skipping.")
        return None

    S_all = pooled.set_index("filename")[negatives]
    truth_all = pooled.set_index("filename")["truth"]

    per_fold: dict[str, dict] = {}
    scores: list[float] = []
    weight_cols: dict[str, list[float]] = {n: [] for n in negatives}
    for _, row in pooled.iterrows():
        fn = row["filename"]
        others = [f for f in S_all.index if f != fn]
        w, meta = learn_weights_from_subscores(
            S_all.loc[others],
            truth_all.loc[others],
            negatives,
            constraint=constraint,
            class_weight=class_weight,
            min_train=min_train,
            warn_imbalance=False,  # reported once by the all-folds fit below
        )
        scores.append(float(sum(w[n] * float(row[n]) for n in negatives)))
        for n in negatives:
            weight_cols[n].append(w[n])
        per_fold[str(row["fold"])] = {
            "held_out": fn,
            "weights": w,
            "n_positive": meta["n_positive"],
            "n_negative": meta["n_negative"],
            "mse": meta["mse"],
        }

    out = pooled[["filename"] + ([c for c in ["file_label"] if c in pooled.columns])
                 + ["truth"] + negatives].copy()
    for n in negatives:
        out[f"w_{n}"] = weight_cols[n]
    out["score"] = np.round(scores, 4)
    out_path = out_dir / "reninness_score_loocv_learned.csv"
    out.to_csv(out_path, index=False)

    # All-folds ("in-sample") fit, for comparison with the out-of-fold numbers.
    w_all, meta_all = learn_weights_from_subscores(
        S_all, truth_all, negatives,
        constraint=constraint, class_weight=class_weight, min_train=min_train,
    )
    summary = {
        "negatives": list(negatives),
        "constraint": constraint,
        "class_weight": class_weight,
        "all_folds_fit": {"weights": w_all, **{k: meta_all[k] for k in
                          ("n_positive", "n_negative", "mse", "auc")}},
        "per_fold_out_of_fold": per_fold,
    }
    try:
        from sklearn.metrics import roc_auc_score

        summary["loocv_auc"] = float(
            roc_auc_score((out["truth"] > 0).astype(int), out["score"])
        )
    except Exception:
        summary["loocv_auc"] = None

    json_path = out_dir / "loocv_learned_weights.json"
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"  Wrote out-of-fold learned scores: {out_path}")
    print(f"  Wrote weights + diagnostics:      {json_path}")
    print(f"  All-folds weights (in-sample):    {w_all}  (auc={meta_all['auc']})")
    print(f"  LOOCV (out-of-fold) AUC:          {summary['loocv_auc']}")
    return out_path


def evaluate_loocv(output_dir: Path, config_name: str, target: str,
                   score_tables: list[tuple[str, Path]],
                   model_tsv: Path | None = None,
                   geniml_eval: bool = True, geniml_eval_workers: int = 10,
                   bin_embed: str | None = None) -> Path:
    """Both accuracy families (+ both clusterability families) for a LOOCV run.

    A single fold holds out one sample, so no fold has two classes and nothing
    can be scored per fold. Everything here is therefore computed on the
    CONCATENATED tables, which is exactly what LOOCV is for: every row was
    scored by a model that never saw it, so these are out-of-sample numbers --
    unlike the in-sample numbers `score --learn-weights` reports on its own
    table.

    Clusterability is computed once on the LOOCV-aggregated model
    (`model_<config>.tsv`), whose label rows are averaged across folds; the
    region rows come from fold 0 (see `aggregate_loocv_artifacts`).

    Writes <output_dir>/loocv_eval_metrics.{csv,json}.
    """
    rows = []
    for label, path in score_tables:
        if not Path(path).exists():
            continue
        try:
            m = _sweep.score_table_metrics(str(path), target=target)
        except Exception as e:
            print(f"  metrics failed for {label}: {e}")
            continue
        if m is None:
            print(f"  {label}: only one class present — skipped")
            continue
        rows.append(dict(config=config_name, table=label, **m))
        print(f"  {label:28s} acc_sign={m['acc_sign']:.4f} bal_acc={m['bal_acc']:.4f} "
              f"prec_rank={m['prec_rank_overall']:.4f} auc={m['auc']:.4f}")

    cl = {}
    if geniml_eval and model_tsv and Path(model_tsv).exists():
        print(f"  geniml.eval on {Path(model_tsv).name} ...", flush=True)
        try:
            cl = _sweep.geniml_eval_metrics(
                str(model_tsv),
                cache_path=str(output_dir / f"{config_name}_base_embed.pt"),
                num_workers=geniml_eval_workers, bin_embed=bin_embed,
            )
            print("  " + "  ".join(f"{k}={v:.4f}" for k, v in cl.items()
                                   if isinstance(v, float)))
        except Exception as e:
            print(f"  geniml.eval failed: {e}")
    if model_tsv and Path(model_tsv).exists():
        try:
            emb, _, _, _ = _sweep.load_model_embeddings(str(model_tsv))
            mk = _sweep.multi_k_clustering(emb)
            cl["sil_k2"] = mk.get(2, {}).get("silhouette")
            cl["sil_k3"] = mk.get(3, {}).get("silhouette")
            cl["best_k"] = max(mk, key=lambda k: mk[k]["silhouette"])
            print(f"  silhouette k=2 {cl['sil_k2']:.4f}  k=3 {cl['sil_k3']:.4f}  "
                  f"best_k={cl['best_k']}")
        except Exception as e:
            print(f"  silhouette failed: {e}")

    df = pd.DataFrame(rows)
    for k, v in cl.items():
        df[k] = v
    out_csv = output_dir / "loocv_eval_metrics.csv"
    df.to_csv(out_csv, index=False)
    with open(output_dir / "loocv_eval_metrics.json", "w") as fh:
        json.dump({"config": config_name, "target": target,
                   "accuracy": rows, "clusterability": cl}, fh, indent=2, default=float)
    print(f"  Wrote {out_csv}")
    return out_csv


def concat_scores(fold_score_paths: list[tuple[str, Path]], out_path: Path) -> None:
    """Concatenate per-fold reninness_score.csv and drop duplicates.

    Each fold's CSV contains the held-out test sample plus one anchor row per
    LABEL (filename == file_label), added by `--include-label-distances`. There
    are as many anchors as the model has labels: 2 for a single `label` column,
    3 for `--label-col Label_1,Label_2` (Renin/Control/Tumoral).

    Anchor scores are NOT identical across folds. Every fold trains its own
    model, so only the structurally-fixed anchors (target = +1, its paired
    negative = -1) repeat; the rest -- e.g. Tumoral under `--method pairwise
    --negatives Control` -- take a different value in every fold, and a plain
    drop_duplicates would leave one anchor row per fold (250+ junk rows). They
    are collapsed here to one row per label carrying the across-fold MEAN, which
    is the same treatment `aggregate_loocv_artifacts` gives the label
    embeddings. (Metrics are unaffected either way: `score_table_metrics` drops
    anchor rows before scoring.)

    Expected result: N held-out samples + one row per label.
    """
    frames = [pd.read_csv(p) for _, p in fold_score_paths]
    combined = pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)

    is_anchor = combined["filename"].astype(str) == combined["file_label"].astype(str)
    anchors, samples = combined[is_anchor], combined[~is_anchor]
    n_anchors = anchors["filename"].nunique()
    if len(anchors) > n_anchors:
        collapsed = (anchors.groupby(["filename", "file_label"], as_index=False)
                     .agg({c: ("mean" if pd.api.types.is_numeric_dtype(anchors[c]) else "first")
                           for c in anchors.columns if c not in ("filename", "file_label")}))
        print(f"  Collapsed {len(anchors)} anchor rows -> {len(collapsed)} "
              f"(mean across folds; anchor scores are per-fold model-dependent)")
        combined = pd.concat([samples, collapsed[combined.columns]],
                             ignore_index=True).reset_index(drop=True)
    combined.to_csv(out_path, index=False)

    n_folds = len(frames)
    expected = n_folds + n_anchors
    actual = len(combined)
    print(f"\nConcatenated {n_folds} folds  {out_path}")
    print(f"  Rows: {actual}   (expected {expected} = {n_folds} samples + "
          f"{n_anchors} label rows)")

    if actual != expected:
        print(f"    Row count mismatch.")
        if actual > expected:
            # Inspect what failed to collapse  most often the Renin/Non_Renin
            # label rows differing across folds, or duplicate sample basenames.
            counts = combined.groupby(["filename", "file_label"]).size().reset_index(name="n")
            dup_anchors = counts[counts["filename"] == counts["file_label"]]
            dup_samples = counts[(counts["filename"] != counts["file_label"]) & (counts["n"] > 1)]
            if (dup_anchors["n"] > 1).any():
                print("    Anchor rows did not collapse  score values likely differ "
                      "across folds (rare; usually means `dl` reference distance "
                      "differs). Affected:")
                print(dup_anchors[dup_anchors["n"] > 1].to_string(index=False))
            if len(dup_samples):
                print("    Sample rows appear with multiple (filename, file_label) "
                      "pairs  most likely an old geniml on the cluster missing the "
                      "meta_preprocessing fix, so file_label='unknown' in some folds. "
                      "Affected:")
                print(dup_samples.to_string(index=False))
        else:
            print(f"    Fewer rows than expected  some folds may have produced no "
                  f"test-sample row, or test basenames collided. Inspect per-fold CSVs.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="Config name from SWEEP_CONFIGS (e.g. d3_e10_m0.1)")
    parser.add_argument("--metadata", required=True, help="Full metadata CSV (one row per sample)")
    parser.add_argument("--data-path", required=True, help="Directory holding BED files referenced in metadata")
    parser.add_argument("--universe", required=True, help="Universe BED file")
    parser.add_argument("--starspace", required=True, help="StarSpace install directory")
    parser.add_argument("--label-col", default="Label_1,Label_2",
                        help="Label column name(s), comma-separated (default: Label_1,Label_2)")
    parser.add_argument(
        "--target", default="Renin",
        help="Positive class name passed to `bedspace score --target` (default: 'Renin')",
    )
    parser.add_argument("--output-dir", required=True, help="LOOCV output directory")
    parser.add_argument(
        "--name-col", default="file_name",
        help="Column in metadata that uniquely identifies a sample (used for fold naming). Default: 'file_name'.",
    )
    parser.add_argument(
        "--method", default="weighted_multi", choices=["pairwise", "weighted_multi"],
        help="Which method's table is copied to reninness_score.csv "
        "(default: weighted_multi). All methods are always written.",
    )
    parser.add_argument(
        "--negatives", default="Control,Tumoral",
        help="Comma-separated negatives for weighted_multi. The first is also "
        "the single negative used by pairwise.",
    )
    parser.add_argument(
        "--line-labels", default="Control,Tumoral",
        help="Two negative labels forming the line for the 'line' method. "
        "Unset -> that method is skipped.",
    )
    parser.add_argument(
        "--learn-weights", action="store_true",
        help="Also produce the pooled OUT-OF-FOLD weight table, as a comparison "
        "against the per-fold train-fit weights. Requires --method "
        "weighted_multi. Writes reninness_score_loocv_learned.csv.",
    )
    parser.add_argument(
        "--truth-column", default=None,
        help="(--learn-weights) Column in --metadata holding the TRUE binary "
        "label (e.g. Renin vs Non-Renin), when it differs from the label the "
        "model is trained on. Default: derive it from the label column(s) the "
        "distance table already carries.",
    )
    parser.add_argument(
        "--truth-positive", default=None,
        help="(--learn-weights) Value in --truth-column marking the positive "
        "class (default: --target).",
    )
    parser.add_argument(
        "--weight-constraint", default="nonneg-scaled", choices=["nonneg-scaled", "simplex"],
        help="(--learn-weights) How the two weights are fit. Both return "
        "non-negative weights summing to 1, keeping the score in [-1, 1]. "
        "'nonneg-scaled' (default): fit with the scale free, then rescale -- "
        "rank-identical to the unconstrained fit, and the more accurate of the "
        "two. 'simplex': impose the sum during the fit.",
    )
    parser.add_argument(
        "--class-weight", default="balanced", choices=["none", "balanced"],
        help="Per-sample weighting of the +1/-1 weight-fitting loss "
        "(default: balanced). Renin samples are rare.",
    )
    parser.add_argument(
        "--skip-eval", action="store_true",
        help="Skip the end-of-run evaluation metrics (both accuracy families "
        "and both clusterability families) written to loocv_eval_metrics.csv.",
    )
    parser.add_argument(
        "--skip-geniml-eval", action="store_true",
        help="Within the evaluation step, skip the geniml.eval region-embedding "
        "tests (CTT/GDST/NPT) and report the silhouette family only.",
    )
    parser.add_argument("--geniml-eval-workers", type=int, default=10,
                        help="Parallel workers for the geniml.eval tests (default: 10)")
    parser.add_argument(
        "--bin-embed", default=None,
        help="Path to a binary-embedding pickle over the same universe; enables "
        "the geniml.eval RCT test. Omitted -> RCT is skipped.",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep per-fold preprocess/model/distance files (default: delete to save space, keep only reninness_score.csv)",
    )
    parser.add_argument(
        "--balance", action="store_true",
        help="Oversample the minority label groups of each FOLD's training "
        "split before StarSpace sees it. Each fold is preprocessed separately, "
        "so this cannot reach the held-out sample. Use it whenever the "
        "sweep/holdout arms being compared against were balanced.",
    )
    parser.add_argument(
        "--balance-strategy", default=_BAL_STRAT, choices=list(_BAL_STRATS),
        help="'labels' (default) equalises LABEL OCCURRENCES, counting a "
        "multi-label document once per label.",
    )
    parser.add_argument("--balance-factors", default=None)
    parser.add_argument("--balance-min-ratio", type=float, default=_BAL_RATIO)
    parser.add_argument("--balance-max-rounds", type=int, default=_BAL_ROUNDS)
    parser.add_argument("--balance-seed", type=int, default=42)
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run every fold from scratch, even if reninness_score.csv already exists. "
        "Default: resume mode  folds whose reninness_score.csv exists are skipped, missing ones are re-run.",
    )
    args = parser.parse_args()

    if args.config not in SWEEP_CONFIGS:
        print(f"Error: config '{args.config}' not found in SWEEP_CONFIGS.", file=sys.stderr)
        print(f"Available: {sorted(SWEEP_CONFIGS.keys())}", file=sys.stderr)
        return 1
    params = SWEEP_CONFIGS[args.config]
    print(f"Config: {args.config}  ({params['description']})")
    print("Weights: fit on each fold's TRAIN split, applied to its held-out sample.")

    meta_df = pd.read_csv(args.metadata)
    if args.name_col not in meta_df.columns:
        print(f"Error: --name-col '{args.name_col}' not in metadata columns {list(meta_df.columns)}", file=sys.stderr)
        return 1
    # --label-col may name several columns (comma-separated, as `bedspace
    # distances -l` accepts) for a multi-label model, e.g. "Label_1,Label_2".
    label_cols = [c.strip() for c in str(args.label_col).split(",") if c.strip()]
    missing_label_cols = [c for c in label_cols if c not in meta_df.columns]
    if missing_label_cols:
        print(f"Error: --label-col {missing_label_cols} not in metadata columns {list(meta_df.columns)}", file=sys.stderr)
        return 1

    if args.learn_weights and args.method != "weighted_multi":
        print("Error: --learn-weights requires --method weighted_multi "
              "(the two learned weights combine a three-label model's "
              "per-negative sub-scores).", file=sys.stderr)
        return 1
    if args.learn_weights and args.truth_column and args.truth_column not in meta_df.columns:
        print(f"Error: --truth-column '{args.truth_column}' not in metadata columns "
              f"{list(meta_df.columns)}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_root = output_dir / "_loocv_artifacts"

    fold_scores: list[tuple[str, Path]] = []
    weight_rows: list[dict] = []
    failed_folds: list[str] = []
    resumed_count = 0
    ran_count = 0
    n = len(meta_df)
    meta_df_reset = meta_df.reset_index(drop=True)

    for i, row in meta_df_reset.iterrows():
        # Derive a filesystem-safe fold name from the held-out sample's name.
        raw_name = str(row[args.name_col])
        held_out = Path(raw_name).stem.replace(".bed", "").replace(".gz", "")

        fold_dir = output_dir / f"fold_{held_out}"
        existing_score = fold_dir / "reninness_score.csv"

        # ---- Resume logic ----
        # If a previous run already produced reninness_score.csv for this fold,
        # skip re-running it (unless --force). This lets you re-invoke the
        # script after fixing failures and only the failed folds get redone.
        if existing_score.exists() and not args.force:
            try:
                _ = pd.read_csv(existing_score, nrows=1)  # sanity check
            except Exception as e:
                print(f"\n{'='*70}\nFOLD {i+1}/{n}: {held_out}  existing score "
                      f"unreadable ({e}); re-running\n{'='*70}")
            else:
                print(f"\nFOLD {i+1}/{n}: {held_out}  already done, skipping "
                      f"(use --force to redo)")
                fold_scores.append((held_out, existing_score))
                wpath = fold_dir / TRAINFIT_WEIGHTS
                if wpath.exists():
                    with open(wpath) as fh:
                        weight_rows.append(
                            dict(fold=held_out, **_flatten_weights(json.load(fh))))
                resumed_count += 1
                continue

        print(f"\n{'='*70}\nFOLD {i+1}/{n}: held out = {held_out}\n{'='*70}")

        # If --force on a previously-attempted fold, wipe the stale fold dir
        # so we don't trip over leftover models / partial preprocess output.
        if args.force and fold_dir.exists():
            shutil.rmtree(fold_dir)

        train_meta = fold_dir / "train_meta.csv"
        test_meta = fold_dir / "test_meta.csv"
        fold_dir.mkdir(parents=True, exist_ok=True)
        write_fold_metadata(meta_df_reset, i, train_meta, test_meta)

        try:
            score_path, winfo = run_fold(
                config_name=args.config,
                params=params,
                held_out_name=held_out,
                train_meta_path=train_meta,
                test_meta_path=test_meta,
                data_path=args.data_path,
                universe=args.universe,
                starspace=args.starspace,
                label_col=args.label_col,
                fold_dir=fold_dir,
                target=args.target,
                method=args.method,
                negatives=args.negatives,
                line_labels=args.line_labels or None,
                weight_constraint=args.weight_constraint,
                class_weight=args.class_weight,
                balance=args.balance,
                balance_strategy=args.balance_strategy,
                balance_factors=args.balance_factors,
                balance_min_ratio=args.balance_min_ratio,
                balance_max_rounds=args.balance_max_rounds,
                balance_seed=args.balance_seed,
            )
        except Exception as e:
            print(f"  FOLD FAILED ({held_out}): {e}", flush=True)
            failed_folds.append(held_out)
            continue

        if score_path is None:
            failed_folds.append(held_out)
            continue

        fold_scores.append((held_out, score_path))
        if winfo:
            weight_rows.append(dict(fold=held_out, **_flatten_weights(winfo)))
        ran_count += 1

        # Stash model.tsv + held-out doc embed outside fold_dir BEFORE the
        # cleanup below wipes them. These feed the LOOCV-wide aggregation at
        # the end so the notebook can build a sample-vs-label embedding.
        cache_fold_artifacts(fold_dir, args.config, artifact_root / held_out)

        # Trim intermediates to keep disk usage reasonable across many folds.
        # Keep reninness_score.csv (and subscores.csv, which the out-of-fold
        # weight fit at the end pools across folds).
        if not args.keep_intermediates:
            keep = {"reninness_score.csv", SUBSCORES_NAME, TRAINFIT_WEIGHTS} | {
                f"reninness_score_{m}.csv" for m in METHODS
            }
            for child in fold_dir.iterdir():
                if child.is_dir():
                    shutil.rmtree(child)
                elif child.name not in keep:
                    child.unlink()

    # ---- Final status report ----
    print(f"\n{'='*70}\nLOOCV SUMMARY\n{'='*70}")
    print(f"  Total folds:      {n}")
    print(f"  Resumed (skipped): {resumed_count}")
    print(f"  Newly run:        {ran_count}")
    print(f"  Failed:           {len(failed_folds)}")
    if failed_folds:
        print(f"  Failed fold names: {failed_folds}")
        print(f"   Re-run this script to retry only the failed folds "
              f"(successful folds will be auto-skipped).")

    if not fold_scores:
        print("\nNo folds produced scores; nothing to concatenate.", file=sys.stderr)
        return 1

    if weight_rows:
        wdf = pd.DataFrame(weight_rows)
        wpath = output_dir / "loocv_trainfit_weights.csv"
        wdf.to_csv(wpath, index=False)
        wcols = [c for c in wdf.columns if c.startswith("w_")]
        print(f"\n  Train-fit weights across {len(wdf)} folds -> {wpath}")
        print(wdf[wcols + [c for c in ("n_positive", "n_negative", "auc") if c in wdf]]
              .describe().loc[["mean", "std", "min", "max"]].round(4).to_string())

    out_path = output_dir / "reninness_score_loocv.csv"
    concat_scores(fold_scores, out_path)

    # Per-method concatenations, so the LOOCV run reports the same four scoring
    # methods 7a and 7c do rather than only the one --method selected.
    per_method: list[tuple[str, Path]] = []
    for m in METHODS:
        paths = [(n, p.parent / f"reninness_score_{m}.csv") for n, p in fold_scores
                 if (p.parent / f"reninness_score_{m}.csv").exists()]
        if not paths:
            continue
        dest = output_dir / f"reninness_score_loocv_{m}.csv"
        concat_scores(paths, dest)
        per_method.append((m, dest))

    if args.learn_weights:
        print(f"\n{'='*70}\nOUT-OF-FOLD WEIGHTS (comparison only)\n{'='*70}")
        fold_subscores = [
            (name, p.parent / SUBSCORES_NAME)
            for name, p in fold_scores
            if (p.parent / SUBSCORES_NAME).exists()
        ]
        n_missing = len(fold_scores) - len(fold_subscores)
        if n_missing:
            print(f"  {n_missing}/{len(fold_scores)} folds have no {SUBSCORES_NAME} "
                  f"(scored before --method weighted_multi was used, or resumed "
                  f"from an older run); re-run those folds with --force.")
        if fold_subscores:
            if args.negatives:
                negatives = [n.strip() for n in args.negatives.split(",") if n.strip()]
            else:
                cols = pd.read_csv(fold_subscores[0][1], nrows=0).columns
                negatives = [
                    c for c in cols if c not in ("filename", "file_label", "truth", "fold")
                ]
            print(f"  Negatives: {negatives}")
            learn_weights_loocv(
                fold_subscores,
                negatives,
                output_dir,
                constraint=args.weight_constraint,
                class_weight=args.class_weight,
            )
        else:
            print("  No sub-score files at all; skipping.")

    print(f"\n{'='*70}\nAGGREGATING EMBEDDINGS FOR PLOTTING\n{'='*70}")
    if artifact_root.exists():
        aggregate_loocv_artifacts(artifact_root, args.config, output_dir)
    else:
        print(
            f"  No artifact dir at {artifact_root} - likely all folds were "
            f"resumed from a prior run that pre-dated artifact caching. "
            f"Re-run with --force to regenerate."
        )

    # Runs whether or not the embedding artifacts exist: the accuracy metrics
    # need only the concatenated score tables. Clusterability needs the
    # aggregated model, and evaluate_loocv skips it when that is absent.
    if not args.skip_eval:
        print(f"\n{'='*70}\nEVALUATION METRICS\n{'='*70}")
        tables = [(m, p) for m, p in per_method] or [("loocv", out_path)]
        learned = output_dir / "reninness_score_loocv_learned.csv"
        if learned.exists():
            # Pooled comparison. `weighted_multi_learned_post` above is the
            # per-fold TRAIN-FIT result, which is what the paper reports.
            tables.append(("weighted_multi_learned_out_of_fold", learned))
        evaluate_loocv(
            output_dir, args.config, args.target, tables,
            model_tsv=output_dir / f"model_{args.config}.tsv",
            geniml_eval=not args.skip_geniml_eval,
            geniml_eval_workers=args.geniml_eval_workers,
            bin_embed=args.bin_embed,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

