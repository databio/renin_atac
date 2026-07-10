#!/usr/bin/env python
"""
7b  Leave-one-out cross-validation (LOOCV) driver for the bedspace pipeline.

For a given param combo from ``SWEEP_CONFIGS`` (loaded from the co-located
``7a_run_reninness_score_in_batch.py``), loops over each sample in the
metadata, holds it out as the test set, runs preprocess  train  distances
 score on the remaining N-1 samples, and concatenates per-fold reninness
scores into a single CSV.

------------------------------------------------------------------
USAGE
------------------------------------------------------------------

    source /home/bingjiexue/Documents/code/reninness_workspace/.venv/bin/activate
    python 7b_run_reninness_score_loo_in_batch.py \\
        --config     d3_e10_m0.1 \\
        --metadata   /path/to/renin_encode_atac.csv \\
        --data-path  /path/to/bed/files/ \\
        --universe   /path/to/universe.bed \\
        --starspace  /path/to/Starspace/ \\
        --label-col  label \\
        --target     Renin \\
        --output-dir /project/shefflab/.../universe_cc/d3_e10_m0.1/loocv/

Output:
    <output-dir>/
        fold_<sample_name>/
            reninness_score.csv      (intermediates removed by default;
                                      pass --keep-intermediates to keep them)
        reninness_score_loocv.csv    concatenated; duplicates dropped, so the
                                     Renin and Non_Renin label rows that appear
                                     in every fold collapse to one row each.
                                     Final shape: (N held-out samples + 2 label
                                     anchor rows) × 3 cols (filename, file_label,
                                     score).
        <config>_starspace_embed.txt LOOCV-wide held-out doc embeddings 
                                     4-row header + (input, embedding) pairs
                                     concatenated across folds. Parsed by the
                                     notebook's `get_sample_embedding`.
        model_<config>.tsv           First-fold header/regions + per-label
                                     means across folds (last 2 rows).
                                     Consumed by `get_label_embedding`.
        _loocv_artifacts/            Per-fold cached model.tsv + held-out doc
                                     embed (survives fold cleanup; inputs to
                                     the aggregation step).

Resume behavior:
    By default the script is resumable. A fold is considered "done" if its
    fold_<name>/reninness_score.csv exists and is readable. Re-running the
    same command will skip done folds and only attempt the failed/missing
    ones. End-of-run summary lists which folds were resumed, run fresh, and
    which failed. Pass --force to wipe every fold dir and re-run from scratch.
"""

from __future__ import annotations

import argparse
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Load SWEEP_CONFIGS from the co-located 7a file. We can't `import` it
# directly because the module name starts with a digit, so use importlib.
_THIS_DIR = Path(__file__).resolve().parent
_SWEEP_FILE = _THIS_DIR / "7a_run_reninness_score_in_batch.py"

if not _SWEEP_FILE.exists():
    print(f"ERROR: Expected sibling config file not found: {_SWEEP_FILE}", file=sys.stderr)
    print("This script reads SWEEP_CONFIGS from 7a_run_reninness_score_in_batch.py.", file=sys.stderr)
    sys.exit(1)

_spec = importlib.util.spec_from_file_location("_sweep_configs_module", _SWEEP_FILE)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
SWEEP_CONFIGS = _mod.SWEEP_CONFIGS


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
) -> Path | None:
    """Run one LOOCV fold; return path to that fold's reninness_score.csv."""
    fold_dir.mkdir(parents=True, exist_ok=True)
    preprocess_dir = fold_dir / "preprocessed"

    # ---- preprocess (only train data is preprocessed; test files are
    #      tokenized on the fly inside `distances`) ----
    run([
        "geniml", "bedspace", "preprocess",
        "-i", data_path,
        "-m", str(train_meta_path),
        "-u", universe,
        "-o", str(preprocess_dir),
        "-l", label_col,
        "--mode", "bulk",
    ])

    # ---- train ----
    train_input = preprocess_dir / "train_input.txt"
    if not train_input.exists():
        raise FileNotFoundError(f"preprocess did not produce {train_input}")

    train_cmd = [
        "geniml", "bedspace", "train",
        "-s", starspace,
        "-i", str(train_input),
        "-o", str(fold_dir),
        "-n", str(params["epoch"]),
        "-d", str(params["dim"]),
        "-l", str(params["lr"]),
        "--neg-search-limit", str(params["negSearchLimit"]),
        "--margin", str(params["margin"]),
        "--min-count", str(params["minCount"]),
    ]
    train_cmd.append("--adagrad" if params.get("adagrad") else "--no-adagrad")
    run(train_cmd)

    # Locate the trained model. `geniml bedspace train` writes a file named
    # like model_<input>_<params>.tsv; we glob to find it.
    candidates = sorted(fold_dir.glob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No model .tsv found in {fold_dir}")
    model_path = candidates[0]
    model_prefix = str(model_path)[:-4]  # strip .tsv

    # ---- distances ----
    run([
        "geniml", "bedspace", "distances",
        "-i", model_prefix,
        "-s", starspace,
        "--metadata-train", str(train_meta_path),
        "--metadata-test", str(test_meta_path),
        "-u", universe,
        "-p", config_name,
        "-f", data_path,
        "-l", label_col,
        "-o", str(fold_dir),
        "--include-label-distances",
    ])

    # ---- score ----
    raw_dist = fold_dir / "raw_cosdist_rl.csv"
    score_out = fold_dir / "reninness_score.csv"
    score_cmd = [
        "geniml", "bedspace", "score",
        "-i", str(raw_dist),
        "-o", str(score_out),
    ]
    if target:
        score_cmd.extend(["--target", target])
    run(score_cmd)

    if not score_out.exists():
        print(f"  WARN: {score_out} was not produced for fold {held_out_name}")
        return None
    return score_out


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
        labels = df.tail(2)
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

    preserved = base_df.iloc[:-2]  # fold-0 header (if any) + region rows
    avg_df = pd.DataFrame(avg_rows, columns=base_df.columns)
    combined = pd.concat([preserved, avg_df], ignore_index=True)
    model_out = output_dir / f"model_{config_name}.tsv"
    combined.to_csv(model_out, sep="\t", index=False, header=False)
    print(
        f"  Wrote LOOCV model: {model_out}  "
        f"({len(preserved) - 1} region rows + {len(avg_df)} averaged label rows; "
        f"averaged across {len(fold_dirs)} folds)"
    )


def concat_scores(fold_score_paths: list[tuple[str, Path]], out_path: Path) -> None:
    """Concatenate per-fold reninness_score.csv and drop duplicates.

    Each fold's CSV contains the held-out test sample plus the two label
    rows (Renin / Non_Renin) that come from `--include-label-distances`.
    The label rows are identical across all folds, so we dedupe after concat.
    Expected result: N held-out samples + 2 label rows.
    """
    frames = [pd.read_csv(p) for _, p in fold_score_paths]
    combined = pd.concat(frames, ignore_index=True).drop_duplicates().reset_index(drop=True)
    combined.to_csv(out_path, index=False)

    n_folds = len(frames)
    expected = n_folds + 2
    actual = len(combined)
    print(f"\nConcatenated {n_folds} folds  {out_path}")
    print(f"  Rows: {actual}   (expected {expected} = {n_folds} samples + 2 label rows)")

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
    parser.add_argument("--label-col", default="label", help="Label column name in metadata (default: 'label')")
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
        "--keep-intermediates", action="store_true",
        help="Keep per-fold preprocess/model/distance files (default: delete to save space, keep only reninness_score.csv)",
    )
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

    meta_df = pd.read_csv(args.metadata)
    if args.name_col not in meta_df.columns:
        print(f"Error: --name-col '{args.name_col}' not in metadata columns {list(meta_df.columns)}", file=sys.stderr)
        return 1
    if args.label_col not in meta_df.columns:
        print(f"Error: --label-col '{args.label_col}' not in metadata columns {list(meta_df.columns)}", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    artifact_root = output_dir / "_loocv_artifacts"

    fold_scores: list[tuple[str, Path]] = []
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
            score_path = run_fold(
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
            )
        except Exception as e:
            print(f"  FOLD FAILED ({held_out}): {e}", flush=True)
            failed_folds.append(held_out)
            continue

        if score_path is None:
            failed_folds.append(held_out)
            continue

        fold_scores.append((held_out, score_path))
        ran_count += 1

        # Stash model.tsv + held-out doc embed outside fold_dir BEFORE the
        # cleanup below wipes them. These feed the LOOCV-wide aggregation at
        # the end so the notebook can build a sample-vs-label embedding.
        cache_fold_artifacts(fold_dir, args.config, artifact_root / held_out)

        # Trim intermediates to keep disk usage reasonable across many folds.
        # Keep only reninness_score.csv in the fold dir.
        if not args.keep_intermediates:
            for child in fold_dir.iterdir():
                if child.is_dir():
                    shutil.rmtree(child)
                elif child.name != "reninness_score.csv":
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

    out_path = output_dir / "reninness_score_loocv.csv"
    concat_scores(fold_scores, out_path)

    print(f"\n{'='*70}\nAGGREGATING EMBEDDINGS FOR PLOTTING\n{'='*70}")
    if artifact_root.exists():
        aggregate_loocv_artifacts(artifact_root, args.config, output_dir)
    else:
        print(
            f"  No artifact dir at {artifact_root}  likely all folds were "
            f"resumed from a prior run that pre-dated artifact caching. "
            f"Re-run with --force to regenerate."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

