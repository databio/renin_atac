#!/usr/bin/env python
"""
7c — Score pre-aggregated single-cell BED files (per-cluster or per-barcode
pseudo-bulks) against already-trained StarSpace models.

Cells are assumed to be already converted into BED files (one BED per cluster
or per barcode), so this script runs the bulk-mode distance + score steps just
like 7a does for true bulk samples — only the model train step is skipped.

For each config in CONFIGS:
    geniml bedspace distances --mode bulk  →  raw_cosdist_rl.csv
    geniml bedspace score                  →  reninness_score.csv

------------------------------------------------------------------
USAGE
------------------------------------------------------------------

    source /home/bingjiexue/Documents/code/reninness_workspace/.venv/bin/activate
    python 7c_run_reninness_score_sc_in_batch.py \\
        --test-meta    /path/to/test_metadata.csv \\
        --output-dir   /project/shefflab/.../universe_cc/ \\
        --result-tag   sc \\
        --bed-dir      /path/to/per_cluster_or_per_barcode_beds/ \\
        --train-meta   /path/to/train_metadata.csv \\
        --models-root  /project/shefflab/.../universe_cc/ \\
        --starspace    /path/to/Starspace/ \\
        --universe     /path/to/universe.bed

Per-config output:
    <output-dir>/<config>/<result-tag>/
        <config>_starspace_embed.txt        test-sample embeddings
        <config>_train_starspace_embed.txt  train-sample embeddings
        raw_cosdist_rl.csv                  sample × label cosine distances
        similarity_score_rl.csv             scaled sample × label similarities
        similarity_score_rr.csv             scaled test × train similarities
        reninness_score.csv                 per-sample reninness score [-1, 1]

Trained-model layout (read-only):
    <models-root>/<config>/model_<config>.tsv

Resume behavior:
    A config is considered done if reninness_score.csv exists in its output
    dir. Re-running the script will skip done configs. Use --force to wipe
    each config's output dir and re-run.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Configs to score: both margins (m0.1 and m1.0) per dimension.
# Note: the d10/m1.0 slot uses `uc_d10_e10_m1.0` — the only d10_*_m1.0 entry
# in SWEEP_CONFIGS (no plain `d10_e10_m1.0` was trained).
CONFIGS = [
    "d3_e10_m0.1",
    "d3_e10_m1.0",
    "d5_e10_m0.1",
    "d5_e10_m1.0",
    "d10_e10_m0.1",
    "uc_d10_e10_m1.0",
    "d50_e10_m0.1",
    "d50_e10_m1.0",
]


def run(cmd: list[str], allow_fail: bool = False) -> int:
    """Stream a subprocess. Raise on non-zero unless `allow_fail`; return rc."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0 and not allow_fail:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {' '.join(cmd)}")
    return res.returncode


def run_config(
    config: str,
    test_meta: Path,
    train_meta: Path,
    bed_dir: str,
    models_root: Path,
    starspace: str,
    universe: str,
    label_col: str,
    target: str,
    out_dir: Path,
    analytic: bool,
) -> Path:
    """Run distances + score for one config; return reninness_score.csv path."""
    model_tsv = models_root / config / f"model_{config}.tsv"
    if not model_tsv.exists():
        raise FileNotFoundError(f"Trained model not found: {model_tsv}")

    out_dir.mkdir(parents=True, exist_ok=True)
    model_prefix = str(model_tsv)[:-4]  # strip .tsv

    # ---- distances (bulk mode) ----
    # In analytic mode (default for this script), the rr step is skipped at
    # the geniml level via --analytic, so the call should succeed even when
    # --bed-dir does not contain the training BED files.
    raw_dist = out_dir / "raw_cosdist_rl.csv"
    cmd = [
        "geniml", "bedspace", "distances",
        "-i", model_prefix,
        "-s", starspace,
        "--metadata-train", str(train_meta),
        "--metadata-test", str(test_meta),
        "-u", universe,
        "-p", config,
        "-f", bed_dir,
        "-l", label_col,
        "-o", str(out_dir),
        "--include-label-distances",
    ]
    if analytic:
        cmd.append("--analytic")
    rc = run(cmd, allow_fail=True)
    if rc != 0 and raw_dist.exists():
        print(f"  WARN: distances failed (rc={rc}) but {raw_dist.name} was produced; "
              f"continuing to `score`. similarity_score_rr.csv will be missing — "
              f"pass --mode analytic to skip the rr step cleanly.", flush=True)
    elif rc != 0:
        raise RuntimeError(f"distances failed (rc={rc}) and no {raw_dist.name} was produced")

    # ---- score ----
    if not raw_dist.exists():
        raise FileNotFoundError(f"distances did not produce {raw_dist}")
    score_out = out_dir / "reninness_score.csv"
    score_cmd = [
        "geniml", "bedspace", "score",
        "-i", str(raw_dist),
        "-o", str(score_out),
    ]
    if target:
        score_cmd.extend(["--target", target])
    run(score_cmd)

    if not score_out.exists():
        raise FileNotFoundError(f"score did not produce {score_out}")
    return score_out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--test-meta", required=True,
        help="Path to test metadata CSV (one row per pseudo-bulk BED file).",
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="Root output directory. Per-config outputs go to "
        "<output-dir>/<config>/<result-tag>/.",
    )
    parser.add_argument(
        "--result-tag", required=True,
        help="Sub-dir under each config dir for this run's outputs (e.g. 'sc', "
        "'pseudobulk'). Lets you co-locate sc / pseudobulk evals with the "
        "trained model without clobbering bulk-eval outputs.",
    )
    parser.add_argument(
        "--bed-dir", required=True,
        help="Directory holding the pseudo-bulk BED files referenced in the test metadata.",
    )
    parser.add_argument(
        "--train-meta", required=True,
        help="Path to the training metadata CSV (needed because bulk-mode distances "
        "also computes test-to-train sample similarities).",
    )
    parser.add_argument(
        "--models-root",
        default="/project/shefflab/brickyard/results_pipeline/gomez_atac/results_pipeline/reninness_score_output_new/universe_cc",
        help="Root containing trained models at <models-root>/<config>/model_<config>.tsv.",
    )
    parser.add_argument(
        "--starspace", required=True,
        help="StarSpace install dir (must contain `embed_doc` binary).",
    )
    parser.add_argument(
        "--universe", required=True,
        help="Universe BED used during training.",
    )
    parser.add_argument(
        "--label-col", default="Group",
        help="Label column name in metadata (default: 'Group').",
    )
    parser.add_argument(
        "--target", default="Renin",
        help="Positive class name passed to `bedspace score --target` (default: 'Renin').",
    )
    parser.add_argument(
        "--mode", choices=["database", "analytic"], default="analytic",
        help="'database' (run full bulk pipeline including test-to-train sample "
        "similarity, requires training BEDs to be in --bed-dir) or 'analytic' "
        "(only sample-to-label distance + score, no rr step). Default: analytic.",
    )
    parser.add_argument(
        "--configs", nargs="+", default=CONFIGS,
        help=f"Config names to run. Default: {CONFIGS}",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Wipe each config's output dir and re-run, even if reninness_score.csv exists.",
    )
    args = parser.parse_args()

    test_meta = Path(args.test_meta).resolve()
    train_meta = Path(args.train_meta).resolve()
    for p, name in [(test_meta, "--test-meta"), (train_meta, "--train-meta")]:
        if not p.exists():
            print(f"ERROR: {name} not found: {p}", file=sys.stderr)
            return 1

    models_root = Path(args.models_root).resolve()
    if not models_root.exists():
        print(f"ERROR: --models-root not found: {models_root}", file=sys.stderr)
        return 1

    output_root = Path(args.output_dir).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    done: list[str] = []
    ran: list[str] = []
    failed: list[tuple[str, str]] = []

    for i, cfg in enumerate(args.configs, 1):
        cfg_out = output_root / cfg / args.result_tag
        score_csv = cfg_out / "reninness_score.csv"

        if score_csv.exists() and not args.force:
            print(f"\nCONFIG {i}/{len(args.configs)}: {cfg} — already done, skipping "
                  f"(use --force to redo)")
            done.append(cfg)
            continue

        if args.force and cfg_out.exists():
            shutil.rmtree(cfg_out)

        print(f"\n{'=' * 70}\nCONFIG {i}/{len(args.configs)}: {cfg}\n{'=' * 70}")
        try:
            run_config(
                config=cfg,
                test_meta=test_meta,
                train_meta=train_meta,
                bed_dir=args.bed_dir,
                models_root=models_root,
                starspace=args.starspace,
                universe=args.universe,
                label_col=args.label_col,
                target=args.target,
                out_dir=cfg_out,
                analytic=(args.mode == "analytic"),
            )
        except Exception as e:
            print(f"  CONFIG FAILED ({cfg}): {e}", flush=True)
            failed.append((cfg, str(e)))
            continue
        ran.append(cfg)

    print(f"\n{'=' * 70}\nSC SCORING SUMMARY\n{'=' * 70}")
    print(f"  Total configs: {len(args.configs)}")
    print(f"  Resumed (skipped): {len(done)}")
    print(f"  Newly run:         {len(ran)}")
    print(f"  Failed:            {len(failed)}")
    if failed:
        for cfg, err in failed:
            print(f"    - {cfg}: {err}")
        print("  → Re-run this script to retry only the failed configs.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

