#!/usr/bin/env python
"""
7e  Score EXTERNAL sample sets -- pseudobulk AND single cell -- against the
    models 7a already trained.

Both sets run through the same code path; `--subdir` and `--metadata-test`
are the only things that differ:

    --subdir pseudobulk  --label-col Group  --metadata-test renin_pseudobulk_data.csv
    --subdir singlecell  --label-col label  --metadata-test renin_singlecell_data_with_label.csv

Nothing is preprocessed and nothing is trained. For each config this runs only
the tail of the 7a pipeline -- `distances` then `score` for every method --
against the frozen model in that config's sweep directory, and writes into a
SUBDIRECTORY of that same config directory:

    sweep/<config>/model_<config>[.tsv]      <- input, untouched
    sweep/<config>/weights_learned_post.json <- input, untouched
    sweep/<config>/<subdir>/                 <- everything this script writes
        raw_cosdist_rl.csv
        reninness_score_<method>.csv
        <config>_starspace_embed.txt

The subdirectory is what keeps the bulk sweep's own score tables intact: they
live at sweep/<config>/reninness_score_*.csv and are never opened here.

The per-part summaries it writes under --summary-dir are the input to
`7_model_config_evaluation.py external`, which aggregates them into
sweep/external_perf/<set>_scored.csv -- the tables figs 2-5 read.

WHY A SEPARATE DRIVER AND NOT 7a --rescore
------------------------------------------
7a's --rescore re-runs distances+score in place, over the metadata the sweep was
built from, and writes to sweep/<config>/ itself. Pointing it at new samples
would overwrite the bulk tables. This driver calls the SAME shared helper 7a,
7b and 7c all call (`_sweep.run_distances_and_score`), so the four drivers
cannot drift on distance flags or on how the scoring methods are invoked -- only
the output directory and the test metadata differ.

METADATA-TRAIN
--------------
`bedspace distances --mode bulk` insists on a readable --metadata-train and
builds a `file_list_train` from it, but that list is used ONLY when --analytic
is off (for the sample-to-sample distance table). The shared helper always
passes --analytic, so the train list is parsed and discarded. The label
EMBEDDINGS come from the model's own .tsv (`get_label_embedding`), not from any
metadata. This script therefore passes the test metadata as both, which keeps
--metadata-train's label columns and file paths trivially valid; passing the
original bulk sheet would fail, because its label columns (Label_1, Label_2) do
not exist in the external sheets and its BED files are in another directory.

WEIGHTS
-------
`weighted_multi_learned_post` needs two weights over the negative labels. By
default this reuses sweep/<config>/weights_learned_post.json -- the weights that
config learned on the 252 bulk training samples -- so the external samples are
scored by a rule fit elsewhere. That is the out-of-sample application and it is
the default for a reason. `--learn-weights-here` instead fits the weights on the
external table itself; that is an IN-SAMPLE number and every row it produces is
flagged `weights_in_sample=True` in the summary CSV.

USAGE
-----
    python3 7e_run_reninness_score_sc_in_batch.py \\
        --sweep-dir    <RUN>/sweep \\
        --subdir       pseudobulk \\
        --metadata-test /home/bx2ur/code/renin_atac/metadata/renin_pseudobulk_data.csv \\
        --data-path    .../pseudobulk_data/cell_type_cluster/ \\
        --universe     .../universe_cc_073126.bed \\
        --starspace    /home/bx2ur/code/reninness_score/tools/Starspace/ \\
        --label-col    Group \\
        --all-configs

Resume: a config is done when its <subdir>/reninness_score_weighted_multi_uniform.csv
is readable. --force redoes every config.

SUPERSEDES the original 7e, which scored single cell only: it hard-coded eight
configs from the retired design, built its own `distances` command instead of
calling the shared helper, and wrote one scoring method rather than all of them.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import pandas as pd

_THIS_DIR = Path(__file__).resolve().parent


def _load(path: Path, name: str):
    """Import a sibling module whose filename starts with a digit."""
    if not path.exists():
        print(f"ERROR: expected sibling file not found: {path}", file=sys.stderr)
        sys.exit(1)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# 7a owns the distances+score driver, the metric definitions and the figure
# helpers. Import it rather than copy any of it.
_sweep = _load(_THIS_DIR / "7a_run_reninness_score_insample_in_batch.py", "_sweep_module")

SWEEP_CONFIGS = _sweep.SWEEP_CONFIGS
SCORING_METHODS = _sweep.SCORING_METHODS


def all_designs() -> list[str]:
    """The 266 configs, in the same order and by the same helpers 7a/7c use."""
    return list(dict.fromkeys(
        _sweep.ofat10_configs()
        + list(_sweep.JOINT_CONFIGS)
        + _sweep.ofat10b_configs()
        + _sweep.ofat10c_configs()
        + _sweep.ofat10d_configs()
    ))


def model_prefix_for(cfg_dir: Path, cfg: str) -> str | None:
    """The `-i` argument for `bedspace distances`: a path with no .tsv suffix.

    Both `model_<cfg>` (the StarSpace BINARY, read by embed_doc) and
    `model_<cfg>.tsv` (the embedding table, read for the label vectors) must be
    present. A directory holding only the .tsv looks complete and is not: that
    is exactly what stranded the 7c holdout arm.
    """
    binary = cfg_dir / f"model_{cfg}"
    table = cfg_dir / f"model_{cfg}.tsv"
    if binary.is_file() and table.is_file():
        return str(binary)
    # Fall back to whatever .tsv is there, but only if its binary exists too.
    for table in sorted(cfg_dir.glob("*.tsv")):
        if table.with_suffix("").is_file():
            return str(table.with_suffix(""))
    return None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sweep-dir", required=True,
                   help="Directory holding one subdirectory per trained config")
    p.add_argument("--subdir", required=True,
                   help="Name of the per-config output subdirectory, e.g. 'pseudobulk'")
    p.add_argument("--metadata-test", required=True,
                   help="Metadata CSV for the external samples")
    p.add_argument("--data-path", required=True, help="Directory holding their BED files")
    p.add_argument("--universe", required=True, help="Universe BED file")
    p.add_argument("--starspace", required=True, help="StarSpace install directory")
    p.add_argument("--label-col", required=True,
                   help="Column in the test metadata carrying the truth label, "
                        "e.g. 'Group'. Comma-separated for multi-label.")
    p.add_argument("--configs", nargs="*", default=None,
                   help="Configs to score. Default: --all-configs")
    p.add_argument("--all-configs", action="store_true",
                   help="Score every config in the 266-config design")
    p.add_argument("--target", default="Renin")
    p.add_argument("--negatives", default="Control,Tumoral",
                   help="MODEL labels used as negatives (not test-set labels)")
    p.add_argument("--line-labels", default="Control,Tumoral")
    p.add_argument("--learn-weights-here", action="store_true",
                   help="Fit weighted_multi weights on the EXTERNAL table (in-sample) "
                        "instead of reusing the config's weights_learned_post.json")
    p.add_argument("--weight-constraint", default="nonneg-scaled",
                   choices=["none", "simplex", "nonneg-scaled"])
    p.add_argument("--class-weight", default="balanced", choices=["none", "balanced"])
    p.add_argument("--summary-prefix", default=None,
                   help="Prefix for the summary CSVs (default: the --subdir name)")
    p.add_argument("--summary-dir", default=None,
                   help="Where the summary CSVs go (default: <sweep-dir>/<subdir>_perf)")
    p.add_argument("--skip-eval", action="store_true")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    sweep_dir = Path(args.sweep_dir)
    if not sweep_dir.is_dir():
        print(f"Error: --sweep-dir not a directory: {sweep_dir}", file=sys.stderr)
        return 1
    for path in (args.metadata_test, args.data_path, args.universe):
        if not Path(path).exists():
            print(f"Error: not found: {path}", file=sys.stderr)
            return 1

    configs = all_designs() if args.all_configs else (args.configs or [])
    if not configs:
        print("Error: pass --all-configs or --configs ...", file=sys.stderr)
        return 1
    unknown = [c for c in configs if c not in SWEEP_CONFIGS]
    if unknown:
        print(f"Error: config(s) not in SWEEP_CONFIGS: {unknown}", file=sys.stderr)
        return 1

    n_test = len(pd.read_csv(args.metadata_test))
    print(f"external set: {args.metadata_test}  ({n_test} samples)")
    print(f"label column: {args.label_col}")
    print(f"configs:      {len(configs)}")
    print(f"writing to:   {sweep_dir}/<config>/{args.subdir}/")
    print(f"weights:      {'learned HERE (in-sample)' if args.learn_weights_here else 'reused from each config weights_learned_post.json'}")

    scored_rows: list[dict] = []
    failed: list[str] = []
    resumed = ran = 0

    for i, cfg in enumerate(configs, 1):
        cfg_dir = sweep_dir / cfg
        out_dir = cfg_dir / args.subdir
        sentinel = out_dir / "reninness_score_weighted_multi_uniform.csv"

        if sentinel.exists() and not args.force:
            try:
                pd.read_csv(sentinel, nrows=1)
            except Exception as e:
                print(f"\nCONFIG {i}/{len(configs)}: {cfg} — existing score unreadable "
                      f"({e}); re-running")
            else:
                print(f"\nCONFIG {i}/{len(configs)}: {cfg} — already done, skipping")
                resumed += 1
                _collect(cfg, out_dir, args, scored_rows)
                continue

        print(f"\n{'='*70}\nCONFIG {i}/{len(configs)}: {cfg}\n{'='*70}", flush=True)
        prefix = model_prefix_for(cfg_dir, cfg)
        if prefix is None:
            print(f"  NO USABLE MODEL in {cfg_dir} (need both the binary and the "
                  f".tsv); skipping", flush=True)
            failed.append(cfg)
            continue

        weights = None
        if not args.learn_weights_here:
            w = cfg_dir / "weights_learned_post.json"
            # "SKIP" makes the shared helper drop weighted_multi_learned_post
            # rather than silently fall back to fitting in-sample.
            weights = str(w) if w.exists() else "SKIP"
            if weights == "SKIP":
                print(f"  WARN: no weights_learned_post.json for {cfg}; "
                      f"weighted_multi_learned_post will be skipped", flush=True)

        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            _sweep.run_distances_and_score(
                model_prefix=prefix, cfg_dir=out_dir, starspace=args.starspace,
                universe=args.universe, data_path=args.data_path,
                label_col=args.label_col,
                # See the module docstring: --analytic makes the train list
                # inert, so the test sheet is a valid stand-in for both.
                metadata_train=args.metadata_test,
                metadata_test=args.metadata_test,
                config_name=cfg, target=args.target, negatives=args.negatives,
                line_labels=args.line_labels or None,
                weights_from=weights,
                weight_constraint=args.weight_constraint,
                class_weight=args.class_weight,
            )
        except Exception as e:
            print(f"  CONFIG FAILED ({cfg}): {e}", flush=True)
            failed.append(cfg)
            continue
        ran += 1
        _collect(cfg, out_dir, args, scored_rows)

    print(f"\n{'='*70}\n7n SUMMARY ({args.subdir})\n{'='*70}")
    print(f"  Configs: {len(configs)}   resumed: {resumed}   ran: {ran}   "
          f"failed: {len(failed)}")
    if failed:
        print(f"  Failed: {failed}")

    if not args.skip_eval and scored_rows:
        prefix = args.summary_prefix or args.subdir
        out = Path(args.summary_dir) if args.summary_dir else sweep_dir / f"{args.subdir}_perf"
        out.mkdir(parents=True, exist_ok=True)
        # Clusterability grades the MODEL, which this script does not touch, so
        # no cluster rows are produced here -- read them from the sweep's own
        # all266_clusterability.csv, which was computed from these same models.
        s_path, _ = _sweep.write_run_summary(out, scored_rows, [], prefix=prefix)
        print(f"\nFigures (clusterability comes from the sweep's own table):")
        print(f"  python3 plot_param_figure.py --scored {s_path} "
              f"--clusterability {sweep_dir}/all266_perf/all266_clusterability.csv "
              f"--prefix {out}/{prefix}")

    return 0 if not failed else 1


def _collect(cfg: str, out_dir: Path, args, scored_rows: list) -> None:
    """Accuracy metrics per scoring method, via 7a's shared implementation.

    silhouette/geniml_eval are off: both grade the MODEL, and the model lives in
    the parent directory and is identical to the one the sweep already graded.
    Recomputing them here would burn hours to reproduce numbers already on disk.
    """
    if args.skip_eval:
        return
    scored, _ = _sweep.collect_config_metrics(
        cfg, out_dir, target=args.target, silhouette=False, geniml_eval=False,
    )
    for row in scored:
        row["dataset"] = args.subdir
        row["weights_in_sample"] = (
            row["method"] == "weighted_multi_learned_post" and args.learn_weights_here
        )
    scored_rows.extend(scored)


if __name__ == "__main__":
    raise SystemExit(main())
