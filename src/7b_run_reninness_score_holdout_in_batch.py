#!/usr/bin/env python
"""
7c — Held-out test-set driver for the bedspace pipeline.

Sibling of 7b (leave-one-out). Where 7b retrains N times holding out one sample
each, this splits the metadata ONCE into train / held-out test, then for each
param combo from ``SWEEP_CONFIGS`` (loaded from the co-located
``7a_run_reninness_score_insample_in_batch.py``) runs preprocess → train → distances →
score on that single split and evaluates on the held-out samples only.

------------------------------------------------------------------
STRATIFIED SPLIT
------------------------------------------------------------------
The test set is sampled WITHIN each unique label group, not from the metadata as
a whole. With ``--label-col Label_1,Label_2`` the groups are the distinct
(Label_1, Label_2) combinations, so a 30% split takes 30% of

    (Renin,   <blank>)     (Renin,   Tumoral)
    (Control, <blank>)     (<blank>, Tumoral)

independently. Sampling the flat list instead would let the small groups —
(Renin, <blank>) is 11 of 246 samples — land entirely in train or entirely in
test by chance, which either removes the class from the test set or removes it
from training. Groups too small to give up a sample at the requested fraction
keep at least one in each side when ``--min-per-group`` allows it, and are
reported.

------------------------------------------------------------------
USAGE
------------------------------------------------------------------

    python 7c_reninness_score_holdout_in_batch.py \\
        --configs    d50_e10_m0.2_neg100 d3_e10_m0.2_neg100 \\
        --metadata   /path/to/renin_annotation_sheet_add_tumoral_l3.csv \\
        --data-path  /path/to/bed/files/ \\
        --universe   /path/to/universe.bed \\
        --starspace  /path/to/Starspace/ \\
        --label-col  Label_1,Label_2 \\
        --target     Renin \\
        --negatives  Control,Tumoral \\
        --test-frac  0.3 \\
        --output-dir /project/shefflab/.../holdout/

------------------------------------------------------------------
LEARNED WEIGHTS
------------------------------------------------------------------
`bedspace score --learn-weights` fits the two weighted_multi weights on the same
table it scores. For a held-out evaluation that would leak, so the default here
(``--learn-weights-from train``) is the two-stage form:

    distances on TRAIN  -> score --learn-weights --weights-out W.json
    distances on TEST   -> score --weights-file W.json

which costs a second `distances` pass per config. ``--learn-weights-from test``
reproduces the in-sample number instead, for comparison; it is not a holdout
result and is labelled as such in the output.

------------------------------------------------------------------
OUTPUT
------------------------------------------------------------------
    <output-dir>/
        holdout_train_meta.csv       the split, written once and reused by
        holdout_test_meta.csv        every config (and by any re-run)
        split_summary.csv            per label group: n_total, n_train, n_test
        <config>/
            raw_cosdist_rl.csv       test-sample distances
            raw_cosdist_rl_train.csv train-sample distances (learned weights)
            weights_learned_post.json
            reninness_score_<method>.csv     one per scoring method
        holdout_scored.csv           tidy config x method accuracy metrics --
                                     BOTH families (bal_acc, prec_rank_*)
        holdout_clusterability.csv   per config, BOTH families (silhouette,
                                     geniml.eval ctt/gdst/npt)

The two summary CSVs use the same schema as the sweep's
``<prefix>_scored.csv`` / ``<prefix>_clusterability.csv``, so the figures come
straight from the existing scripts:

    python3 plot_param_figure.py  --scored holdout_scored.csv \\
        --clusterability holdout_clusterability.csv --prefix holdout
    python3 plot_param_scatter.py --scored holdout_scored.csv \\
        --clusterability holdout_clusterability.csv --prefix holdout_scatter

NOTE: the LINE figures lay their panels out from the ofat10 level table, so they
only fill in for configs that belong to that design; other configs simply leave
gaps. The SCATTER figures take whatever configs are present, so they are the
ones to use for an arbitrary --configs list.

Resume behavior:
    A config is "done" if its reninness_score_pairwise.csv exists and is
    readable; re-running skips those. --force wipes and redoes every config.
    The split itself is written once and REUSED on re-runs unless --force-split
    is given, so resumed and fresh configs are always evaluated on the same
    held-out samples.
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

# The balancing itself happens inside `geniml bedspace train` (see run_config).
# These are imported rather than restated so this parser can never offer a
# strategy geniml does not implement, or refuse one it does. A hand-copied
# choices list is exactly what failed job 19006294: geniml had gained 'labels'
# and this file still listed only ('to-majority', 'factors'), so the run died at
# argument parsing 40 s in.
try:
    from geniml.bedspace.balance import (
        DEFAULT_MAX_ROUNDS as _BAL_MAX_ROUNDS,
        DEFAULT_MIN_RATIO as _BAL_MIN_RATIO,
        DEFAULT_STRATEGY as _BAL_STRATEGY,
        STRATEGIES as _BAL_STRATEGIES,
    )
except Exception:  # geniml too old to balance; --balance will fail loudly there
    _BAL_STRATEGIES = ("labels", "to-majority", "factors")
    _BAL_STRATEGY, _BAL_MIN_RATIO, _BAL_MAX_ROUNDS = "labels", 2.0, 1

# Load SWEEP_CONFIGS (and the shared evaluation metrics) from the co-located 7a
# file. We can't `import` it directly because the module name starts with a
# digit, so use importlib -- same approach as 7b.
_THIS_DIR = Path(__file__).resolve().parent
_SWEEP_FILE = _THIS_DIR / "7a_run_reninness_score_insample_in_batch.py"

if not _SWEEP_FILE.exists():
    print(f"ERROR: Expected sibling config file not found: {_SWEEP_FILE}", file=sys.stderr)
    print("This script reads SWEEP_CONFIGS from 7a_run_reninness_score_insample_in_batch.py.",
          file=sys.stderr)
    sys.exit(1)

_spec = importlib.util.spec_from_file_location("_sweep_configs_module", _SWEEP_FILE)
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
SWEEP_CONFIGS = _mod.SWEEP_CONFIGS
_sweep = _mod

TRAIN_META = "holdout_train_meta.csv"
TEST_META = "holdout_test_meta.csv"

# Scoring methods and the metric/figure helpers all come from 7a, so 7a, 7b and
# 7c cannot drift on what they compute or how the figures are drawn.
METHODS = _sweep.SCORING_METHODS


def run(cmd: list[str]) -> None:
    """Run a subprocess, stream output, raise on failure."""
    print(f"  $ {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {' '.join(cmd)}")


# ---------------------------------------------------------------- split
def stratified_holdout_split(
    meta_df: pd.DataFrame,
    label_cols: list[str],
    test_frac: float,
    seed: int,
    min_per_group: int = 1,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Split by unique label-column COMBINATION, not across the flat list.

    Each distinct tuple of `label_cols` values (blanks included, so (Renin, '')
    and (Renin, 'Tumoral') are different groups) gives up `test_frac` of its own
    rows. Rounding is by `round`, then clamped so that a group with at least
    2*min_per_group rows keeps at least `min_per_group` on each side -- without
    that, the 11-sample (Renin, '') group can round to zero test samples.

    Returns (train_df, test_df, summary_df).
    """
    rng = np.random.default_rng(seed)
    key = (
        meta_df[label_cols]
        .fillna("")
        .astype(str)
        .apply(lambda r: " | ".join(v.strip() for v in r), axis=1)
    )
    test_idx: list = []
    summary = []
    for group, idx in key.groupby(key).groups.items():
        idx = np.array(list(idx))
        n = len(idx)
        n_test = int(round(n * test_frac))
        if n >= 2 * min_per_group:
            n_test = min(max(n_test, min_per_group), n - min_per_group)
        else:
            # Too small to split without emptying a side; keep it all in train
            # so the class still trains, and say so.
            n_test = 0
        picked = rng.choice(idx, size=n_test, replace=False) if n_test else np.array([], dtype=int)
        test_idx.extend(picked.tolist())
        summary.append(dict(label_group=group, n_total=n, n_train=n - n_test,
                            n_test=n_test,
                            frac_test=(n_test / n) if n else np.nan))
    test_mask = meta_df.index.isin(test_idx)
    return (meta_df[~test_mask].copy(), meta_df[test_mask].copy(),
            pd.DataFrame(summary).sort_values("label_group").reset_index(drop=True))


def write_split(output_dir: Path, meta_df: pd.DataFrame, label_cols: list[str],
                test_frac: float, seed: int, min_per_group: int,
                force: bool) -> tuple[Path, Path]:
    """Create the split once and reuse it, so every config sees the same test set."""
    train_path, test_path = output_dir / TRAIN_META, output_dir / TEST_META
    if train_path.exists() and test_path.exists() and not force:
        tr, te = pd.read_csv(train_path), pd.read_csv(test_path)
        print(f"Reusing existing split: {len(tr)} train / {len(te)} test "
              f"({train_path.name}). Pass --force-split to resample.")
        return train_path, test_path

    train_df, test_df, summary = stratified_holdout_split(
        meta_df.reset_index(drop=True), label_cols, test_frac, seed, min_per_group
    )
    train_df.to_csv(train_path, index=False)
    test_df.to_csv(test_path, index=False)
    summary.to_csv(output_dir / "split_summary.csv", index=False)
    print(f"\nStratified split by {label_cols} (test_frac={test_frac}, seed={seed}):")
    print(summary.to_string(index=False))
    print(f"  TOTAL  train={len(train_df)}  test={len(test_df)}")
    dropped = summary[summary.n_test == 0]
    if len(dropped):
        print(f"  NOTE: {len(dropped)} group(s) too small to split "
              f"(< {2 * min_per_group} rows) kept entirely in train: "
              f"{dropped.label_group.tolist()}")
    return train_path, test_path


# ---------------------------------------------------------------- one config
def run_config(config_name: str, params: dict, train_meta: Path, test_meta: Path,
               data_path: str, universe: str, starspace: str, label_col: str,
               cfg_dir: Path, target: str, negatives: str | None,
               line_labels: str | None, learn_weights_from: str,
               weight_constraint: str, class_weight: str,
               truth_column: str | None, truth_positive: str | None,
               name_col: str, balance: bool = False,
               balance_strategy: str = "labels",
               balance_factors: str | None = None,
               balance_min_ratio: float = 2.0,
               balance_max_rounds: int = 1,
               balance_seed: int = 42,
               fit_frac: float = 0.1) -> dict[str, Path]:
    """preprocess -> train -> distances -> score (4 methods) for one config.

    The distances+score stage is 7a's `run_distances_and_score`, so the holdout
    driver emits exactly the files the sweep driver does.
    """
    cfg_dir.mkdir(parents=True, exist_ok=True)
    preprocess_dir = cfg_dir / "preprocessed"

    # ---- preprocess (train only; test BEDs are tokenized inside `distances`) ----
    run(["geniml", "bedspace", "preprocess",
         "-i", data_path, "-m", str(train_meta), "-u", universe,
         "-o", str(preprocess_dir), "-l", label_col, "--mode", "bulk"])

    train_input = preprocess_dir / "train_input.txt"
    if not train_input.exists():
        raise FileNotFoundError(f"preprocess did not produce {train_input}")

    # ---- train ----
    train_cmd = [
        "geniml", "bedspace", "train",
        "-s", starspace, "-i", str(train_input), "-o", str(cfg_dir),
        "-n", str(params["epoch"]), "-d", str(params["dim"]), "-l", str(params["lr"]),
        "--neg-search-limit", str(params["negSearchLimit"]),
        "--margin", str(params["margin"]), "--min-count", str(params["minCount"]),
    ]
    train_cmd.append("--adagrad" if params.get("adagrad") else "--no-adagrad")
    # CLASS BALANCE. `train_input` above was preprocessed from
    # holdout_train_meta.csv ALONE, so it holds only the 70% training split --
    # the held-out BEDs are tokenized later, inside `distances`. Oversampling
    # inside `bedspace train` therefore touches the training split and nothing
    # else, which is the property that makes this safe to do per config.
    # geniml's train writes train_input_balanced.txt + balance_receipt.json into
    # cfg_dir and trains from the balanced file; if the split turns out to be
    # within --balance-min-ratio it trains from the original and says so.
    if balance:
        train_cmd.append("--balance")
        train_cmd += ["--balance-strategy", balance_strategy,
                      "--balance-min-ratio", str(balance_min_ratio),
                      "--balance-max-rounds", str(balance_max_rounds),
                      "--balance-seed", str(balance_seed)]
        if balance_factors:
            train_cmd += ["--balance-factors", balance_factors]
    run(train_cmd)

    candidates = sorted(cfg_dir.glob("*.tsv"))
    if not candidates:
        raise FileNotFoundError(f"No model .tsv found in {cfg_dir}")
    model_prefix = str(candidates[0])[:-4]

    # ---- region-embedding t-SNE / UMAP ------------------------------------
    # The holdout arm never drew these: only 7a called evaluate_model, so a
    # holdout config produced scores and no picture of the embedding they came
    # from. That made the two arms impossible to eyeball side by side.
    #
    # `fit_frac` (default 0.1) fits the projections on a FRACTION of the
    # vocabulary. It is not an optimisation to skip lightly: the full 422k-region
    # fit did not finish one config's t-SNE in 10 minutes, which is longer than a
    # whole config otherwise takes. At 10% (~42k regions) it costs ~2 min.
    #
    # NOTHING IN THE ANALYSIS READS THESE. sil_k2 comes from multi_k_clustering
    # on the FULL embeddings, ctt/gdst from geniml.eval on the FULL embeddings,
    # nl_auprc from raw_cosdist_rl.csv. The fraction changes the pictures only.
    if fit_frac:
        try:
            _sweep.evaluate_model(
                str(candidates[0]), output_name=config_name,
                output_dir=str(cfg_dir), deseq_path=None,
                fit_max_regions=fit_frac,
            )
        except Exception as e:
            print(f"    region-embedding plots failed: {e}", flush=True)

    # ---- weights fit on TRAIN (the honest holdout form) ----
    # A second distances pass over the training samples, scored with
    # --learn-weights, produces weights.json; the held-out scoring below then
    # consumes it with --weights-file so no held-out sample enters the fit.
    # `--learn-weights-from test` skips this and fits on the held-out table
    # itself, which is IN-SAMPLE and labelled as such in the outputs.
    weights_json = None
    if learn_weights_from == "train":
        train_dir = cfg_dir / "_train_distances"
        try:
            _sweep.run_distances_and_score(
                model_prefix=model_prefix, cfg_dir=train_dir, starspace=starspace,
                universe=universe, data_path=data_path, label_col=label_col,
                metadata_train=str(train_meta), metadata_test=str(train_meta),
                config_name=f"{config_name}_train", target=target,
                negatives=negatives, line_labels=line_labels,
                weight_constraint=weight_constraint, class_weight=class_weight,
            )
        except Exception as e:
            # A degenerate model can fail the weight fit outright (NNLS returns
            # all-zero weights). That must not cost this config its other three
            # scorings on the held-out samples, so carry on without a
            # weights.json -- the learned method is then skipped below.
            print(f"  WARN: train-side weight fit failed for {config_name}: {e}")
        weights_json = train_dir / "weights_learned_post.json"
        if weights_json.exists():
            shutil.copy2(weights_json, cfg_dir / "weights_learned_post_trainfit.json")
        else:
            # Deliberately NOT falling back to learning on the held-out table:
            # that would be an in-sample number reported in a holdout run, and
            # `weights_in_sample` would label it False because the run was
            # launched with --learn-weights-from train. Skipping the method
            # leaves an honest gap in the metric CSV instead.
            print(f"  WARN: no weights.json from the train pass for "
                  f"{config_name}; SKIPPING weighted_multi_learned_post for "
                  f"this config (the other methods still run)")
            weights_json = "SKIP"
        train_dist = train_dir / "raw_cosdist_rl.csv"
        if train_dist.exists():
            shutil.copy2(train_dist, cfg_dir / "raw_cosdist_rl_train.csv")

    # ---- distances + score on the HELD-OUT samples ----
    return _sweep.run_distances_and_score(
        model_prefix=model_prefix, cfg_dir=cfg_dir, starspace=starspace,
        universe=universe, data_path=data_path, label_col=label_col,
        metadata_train=str(train_meta), metadata_test=str(test_meta),
        config_name=config_name, target=target, negatives=negatives,
        line_labels=line_labels, weights_from=weights_json,
        weight_constraint=weight_constraint, class_weight=class_weight,
    )


def prune_intermediates(cfg_dir: Path) -> None:
    """Drop a config's working files, keeping what later steps actually read.

    Two things are KEPT that were not, before 2026-08-25, and the reason is the
    2026-08-22 `line` fix: neither could be repaired here without retraining all
    120 configs, because this function had deleted the inputs.

      starspace_trained_model  -- the StarSpace BINARY. `bedspace distances`
          shells out to `embed_doc <model> <docs>`, which reads the binary, NOT
          the .tsv. The old rule kept `*.tsv` only, so every holdout config had
          a model that looked present and could not actually be re-run:
          `embed_doc` exits 1. ~110 MB/config. This is the file that makes a
          re-scoring of this arm possible at all.
      raw_cosdist_rl.csv       -- the distance table, sole input to `bedspace
          score`. ~130 KB/config, and it makes a pure scoring-method fix a
          `score` re-run with no model needed.

    The sweep arm kept both and was repaired in one pass (job 18808473); this
    arm kept neither and needs a full retrain. Do not "tidy" either back out.
    """
    # balance_receipt.json / training_params.json are a few hundred bytes each
    # and are the ONLY record of whether a model was trained on a resampled
    # corpus and with which duplication factors. A model whose provenance is
    # gone cannot be compared to anything.
    keep_exact = {"raw_cosdist_rl.csv", "balance_receipt.json",
                  "training_params.json"}
    # the region-embedding pictures, and the stats behind them
    keep_suffix = ("_tsne.png", "_umap.png", "_report.txt", "_stats.json",
                   "_silhouette_vs_k.png", "_elbow.png")
    for child in cfg_dir.iterdir():
        if child.is_dir():
            shutil.rmtree(child)
        elif not (child.name in keep_exact
                  or child.name.endswith(keep_suffix)
                  or child.name.startswith("starspace_trained_model")
                  or child.name.startswith("reninness_score_")
                  or child.name.startswith("weights_")
                  or child.name.startswith("subscores")
                  or child.name.endswith(".tsv")):
            child.unlink()


# ---------------------------------------------------------------- evaluation
def evaluate(output_dir: Path, scored_rows: list[dict], cluster_rows: list[dict],
             arm_figures: bool = True, prefix: str = "holdout") -> None:
    """Write the two tidy CSVs the figure scripts consume, then draw the figures.

    `prefix` MUST differ between concurrent jobs sharing one --output-dir. The
    144-config grid is split into four jobs by dim; with the prefix hard-coded
    all four wrote <output_dir>/holdout_scored.csv and the last to finish left a
    36-row table where 144 were expected -- silently, since the file exists and
    parses. The four are concatenated afterwards.
    """
    s_path, c_path = _sweep.write_run_summary(
        output_dir, scored_rows, cluster_rows, prefix=prefix
    )
    s, c = pd.DataFrame(scored_rows), pd.DataFrame(cluster_rows)
    if len(s):
        show = [x for x in ("config", "method", "n_pos", "n_neg", "bal_acc",
                            "prec_rank_overall", "prec_rank_renin", "auc")
                if x in s.columns]
        print("\n" + s[show].round(4).to_string(index=False))
    if len(c):
        show = [x for x in ("config", "sil_k2", "best_k", "ctt", "gdst",
                            "npt_snpr_mean", "rct") if x in c.columns]
        print("\n" + c[show].round(4).to_string(index=False))
    if not arm_figures:
        print("\n--no-arm-figures: summary CSVs written, per-arm parameter "
              "figures skipped")
    elif s.config.nunique() >= 2:
        print(f"\n{'='*70}\nFIGURES\n{'='*70}")
        _sweep.make_figures(s_path, c_path, output_dir / "holdout")
    else:
        print("\nOnly one config — the parameter figures need >= 2; skipping.")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--configs", nargs="+", required=True,
                        help="Config name(s) from SWEEP_CONFIGS")
    parser.add_argument("--metadata", required=True, help="Full metadata CSV")
    parser.add_argument("--data-path", required=True,
                        help="Directory holding BED files referenced in metadata")
    parser.add_argument("--universe", required=True, help="Universe BED file")
    parser.add_argument("--starspace", required=True, help="StarSpace install directory")
    parser.add_argument("--label-col", default="Label_1,Label_2",
                        help="Label column name(s), comma-separated. The unique "
                             "COMBINATION of these defines the strata for the "
                             "holdout split. Default: 'Label_1,Label_2'.")
    parser.add_argument("--target", default="Renin", help="Positive class (default: Renin)")
    parser.add_argument("--negatives", default="Control,Tumoral",
                        help="Comma-separated negatives for weighted_multi "
                             "(default: 'Control,Tumoral'). The first is also "
                             "the single negative used by --method pairwise.")
    parser.add_argument("--line-labels", default="Control,Tumoral",
                        help="Two negative labels forming the line for "
                             "--method line. Empty string disables it.")
    parser.add_argument("--output-dir", required=True, help="Holdout output directory")
    parser.add_argument("--name-col", default="file_name",
                        help="Column uniquely identifying a sample (default: file_name)")

    parser.add_argument("--test-frac", type=float, default=0.3,
                        help="Fraction held out WITHIN each label group (default: 0.3)")
    parser.add_argument("--seed", type=int, default=42, help="Split seed (default: 42)")
    parser.add_argument("--min-per-group", type=int, default=1,
                        help="Minimum rows kept on each side of the split for a "
                             "group large enough to give one up (default: 1)")
    parser.add_argument("--force-split", action="store_true",
                        help="Resample the split even if one already exists. "
                             "Invalidates comparison with previously-run configs.")

    parser.add_argument("--learn-weights-from", choices=["train", "test"], default="train",
                        help="Where the two weighted_multi weights are fit. "
                             "'train' (default): honest holdout -- a second "
                             "distances pass on the training samples, weights "
                             "applied to the held-out ones. 'test': fit on the "
                             "held-out table itself (IN-SAMPLE, for comparison).")
    parser.add_argument("--weight-constraint", default="nonneg-scaled",
                        choices=["nonneg-scaled", "simplex"])
    parser.add_argument("--class-weight", default="balanced", choices=["none", "balanced"])
    parser.add_argument("--truth-column", default=None,
                        help="Column holding the TRUE binary label when it "
                             "differs from the training labels. Default: derive "
                             "from the label column(s) in the distance table.")
    parser.add_argument("--truth-positive", default=None,
                        help="Value in --truth-column marking the positive class "
                             "(default: --target)")

    parser.add_argument("--skip-eval", action="store_true",
                        help="Skip the evaluation metrics / summary CSVs")
    parser.add_argument("--skip-geniml-eval", action="store_true",
                        help="Skip the geniml.eval tests (CTT/GDST/NPT)")
    parser.add_argument("--skip-silhouette", action="store_true",
                        help="Skip the multi-k silhouette (it loads the full model)")
    parser.add_argument("--geniml-eval-workers", type=int, default=10)
    parser.add_argument("--bin-embed", default=None,
                        help="Binary-embedding pickle enabling the RCT test")

    # ---- class-imbalance oversampling of the TRAINING SPLIT ---------------
    # Delegated to `geniml bedspace train`, which is the only stage that reads
    # the training file. `preprocess` here runs against holdout_train_meta.csv
    # alone, so the file it balances contains only the 70% split; the held-out
    # samples are tokenized inside `distances` and never reach it.
    parser.add_argument("--balance", action="store_true",
                        help="Oversample the minority label groups of each "
                             "config's TRAINING split before StarSpace sees it. "
                             "Off by default: it resamples the corpus, so a "
                             "balanced arm is not comparable to an unbalanced "
                             "one.")
    parser.add_argument("--balance-strategy", default=_BAL_STRATEGY,
                        choices=list(_BAL_STRATEGIES),
                        help="'labels' (default) equalises LABEL OCCURRENCES: "
                             "it duplicates every document carrying the "
                             "scarcest label until that label reaches the "
                             "commonest one's count, counting a multi-label "
                             "document once per label. 'to-majority' works on "
                             "label-SET groups instead. 'factors' applies "
                             "--balance-factors verbatim -- use it to match a "
                             "sweep arm that used explicit factors. The choices "
                             "come from geniml.bedspace.balance.STRATEGIES; "
                             "this parser only forwards them.")
    parser.add_argument("--balance-max-rounds", type=int, default=_BAL_MAX_ROUNDS,
                        help="'labels' only: passes of the rule (default 1). "
                             "Iterating oscillates -- duplicating a compound "
                             "document raises two labels at once.")
    parser.add_argument("--balance-factors", default=None,
                        help="For --balance-strategy factors: 'GROUP=N,...', "
                             "e.g. 'Renin=15,Renin|Tumoral=15'. The group key "
                             "is the '|'-joined label SET (a two-label document "
                             "is its own group), order-insensitive.")
    parser.add_argument("--balance-min-ratio", type=float, default=_BAL_MIN_RATIO,
                        help="Detection threshold on largest/smallest group; "
                             "below it the split is left untouched (default 2.0)")
    parser.add_argument("--balance-seed", type=int, default=42)
    parser.add_argument("--summary-prefix", default="holdout",
                        help="Basename for this job's two summary CSVs "
                             "(default 'holdout'). Concurrent jobs sharing one "
                             "--output-dir MUST each pass a distinct value or "
                             "they overwrite each other's table.")
    parser.add_argument("--no-arm-figures", dest="arm_figures",
                        action="store_false", default=True,
                        help="Write the run-level summary CSVs but SKIP the old per-arm parameter figures. Those are superseded by the five deliverable figures (7z_upsampled_final_figures.py); drawing them here just refills the arm directory with plots nobody reads.")
    parser.add_argument("--fit-frac", type=float, default=0.1,
                        help="Fraction of each config's vocabulary the region "
                             "t-SNE / UMAP are FIT on (default 0.1, ~42k of "
                             "422k). 0 disables the plots. They are pictures "
                             "only -- silhouette, CTT/GDST and nl_auprc are all "
                             "computed on the FULL embeddings.")
    parser.add_argument("--keep-intermediates", action="store_true",
                        help="Keep per-config preprocess/model/distance files")
    parser.add_argument("--force", action="store_true",
                        help="Re-run every config from scratch. Default: resume — "
                             "configs with a readable reninness_score_pairwise.csv "
                             "are skipped.")
    parser.add_argument("--rescore", action="store_true",
                        help="Reuse the model .tsv already in each config dir and "
                             "re-run distances + score against it, instead of "
                             "skipping the config. Configs with no model on disk "
                             "still train normally, so one --rescore run both "
                             "repairs existing configs and fills in new ones. "
                             "Mirrors 7a's --rescore.")
    args = parser.parse_args()

    unknown = [c for c in args.configs if c not in SWEEP_CONFIGS]
    if unknown:
        print(f"Error: config(s) {unknown} not in SWEEP_CONFIGS.", file=sys.stderr)
        print(f"Available: {sorted(SWEEP_CONFIGS.keys())}", file=sys.stderr)
        return 1

    meta_df = pd.read_csv(args.metadata)
    label_cols = [c.strip() for c in str(args.label_col).split(",") if c.strip()]
    missing = [c for c in label_cols + [args.name_col] if c not in meta_df.columns]
    if missing:
        print(f"Error: column(s) {missing} not in metadata {list(meta_df.columns)}",
              file=sys.stderr)
        return 1
    if args.truth_column and args.truth_column not in meta_df.columns:
        print(f"Error: --truth-column '{args.truth_column}' not in metadata "
              f"{list(meta_df.columns)}", file=sys.stderr)
        return 1
    if not 0 < args.test_frac < 1:
        print("Error: --test-frac must be strictly between 0 and 1.", file=sys.stderr)
        return 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_meta, test_meta = write_split(
        output_dir, meta_df, label_cols, args.test_frac, args.seed,
        args.min_per_group, args.force_split,
    )

    scored_rows: list[dict] = []
    cluster_rows: list[dict] = []
    failed: list[str] = []
    resumed = ran = 0

    for i, cfg in enumerate(args.configs, 1):
        cfg_dir = output_dir / cfg

        # --rescore: the model on disk is fine, the SCORING changed. Re-run
        # distances + score against it rather than skipping the config (plain
        # resume) or retraining it (--force, which would replace the models
        # every published holdout number was computed from).
        #
        # The train-side weight fit is deliberately NOT redone. Its output,
        # weights_learned_post_trainfit.json, was written by the original run
        # from a distance table the scoring method never touches, so reusing it
        # keeps weighted_multi_learned_post bit-identical and saves the second
        # distances pass over the training split -- half the cost of the config.
        # A config with no model on disk falls through to the normal path and
        # trains, so ONE --rescore run repairs the existing configs and fills in
        # the ones this arm never covered.
        existing_model = (next(iter(sorted(cfg_dir.glob("*.tsv"))), None)
                          if cfg_dir.is_dir() else None)
        if args.rescore and existing_model is not None and not args.force:
            print(f"\n{'='*70}\nCONFIG {i}/{len(args.configs)}: {cfg} \u2014 model on "
                  f"disk, re-running distances + score (no retraining)\n{'='*70}")
            wt = cfg_dir / "weights_learned_post_trainfit.json"
            try:
                _sweep.run_distances_and_score(
                    model_prefix=str(existing_model)[:-4],
                    cfg_dir=cfg_dir, starspace=args.starspace,
                    universe=args.universe, data_path=args.data_path,
                    label_col=args.label_col, metadata_train=str(train_meta),
                    metadata_test=str(test_meta), config_name=cfg,
                    target=args.target, negatives=args.negatives,
                    line_labels=args.line_labels or None,
                    weights_from=(str(wt) if wt.exists() else "SKIP"),
                    weight_constraint=args.weight_constraint,
                    class_weight=args.class_weight,
                )
            except Exception as e:
                print(f"  RESCORE FAILED ({cfg}): {e}", flush=True)
                failed.append(cfg)
                continue
            ran += 1
            _collect(cfg, cfg_dir, args, scored_rows, cluster_rows)
            if not args.keep_intermediates:
                prune_intermediates(cfg_dir)
            continue

        sentinel = cfg_dir / "reninness_score_pairwise.csv"
        if sentinel.exists() and not args.force:
            try:
                pd.read_csv(sentinel, nrows=1)
            except Exception as e:
                print(f"\nCONFIG {i}/{len(args.configs)}: {cfg} — existing score "
                      f"unreadable ({e}); re-running")
            else:
                print(f"\nCONFIG {i}/{len(args.configs)}: {cfg} — already done, "
                      f"skipping (use --force to redo)")
                resumed += 1
                # Still evaluated: the metrics are cheap next to retraining, and
                # a resumed run must appear in the summary CSVs like any other.
                _collect(cfg, cfg_dir, args, scored_rows, cluster_rows)
                continue

        print(f"\n{'='*70}\nCONFIG {i}/{len(args.configs)}: {cfg}\n{'='*70}")
        if args.force and cfg_dir.exists():
            shutil.rmtree(cfg_dir)
        try:
            run_config(
                config_name=cfg, params=SWEEP_CONFIGS[cfg],
                train_meta=train_meta, test_meta=test_meta,
                data_path=args.data_path, universe=args.universe,
                starspace=args.starspace, label_col=args.label_col,
                cfg_dir=cfg_dir, target=args.target, negatives=args.negatives,
                line_labels=args.line_labels or None,
                learn_weights_from=args.learn_weights_from,
                weight_constraint=args.weight_constraint,
                class_weight=args.class_weight,
                truth_column=args.truth_column, truth_positive=args.truth_positive,
                name_col=args.name_col,
                balance=args.balance, balance_strategy=args.balance_strategy,
                balance_factors=args.balance_factors,
                balance_min_ratio=args.balance_min_ratio,
                balance_max_rounds=args.balance_max_rounds,
                balance_seed=args.balance_seed,
                fit_frac=args.fit_frac,
            )
        except Exception as e:
            print(f"  CONFIG FAILED ({cfg}): {e}", flush=True)
            failed.append(cfg)
            continue
        ran += 1
        _collect(cfg, cfg_dir, args, scored_rows, cluster_rows)

        if not args.keep_intermediates:
            prune_intermediates(cfg_dir)

    print(f"\n{'='*70}\nHOLDOUT SUMMARY\n{'='*70}")
    print(f"  Configs:  {len(args.configs)}   resumed: {resumed}   ran: {ran}   "
          f"failed: {len(failed)}")
    if failed:
        print(f"  Failed: {failed}")
        print("   Re-run this script to retry only the failed configs.")

    if not args.skip_eval:
        evaluate(output_dir, scored_rows, cluster_rows, args.arm_figures,
                 args.summary_prefix)
        print("\nFigures:")
        print(f"  python3 plot_param_figure.py  --scored {output_dir/'holdout_scored.csv'} "
              f"--clusterability {output_dir/'holdout_clusterability.csv'} "
              f"--prefix {output_dir/'holdout'}")
        print(f"  python3 plot_param_scatter.py --scored {output_dir/'holdout_scored.csv'} "
              f"--clusterability {output_dir/'holdout_clusterability.csv'} "
              f"--prefix {output_dir/'holdout_scatter'}")
    return 0 if not failed else 1


def _collect(cfg: str, cfg_dir: Path, args, scored_rows: list, cluster_rows: list) -> None:
    """Accuracy metrics per scoring method + clusterability metrics per config.

    Delegates to the shared implementation in 7a so the holdout, LOOCV and sweep
    drivers all report the identical metric set.
    """
    if args.skip_eval:
        return
    scored, crow = _sweep.collect_config_metrics(
        cfg, cfg_dir, target=args.target,
        silhouette=not args.skip_silhouette,
        geniml_eval=not args.skip_geniml_eval,
        geniml_eval_workers=args.geniml_eval_workers,
        bin_embed=args.bin_embed,
    )
    for row in scored:
        # The in-sample variant is a comparison point, not a holdout result.
        row["weights_in_sample"] = (row["method"] == "weighted_multi_learned_post"
                                    and args.learn_weights_from == "test")
    scored_rows.extend(scored)
    if len(crow) > 1:
        cluster_rows.append(crow)


if __name__ == "__main__":
    raise SystemExit(main())
