#!/usr/bin/env python3
"""
7  Model-config evaluation: region metrics, the model axis, SUMMARY_TABLE, and
   the figures_final set.

One file for the whole model-selection pipeline. Absorbs what used to be
7q (region geometry), 7u + 7i (nearest label), 7z_upsampled_summary_table
(SUMMARY_TABLE + main effects) and 7z_upsampled_final_figures (fig1-fig5).

Stages, in order
  geom      per-config region geometry: sil_k2, dip_axis, pct_lt_025, eff_dim
  rand      dimension-matched random baseline for those
  rcs       reconstruction test; the long pole
  nlauc     model axis: AUC/AUPRC of the nearest-label margin
  nlhard    nearest-label hard assignment
  merge     joins them; derives dip_ratio and pct_excess
  summary   SUMMARY_TABLE + MAIN_EFFECTS + MAIN_EFFECT_SPREAD
  external  pseudobulk + single-cell metrics on Renin vs Non_Renin
  figures   figures_final/fig1-fig5 + score_distributions.csv

Why each stage
  - `rand` is not optional: pct_lt_025 is not comparable across dimensions
    (random data alone gives 7.97% at d=2 and 0.00% at d>=15), so the gate is
    on pct_excess = pct_lt_025 - pct_lt_025_rand.
  - `rcs` is the only test that detects an UNTRAINED model.
  - `nlauc` is the model axis every figure uses: weight-free, threshold-free.
    margin(s) = min_{L != Renin} d(s,L) - d(s,Renin); sign(margin) reproduces
    the hard nearest-label call exactly.

The gate
      ctt > 0.5              AND        pct_excess <= 0.1
      untrained detector                degeneracy / collapse detector
  - Neither gate can be dropped. pct_excess cannot see an untrained model
    (a random projection is maximally spread); ctt cannot see a collapsed one
    (it rates collapse highest).
  - CTT is a GATE and never a score: among gated configs it carries no
    continuous quality signal.
  - pct_excess <= 0.1 is eff_n_distinct >= 1000, an absolute capacity
    requirement, not a property of the config set.
  - notes/REGION_METRICS_2026-08-26.md is the methods source for all of this.

Marker convention, every figure
  coloured      passes both gates
  grey FILLED   collapsed  (pct_excess > cut)
  grey OUTLINE  untrained  (ctt <= cut)
  a config failing BOTH is drawn as collapsed, the more specific pathology

Axes
  - AUPRC axes are refitted to the data and SHARED across every panel of a
    figure, so panels compare by eye. `--auprc-lo/--auprc-hi` pin them.
  - Chance is the PREVALENCE (0.091 sweep / 0.093 holdout / 0.081 single cell),
    never 0.5. Stated in each panel title; no reference lines are drawn.
  - Never rank on plain accuracy: with 23 positives in 252, all-negative scores
    0.909.

RCS needs an UNBALANCED corpus
  A balanced train_input.txt repeats whole documents (413 lines for 252 files);
  duplicate columns inflate RCS. Pass --train-input train_input_original.txt on
  a class-balanced arm.

Usage
  # everything, in order
  python3 7_model_config_evaluation.py all --run <RUN> \\
      --configs-file <RUN>/sweep/upsgrid_perf/upsgrid_configs.txt \\
      --train-input <RUN>/sweep/preprocessed/train_input_original.txt

  # one stage
  python3 7_model_config_evaluation.py geom    --run <RUN> --pool 16
  python3 7_model_config_evaluation.py summary --run <RUN>
  python3 7_model_config_evaluation.py figures --run <RUN> --figures 1,2

  # the holdout arm
  python3 7_model_config_evaluation.py all --run <RUN> --arm-dir holdout \\
      --model-glob starspace_trained_model.tsv

Outputs, all under <perf-dir> (default <run>/<arm-dir>/upsgrid_perf)
  <prefix>_region_geometry.csv          geom
  <prefix>_region_geometry_RANDOM.csv   rand
  <prefix>_rcs.csv                      rcs
  <prefix>_nl_auc.csv                   nlauc
  <prefix>_nearest.csv                  nlhard
  <prefix>_REGION_METRICS.csv           merge
  <prefix>_SUMMARY_TABLE.csv            summary
  <prefix>_MAIN_EFFECTS.csv             summary
  <prefix>_MAIN_EFFECT_SPREAD.csv       summary
  <prefix>_TWOWAY_CELLS.csv             summary --interactions
  <run>/sweep/external_perf/*_scored.csv   external
  <run>/figures_final/fig{1..5}*        figures
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import glob
import importlib.util
import multiprocessing as mp
import sys
import traceback

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

SRC = os.path.dirname(os.path.abspath(__file__))
SCATTER = "/home/bx2ur/code/geniml_dev/geniml/bedspace/tests/plot_param_scatter.py"

# ---- region geometry -------------------------------------------------------
N_FAST = 4000        # regions for silhouette / dip / pct  (pdist is O(n^2))
N_RCS = 60000        # regions for the reconstruction test
SEED = 42
THRESHOLDS = (0.10, 0.25, 0.40)
N_NULL = 5
MODEL_GLOB = "model_*.tsv"   # sweep arm; holdout arm: starspace_trained_model.tsv

# ---- nearest label ---------------------------------------------------------
LABELS = ["Renin", "Control", "Tumoral"]
TARGET = "Renin"


# =========================================================================== #
# PART 1 -- region geometry (fast / rand / rcs arms)
# =========================================================================== #

def _glob_init(model_glob):
    """Set MODEL_GLOB inside a worker.

    The pool start method is forkserver (see the bottom of this file), so each
    worker RE-IMPORTS this module and gets the default above -- the
    `global MODEL_GLOB` assignment in main() applies to the parent only. Without
    this initializer, --model-glob silently works for discovering config dirs in
    the parent and is ignored by every worker, which on the holdout arm (models
    named starspace_trained_model.tsv) made all 266 configs fail with
    `IndexError: list index out of range` from load_regions and still wrote a
    CSV -- one column, no metrics, exit status 0.
    """
    global MODEL_GLOB
    MODEL_GLOB = model_glob


def load_regions(cfg_dir):
    """(names, embedding) for the region rows of a config's StarSpace .tsv."""
    p = sorted(glob.glob(os.path.join(cfg_dir, MODEL_GLOB)))[0]
    df = pd.read_csv(p, sep="\t", header=None, engine="c")
    nm = df[0].astype(str)
    keep = ~nm.str.startswith("__label__")
    R = df.loc[keep].iloc[:, 1:].to_numpy(np.float64)
    ok = np.linalg.norm(R, axis=1) > 0
    return nm[keep].to_numpy()[ok], R[ok]


def l2(X):
    return X / np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)


# --------------------------------------------------------------------------- #
# fast arm
# --------------------------------------------------------------------------- #
def _dip_axis(X, diptest, KMeans):
    """dip of the 1-D projection onto the k=2 centroid-difference direction."""
    lb = KMeans(2, n_init=10, random_state=0).fit_predict(X)
    a, b = X[lb == 0], X[lb == 1]
    if len(a) < 50 or len(b) < 50:
        return np.nan, lb
    w = a.mean(0) - b.mean(0)
    n = np.linalg.norm(w)
    if n < 1e-12:
        return np.nan, lb
    return float(diptest.dipstat(np.sort(X @ (w / n)))), lb


def fast_one(cfg_dir):
    import diptest
    from scipy.spatial.distance import pdist
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    cfg = os.path.basename(cfg_dir)
    try:
        names, R = load_regions(cfg_dir)
        X = l2(R)
        X = np.ascontiguousarray(X[np.random.RandomState(SEED).permutation(len(X))[:N_FAST]])
        row = dict(config=cfg, d=X.shape[1], n_regions=len(R), n_used=len(X))

        d = pdist(X)
        row["dip_dist"] = float(diptest.dipstat(d))
        for t in THRESHOLDS:
            row[f"pct_lt_{str(t).replace('.', '')}"] = float(100 * (d < t).mean())
        f = (d < 0.25).mean()
        row["eff_n_distinct"] = float(1.0 / f) if f > 1e-9 else np.inf
        row["median_dist"] = float(np.median(d))

        da, lb = _dip_axis(X, diptest, KMeans)
        row["dip_axis"] = da
        row["sil_k2"] = float(silhouette_score(X, lb, sample_size=min(4000, len(X)),
                                               random_state=0)) if len(set(lb)) > 1 else np.nan

        # matched UNIMODAL null for dip_axis: same covariance, same n, same
        # normalisation, same k-means-then-project procedure.
        rs = np.random.RandomState(0)
        mu, C = X.mean(0), np.cov(X, rowvar=False)
        C = np.atleast_2d(C) + np.eye(X.shape[1]) * 1e-10
        L = np.linalg.cholesky(C)
        nulls = []
        for _ in range(N_NULL):
            G = l2(mu + rs.randn(len(X), X.shape[1]) @ L.T)
            v, _ = _dip_axis(G, diptest, KMeans)
            if np.isfinite(v):
                nulls.append(v)
        if nulls:
            row["dip_axis_null"] = float(np.mean(nulls))
            row["dip_axis_null_sd"] = float(np.std(nulls))
            row["dip_axis_z"] = float((da - np.mean(nulls)) / max(np.std(nulls), 1e-12))
        # effective dimensionality (participation ratio) -- ABSOLUTE, not the
        # fraction: the fraction does not separate (oc 0.060 vs ob 0.034).
        ev = np.linalg.eigvalsh(np.cov(X - X.mean(0), rowvar=False))
        ev = ev[ev > 0]
        row["eff_dim"] = float(ev.sum() ** 2 / (ev ** 2).sum())
        return row, None
    except Exception:
        return {"config": cfg}, traceback.format_exc(limit=3)


# --------------------------------------------------------------------------- #
# rcs arm
# --------------------------------------------------------------------------- #
_Y = _SEL = None


def _rcs_init(Y, sel, model_glob=None):
    if model_glob is not None:
        _glob_init(model_glob)
    global _Y, _SEL
    _Y, _SEL = Y, sel


def rcs_one(cfg_dir):
    import warnings
    warnings.filterwarnings("ignore")
    import sklearn.neural_network as nn
    from sklearn.compose import TransformedTargetRegressor
    from sklearn.model_selection import KFold, cross_validate
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    cfg = os.path.basename(cfg_dir)
    try:
        names, R = load_regions(cfg_dir)
        X = l2(R)
        pos = {r: i for i, r in enumerate(names)}
        idx = [(pos.get(r), j) for j, r in enumerate(_SEL)]
        keep = [(i, j) for i, j in idx if i is not None]
        if len(keep) < 2000:
            return {"config": cfg, "rcs": np.nan, "n_used": len(keep)}, None
        Xs = X[[i for i, _ in keep]]
        Ys = _Y[[j for _, j in keep]]

        def score(Xa):
            reg = nn.MLPRegressor(hidden_layer_sizes=(200,), activation="relu",
                                  solver="adam", alpha=1e-4, learning_rate_init=1e-3,
                                  max_iter=800, shuffle=True, random_state=SEED,
                                  tol=1e-4, early_stopping=True,
                                  validation_fraction=0.1, n_iter_no_change=15)
            m = TransformedTargetRegressor(
                regressor=make_pipeline(StandardScaler(), reg), transformer=StandardScaler())
            kf = KFold(n_splits=5, shuffle=True, random_state=SEED)
            o = cross_validate(m, Xa, Ys, cv=kf, n_jobs=1, return_estimator=True)
            its = [e.regressor_.named_steps["mlpregressor"].n_iter_ for e in o["estimator"]]
            return float(o["test_score"].mean()), max(its)

        r, it = score(Xs)
        perm, _ = score(Xs[np.random.RandomState(0).permutation(len(Xs))])
        return dict(config=cfg, d=Xs.shape[1], n_used=len(Xs), rcs=r,
                    rcs_perm=perm, rcs_max_iter=it, rcs_converged=bool(it < 800)), None
    except Exception:
        return {"config": cfg}, traceback.format_exc(limit=3)



# --------------------------------------------------------------------------- #
# rand arm -- per-config RANDOM-EMBEDDING baseline
# --------------------------------------------------------------------------- #
# A random embedding of the SAME SHAPE as the real one (n_used x d), drawn iid
# gaussian and L2-normalised = uniform on the unit sphere: the "no structure at
# all" floor for every geometry column.
#
# It depends only on (n_used, d), never on the config's scale, because
#   * the geometry metrics are computed on L2-normalised vectors, and
#   * CTT is scale-invariant -- its uniform-probe box is built from the data's
#     own 5-95% quantiles, so multiplying the embedding by any constant leaves
#     y/(x+y) unchanged (verified on synthetics: collapsing 500x moved it <1%).
# So the baseline is computed from (config, d, n_used) in the geometry table and
# needs no re-read of the 266 model .tsv files.
def _ctt_raw(data, seed, num_data=10000):
    """geniml/eval/ctt.py inlined: bare ratio y/(x+y); uniform ~ 0.5."""
    from sklearn.neighbors import NearestNeighbors
    np.random.seed(seed)
    n0, dim = data.shape
    num = min(num_data, n0)
    d = data[np.random.choice(n0, num)]
    ns = max(int(num * 0.1), 10)
    ds = d[np.random.choice(num, ns)]
    neigh = NearestNeighbors(n_neighbors=2, n_jobs=1).fit(d)
    sd, _ = neigh.kneighbors(ds)
    mx = np.quantile(d, 0.95, axis=0)
    mn = np.quantile(d, 0.05, axis=0)
    rd, _ = neigh.kneighbors(np.random.uniform(mn, mx, (ns, dim)), n_neighbors=1)
    x, y = np.sum(sd[:, 1] ** 2), np.sum(rd[:, 0] ** 2)
    return float(y / (x + y)) if (x + y) > 0 else np.nan


def rand_one(spec):
    import diptest
    from scipy.spatial.distance import pdist
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    cfg, d, n_used, n_seeds = spec
    try:
        acc = {}
        for s in range(n_seeds):
            rs = np.random.RandomState(1000 + s)
            G = rs.randn(int(n_used), int(d))
            X = l2(G)
            dd = pdist(X)
            lb = KMeans(2, n_init=10, random_state=0).fit_predict(X)
            da, _ = _dip_axis(X, diptest, KMeans)
            ev = np.linalg.eigvalsh(np.cov(X - X.mean(0), rowvar=False))
            ev = ev[ev > 0]
            f = (dd < 0.25).mean()
            for k, v in (("sil_k2_rand", silhouette_score(X, lb, sample_size=min(4000, len(X)),
                                                          random_state=0)),
                         ("ctt_rand", _ctt_raw(G, seed=42 + s)),
                         ("dip_dist_rand", diptest.dipstat(dd)),
                         ("dip_axis_rand", da),
                         ("pct_lt_025_rand", 100 * f),
                         ("median_dist_rand", np.median(dd)),
                         ("eff_dim_rand", ev.sum() ** 2 / (ev ** 2).sum())):
                acc.setdefault(k, []).append(float(v))
        row = {"config": cfg}
        for k, v in acc.items():
            row[k] = float(np.mean(v))
        row["eff_n_distinct_rand"] = (1.0 / (row["pct_lt_025_rand"] / 100)
                                      if row["pct_lt_025_rand"] > 1e-7 else np.inf)
        return row, None
    except Exception:
        return {"config": cfg}, traceback.format_exc(limit=3)


def build_occurrence(sweep_dir, train_input=None):
    """region x training-FILE binary matrix from the shared preprocess corpus.

    ONE COLUMN PER LINE, and each line is one training BED file -- so the file
    passed here must hold each training sample EXACTLY ONCE.

    That is not automatic since class balancing landed (2026-08-28). A balanced
    `train_input.txt` repeats whole documents -- the upsampled sweep arm has 413
    lines for 252 distinct files -- and feeding it here would build duplicate
    columns, making every region's occurrence vector artificially easy to
    reconstruct and inflating RCS for reasons that have nothing to do with the
    embedding. Pass `train_input` pointing at the UNBALANCED corpus
    (`train_input_original.txt` for a balanced sweep arm); the default keeps the
    historical behaviour for arms that were never balanced.
    """
    path = train_input or os.path.join(sweep_dir, "preprocessed", "train_input.txt")
    docs = []
    with open(path) as f:
        for line in f:
            docs.append({t for t in line.split() if not t.startswith("__label__")})
    universe = sorted(set().union(*docs))
    rs = np.random.RandomState(SEED)
    sel = [universe[i] for i in rs.permutation(len(universe))[:N_RCS]]
    idx = {r: i for i, r in enumerate(sel)}
    Y = np.zeros((len(sel), len(docs)), np.float32)
    for j, dd in enumerate(docs):
        for r in dd:
            i = idx.get(r)
            if i is not None:
                Y[i, j] = 1.0
    print(f"occurrence matrix {Y.shape} from {len(docs)} files "
          f"({path}), universe {len(universe):,}", flush=True)
    if len(docs) != len({frozenset(d) for d in docs}):
        print(f"  WARNING: {len(docs) - len({frozenset(d) for d in docs})} "
              f"DUPLICATE documents in {path}. This looks like a BALANCED "
              f"train_input; RCS needs each training file exactly once. Pass "
              f"--train-input pointing at the unbalanced corpus.", flush=True)
    return Y, sel



# =========================================================================== #
# PART 2 -- nearest label: the model axis
# =========================================================================== #

def nl_scores(dist_csv, labels=LABELS, target=TARGET):
    r = pd.read_csv(dist_csv)
    r["filename"] = r.filename.apply(lambda x: str(x).rsplit("/", 1)[-1])
    w = r.pivot_table(index="filename", columns="search_term", values="score",
                      aggfunc="first")
    fl = r.drop_duplicates("filename").set_index("filename")["file_label"]
    w = w.join(fl)
    w = w[~w.index.isin(set(labels))]          # drop the label anchor rows
    cols = [c for c in labels if c in w.columns]
    if target not in cols or len(cols) < 2 or w.empty:
        return None
    others = [c for c in cols if c != target]
    margin = w[others].min(axis=1) - w[target]   # >0 => nearest label IS Renin
    y = (w["file_label"].astype(str) == target).astype(int).to_numpy()
    s = margin.to_numpy(float)
    if y.sum() == 0 or y.sum() == len(y):
        return None
    o = np.argsort(-s, kind="stable"); ys = y[o]
    tp = np.cumsum(ys)
    prec = tp / np.arange(1, len(ys) + 1)
    auprc = float((prec * ys).sum() / ys.sum())
    rk = pd.Series(s).rank().to_numpy()
    npos, nneg = int(y.sum()), int((y == 0).sum())
    auc = float((rk[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))
    # the hard-assignment numbers, recomputed here so both live in one table
    pred = (s > 0)
    sens = float(pred[y == 1].mean()); spec = float((~pred[y == 0]).mean())
    return dict(n=len(y), n_pos=npos, prevalence=npos / len(y),
                nl_auc=auc, nl_auprc=auprc,
                nl_bal_acc=(sens + spec) / 2, nl_sens=sens, nl_spec=spec,
                nl_margin_median_pos=float(np.median(s[y == 1])),
                nl_margin_median_neg=float(np.median(s[y == 0])))



def nearest_label_metrics(dist_csv, labels, target="Renin"):
    """Nearest-label classification straight from the distance table.

    Assign each sample to the label vector it is closest to (min cosine
    distance), then score that assignment against file_label. This is a
    MODEL-quality measure: no score formula, no weights, no threshold.
    """
    r = pd.read_csv(dist_csv)
    r["filename"] = r.filename.apply(lambda x: str(x).rsplit("/", 1)[-1])
    w = r.pivot_table(index="filename", columns="search_term", values="score",
                      aggfunc="first")
    fl = r.drop_duplicates("filename").set_index("filename")["file_label"]
    w = w.join(fl)
    w = w[~w.index.isin(set(labels))]          # drop the label anchor rows
    cols = [c for c in labels if c in w.columns]
    if len(cols) < 2 or w.empty:
        return None
    pred_is_target = (w[cols].idxmin(axis=1) == target).to_numpy()
    y = (w["file_label"].astype(str) == target).to_numpy()
    if y.sum() == 0 or y.sum() == len(y):
        return None
    sens = pred_is_target[y].mean()
    spec = (~pred_is_target[~y]).mean()
    return dict(nl_acc=float((pred_is_target == y).mean()),
                nl_bal_acc=float((sens + spec) / 2),
                nl_sens=float(sens), nl_spec=float(spec),
                n=len(w), n_pos=int(y.sum()))



# =========================================================================== #
# PART 3 -- SUMMARY_TABLE and parameter main effects
# =========================================================================== #

SUMMARY_COLS = ["row_type", "config", "d", "sil_k2", "sil_k2_multik", "sil_k3",
                "ctt", "rcs", "dip_dist", "dip_axis", "dip_axis_null",
                "dip_ratio", "pct_lt_025", "pct_excess", "eff_n_distinct",
                "nl_bal_acc", "nl_auprc"]

# Two DIFFERENT silhouettes live here; they are not a k-series.
#   sil_k2         geom arm -- 4,000 L2-NORMALISED regions, labels from the dip
#                  axis, silhouette_score with its DEFAULT euclidean metric.
#   sil_k2_multik  7a multi_k_clustering -- 5,000 RAW regions, plain
#   sil_k3         KMeans(n_init=10), silhouette_score(metric="cosine").
# Only sil_k2_multik vs sil_k3 is a like-for-like pair.

AXES = [("dim", "dim"), ("epoch", "epoch"), ("margin", "margin"), ("neg", "neg")]
MODEL_METRICS = ["nl_auprc", "nl_auc", "nl_bal_acc"]
SCORE_METRICS = ["auprc", "auc", "bal_acc", "prec_rank_overall"]



def build_summary(metrics, nl_auc, rand, cluster=None):
    t = metrics.copy()
    if cluster is not None:
        keep = ["config"] + [c for c in ("sil_k2", "sil_k3") if c in cluster.columns]
        cl = cluster[keep].rename(columns={"sil_k2": "sil_k2_multik"})
        t = t.merge(cl, on="config", how="left")
    if nl_auc is not None and "nl_auprc" in nl_auc.columns:
        keep = ["config"] + [c for c in ("nl_auprc", "nl_auc") if c in nl_auc.columns]
        t = t.drop(columns=[c for c in ("nl_auprc", "nl_auc") if c in t.columns],
                   errors="ignore").merge(nl_auc[keep], on="config", how="left")
    t.insert(0, "row_type", "CONFIG")

    rows = [t.reindex(columns=SUMMARY_COLS)]
    # one control row per distinct d, averaged over the rand arm
    if rand is not None and "d" in metrics.columns:
        r = metrics[["config", "d"]].merge(rand, on="config", how="inner")
        ren = {c + "_rand": c for c in ("sil_k2", "ctt", "dip_dist", "dip_axis",
                                        "pct_lt_025", "eff_n_distinct")}
        r = r.rename(columns={k: v for k, v in ren.items() if k in r.columns})
        agg = {c: "mean" for c in SUMMARY_COLS
               if c in r.columns and c not in ("row_type", "config", "d")}
        if agg:
            ctrl = r.groupby("d", as_index=False).agg(agg)
            ctrl["row_type"] = "RANDOM_CONTROL"
            ctrl["config"] = ctrl.d.map(lambda x: f"__RANDOM__d{int(x)}")
            rows.append(ctrl.reindex(columns=SUMMARY_COLS))
    return pd.concat(rows, ignore_index=True)


def main_effects(cfgtab, scored):
    out = []
    for col, label in AXES:
        if col not in cfgtab.columns:
            continue
        for lv, g in cfgtab.groupby(col):
            row = {"axis": label, "level": lv, "n_configs": len(g), "method": ""}
            for m in MODEL_METRICS:
                if m in g.columns:
                    row[m] = float(np.nanmedian(g[m]))
            out.append(row)
    if scored is not None:
        for method, sm in scored.groupby("method"):
            for col, label in AXES:
                if col not in sm.columns:
                    continue
                for lv, g in sm.groupby(col):
                    row = {"axis": label, "level": lv, "n_configs": g.config.nunique(),
                           "method": method}
                    for m in SCORE_METRICS:
                        if m in g.columns:
                            row[m] = float(np.nanmedian(g[m]))
                    out.append(row)
    d = pd.DataFrame(out)
    if d.empty:
        return d, pd.DataFrame()

    # spread = max(level median) - min(level median), the effect size per axis
    sp = []
    for (axis, method), g in d.groupby(["axis", "method"], dropna=False):
        row = {"axis": axis, "method": method, "n_levels": len(g)}
        for m in MODEL_METRICS + SCORE_METRICS:
            if m in g.columns and g[m].notna().any():
                row[f"spread_{m}"] = float(g[m].max() - g[m].min())
        sp.append(row)
    return d, pd.DataFrame(sp)



# =========================================================================== #
# PART 4 -- figures_final (fig1 .. fig5)
# =========================================================================== #

MM = 1 / 25.4
SURFACE, INK, MUTED, GRID = "#fcfcfb", "#1a1a19", "#55524e", "#d8d7d2"
EXCL_GREY = "#9a9691"

# (key, per-config subdirectory relative to the arm root, label)
METHODS = [("weighted_multi_learned_post", "weighted_multi (learned, post)", "#2a78d6"),
           ("pairwise", "pairwise", "#eb6834"),
           ("weighted_multi_uniform", "weighted_multi (uniform)", "#1baf7a"),
           ("line", "line", "#eda100"),
           ("projection", "signed projection", "#7d5ba6")]

MODEL_YCOL = "nl_auprc"
MODEL_YLAB = "model AUPRC (nearest-label AUPRC)"
SCORE_XLAB = "score AUPRC"


def load_scatter():
    spec = importlib.util.spec_from_file_location("pps", SCATTER)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------- arms
def arm_paths(run, arm):
    """Where each arm's four tables live. One place, so a figure cannot drift."""
    if arm in ("sweep", "holdout"):
        p = os.path.join(run, arm, "upsgrid_perf")
        return dict(perf=p,
                    metrics=os.path.join(p, "upsgrid_REGION_METRICS.csv"),
                    nl=os.path.join(p, "upsgrid_nl_auc.csv"),
                    scored=os.path.join(p, "upsgrid_scored.csv"))
    # external arms are scored with the sweep's models
    p = os.path.join(run, "sweep", "external_perf")
    return dict(perf=p, metrics=None, nl=None,
                scored=os.path.join(p, f"{arm}_scored.csv"))


ARM_LABEL = {
    "sweep": "in-sample sweep",
    "holdout": "holdout",
    "pseudobulk": "pseudobulk",
    "singlecell": "single cell",
}
ARMS = ["sweep", "holdout", "pseudobulk", "singlecell"]


def load_arm(run, arm, ctt_min, pct_max, sweep_gate=None):
    """One arm's config table: params, gate flags, model axis, score metrics.

    Returns None when the arm has no scored table -- that is the normal state
    for an external arm a run has not scored yet, and the caller draws an empty
    panel for it rather than dropping it.
    """
    P = arm_paths(run, arm)
    if not P["scored"] or not os.path.exists(P["scored"]):
        return None
    sc = pd.read_csv(P["scored"])
    if not len(sc):
        return None

    # ---- model axis + gate -------------------------------------------------
    if P["nl"] and os.path.exists(P["nl"]):
        nl = pd.read_csv(P["nl"])
        keep = ["config"] + [c for c in ("nl_auprc", "nl_auc", "nl_bal_acc",
                                         "n", "n_pos", "prevalence")
                             if c in nl.columns]
        sc = sc.drop(columns=[c for c in keep if c != "config" and c in sc.columns],
                     errors="ignore").merge(nl[keep], on="config", how="left")

    if P["metrics"] and os.path.exists(P["metrics"]):
        mt = pd.read_csv(P["metrics"])
        cols = ["config"] + [c for c in ("ctt", "pct_excess", "dip_ratio", "rcs", "d")
                             if c in mt.columns]
        gate = mt[cols].copy()
        cut = pct_max
        if cut is None and "pct_excess" in gate:
            cut = knee_cut(gate.pct_excess)
            print(f"    {arm}: pct_excess cut from the data (knee) = "
                  f"{cut:.4f}" if cut is not None else
                  f"    {arm}: knee undetermined; no degeneracy gate")
        gate["passes_ctt"] = (gate.ctt > ctt_min) if ctt_min is not None else True
        gate["passes_pct"] = ((gate.pct_excess <= cut)
                              if ("pct_excess" in gate and cut is not None) else True)
        gate["passes"] = gate.passes_ctt & gate.passes_pct
        gate["pct_cut"] = cut if cut is not None else np.nan
    elif sweep_gate is not None:
        # external arm: the SWEEP's models were used, so the SWEEP's gate applies
        gate = sweep_gate.copy()
    else:
        gate = pd.DataFrame({"config": sc.config.unique()})
        gate["passes_ctt"] = gate["passes_pct"] = gate["passes"] = True

    sc = sc.merge(gate, on="config", how="left")
    for c in ("passes", "passes_ctt", "passes_pct"):
        if c in sc.columns:
            sc[c] = sc[c].fillna(False).astype(bool)
    sc["arm"] = arm
    return sc


def knee_cut(values, trim=0.025):
    """The pct_excess cut, found from the data: where the sorted curve turns up.

    NOT USED BY DEFAULT, AND THAT IS DELIBERATE. A knee is a property of the
    DESIGN, not of degeneracy: sweep only margins 0.6-0.8 and most configs
    collapse so the knee moves up; sweep only margin 0.2 and it is meaningless.
    Letting it define "collapsed" makes the threshold depend on which configs
    were chosen -- and the configs are what is being evaluated. It also gave the
    two arms different definitions of the same word (0.2433 vs 0.1960), which is
    a fact about two config populations, not about the embeddings.

    The shipped gate is an absolute capacity requirement instead:
    eff_n_distinct >= 1000, i.e. pct_excess <= 0.1. Kept here because "where
    does this run's distribution bend" is still a reasonable diagnostic to ask.

    THE METHOD is Kneedle: sort the values, normalise rank and value to [0, 1],
    and take the point furthest BELOW the chord joining the first and last --
    the point of maximum curvature, i.e. where the curve stops creeping and
    starts climbing. That is literally "where it sharply increases".

    WHY NOT THE LARGEST GAP, which is the obvious first idea: this distribution
    spans 0.001 to 1.9 and is heavy-tailed, so the biggest raw jumps all sit in
    the far tail (1.61 -> 1.82) and would put the cut at ~1.6, passing 142/144.
    The largest RATIO jump fails the other way and is unstable across arms
    (19/144 on the sweep vs 72/144 on the holdout). Kneedle on the raw sorted
    curve is the one that is both principled and stable.

    WHY TRIM THE TOP 2.5% FIRST. Kneedle's chord is anchored on the extremes, so
    one runaway config drags it. Measured here: untrimmed, the holdout cut is
    0.4565 (118/144); trimming the two most extreme configs it is 0.1960
    (95/144) and stays there at 5% and 10% trim. The sweep is insensitive
    (0.2433 at every trim). The trimmed values are also far closer to each
    other, which is what you want from a rule applied per arm.
    """
    v = np.sort(np.asarray([x for x in np.ravel(values) if np.isfinite(x)],
                           dtype=float))
    if len(v) < 10:
        return None
    core = v[v <= np.quantile(v, 1.0 - trim)]
    if len(core) < 10 or core[-1] <= core[0]:
        return None
    x = np.arange(len(core), dtype=float)
    xn = (x - x[0]) / (x[-1] - x[0])
    yn = (core - core[0]) / (core[-1] - core[0])
    return float(core[int(np.argmax(xn - yn))])


def excluded_split(d):
    """(collapsed, untrained). A config failing both counts as collapsed."""
    if "passes_pct" not in d.columns:
        return d.iloc[:0], d.iloc[:0]
    collapsed = d[~d.passes_pct]
    untrained = d[d.passes_pct & ~d.passes_ctt]
    return collapsed, untrained


def draw_excluded(ax, d, xcol, ycol, s=13):
    """Grey context markers: FILLED collapsed, OUTLINE untrained. Below the data."""
    col, unt = excluded_split(d)
    if len(col):
        ax.scatter(col[xcol], col[ycol], s=s, c=EXCL_GREY, linewidths=0,
                   zorder=2, label="_nolegend_")
    if len(unt):
        ax.scatter(unt[xcol], unt[ycol], s=s, facecolors="none",
                   edgecolors=EXCL_GREY, linewidths=0.7, zorder=2,
                   label="_nolegend_")


# fig1 marker convention. Shape carries the GATE, colour carries the parameter --
# so a collapsed config still shows which dim/epoch/margin/neg produced it.
# Greying them out discarded exactly the information the figure is colouring by.
GATE_MARKER = {"pass": ("o", "passes gate"),
               "collapsed": ("X", "collapsed (pct_excess > cut)"),
               "untrained": ("^", "untrained (ctt <= cut)")}


def fig1_panel(ax, groups, col, title, use_log, xcol, ycol, show_y, show_x):
    """One scatter, coloured by parameter `col`, shaped by gate status.

    ONE norm across every group, computed over all configs. Drawing the groups
    with separate scatter calls would otherwise give each its own colour scale
    and make a collapsed dim=150 the same colour as a passing dim=25.
    """
    from matplotlib.colors import LogNorm, Normalize
    allv = pd.to_numeric(pd.concat([g[col] for g in groups.values() if len(g)]),
                         errors="coerce").to_numpy(dtype=float)
    good = np.isfinite(allv)
    norm = None
    if good.any():
        lo_, hi_ = np.nanmin(allv[good]), np.nanmax(allv[good])
        norm = (LogNorm(vmin=lo_, vmax=hi_)
                if use_log and lo_ > 0 and hi_ > 0 and hi_ / lo_ > 10
                else Normalize(vmin=lo_, vmax=hi_))
    sc = None
    for key in ("pass", "collapsed", "untrained"):
        g = groups.get(key)
        if g is None or not len(g):
            continue
        mk, _lab = GATE_MARKER[key]
        v = pd.to_numeric(g[col], errors="coerce").to_numpy(dtype=float)
        sc = ax.scatter(g[xcol].to_numpy(float), g[ycol].to_numpy(float),
                        s=26 if key != "pass" else 16, c=v, cmap="viridis",
                        norm=norm, marker=mk, edgecolors="white",
                        linewidths=0.35, zorder=3 if key == "pass" else 2)
    if sc is not None:
        cb = ax.figure.colorbar(sc, ax=ax, pad=0.02, fraction=0.035, aspect=28)
        cb.ax.tick_params(labelsize=5.2, length=2, color=MUTED)
        cb.outline.set_visible(False)
    ax.set_title(title, pad=3)
    ax.grid(color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelleft=show_y, labelbottom=show_x, labelsize=6)


def gate_legend(fig, ax, n_col, n_unt):
    handles = [
        Line2D([], [], marker="o", ls="", ms=4.5, mfc=EXCL_GREY, mec=EXCL_GREY,
               label=f"collapsed ({n_col})"),
        Line2D([], [], marker="o", ls="", ms=4.5, mfc="none", mec=EXCL_GREY,
               label=f"untrained, ctt ≤ gate ({n_unt})"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=5.4, loc="lower right",
              handletextpad=0.3, borderpad=0.2, labelspacing=0.3)


def shared_range(values, pad_frac=0.04, lo=None, hi=None):
    """One axis range for every panel of a figure.

    Refitted to the data because every AUPRC here sits high in a narrow band and
    a 0..1 axis would spend most of its length empty -- but SHARED, so panels
    stay comparable by eye. That combination is the whole point; refitting each
    panel independently would make them prettier and incomparable.
    """
    v = np.asarray([x for x in np.ravel(values) if np.isfinite(x)], dtype=float)
    if lo is None and hi is None and not len(v):
        return (0.0, 1.0)
    a = float(lo) if lo is not None else float(v.min())
    b = float(hi) if hi is not None else float(v.max())
    if b <= a:
        a, b = a - 0.01, b + 0.01
    # Pad even when both ends are PINNED. A point sitting exactly on a pinned
    # limit is half-clipped by the axis and reads as missing -- and the values
    # here pile up on their ceiling (45/68 configs at nl_auprc = 1.000), so the
    # clipped marks are the ones that matter most.
    pad = (b - a) * pad_frac
    return (a - pad, b + pad)


def style(ax, show_y=True, show_x=True):
    ax.grid(color=GRID, lw=0.4, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(labelsize=6, labelleft=show_y, labelbottom=show_x)


def prevalence_note(d):
    if "prevalence" in d.columns and d.prevalence.notna().any():
        return f"chance = prevalence {float(d.prevalence.dropna().iloc[0]):.3f}"
    if {"n_pos", "n_samples"} <= set(d.columns) and d.n_pos.notna().any():
        return (f"chance = prevalence "
                f"{float(d.n_pos.iloc[0]) / float(d.n_samples.iloc[0]):.3f}")
    return "chance = prevalence"


def empty_panel(ax, label, why="no data for this run"):
    """A visible placeholder, NOT a removed panel.

    Ticks are hidden with tick_params, never with set_xticks([]) -- on a figure
    built with sharex/sharey that CLEARS THE SHARED AXIS, so one empty panel
    silently strips the tick labels off every populated panel in its row or
    column. That is how the holdout row lost its landmark labels the first time.
    """
    ax.text(0.5, 0.5, f"{label}\n{why}", ha="center", va="center",
            fontsize=6.2, color=MUTED, transform=ax.transAxes)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def label_bottom_populated(ax, populated):
    """Put the x tick labels on the lowest POPULATED row of each column.

    `show_x=(ri == nrow-1)` is wrong whenever the last row is an arm this run
    has no data for: the labels land on an empty placeholder and every real
    panel above it loses them. `populated` is a boolean grid of the same shape.
    """
    nrow, ncol = len(populated), len(populated[0])
    for ci in range(ncol):
        rows = [ri for ri in range(nrow) if populated[ri][ci]]
        if not rows:
            continue
        for ri in rows:
            ax[ri][ci].tick_params(labelbottom=(ri == rows[-1]))


def save(fig, out_base):
    for ext in ("png", "pdf"):
        fig.savefig(f"{out_base}.{ext}", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {out_base}.png / .pdf")


# ============================================================== FIG 1
def fig1(run, arms, outdir, ctt_min, pct_max, top_n, lo, hi, arm_name=None):
    """CTT vs model AUPRC, one ROW per arm, one COLUMN per swept parameter.

    BOTH ARMS IN ONE FIGURE. fig1 used to be arm-specific -- one embedding, one
    CTT, one model axis -- which forced a separate figure set per arm. Stacking
    the arms as rows puts the in-sample and held-out views side by side on ONE
    shared pair of axes, which is the comparison anyone reading this actually
    wants, and removes the duplicate output directories.

    EXPECT A FLAT CLOUD among the gated configs. CTT is a gate, not a score:
    once the untrained and collapsed configs are removed it carries no
    continuous quality signal. A flat panel is the result, not a plotting bug.

    ONLY PARAMETERS THAT VARY GET A COLUMN. This design holds lr, minCount and
    adagrad constant, so those panels would be one colour with a colourbar
    spanning floating-point noise (0.00090-0.00110 for a constant lr = 0.001),
    inviting a reader to see a gradient that cannot exist.

    Shape carries the GATE, colour the parameter, on ONE scale across all three
    groups -- so a collapsed dim=150 and a passing dim=150 are the same colour.
    """
    pps = load_scatter()
    rows = [a for a in ("sweep", "holdout")
            if arms.get(a) is not None and MODEL_YCOL in arms[a].columns]
    if not rows:
        print("  fig1: no arm with a REGION_METRICS/nl table; skipping"); return

    per = {}
    for a in rows:
        cfg = arms[a].drop_duplicates("config").dropna(subset=["ctt", MODEL_YCOL])
        per[a] = cfg
    allcfg = pd.concat(per.values(), ignore_index=True)

    axes_all = [(c, t, lg) for c, t, lg in pps.AXES if c in allcfg.columns]
    varying = [(c, t, lg) for c, t, lg in axes_all
               if allcfg[c].nunique(dropna=True) > 1]
    dropped = [t for c, t, _ in axes_all if allcfg[c].nunique(dropna=True) <= 1]

    # one shared pair of axes across every panel and both arms
    yr = shared_range(allcfg[MODEL_YCOL], lo=lo, hi=hi)
    xr = shared_range(allcfg["ctt"])

    ncol = len(varying) + 1                       # +1 for the per-arm summary
    fig, ax = plt.subplots(len(rows), ncol,
                           figsize=(46 * ncol * MM, 56 * len(rows) * MM),
                           sharex=False, sharey=True, squeeze=False)
    for ri, a in enumerate(rows):
        cfg = per[a]
        passing = cfg[cfg.passes]
        collapsed, untrained = excluded_split(cfg)
        groups = {"pass": passing, "collapsed": collapsed, "untrained": untrained}
        pool = cfg[cfg.passes_pct]
        top = pps.select_passing(pool, "ctt", MODEL_YCOL, -np.inf, -np.inf, top_n)

        for ci, (col, title, use_log) in enumerate(varying):
            axx = ax[ri][ci]
            fig1_panel(axx, groups, col, title if ri == 0 else "", use_log,
                       "ctt", MODEL_YCOL, show_y=(ci == 0),
                       show_x=(ri == len(rows) - 1))
            if len(top):
                axx.scatter(top["ctt"], top[MODEL_YCOL], s=48, facecolors="none",
                            edgecolors=INK, linewidths=0.7, zorder=5)
                for k, (_, r) in enumerate(top.iterrows()):
                    axx.annotate(r["tag"], (r["ctt"], r[MODEL_YCOL]),
                                 textcoords="offset points",
                                 xytext=(4, 3 if k % 2 == 0 else -7),
                                 fontsize=5.0, color=INK, zorder=6,
                                 annotation_clip=True)
            axx.set_xlim(*xr); axx.set_ylim(*yr)
            if ci == 0:
                axx.set_ylabel(f"{ARM_LABEL[a]}\n{prevalence_note(arms[a])}",
                               fontsize=6.2, color=MUTED)

        rho, pv = (spearmanr(passing["ctt"], passing[MODEL_YCOL])
                   if len(passing) > 2 else (np.nan, np.nan))
        col_n, unt_n = len(collapsed), len(untrained)
        ymax = float(cfg[MODEL_YCOL].max())
        n_sat = int((passing[MODEL_YCOL] >= ymax - 1e-4).sum())
        n_vis = passing.round({"ctt": 3, MODEL_YCOL: 3}).drop_duplicates(
            ["ctt", MODEL_YCOL]).shape[0]
        txt = (f"{ARM_LABEL[a]}\n"
               f"{len(passing)} of {len(cfg)} pass the gate\n"
               f"ctt > {ctt_min:g} and pct_excess ≤ {pct_max:g}\n"
               f"(= eff_n_distinct ≥ {int(round(100 / pct_max))})\n"
               f"Spearman ρ = {rho:.2f} gated (p = {pv:.1e})\n"
               f"● passes ({len(passing)})  ✕ collapsed ({col_n})"
               f"  ▲ untrained ({unt_n})\n"
               f"{n_sat}/{len(passing)} at the {ymax:.3f} ceiling"
               f" ({n_vis} distinct positions)\n\n"
               f"top {len(top)} by model AUPRC AND CTT:\n"
               + "\n".join(f"{r.tag}. {r.config}" for r in top.itertuples()))
        ax[ri][-1].text(0.5, 0.99, txt, ha="center", va="top", fontsize=5.0,
                        color=MUTED, transform=ax[ri][-1].transAxes,
                        linespacing=1.45)
        ax[ri][-1].axis("off")
        if len(top):
            top_cols = ["config", "ctt", MODEL_YCOL] + [c for c, _, _ in varying]
            top.assign(arm=a)[[c for c in top_cols if c in top.columns] + ["arm"]] \
               .to_csv(os.path.join(outdir,
                       f"fig1_ctt_vs_nearestlabel_AUPRC_top_{a}.csv"), index=False)

    fig.supxlabel("CTT (geniml.eval cluster tendency)   [all configs drawn; "
                  "shape = gate status, colour = the parameter above]"
                  + (f"   ·   held constant: {', '.join(d.split()[0] for d in dropped)}"
                     if dropped else ""),
                  fontsize=7.5, y=0.015)
    fig.supylabel(MODEL_YLAB, fontsize=7.5, x=0.005)
    fig.suptitle("CTT vs model AUPRC — in-sample sweep and holdout, "
                 "one shared axis pair", fontsize=8.5, y=0.997)
    fig.tight_layout(rect=(0.02, 0.025, 1, 0.96))
    save(fig, os.path.join(outdir, "fig1_ctt_vs_nearestlabel_AUPRC"))


# ============================================================== FIG 2
def fig2(run, arms, outdir, lo, hi):
    """Score AUPRC vs model AUPRC, one panel per scoring method x arm.

    ALL configs are drawn; the excluded ones in grey. They are context, not
    data: a collapsed embedding can score well at sample level for reasons
    unrelated to having learned anything, so the correlation printed on each
    panel is computed on the GATED subset only.
    """
    rows = [a for a in ARMS]
    sweep = arms.get("sweep")
    # The external arms have no embedding of their own -- they are the SWEEP's
    # models applied to another test set -- so their model axis IS the sweep's
    # nl_auprc, joined by config. Without this join they have no y value and the
    # panels were drawn empty.
    plotted = {}
    for a in rows:
        d = arms.get(a)
        if d is None:
            continue
        d = d.copy()
        if MODEL_YCOL not in d.columns and sweep is not None:
            d = d.merge(sweep.drop_duplicates("config")[["config", MODEL_YCOL]],
                        on="config", how="left")
        plotted[a] = d
    fig, ax = plt.subplots(len(rows), len(METHODS),
                           figsize=(40 * len(METHODS) * MM, 42 * len(rows) * MM),
                           squeeze=False)
    # Range over what is ACTUALLY DRAWN. Computing it over every loaded arm put
    # pseudobulk's 0.22 on the axis while its panels were empty, so the figure
    # ran from 0.25 with nothing below 0.5.
    xs = [d["auprc"].to_numpy(float) for d in plotted.values() if "auprc" in d]
    ys = [d[MODEL_YCOL].to_numpy(float) for d in plotted.values()
          if MODEL_YCOL in d.columns]
    xr = shared_range(np.concatenate(xs) if xs else [0, 1], lo=lo, hi=hi)
    yr = shared_range(np.concatenate(ys) if ys else [0, 1], lo=lo, hi=hi)

    populated = [[False] * len(METHODS) for _ in rows]
    for ri, a in enumerate(rows):
        d = plotted.get(a)
        for ci, (meth, mlab, colr) in enumerate(METHODS):
            axx = ax[ri][ci]
            if d is None or MODEL_YCOL not in d.columns:
                empty_panel(axx, f"{ARM_LABEL[a]}\n{mlab}")
                continue
            sub = d[d.method == meth].dropna(subset=["auprc", MODEL_YCOL])
            if len(sub) < 2:
                empty_panel(axx, f"{ARM_LABEL[a]}\n{mlab}", "too few configs")
                continue
            draw_excluded(axx, sub, "auprc", MODEL_YCOL, s=12)
            g = sub[sub.passes]
            axx.scatter(g["auprc"], g[MODEL_YCOL], s=12, c=colr,
                        edgecolors="white", linewidths=0.3, zorder=3)
            rho, _ = spearmanr(g["auprc"], g[MODEL_YCOL]) if len(g) > 2 else (np.nan, 0)
            # A PINNED AXIS CAN HIDE POINTS. With --auprc-lo/--auprc-hi the range
            # no longer follows the data, so say per panel how many fall outside
            # rather than let a reader infer the panel is sparse.
            out = int(((sub["auprc"] < xr[0]) | (sub["auprc"] > xr[1]) |
                       (sub[MODEL_YCOL] < yr[0]) | (sub[MODEL_YCOL] > yr[1])).sum())
            ttl = f"ρ = {rho:.2f} gated (n={len(g)})"
            if out:
                ttl += f"\n{out}/{len(sub)} outside the axis"
            axx.set_title((f"{mlab}\n" if ri == 0 else "") + ttl,
                          fontsize=6.0, pad=3)
            axx.set_xlim(*xr); axx.set_ylim(*yr)
            style(axx, show_y=(ci == 0), show_x=True)
            populated[ri][ci] = True
            if ci == 0:
                # arm name only -- the metric name is on the figure once, or it
                # prints four times and collides with the row labels
                axx.set_ylabel(ARM_LABEL[a], fontsize=6.6, color=MUTED)
    label_bottom_populated(ax, populated)
    col_n = unt_n = 0
    if arms.get("sweep") is not None:
        c, u = excluded_split(arms["sweep"].drop_duplicates("config"))
        col_n, unt_n = len(c), len(u)
    gate_legend(fig, ax[0][-1], col_n, unt_n)
    fig.supylabel(MODEL_YLAB, fontsize=7.5, x=0.004)
    fig.supxlabel(SCORE_XLAB + "   [all configs drawn; ρ computed on the gated "
                  "subset; one shared axis pair across every panel]",
                  fontsize=7.5, y=0.01)
    fig.suptitle("Score AUPRC vs model AUPRC — every scoring method, every arm",
                 fontsize=8.5, y=0.997)
    fig.tight_layout(rect=(0.025, 0.02, 1, 0.97))
    save(fig, os.path.join(outdir, "fig2_score_AUPRC_vs_nearestlabel_AUPRC"))


# ============================================================== FIG 3
def fig3(run, arms, outdir, lo, hi):
    """Method comparison on the GATED configs only, one panel per arm.

    Gated only, unlike fig2. This figure is a claim about which scoring method
    to prefer, and a collapsed or untrained embedding scores well at sample
    level for reasons unrelated to having learned anything -- including them
    would let a pathology vote on the ranking.

    The external arms carry no embedding of their own: they are the SWEEP's
    models applied to new test samples, so their model axis IS the sweep's
    nl_auprc, joined by config.
    """
    sweep = arms.get("sweep")
    fig, ax = plt.subplots(1, len(ARMS), figsize=(52 * len(ARMS) * MM, 62 * MM),
                           squeeze=False)
    ax = ax[0]
    xs, ys = [], []
    for a in ARMS:
        d = arms.get(a)
        if d is None:
            continue
        g = d[d.passes] if "passes" in d.columns else d
        if "auprc" in g:
            xs.append(g["auprc"].to_numpy(float))
        if MODEL_YCOL in g:
            ys.append(g[MODEL_YCOL].to_numpy(float))
    xr = shared_range(np.concatenate(xs) if xs else [0, 1], lo=lo, hi=hi)
    yr = shared_range(np.concatenate(ys) if ys else [0, 1], lo=lo, hi=hi)

    for i, a in enumerate(ARMS):
        d = arms.get(a)
        axx = ax[i]
        if d is None:
            empty_panel(axx, ARM_LABEL[a]); continue
        g = d[d.passes].copy() if "passes" in d.columns else d.copy()
        if MODEL_YCOL not in g.columns and sweep is not None:
            g = g.merge(sweep.drop_duplicates("config")[["config", MODEL_YCOL]],
                        on="config", how="left")
        g = g.dropna(subset=["auprc", MODEL_YCOL])
        if not len(g):
            empty_panel(axx, ARM_LABEL[a], "no gated configs"); continue
        for meth, mlab, colr in METHODS:
            s = g[g.method == meth]
            if not len(s):
                continue
            axx.scatter(s["auprc"], s[MODEL_YCOL], s=11, c=colr, alpha=0.75,
                        edgecolors="white", linewidths=0.25, zorder=3)
            axx.scatter([s["auprc"].median()], [s[MODEL_YCOL].median()], s=70,
                        marker="D", c=colr, edgecolors=INK, linewidths=0.7,
                        zorder=5)
        axx.set_xlim(*xr); axx.set_ylim(*yr)
        style(axx, show_y=(i == 0), show_x=True)
        axx.set_title(f"{ARM_LABEL[a]}\n{g.config.nunique()} gated configs · "
                      f"{prevalence_note(d)}", fontsize=6.4, pad=3)
        axx.set_xlabel(SCORE_XLAB, fontsize=6.4, color=MUTED)
        if i == 0:
            axx.set_ylabel(MODEL_YLAB, fontsize=6.4, color=MUTED)
    # Explicit handles, not collected from the artists: the method legend must
    # list every method even when the panel it would have been collected from is
    # an arm this run has no data for.
    handles = [Line2D([], [], marker="o", ls="", ms=4.5, mfc=c, mec="white",
                      label=l) for _k, l, c in METHODS]
    ax[0].legend(handles=handles, frameon=False, fontsize=5.6, loc="lower right",
                 handletextpad=0.3, borderpad=0.2, labelspacing=0.3)
    fig.suptitle("Scoring method comparison — gated configs only "
                 "(◆ = per-method median)", fontsize=8.5, y=0.99)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.95))
    save(fig, os.path.join(outdir, "fig3_scoring_method_comparion_AUPRC"))


# ================================================= score distributions (5 + 4)
def score_table(run, arms, cache, rebuild=False):
    """Per config x method x arm: medians, gap, extremes, and the rank boundary.

    THE RANK BOUNDARY is the cut the class sizes imply: sort the scores
    descending and take the midpoint between the x-th and (x+1)-th, where x is
    the number of true Renin. If the ordering were perfect every sample above it
    is Renin and every one below is not; if the score were also calibrated it
    would sit at 0. Its spread across configs is what says whether a single
    fixed threshold could ever work for this method.
    """
    if os.path.exists(cache) and not rebuild:
        return pd.read_csv(cache)
    rows = []
    for a in ARMS:
        d = arms.get(a)
        if d is None:
            continue
        for cfg in sorted(d.config.unique()):
            for meth, _mlab, _c in METHODS:
                # the external arms are the SWEEP's models applied to another
                # test set, so their score tables live under the sweep config
                f = (os.path.join(run, a, cfg, f"reninness_score_{meth}.csv")
                     if a in ("sweep", "holdout") else
                     os.path.join(run, "sweep", cfg, a, f"reninness_score_{meth}.csv"))
                if not os.path.exists(f):
                    continue
                try:
                    t = pd.read_csv(f)
                    # drop the label-to-label anchor rows, per the convention in
                    # 7a.score_table_metrics -- they are not samples
                    t = t[t["filename"].astype(str) != t["file_label"].astype(str)]
                    if a in ("pseudobulk", "singlecell"):
                        # Renin vs Non_Renin ONLY. `Recruited` is an intermediate
                        # level, expected above Non_Renin and below Renin, so it
                        # belongs to neither side of a binary cut -- scoring it as
                        # a negative penalises the models that order it correctly.
                        t = t[t["file_label"].astype(str).isin(["Renin", "Non_Renin"])]
                    y = (t["file_label"].astype(str) == "Renin").to_numpy()
                    s = t["score"].to_numpy(float)
                except Exception:
                    continue
                if not len(s) or y.sum() == 0 or y.sum() == len(y):
                    continue
                x = int(y.sum())
                srt = np.sort(s)[::-1]
                boundary = float(srt[x - 1] + srt[x]) / 2.0 if len(srt) > x else float(srt[-1])
                rows.append(dict(
                    arm=a, config=cfg, method=meth, n=len(s), n_pos=x,
                    med_pos=float(np.median(s[y])), med_neg=float(np.median(s[~y])),
                    gap=float(np.median(s[y]) - np.median(s[~y])),
                    max_pos=float(s[y].max()), min_neg=float(s[~y].min()),
                    score_min=float(s.min()), score_max=float(s.max()),
                    rank_boundary=boundary))
    out = pd.DataFrame(rows)
    if len(out):
        out.to_csv(cache, index=False)
        print(f"  wrote {cache}  ({len(out)} arm x config x method rows)")
    return out


# ============================================================== FIG 4
def fig4(run, arms, dist, outdir, lo, hi):
    """Score AUPRC vs the MEDIAN CLASS GAP, in-sample and held out.

    The gap is median(Renin) - median(non-Renin): a SCALE, not a decision, so it
    says whether the classes are pushed apart rather than merely put on the right
    side of a threshold. A method can rank well (high AUPRC) while leaving the
    classes almost superimposed, and this figure is where that shows.

    AXES. y is score AUPRC in both rows and shares ONE range. x is the gap, fixed
    and symmetric so the rows are directly comparable rather than each scaled to
    its own spread.

    THE FIXED RANGE IS [-2, 2], NOT [-1, 1]. Each median lies in [-1, 1] because
    the scores are clipped there, so their DIFFERENCE lies in [-2, 2] -- a method
    that puts the Renin median at +1 and the non-Renin median at -1 has a gap of
    2. Measured here, 613 of 1440 config x method gaps exceed 1, and `pairwise`
    never falls below 1.25, so an x limit of [-1, 1] renders its panels
    completely empty. [-2, 2] is the natural bound and hides nothing.

    The AUC-vs-AUPRC rows the earlier version carried were dropped: both are
    threshold-free summaries of the same ranking, so they agree almost by
    construction and said nothing the gap rows do not.
    """
    # ALL FOUR ARMS as rows. An earlier revision cut this to two when the
    # AUC-vs-AUPRC rows were dropped -- but what was asked for was removing that
    # CONTENT, not the external arms, which belong here exactly as they do in
    # fig2 and fig5.
    blocks = list(ARMS)
    fig, ax = plt.subplots(len(blocks), len(METHODS),
                           figsize=(50 * len(METHODS) * MM, 46 * len(blocks) * MM),
                           squeeze=False)
    yv = [arms[a]["auprc"].to_numpy(float) for a in blocks
          if arms.get(a) is not None and "auprc" in arms[a].columns]
    yr = shared_range(np.concatenate(yv) if yv else [0, 1], lo=lo, hi=hi)
    # each median is in [-1, 1], so the gap is in [-2, 2]; see the docstring
    xr = (-2.05, 2.05)

    populated = [[False] * len(METHODS) for _ in blocks]
    for ri, a in enumerate(blocks):
        d = arms.get(a)
        for ci, (meth, mlab, colr) in enumerate(METHODS):
            axx = ax[ri][ci]
            if d is None:
                empty_panel(axx, f"{ARM_LABEL[a]}\n{mlab}"); continue
            g = dist[(dist.arm == a) & (dist.method == meth)][["config", "gap"]] \
                if len(dist) else pd.DataFrame(columns=["config", "gap"])
            sub = d[d.method == meth].merge(g, on="config", how="left") \
                   .dropna(subset=["gap", "auprc"])
            if len(sub) < 2:
                empty_panel(axx, f"{ARM_LABEL[a]}\n{mlab}", "too few configs"); continue
            draw_excluded(axx, sub, "gap", "auprc", s=12)
            gg = sub[sub.passes]
            axx.scatter(gg["gap"], gg["auprc"], s=12, c=colr, edgecolors="white",
                        linewidths=0.3, zorder=3)
            rho, _ = spearmanr(gg["gap"], gg["auprc"]) if len(gg) > 2 else (np.nan, 0)
            axx.set_xlim(*xr); axx.set_ylim(*yr)
            axx.axvline(0.0, color=GRID, lw=0.8, zorder=1)
            style(axx, show_y=(ci == 0), show_x=True)
            populated[ri][ci] = True
            if ri == 0:
                axx.set_title(mlab, fontsize=6.4, pad=3)
            axx.text(0.03, 0.05, f"ρ = {rho:.2f}", transform=axx.transAxes,
                     fontsize=5.4, color=MUTED)
            if ci == 0:
                axx.set_ylabel(f"{ARM_LABEL[a]}\nscore AUPRC", fontsize=6.2,
                               color=MUTED)
    label_bottom_populated(ax, populated)
    col_n = unt_n = 0
    if arms.get("sweep") is not None:
        c, u = excluded_split(arms["sweep"].drop_duplicates("config"))
        col_n, unt_n = len(c), len(u)
    gate_legend(fig, ax[0][-1], col_n, unt_n)
    fig.supxlabel("median class gap   median(Renin) \u2212 median(non-Renin)   "
                  "[fixed to [-2, 2], the range a difference of two clipped "
                  "medians can take]", fontsize=7.5, y=0.01)
    fig.suptitle("Score AUPRC vs the median class gap \u2014 y shared, "
                 "x fixed to [-2, 2]", fontsize=8.5, y=0.995)
    fig.tight_layout(rect=(0.01, 0.02, 1, 0.96))
    save(fig, os.path.join(outdir, "fig4_score_matrix_agreement"))


# ============================================================== FIG 5
def fig5(run, arms, dist, outdir, ctt_min=0.5, pct_max=0.1):
    """Where the scores sit: class gap, clipping range, and the implied cut.

    Layout follows the published Fig5. Rows = arm, columns = scoring method, and
    within a panel five landmarks along the score axis:

        min                 the floor the method reaches
        median non-Renin
        RANK BOUNDARY       midpoint of the x-th and (x+1)-th highest score,
                            x = the number of true Renin -- the cut the CLASS
                            SIZES imply. Perfect ordering puts every sample above
                            it in the positive class; perfect calibration would
                            also put it at 0.
        median Renin
        max

    EACH LANDMARK IS DRAWN TWICE: pale grey over ALL configs, coloured over the
    GATED ones. The pair is the point -- it shows whether the gate is doing
    anything to the distribution, which a single box cannot.

    BOXES ARE QUARTILES OVER CONFIGS, whiskers the full range. So the spread is
    across models, not across samples: a wide rank-boundary box means no single
    fixed threshold transfers between configs, which is the practical question
    about a score.

    The red line at 0 is the FIXED DECISION CUT the scores are built around, not
    a chance level. Where the rank boundary sits away from it, thresholding at
    zero is miscalibrated even when the ranking is good.
    """
    QUANTS = [("score_min", "min"), ("med_neg", "median\nnon-Renin"),
              ("rank_boundary", "rank\nboundary"), ("med_pos", "median\nRenin"),
              ("score_max", "max")]
    if not len(dist):
        print("  fig5: no score-distribution table; skipping"); return

    gate = {}
    for a in ARMS:
        d = arms.get(a)
        if d is not None and "passes" in d.columns:
            gate[a] = set(d[d.passes].config.unique())

    fig, ax = plt.subplots(len(ARMS), len(METHODS),
                           # 62mm was too narrow for the two-line landmark
                           # labels: "median non-Renin" and "rank boundary" ran
                           # together into "non-Reninboundary".
                           figsize=(78 * len(METHODS) * MM, 54 * len(ARMS) * MM),
                           squeeze=False, sharex=True, sharey=True)
    populated = [[False] * len(METHODS) for _ in ARMS]
    for ri, a in enumerate(ARMS):
        sub_arm = dist[dist.arm == a]
        for ci, (meth, mlab, colr) in enumerate(METHODS):
            axx = ax[ri][ci]
            s_all = sub_arm[sub_arm.method == meth]
            if not len(s_all):
                empty_panel(axx, f"{ARM_LABEL[a]}\n{mlab}"); continue
            s_gate = s_all[s_all.config.isin(gate.get(a, set()))]
            axx.axhline(0.0, color="#c0392b", lw=1.0, ls=(0, (4, 3)), zorder=1)
            for qi, (col, _lab) in enumerate(QUANTS):
                for k, (src, face, edge) in enumerate(
                        ((s_all, "#e2e0dc", INK), (s_gate, colr, INK))):
                    v = src[col].dropna().to_numpy(float)
                    if not len(v):
                        continue
                    axx.boxplot([v], positions=[qi + (-0.19 if k == 0 else 0.19)],
                                widths=0.34, whis=(0, 100), showfliers=False,
                                patch_artist=True, medianprops=dict(color=edge, lw=1.1),
                                boxprops=dict(facecolor=face, edgecolor=edge, lw=0.6),
                                whiskerprops=dict(color=edge, lw=0.6),
                                capprops=dict(color=edge, lw=0.6), zorder=3)
            axx.set_ylim(-1.09, 1.09)   # scores are clipped to [-1, 1]; pad the rails
            axx.set_xticks(range(len(QUANTS)))
            axx.set_xticklabels([q[1] for q in QUANTS], fontsize=5.6,
                                linespacing=1.2)
            axx.set_xlim(-0.6, len(QUANTS) - 0.4)
            g_all, g_gate = s_all.gap.median(), (s_gate.gap.median() if len(s_gate) else np.nan)
            rb = s_gate.rank_boundary if len(s_gate) else s_all.rank_boundary
            axx.set_title((f"{mlab}\n" if ri == 0 else "")
                          + f"median gap {g_all:.2f} all · {g_gate:.2f} gated\n"
                            f"rank boundary {rb.min():+.2f} to {rb.max():+.2f} (gated)",
                          fontsize=5.6, pad=3, linespacing=1.3)
            style(axx, show_y=(ci == 0), show_x=True)
            populated[ri][ci] = True
            if ci == 0:
                n = int(s_all.n.iloc[0]); npos = int(s_all.n_pos.iloc[0])
                extra = ("\n(Recruited dropped)" if a in ("pseudobulk", "singlecell")
                         else "")
                axx.set_ylabel(f"{ARM_LABEL[a]}\n{n:,} scored, {npos} Renin{extra}"
                               f"\n\nreninness score", fontsize=5.8, color=MUTED)
    label_bottom_populated(ax, populated)
    handles = [
        Line2D([], [], marker="s", ls="", ms=6, mfc="#e2e0dc", mec=INK,
               label=f"all {int(dist.config.nunique())} configs"),
        Line2D([], [], marker="s", ls="", ms=6, mfc=MUTED, mec=INK,
               # state the ACTUAL cut in use -- this label read
               # "data-derived knee" for a while after the gate had been
               # restored to the fixed capacity requirement, which is exactly
               # the kind of stale caption nobody re-reads
               label=(f"passes gate (ctt > {ctt_min:g}, pct_excess ≤ "
                      f"{pct_max:g} = eff_n_distinct ≥ "
                      f"{int(round(100 / pct_max))})")),
        Line2D([], [], color="#c0392b", ls=(0, (4, 3)), lw=1.0,
               label="score = 0, the fixed decision cut"),
    ]
    fig.legend(handles=handles, frameon=False, fontsize=6.2, loc="upper right",
               bbox_to_anchor=(0.999, 0.999))
    fig.suptitle(f"Score distributions across all {int(dist.config.nunique())} "
                 f"configs — class gap, clipping range, and the cut the class "
                 f"sizes imply\nboxes are quartiles over CONFIGS, whiskers the "
                 f"full range; rank boundary = mean of the x-th and (x+1)-th "
                 f"highest score, x = number of true Renin",
                 fontsize=7.5, y=0.998, linespacing=1.4)
    fig.tight_layout(rect=(0.01, 0.01, 1, 0.955))
    save(fig, os.path.join(outdir, "fig5_score_distribution"))


# ============================================================== main


# =========================================================================== #
# PART 5 -- external arms (pseudobulk, single cell)
# =========================================================================== #
# Scored with the SWEEP's models: they add test samples, not embeddings.
# Metrics are the Renin vs Non_Renin contrast ONLY -- the intermediate
# `Recruited` level is dropped, since it belongs to neither side of a
# binary cut.
SETS = ("pseudobulk", "singlecell")
META_DIR = "/home/bx2ur/code/renin_atac/metadata"
PB_META = os.path.join(META_DIR, "renin_pseudobulk_data.csv")
SC_META = os.path.join(META_DIR, "renin_singlecell_data_with_label.csv")

_PB_LABELS = None


def ext_pseudobulk_labels():
    """filename -> Renin / Recruited / Non_Renin for the pseudobulk groups.

    WHY THIS EXISTS. `renin_pseudobulk_data.csv` has a `Group` column carrying
    only Renin and Non_Renin, and 7r was run with `--label-col Group`, so its
    score tables call all 16 non-Renin groups `Non_Renin`. But pseudobulk is the
    aggregated form of the SAME single-cell study, and the single-cell sheet
    maps each cell type to a THREE-level scheme in which seven of those groups
    are `Recruited` -- an intermediate level that belongs to neither side of a
    binary cut.

    Remapping by `cell_type` through the single-cell sheet recovers the
    published composition exactly: 18 groups, 2 Renin / 7 Recruited /
    9 Non_Renin, with no unmapped cell type. Without it the Renin-vs-Non_Renin
    metric silently scores the 7 Recruited as negatives -- 2 vs 16 rather than
    2 vs 9 -- which is the opposite of what dropping `Recruited` is meant to do.

    The 7: PC, POD_PC_MC, early_MC, early_SMC, late_MC, late_SMC,
    renin_precursors.
    """
    global _PB_LABELS
    if _PB_LABELS is None:
        sc = pd.read_csv(SC_META)
        pb = pd.read_csv(PB_META)
        ct2lab = sc.drop_duplicates("Group").set_index("Group")["label"].to_dict()
        pb["lab"] = pb["cell_type"].map(ct2lab)
        missing = pb[pb.lab.isna()]["cell_type"].tolist()
        if missing:
            print(f"WARNING: {len(missing)} pseudobulk cell types have no "
                  f"single-cell label: {missing}", file=sys.stderr)
        _PB_LABELS = dict(zip(pb["file_name"], pb["lab"]))
    return _PB_LABELS
EXT_METHODS = ("pairwise", "weighted_multi_uniform", "weighted_multi_learned_post",
           "projection", "line")
EXT_POS, EXT_NEG = "Renin", "Non_Renin"


def ext_rvsn_metrics(path, setkey=None):
    """Score-side metrics on the Renin vs Non_Renin contrast only."""
    from sklearn.metrics import average_precision_score, roc_auc_score
    d = pd.read_csv(path)
    # anchor rows (filename == file_label) are label-to-label distances, not samples
    d = d[d["filename"].astype(str) != d["file_label"].astype(str)]
    if setkey == "pseudobulk":
        # its own file_label collapses Recruited into Non_Renin -- see
        # ext_pseudobulk_labels()
        lab = ext_pseudobulk_labels()
        d = d.assign(file_label=d["filename"].astype(str).map(
            lambda x: lab.get(os.path.basename(x), None)))
        d = d.dropna(subset=["file_label"])
    d = d[d["file_label"].astype(str).isin([EXT_POS, EXT_NEG])]      # drop Recruited
    if not len(d):
        return None
    y = (d["file_label"].astype(str) == EXT_POS).to_numpy().astype(int)
    s = d["score"].to_numpy(float)
    if y.sum() == 0 or y.sum() == len(y):
        return None
    tp = int(((y == 1) & (s > 0)).sum()); fn = int(((y == 1) & (s <= 0)).sum())
    tn = int(((y == 0) & (s < 0)).sum()); fp = int(((y == 0) & (s >= 0)).sum())
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    # precision-at-rank: the top n_pos by score should be the n_pos true Renin
    o = np.argsort(-s, kind="stable")
    npos = int(y.sum())
    return dict(n_samples=len(d), n_pos=npos, n_neg=int((y == 0).sum()),
                prevalence=npos / len(d),
                bal_acc=float((sens + spec) / 2), sens=sens, spec=spec,
                auc=float(roc_auc_score(y, s)),
                auprc=float(average_precision_score(y, s)),
                prec_rank_overall=float(y[o][:npos].mean()))


def _ext_one(task):
    run, setkey, cfg = task
    out = []
    for meth in EXT_METHODS:
        f = os.path.join(run, "sweep", cfg, setkey, f"reninness_score_{meth}.csv")
        if not os.path.exists(f):
            continue
        try:
            m = ext_rvsn_metrics(f, setkey)
        except Exception:
            continue
        if m:
            out.append(dict(config=cfg, method=meth, **m))
    return out


# =========================================================================== #
# STAGE DRIVERS
# =========================================================================== #
def _want_configs(path):
    """Config allow-list, one name per line.

    Without it the arms glob EVERY directory under the arm dir, which silently
    pulls a superseded design sharing that directory into the table -- the
    class-balanced run keeps 36 configs from a previous grid beside its 144.
    """
    if not path:
        return None
    with open(path) as fh:
        return {ln.strip() for ln in fh if ln.strip()}


def _config_dirs(arm_dir, want, model_glob):
    dirs = [d for d in sorted(glob.glob(os.path.join(arm_dir, "*")))
            if os.path.isdir(d) and glob.glob(os.path.join(d, model_glob))]
    if want is not None:
        before = len(dirs)
        dirs = [d for d in dirs if os.path.basename(d) in want]
        print(f"  allow-list: {len(dirs)} of {before} config dirs kept "
              f"({len(want)} named)", flush=True)
        missing = want - {os.path.basename(d) for d in dirs}
        if missing:
            print(f"  WARNING: {len(missing)} named configs have no model on "
                  f"disk, e.g. {sorted(missing)[:3]}", flush=True)
    return dirs


def stage_geom(a, P, arm_dir, want):
    out = os.path.join(P, f"{a.prefix}_region_geometry.csv")
    dirs = _config_dirs(arm_dir, want, a.model_glob)
    print(f"  {len(dirs)} configs, pool={a.pool}", flush=True)
    rows = []
    with mp.Pool(a.pool, initializer=_glob_init, initargs=(a.model_glob,)) as pool:
        for i, (row, err) in enumerate(pool.imap_unordered(fast_one, dirs), 1):
            rows.append(row)
            print(f"  {i}/{len(dirs)} {row['config'][:44]:44s} "
                  f"dip_axis={row.get('dip_axis', float('nan')):.5f} "
                  f"pct<0.25={row.get('pct_lt_025', float('nan')):.1f}", flush=True)
            if err:
                print(err, flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"  wrote {out}")


def stage_rand(a, P, arm_dir, want):
    geom = os.path.join(P, f"{a.prefix}_region_geometry.csv")
    out = os.path.join(P, f"{a.prefix}_region_geometry_RANDOM.csv")
    if not os.path.exists(geom):
        raise SystemExit(f"missing {geom} -- run the geom stage first")
    g = pd.read_csv(geom)
    specs = [(r.config, r.d, r.n_used, a.n_seeds) for r in g.itertuples()
             if np.isfinite(r.d) and np.isfinite(r.n_used)]
    print(f"  {len(specs)} configs, seeds={a.n_seeds}, pool={a.pool}", flush=True)
    rows = []
    with mp.Pool(a.pool) as pool:
        for i, (row, err) in enumerate(pool.imap_unordered(rand_one, specs), 1):
            rows.append(row)
            if i % 25 == 0 or i == len(specs):
                print(f"  {i}/{len(specs)}", flush=True)
            if err:
                print(err, flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"  wrote {out}")


def stage_rcs(a, P, arm_dir, want):
    out = os.path.join(P, f"{a.prefix}_rcs.csv")
    dirs = _config_dirs(arm_dir, want, a.model_glob)
    Y, sel = build_occurrence(arm_dir, a.train_input)
    print(f"  {len(dirs)} configs, pool={a.rcs_pool}", flush=True)
    rows = []
    with mp.Pool(a.rcs_pool, initializer=_rcs_init,
                 initargs=(Y, sel, a.model_glob)) as pool:
        for i, (row, err) in enumerate(pool.imap_unordered(rcs_one, dirs), 1):
            rows.append(row)
            print(f"  {i}/{len(dirs)} {row['config'][:44]:44s} "
                  f"RCS={row.get('rcs', float('nan')):+.4f}", flush=True)
            if err:
                print(err, flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"  wrote {out}")


def stage_nlauc(a, P, arm_dir, want):
    out = os.path.join(P, f"{a.prefix}_nl_auc.csv")
    labels = [x for x in a.labels.split(",") if x]
    rows = []
    for f in sorted(glob.glob(os.path.join(arm_dir, "*", "raw_cosdist_rl.csv"))):
        cfg = os.path.basename(os.path.dirname(f))
        if want is not None and cfg not in want:
            continue
        m = nl_scores(f, labels=labels, target=a.target)
        if m:
            rows.append(dict(config=cfg, **m))
    d = pd.DataFrame(rows)
    d.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(d)} configs)")
    if len(d):
        print(f"  prevalence {d.prevalence.iloc[0]:.3f} | "
              f"nl_bal_acc == 1.000 for {(d.nl_bal_acc >= 0.999).sum()} configs | "
              f"nl_auprc median {d.nl_auprc.median():.3f}")


def stage_nlhard(a, P, arm_dir, want):
    out = os.path.join(P, f"{a.prefix}_nearest.csv")
    labels = [x for x in a.labels.split(",") if x]
    rows = []
    for f in sorted(glob.glob(os.path.join(arm_dir, "*", "raw_cosdist_rl.csv"))):
        cfg = os.path.basename(os.path.dirname(f))
        if want is not None and cfg not in want:
            continue
        m = nearest_label_metrics(f, labels, target=a.target)
        if m is None:
            continue
        row = dict(config=cfg, **m)
        # the anchor's own score row, carried so the table matches the
        # published nearest_*.csv schema
        sp = os.path.join(os.path.dirname(f), a.score_file)
        if a.anchor and os.path.exists(sp):
            s = pd.read_csv(sp)
            v = s[s.filename.astype(str) == a.anchor]["score"]
            row["anchor_score"] = float(v.iloc[0]) if len(v) else np.nan
        rows.append(row)
    d = pd.DataFrame(rows)
    d.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(d)} configs)")


def stage_merge(a, P, arm_dir, want):
    geom = os.path.join(P, f"{a.prefix}_region_geometry.csv")
    if not os.path.exists(geom):
        raise SystemExit(f"missing {geom} -- run the geom stage first")
    t = pd.read_csv(geom)
    extras = (
        (f"{a.prefix}_rcs.csv", ["config", "rcs", "rcs_perm", "rcs_converged"]),
        (f"{a.prefix}_clusterability.csv", ["config", "ctt", "gdst"]),
        (f"{a.prefix}_nearest.csv", ["config", "nl_bal_acc", "nl_acc"]),
        (f"{a.prefix}_region_geometry_RANDOM.csv",
         ["config", "sil_k2_rand", "ctt_rand", "dip_dist_rand", "dip_axis_rand",
          "pct_lt_025_rand", "eff_dim_rand", "eff_n_distinct_rand"]),
    )
    for name, cols in extras:
        p = os.path.join(P, name)
        if not os.path.exists(p):
            print(f"  note: {name} absent; its columns will be missing")
            continue
        e = pd.read_csv(p)
        t = t.merge(e[[c for c in cols if c in e.columns]], on="config", how="left")

    # Derived columns, persisted so selection never recomputes them.
    # dip_ratio: dip on the k=2 separating axis relative to that config's own
    # matched-unimodal null. The raw dip_axis is optimistically biased -- kmeans
    # chooses the axis on the same data -- so the ratio is what to rank on.
    if {"dip_axis", "dip_axis_null"} <= set(t.columns):
        t["dip_ratio"] = t.dip_axis / t.dip_axis_null.replace(0, np.nan)
    # pct_excess: degeneracy ABOVE the dimension-matched random floor.
    if {"pct_lt_025", "pct_lt_025_rand"} <= set(t.columns):
        t["pct_excess"] = t.pct_lt_025 - t.pct_lt_025_rand
    if {"ctt", "ctt_rand"} <= set(t.columns):
        t["ctt_excess"] = t.ctt - t.ctt_rand

    order = [c for c in ["config", "d", "n_regions", "sil_k2", "ctt", "rcs",
                         "dip_dist", "dip_axis", "dip_axis_null", "dip_ratio",
                         "pct_lt_025", "pct_excess", "eff_n_distinct", "eff_dim",
                         "nl_bal_acc"] if c in t.columns]
    t = t[order + [c for c in t.columns if c not in order]]
    out = os.path.join(P, f"{a.prefix}_REGION_METRICS.csv")
    t.to_csv(out, index=False)
    print(f"  wrote {out}  ({len(t)} configs, {len(t.columns)} columns)")


def stage_summary(a, P, arm_dir, want):
    def rd(name):
        p = os.path.join(P, name)
        return pd.read_csv(p) if os.path.exists(p) else None

    metrics = rd(f"{a.prefix}_REGION_METRICS.csv")
    if metrics is None:
        raise SystemExit(f"missing {P}/{a.prefix}_REGION_METRICS.csv -- run merge first")
    nl_auc = rd(f"{a.prefix}_nl_auc.csv")
    rand = rd(f"{a.prefix}_region_geometry_RANDOM.csv")
    scored = rd(f"{a.prefix}_scored.csv")
    cluster = rd(f"{a.prefix}_clusterability.csv")

    summary = build_summary(metrics, nl_auc, rand, cluster)
    out = os.path.join(P, f"{a.prefix}_SUMMARY_TABLE.csv")
    summary.to_csv(out, index=False)
    n_cfg = (summary.row_type == "CONFIG").sum()
    n_ctl = (summary.row_type == "RANDOM_CONTROL").sum()
    print(f"  wrote {out}  ({n_cfg} configs + {n_ctl} random controls, "
          f"{len(summary.columns)} columns)")

    cfgtab = summary[summary.row_type == "CONFIG"].copy()
    if scored is not None:
        params = [c for c in ("dim", "epoch", "margin", "neg", "mc", "lr", "ada")
                  if c in scored.columns]
        cfgtab = cfgtab.merge(scored[["config"] + params].drop_duplicates("config"),
                              on="config", how="left")
        if nl_auc is not None:
            scored = scored.merge(
                nl_auc[["config"] + [c for c in MODEL_METRICS if c in nl_auc.columns]],
                on="config", how="left", suffixes=("", "_nl"))
    if nl_auc is not None and "nl_auc" in nl_auc.columns:
        cfgtab = cfgtab.drop(columns=["nl_auc"], errors="ignore").merge(
            nl_auc[["config", "nl_auc"]], on="config", how="left")

    eff, spread = main_effects(cfgtab, scored)
    if eff.empty:
        return
    p1 = os.path.join(P, f"{a.prefix}_MAIN_EFFECTS.csv")
    p2 = os.path.join(P, f"{a.prefix}_MAIN_EFFECT_SPREAD.csv")
    eff.to_csv(p1, index=False)
    spread.to_csv(p2, index=False)
    print(f"  wrote {p1}  ({len(eff)} axis x level x method rows)")
    print(f"  wrote {p2}")
    model = spread[spread.method.fillna("") == ""]
    if "spread_nl_auprc" in model.columns:
        print("\n  MODEL AXIS (nl_auprc) -- max minus min of the level medians:")
        for r in model.sort_values("spread_nl_auprc", ascending=False).itertuples():
            print(f"    {r.axis:<8} spread {r.spread_nl_auprc:.4f}  "
                  f"over {r.n_levels} levels")

    if a.interactions and scored is not None:
        cells = []
        for i, (c1, _) in enumerate(AXES):
            for c2, _ in AXES[i + 1:]:
                if c1 not in cfgtab.columns or c2 not in cfgtab.columns:
                    continue
                for (l1, l2), g in cfgtab.groupby([c1, c2]):
                    row = {"axis_1": c1, "level_1": l1, "axis_2": c2,
                           "level_2": l2, "n_configs": len(g)}
                    for m in MODEL_METRICS:
                        if m in g.columns:
                            row[m] = float(np.nanmedian(g[m]))
                    cells.append(row)
        if cells:
            p3 = os.path.join(P, f"{a.prefix}_TWOWAY_CELLS.csv")
            pd.DataFrame(cells).to_csv(p3, index=False)
            print(f"  wrote {p3}  ({len(cells)} two-way cells)")


def stage_figures(a, P, arm_dir, want):
    outdir = a.outdir or os.path.join(a.run, "figures_final")
    os.makedirs(outdir, exist_ok=True)
    ctt_min = a.ctt_min if a.ctt_min >= 0 else None
    keep = {s.strip() for s in a.figures.split(",") if s.strip()}
    print(f"  outdir : {outdir}")
    print(f"  gate   : ctt > {ctt_min} and pct_excess <= {a.pct_max}")

    arms, sweep_gate = {}, None
    for arm in ARMS:
        d = load_arm(a.run, arm, ctt_min, a.pct_max, sweep_gate)
        arms[arm] = d
        if arm == "sweep" and d is not None:
            cols = [c for c in ("config", "ctt", "pct_excess", "dip_ratio", "rcs",
                                "d", "passes", "passes_ctt", "passes_pct")
                    if c in d.columns]
            sweep_gate = d[cols].drop_duplicates("config")
        print(f"    {arm:<11} " + ("absent" if d is None
                                   else f"{d.config.nunique()} configs"))
    if all(v is None for v in arms.values()):
        raise SystemExit("nothing to plot: no arm has a scored table yet")

    cache = os.path.join(outdir, "score_distributions.csv")
    dist = pd.DataFrame()
    if keep & {"4", "5"}:
        print("\n  building the score-distribution table ...")
        dist = score_table(a.run, arms, cache, a.rebuild_distributions)

    if "1" in keep:
        fig1(a.run, arms, outdir, ctt_min, a.pct_max, a.top_n, a.auprc_lo, a.auprc_hi)
    if "2" in keep:
        fig2(a.run, arms, outdir, a.auprc_lo, a.auprc_hi)
    # --auprc-lo/hi pins fig2 ONLY. fig3/fig4 carry the external arms, whose
    # score AUPRC runs down to 0.22 (pseudobulk); pinning them to fig2's floor
    # would cut those rows off. Both still share one range within their figure.
    if "3" in keep:
        fig3(a.run, arms, outdir, None, None)
    if "4" in keep:
        fig4(a.run, arms, dist, outdir, None, None)
    if "5" in keep:
        fig5(a.run, arms, dist, outdir, ctt_min, a.pct_max)


def stage_external(a, P, arm_dir, want):
    """Pseudobulk + single-cell metrics on the Renin vs Non_Renin contrast.

    Reads <run>/external/_parts/<set>_part*_scored.csv (written by 7r) and the
    per-config score tables under <run>/sweep/<cfg>/<set>/, and writes
    <run>/sweep/external_perf/<set>_scored.csv -- the tables figs 2-5 read for
    the external rows.
    """
    parts_dir = os.path.join(a.run, "external", "_parts")
    out_dir = os.path.join(a.run, "sweep", "external_perf")
    os.makedirs(out_dir, exist_ok=True)

    for setkey in SETS:
        parts = sorted(glob.glob(os.path.join(parts_dir, f"{setkey}_part*_scored.csv")))
        if not parts:
            print(f"  {setkey}: no parts in {parts_dir}; skipping")
            continue
        base = pd.concat([pd.read_csv(f) for f in parts], ignore_index=True)
        base = base.drop_duplicates(subset=["config", "method"])
        cfgs = sorted(base.config.unique())
        print(f"  {setkey}: {len(parts)} parts -> {len(cfgs)} configs, {len(base)} rows")

        with mp.Pool(a.pool) as pool:
            rows = [r for chunk in pool.imap_unordered(
                        _ext_one, [(a.run, setkey, c) for c in cfgs]) for r in chunk]
        rv = pd.DataFrame(rows)
        if not len(rv):
            print(f"  {setkey}: no per-config score tables; writing the "
                  f"all-classes metrics unchanged")
            base.to_csv(os.path.join(out_dir, f"{setkey}_scored.csv"), index=False)
            continue

        keep_all = base[["config", "method"] + [c for c in ("auprc", "auc", "bal_acc")
                                                if c in base.columns]].rename(
            columns={"auprc": "auprc_all_classes", "auc": "auc_all_classes",
                     "bal_acc": "bal_acc_all_classes"})
        params = [c for c in ("dim", "epoch", "neg", "lr", "mc", "margin", "ada",
                              "in_sweep") if c in base.columns]
        merged = (rv.merge(base[["config"] + params].drop_duplicates("config"),
                           on="config", how="left")
                    .merge(keep_all, on=["config", "method"], how="left"))
        p_out = os.path.join(out_dir, f"{setkey}_scored.csv")
        merged.to_csv(p_out, index=False)
        n = merged.iloc[0]
        print(f"  wrote {p_out}: {merged.config.nunique()} configs x "
              f"{merged.method.nunique()} methods = {len(merged)} rows")
        print(f"    Renin vs Non_Renin: n={int(n.n_samples)} "
              f"({int(n.n_pos)} Renin, prevalence {n.prevalence:.3f}); "
              f"Recruited and anchor rows dropped")


STAGES = [("geom", stage_geom), ("rand", stage_rand), ("rcs", stage_rcs),
          ("nlauc", stage_nlauc), ("nlhard", stage_nlhard),
          ("merge", stage_merge), ("summary", stage_summary),
          ("external", stage_external), ("figures", stage_figures)]
STAGE_NAMES = [n for n, _ in STAGES]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=STAGE_NAMES + ["all"],
                    help="'all' runs every stage in order")
    ap.add_argument("--run", required=True, help="the run directory")
    ap.add_argument("--arm-dir", default="sweep",
                    help="arm under --run holding the per-config dirs (default: sweep)")
    ap.add_argument("--perf-dir", default=None,
                    help="default <run>/<arm-dir>/upsgrid_perf")
    ap.add_argument("--prefix", default="upsgrid")
    ap.add_argument("--configs-file", default=None,
                    help="config allow-list, one name per line; without it every "
                         "directory under the arm dir is globbed")
    ap.add_argument("--model-glob", default="model_*.tsv",
                    help="sweep arm: model_*.tsv (default); "
                         "holdout arm: starspace_trained_model.tsv")
    ap.add_argument("--train-input", default=None,
                    help="rcs: the corpus whose lines are the training files "
                         "(default <arm-dir>/preprocessed/train_input.txt). Must "
                         "list each training sample EXACTLY ONCE -- point it at "
                         "train_input_original.txt on a class-balanced arm.")
    ap.add_argument("--pool", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)))
    ap.add_argument("--rcs-pool", type=int, default=12,
                    help="rcs runs a 5-fold MLP per config; keep this modest")
    ap.add_argument("--n-seeds", type=int, default=3, help="rand: seeds per config")
    ap.add_argument("--labels", default=",".join(LABELS))
    ap.add_argument("--target", default=TARGET)
    ap.add_argument("--anchor", default="Control",
                    help="nlhard: label whose score row is carried as "
                         "anchor_score; '' to omit")
    ap.add_argument("--score-file",
                    default="reninness_score_weighted_multi_learned_post.csv",
                    help="nlhard: per-config score table the anchor is read from")
    ap.add_argument("--interactions", action="store_true",
                    help="summary: also write the two-way cell medians (3-6 "
                         "configs per cell; directional, not estimates)")
    ap.add_argument("--outdir", default=None, help="figures: default <run>/figures_final")
    ap.add_argument("--ctt-min", type=float, default=0.5,
                    help="untrained-embedding gate (default 0.5); negative disables")
    ap.add_argument("--pct-max", type=float, default=0.1,
                    help="degeneracy gate on pct_excess (default 0.1 = "
                         "eff_n_distinct >= 1000)")
    ap.add_argument("--top-n", type=int, default=5,
                    help="fig1: configs circled, by model AUPRC AND CTT")
    ap.add_argument("--auprc-lo", type=float, default=None,
                    help="pin the shared AUPRC axis low end (default: fit the data)")
    ap.add_argument("--auprc-hi", type=float, default=None)
    ap.add_argument("--figures", default="1,2,3,4,5",
                    help="figures: comma-separated subset, e.g. '1,2'")
    ap.add_argument("--rebuild-distributions", action="store_true",
                    help="figures: recompute the score-distribution table "
                         "instead of reading the cache")
    a = ap.parse_args()

    arm_dir = os.path.join(a.run, a.arm_dir)
    P = a.perf_dir or os.path.join(arm_dir, f"{a.prefix}_perf")
    os.makedirs(P, exist_ok=True)
    want = _want_configs(a.configs_file)

    global MODEL_GLOB
    MODEL_GLOB = a.model_glob

    todo = STAGE_NAMES if a.stage == "all" else [a.stage]
    print(f"run      : {a.run}\narm      : {arm_dir}\nperf-dir : {P}\n"
          f"stages   : {', '.join(todo)}\n")
    for name in todo:
        fn = dict(STAGES)[name]
        print(f"{'='*70}\n{name}\n{'='*70}")
        fn(a, P, arm_dir, want)
        print()
    return 0


if __name__ == "__main__":
    mp.set_start_method("forkserver")
    raise SystemExit(main())
