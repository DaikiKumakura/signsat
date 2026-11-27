"""Mechanism-hypothesis test on GSE329052 (MEK-inhibitor organoids): does the known drug action explain more of the
observed TF changes than random hypotheses of the same reach?

Hypothesis h = a set of intervened root nodes with imposed states (the drug action). For one patient and one observed TF f
(with observed state y_f), h explains f if SignSAT finds a simple path from a root of h to f whose signs give y_f
(roots = the nodes of h, states fixed: `SampleModel(..., fixed_states=...)`, single-TF satisfiability). n_explained(h, patient)
= number of the k observed TFs explained one by one.
 * MEKi data: h = MAP2K1 = MAP2K2 = -1 (trametinib) versus the wrong direction (+1) and versus RANDOM_PAIRS random node pairs
   (both set to -1) whose number of reachable usable TFs is within a factor 2 of MEK's.
 * BCLxLi data (negative control): the same MEK hypothesis should not rank high.
 * Two TF-activity variants: 'z' (z-scored expression, as in the earlier scripts) and 'lfc' (paired log2 changes, no z-score).
Empirical p = share of random hypotheses with n_explained >= that of MEK(-1), per patient (and summed over patients).
Run: python -X utf8 validation/gse329052_hypothesis_test.py [k]
"""
from collections import deque
import random
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

import signsat as ss
from signsat.encoding import SampleModel
from signsat.graph import Graph
from signsat.sources import signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
K = int(sys.argv[1]) if len(sys.argv) > 1 else 15
RANDOM_PAIRS = 2000
SEED = 20261127
OUT = "validation/gse329052_hypothesis"
LIB = {"UPL25-007383": ("p1", "ctrl"), "UPL25-007384": ("p1", "MEKi"), "UPL25-007385": ("p1", "BCLxLi"),
       "UPL25-007386": ("p2", "ctrl"), "UPL25-007387": ("p2", "MEKi"), "UPL25-007388": ("p2", "BCLxLi"),
       "UPL25-007389": ("p3", "ctrl"), "UPL25-007390": ("p3", "MEKi"), "UPL25-007391": ("p3", "BCLxLi"),
       "UPL25-007392": ("p4", "ctrl"), "UPL25-007393": ("p4", "MEKi"), "UPL25-007394": ("p4", "BCLxLi")}
PATS = ("p1", "p2", "p3", "p4")

d = pd.read_csv("data/raw/gse329052/GSE329052_all.normalised_counts.tsv.gz", sep="\t")
d["gene_id"] = d.gene_id.str.split(".").str[0]
d = d[d.drop(columns="gene_id").sum(axis=1) >= 20]
m = pd.read_csv("data/raw/gse329052/ens2sym.csv")
lg = np.log2(d.merge(m, on="gene_id").drop(columns="gene_id").groupby("symbol").sum() + 1)
z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0).dropna()
ct = pd.read_csv("data/raw/collectri.tsv", sep="\t", dtype=str)
ct = ct[(ct.is_stimulation == "True") ^ (ct.is_inhibition == "True")]
ct["w"] = np.where(ct.is_stimulation == "True", 1.0, -1.0)
ct = ct.groupby(["source_genesymbol", "target_genesymbol"], as_index=False)["w"].mean()
col = {(p, c): k for k, (p, c) in LIB.items()}


def tf_delta(variant: str, drug: str) -> pd.DataFrame:
    """TF x patient change (drug minus the patient's control)."""
    rows = {}
    for tf, g in ct.groupby("source_genesymbol"):
        if variant == "z":
            g = g[g.target_genesymbol.isin(z.index)]
            if len(g) < 8:
                continue
            sc = z.loc[g.target_genesymbol].mul(g.w.to_numpy(), axis=0).mean(axis=0)
            rows[tf] = {p: sc[col[(p, drug)]] - sc[col[(p, "ctrl")]] for p in PATS}
        else:
            g = g[g.target_genesymbol.isin(lg.index)]
            if len(g) < 8:
                continue
            rows[tf] = {p: float((lg.loc[g.target_genesymbol, col[(p, drug)]] - lg.loc[g.target_genesymbol, col[(p, "ctrl")]]).mul(g.w.to_numpy()).mean())
                        for p in PATS}
    return pd.DataFrame(rows).T[list(PATS)]


df, _ = signor_pkn("data/raw/signor_release.tsv")
P = ss.read_pkn(df)
G = Graph.from_pkn(P)
roots_all = set(G.structural_roots())
fwd = G.forward_reachable(sorted(roots_all))
usable_set = {g for g in G.nodes if fwd[G.index[g]] and G.index[g] not in roots_all}
rng = random.Random(SEED)


def reach_states(root_names, sts):
    """All (node index, state) pairs reachable from the intervened roots when every edge maps state z to sign * z.
    For a single observed TF this is exactly 'a path with the right sign exists' as a walk; it agreed with the SAT
    check (one simple path, roots fixed) on 120 of 120 (TF, state, hypothesis) cases (see docs/validation.md), so it is
    used for the large null; the SAT version stays in signsat.infer / SampleModel."""
    seen, dq = set(), deque()
    for r, st in zip(root_names, sts):
        x0 = (G.index[r], st)
        seen.add(x0)
        dq.append(x0)
    while dq:
        u, zu = dq.popleft()
        for k in G.out_edges[u]:
            v = (int(G.dst[k]), int(G.sign[k]) * zu)
            if v not in seen:
                seen.add(v)
                dq.append(v)
    return seen


def count(seen, obs):
    return sum((G.index[tf], y) in seen for tf, y in obs.items())


def reach_size(node: str, usable: list[str]) -> int:
    r = G.forward_reachable([G.index[node]])
    return int(sum(r[G.index[t]] for t in usable))


# ---- 1. observed TF sets for every (variant, drug)
sets = {}
for variant in ("z", "lfc"):
    for drug in ("MEKi", "BCLxLi"):
        delta = tf_delta(variant, drug)
        tt = stats.ttest_1samp(delta, 0, axis=1)
        rank = pd.Series(np.abs(tt.statistic), index=delta.index).dropna().sort_values(ascending=False)
        top = [g for g in rank.index if g in usable_set][:K]
        sets[(variant, drug)] = (delta, top)
all_tfs = sorted({t for _, top in sets.values() for t in top})
print("distinct observed TFs over all conditions:", len(all_tfs), flush=True)
# ---- 2. hypotheses: MEK (both directions), reach-matched random pairs, every single node
mek_reach = np.mean([reach_size("MAP2K1", all_tfs), reach_size("MAP2K2", all_tfs)])
pool = [g for g in G.nodes if g not in all_tfs and g not in ("MAP2K1", "MAP2K2") and len(G.out_edges[G.index[g]]) >= 2]
reach_cache = {g: reach_size(g, all_tfs) for g in pool}
pairs, tries = [], 0
while len(pairs) < RANDOM_PAIRS and tries < 200000:
    tries += 1
    a, b = rng.sample(pool, 2)
    if 0.5 * mek_reach <= (reach_cache[a] + reach_cache[b]) / 2 <= 2 * mek_reach:
        pairs.append((a, b))
print(f"MEK reach over those TFs {mek_reach:.1f}; random pairs {len(pairs)}", flush=True)
mek_m, mek_p = reach_states(["MAP2K1", "MAP2K2"], (-1, -1)), reach_states(["MAP2K1", "MAP2K2"], (1, 1))
null_seen = [reach_states([a, b], (-1, -1)) for a, b in pairs]
single = {g: reach_states([g], (-1,)) for g in G.nodes}
rows = []
for (variant, drug), (delta, top) in sets.items():
    tot = {"MEK-1": 0, "MEK+1": 0}
    nulls_sum = np.zeros(len(pairs))
    for p in PATS:
        obs = {g: int(np.sign(delta.loc[g, p])) for g in top if delta.loc[g, p] != 0}
        n_minus, n_plus = count(mek_m, obs), count(mek_p, obs)
        null = np.array([count(t, obs) for t in null_seen])
        nulls_sum += null
        tot["MEK-1"] += n_minus
        tot["MEK+1"] += n_plus
        sc = {g: count(sn, obs) for g, sn in single.items()}
        vals = np.array(list(sc.values()))
        rank = lambda g: float((vals > sc[g]).sum() + 0.5 * ((vals == sc[g]).sum() - 1) + 1)
        rows.append({"variant": variant, "drug": drug, "patient": p, "observed_TFs": len(obs), "MEK-1": n_minus, "MEK+1": n_plus,
                     "null_median": float(np.median(null)), "null_max": int(null.max()), "p_ge_MEK-1": float((null >= n_minus).mean()),
                     "p_ge_MEK+1": float((null >= n_plus).mean()), "null_n": len(null),
                     "single_best": int(vals.max()), "rank_MAP2K1(-1)": rank("MAP2K1"), "rank_MAP2K2(-1)": rank("MAP2K2"),
                     "rank_BCL2L1(-1)": rank("BCL2L1") if "BCL2L1" in sc else None, "n_single_hypotheses": len(vals)})
        print(rows[-1], flush=True)
    rows.append({"variant": variant, "drug": drug, "patient": "all 4 (summed)", "MEK-1": tot["MEK-1"], "MEK+1": tot["MEK+1"],
                 "null_median": float(np.median(nulls_sum)), "null_max": int(nulls_sum.max()),
                 "p_ge_MEK-1": float((nulls_sum >= tot["MEK-1"]).mean()), "p_ge_MEK+1": float((nulls_sum >= tot["MEK+1"]).mean()), "null_n": len(pairs)})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}_k{K}.csv", index=False)
