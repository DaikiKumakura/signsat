"""Replication on GSE325686 (A375, DMSO vs vemurafenib), following docs/replication_protocol.md exactly.

Run: python -X utf8 validation/replication_a375_vemurafenib.py
"""
import random
import sys
from collections import deque

import numpy as np
import pandas as pd
from scipy import stats

import signsat as ss
from signsat.graph import Graph
from signsat.sources import signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
SEED = 20261128
K = 15
N_NULL = 2000
PRIMARY_L = 3
OUT = "validation/replication_a375"

# ---------------------------------------------------------------------------- data (protocol 1)
c = pd.read_csv("data/raw/gse325686/readcount.xls.gz", sep="\t", low_memory=False)
c = c[~c.Symbol.astype(str).str.startswith("MSTRG")].drop(columns="GeneID").groupby("Symbol").sum()
c = c[c.sum(axis=1) >= 20]
cpm = c / c.sum(axis=0) * 1e6
lg = np.log2(cpm + 1)
z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0).dropna()
UNITS = {"WT-r1": ("WT-IBRAF-1", "WT-1"), "WT-r2": ("WT-IBRAF-2", "WT-2"), "KO-r1": ("KO-IBRAF-1", "KO-1"), "KO-r2": ("KO-IBRAF-2", "KO-2")}

# ---------------------------------------------------------------------------- sanity gate (protocol 2)
chk = {g: round(float(np.mean([lg.loc[g, a] - lg.loc[g, b] for a, b in UNITS.values()])), 2)
       for g in ("DUSP6", "SPRY4", "ETV4", "ETV5", "DUSP4") if g in lg.index}
print("sanity gate, mean paired log2 change under vemurafenib:", chk, flush=True)
if sum(v < 0 for v in chk.values()) < 4:
    print("GATE FAILED: dataset unusable per protocol; stopping.")
    sys.exit(0)

# ---------------------------------------------------------------------------- TF scores (protocol 3)
ct = pd.read_csv("data/raw/collectri.tsv", sep="\t", dtype=str)
ct = ct[(ct.is_stimulation == "True") ^ (ct.is_inhibition == "True")]
ct["w"] = np.where(ct.is_stimulation == "True", 1.0, -1.0)
ct = ct.groupby(["source_genesymbol", "target_genesymbol"], as_index=False)["w"].mean()


def tf_delta(variant):
    rows = {}
    for tf, g in ct.groupby("source_genesymbol"):
        if variant == "z":
            g = g[g.target_genesymbol.isin(z.index)]
            if len(g) < 8:
                continue
            sc = z.loc[g.target_genesymbol].mul(g.w.to_numpy(), axis=0).mean(axis=0)
            rows[tf] = {u: sc[a] - sc[b] for u, (a, b) in UNITS.items()}
        else:
            g = g[g.target_genesymbol.isin(lg.index)]
            if len(g) < 8:
                continue
            rows[tf] = {u: float((lg.loc[g.target_genesymbol, a] - lg.loc[g.target_genesymbol, b]).mul(g.w.to_numpy()).mean())
                        for u, (a, b) in UNITS.items()}
    return pd.DataFrame(rows).T[list(UNITS)]


df, _ = signor_pkn("data/raw/signor_release.tsv")
P = ss.read_pkn(df)
G = Graph.from_pkn(P)
roots_all = set(G.structural_roots())
fwd = G.forward_reachable(sorted(roots_all))
usable = {g for g in G.nodes if fwd[G.index[g]] and G.index[g] not in roots_all}
rng = random.Random(SEED)


def reach_depth(names, sts, L):
    seen = {}
    dq = deque()
    for r, st in zip(names, sts):
        seen[(G.index[r], st)] = 0
        dq.append((G.index[r], st))
    while dq:
        u, zu = dq.popleft()
        d = seen[(u, zu)]
        if L is not None and d >= L:
            continue
        for k in G.out_edges[u]:
            v = (int(G.dst[k]), int(G.sign[k]) * zu)
            if v not in seen:
                seen[v] = d + 1
                dq.append(v)
    return set(seen)


def cnt(seen, obs):
    return sum((G.index[t], y) in seen for t, y in obs.items())


def unlimited_reach(node, tfs):
    r = G.forward_reachable([G.index[node]])
    return int(sum(r[G.index[t]] for t in tfs))


rows = []
for variant in ("lfc", "z"):
    delta = tf_delta(variant)
    tt = stats.ttest_1samp(delta, 0, axis=1)
    rank = pd.Series(np.abs(tt.statistic), index=delta.index).dropna().sort_values(ascending=False)
    top = [g for g in rank.index if g in usable][:K]
    obs_u = {u: {g: int(np.sign(delta.loc[g, u])) for g in top if delta.loc[g, u] != 0} for u in UNITS}
    braf_reach = unlimited_reach("BRAF", top)
    pool = [g for g in G.nodes if g not in top and g != "BRAF" and len(G.out_edges[G.index[g]]) >= 2]
    rc = {g: unlimited_reach(g, top) for g in pool}
    matched = [g for g in pool if 0.5 * braf_reach <= rc[g] <= 2 * braf_reach]
    nodes_null = [rng.choice(matched) for _ in range(N_NULL)]
    print(f"\n== variant {variant}: top-{K} TFs {top[:10]}...; BRAF reach {braf_reach}; matched pool {len(matched)}", flush=True)
    for L in (PRIMARY_L, 2, 4):
        b_m, b_p = reach_depth(["BRAF"], (-1,), L), reach_depth(["BRAF"], (1,), L)
        m_m = reach_depth(["MAP2K1", "MAP2K2"], (-1, -1), L)
        null = {n: reach_depth([n], (-1,), L) for n in set(nodes_null)}
        single = {g: reach_depth([g], (-1,), L) for g in G.nodes}
        for scope, units in (("pooled (4 units)", list(UNITS)), ("WT", ["WT-r1", "WT-r2"]), ("KO", ["KO-r1", "KO-r2"])):
            tot_m = sum(cnt(b_m, obs_u[u]) for u in units)
            tot_p = sum(cnt(b_p, obs_u[u]) for u in units)
            tot_mek = sum(cnt(m_m, obs_u[u]) for u in units)
            nl = np.array([sum(cnt(null[n], obs_u[u]) for u in units) for n in nodes_null])
            vals = np.array([sum(cnt(s_, obs_u[u]) for u in units) for s_ in single.values()])
            rows.append({"variant": variant, "L": L, "scope": scope, "observed_total": sum(len(obs_u[u]) for u in units),
                         "BRAF-1": tot_m, "BRAF+1": tot_p, "MEK-1": tot_mek, "null_median": float(np.median(nl)), "null_max": int(nl.max()),
                         "p_ge_BRAF-1": float((nl >= tot_m).mean()), "p_ge_BRAF+1": float((nl >= tot_p).mean()),
                         "rank_BRAF-1_among_single_nodes": int((vals > tot_m).sum() + 1)})
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT + "_results.csv", index=False)
prim = [r for r in rows if r["variant"] == "lfc" and r["L"] == PRIMARY_L and r["scope"].startswith("pooled")][0]
ok = prim["p_ge_BRAF-1"] < 0.05 and prim["p_ge_BRAF+1"] >= 0.05
print("\nPRIMARY (protocol 6): p(BRAF-1) =", prim["p_ge_BRAF-1"], "p(BRAF+1) =", prim["p_ge_BRAF+1"], "->", "criterion MET" if ok else "criterion NOT met")
