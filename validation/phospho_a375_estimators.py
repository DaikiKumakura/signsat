"""Second estimators following docs/phospho_protocol_2.md (PXD013923, A375, 30 min, Table S1): KSEA z (primary), kinase-filtered mean, rank-sum.

Run: python -X utf8 validation/phospho_a375_estimators.py
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
SEED = 20261129
K = 15
N_NULL = 2000
PRIMARY_L = 3
OUT = "validation/phospho_a375_estimators"
DRUGS = {"Dabrafenib": (["BRAF"], (-1,)), "Trametinib": (["MAP2K1", "MAP2K2"], (-1, -1)), "SCH772984": (["MAPK1", "MAPK3"], (-1, -1))}
OWN = {"Dabrafenib": {"BRAF"}, "Trametinib": {"MAP2K1", "MAP2K2"}, "SCH772984": {"MAPK1", "MAPK3"}}

# ---------------------------------------------------------------- data (protocol 1)
raw = pd.read_csv("data/raw/phos_a375/tableS1.csv", header=None, low_memory=False)  # converted once from the xlsx (MOESM2) with pandas+openpyxl
body = raw.iloc[2:].reset_index(drop=True)
hdr = {"AA": 1, "Gene": 4, "Pos": 6}
drug_cols = {"Dabrafenib": [7, 8, 9, 10], "Trametinib": [12, 13, 14, 15], "SCH772984": [17, 18, 19, 20]}
sites = pd.DataFrame({"gene": body[4].astype(str).str.split().str[0], "aa": body[1], "pos": pd.to_numeric(body[6], errors="coerce")})
for d, cols in drug_cols.items():
    for r, c in enumerate(cols[:3], 1):
        sites[f"{d}_r{r}"] = pd.to_numeric(body[c], errors="coerce").replace(0, np.nan)
    sites[f"{d}_avg"] = pd.to_numeric(body[cols[3]], errors="coerce").replace(0, np.nan)
print("sites in table:", len(sites), flush=True)

# ---------------------------------------------------------------- sanity gate (protocol 2)
act = [("MAPK1", "T", 185), ("MAPK1", "Y", 187), ("MAPK3", "T", 202), ("MAPK3", "Y", 204)]
gate_ok = True
for d in DRUGS:
    vals = []
    for g, a, p in act:
        m = sites[(sites.gene == g) & (sites.aa == a) & (sites.pos == p)][f"{d}_avg"].dropna()
        if len(m):
            vals.append(float(m.iloc[0]))
    print(f"gate {d}: n_sites={len(vals)} values={vals}", flush=True)
    if len(vals) < 3 or np.median(vals) >= 0:
        gate_ok = False
if not gate_ok:
    print("GATE FAILED: data unusable per protocol; stopping.")
    sys.exit(0)

# ---------------------------------------------------------------- kinase scores (protocol 3)
es = pd.read_csv("data/raw/omnipath_enzsub.tsv", sep="\t")
es = es[es.modification == "phosphorylation"][["enzyme_genesymbol", "substrate_genesymbol", "residue_type", "residue_offset"]].drop_duplicates()
mk = sites.merge(es, left_on=["gene", "aa", "pos"], right_on=["substrate_genesymbol", "residue_type", "residue_offset"], how="inner")
print(f"table sites matched to a kinase-substrate pair: {sites[['gene','aa','pos']].drop_duplicates().merge(es, left_on=['gene','aa','pos'], right_on=['substrate_genesymbol','residue_type','residue_offset']).drop_duplicates(['gene','aa','pos']).shape[0]} of {len(sites)} sites; pairs {len(mk)}", flush=True)


KIN = set(pd.read_csv("data/raw/kinasecom.tsv", sep="	").genesymbol)
mk = mk[mk.enzyme_genesymbol.isin(KIN)]
print("kinase-substrate pairs after kinase.com filter:", len(mk), "kinases:", mk.enzyme_genesymbol.nunique(), flush=True)


def kinase_scores(d, method):
    out = {}
    allv = {r: sites[f"{d}_r{r}"].dropna().to_numpy() for r in (1, 2, 3)}
    for k, g in mk.groupby("enzyme_genesymbol"):
        row = {}
        for r in (1, 2, 3):
            v = g[f"{d}_r{r}"].dropna().to_numpy()
            m = len(v)
            if m < 3:
                row[r] = np.nan
            elif method == "C":
                row[r] = float(v.mean())
            elif method == "B":
                row[r] = float((v.mean() - allv[r].mean()) * np.sqrt(m) / allv[r].std(ddof=1))
            else:
                row[r] = float(stats.mannwhitneyu(v, allv[r], alternative="two-sided").statistic)
                n1, n2 = m, len(allv[r])
                row[r] = (row[r] - n1 * n2 / 2) / np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        if not np.isnan(list(row.values())).any():
            out[k] = row
    return pd.DataFrame(out).T[[1, 2, 3]]


# ---------------------------------------------------------------- PKN and helpers
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
        dd = seen[(u, zu)]
        if L is not None and dd >= L:
            continue
        for k in G.out_edges[u]:
            v = (int(G.dst[k]), int(G.sign[k]) * zu)
            if v not in seen:
                seen[v] = dd + 1
                dq.append(v)
    return set(seen)


def cnt(seen, obs):
    return sum((G.index[t], y) in seen for t, y in obs.items())


def reach_n(node, tfs):
    r = G.forward_reachable([G.index[node]])
    return int(sum(r[G.index[t]] for t in tfs))


rows = []
for method, drug in [(m_, d_) for m_ in ("B", "C", "D") for d_ in DRUGS]:
    names, sts = DRUGS[drug]
    sc = kinase_scores(drug, method)
    tt = stats.ttest_1samp(sc, 0, axis=1)
    rank = pd.Series(np.abs(tt.statistic), index=sc.index).dropna().sort_values(ascending=False)
    top = [g for g in rank.index if g in usable and g not in OWN[drug]][:K]
    obs_u = {r: {g: int(np.sign(sc.loc[g, r])) for g in top if sc.loc[g, r] != 0} for r in (1, 2, 3)}
    h_reach = float(np.mean([reach_n(n, top) for n in names]))
    pool = [g for g in G.nodes if g not in top and g not in OWN[drug] and len(G.out_edges[G.index[g]]) >= 2]
    rc = {g: reach_n(g, top) for g in pool}
    if len(names) == 1:
        nulls = [(rng.choice([g for g in pool if 0.5 * h_reach <= rc[g] <= 2 * h_reach]),) for _ in range(N_NULL)]
    else:
        nulls, tries = [], 0
        while len(nulls) < N_NULL and tries < 400000:
            tries += 1
            a, b = rng.sample(pool, 2)
            if 0.5 * h_reach <= (rc[a] + rc[b]) / 2 <= 2 * h_reach:
                nulls.append((a, b))
    print(f"\n== {drug}: kinases defined {len(sc)}, observed {len(top)} {top}; hypothesis reach {h_reach:.1f}; null {len(nulls)}", flush=True)
    sign_dir = {"-1": -1, "+1": 1}
    for L in (PRIMARY_L, 2, 4, None):
        hm = reach_depth(names, sts, L)
        hp = reach_depth(names, tuple(-s for s in sts), L)
        nul = {n: reach_depth(list(n), tuple(-1 for _ in n), L) for n in set(nulls)}
        single = {g: reach_depth([g], (-1,), L) for g in G.nodes}
        tot_m = sum(cnt(hm, obs_u[r]) for r in obs_u)
        tot_p = sum(cnt(hp, obs_u[r]) for r in obs_u)
        nl = np.array([sum(cnt(nul[n], obs_u[r]) for r in obs_u) for n in nulls])
        vals = np.array([sum(cnt(s_, obs_u[r]) for r in obs_u) for s_ in single.values()])
        rows.append({"method": method, "drug": drug, "L": "inf" if L is None else L, "observed_total": sum(len(v) for v in obs_u.values()),
                     "H(-1)": tot_m, "H(+1)": tot_p, "null_median": float(np.median(nl)), "null_max": int(nl.max()), "null_n": len(nl),
                     "p_ge_H(-1)": float((nl >= tot_m).mean()), "p_ge_H(+1)": float((nl >= tot_p).mean()),
                     "rank_best_single_node": int((vals > tot_m).sum() + 1) if len(names) == 1 else None})
        print(rows[-1], flush=True)
res = pd.DataFrame(rows)
res.to_csv(OUT + "_results.csv", index=False)
prim = res[(res.method == "B") & (res.drug == "Trametinib") & (res.L == PRIMARY_L)].iloc[0]
ok = prim["p_ge_H(-1)"] < 0.05 and prim["p_ge_H(+1)"] >= 0.05
print("\nPRIMARY (protocol 2, estimator B, trametinib L=3): p(MEK-1) =", prim["p_ge_H(-1)"], "p(MEK+1) =", prim["p_ge_H(+1)"], "->", "criterion MET" if ok else "criterion NOT met")
for m_ in ("B", "C", "D"):
  for d in DRUGS:
    r = res[(res.method == m_) & (res.drug == d) & (res.L == PRIMARY_L)].iloc[0]
    print(m_, d, "L=3 met:", bool(r["p_ge_H(-1)"] < 0.05 and r["p_ge_H(+1)"] >= 0.05))
