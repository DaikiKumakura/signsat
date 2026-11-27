"""MEK-inhibitor organoids (GSE329052): how many of the strongest TF changes can be explained at once, and do the
explanations route through MEK?

Uses the same data and TF scoring as gse329052_organoids.py (imported logic is repeated here to keep each script
self-contained). For k = 5..40 (top-k TFs by paired |t|, MEKi vs control, 4 patients):
  * real data: SAT / UNSAT per patient, minimal conflicting TF set for UNSAT;
  * null: the same k TFs with the signs of the same patient shuffled across TFs (20 shuffles per patient): SAT rate;
  * for SAT patients: number of observed TFs with a consistent explanation through MEK (MAP2K1/2 = -1) and through
    the wrong direction (+1), 60 s limit, edges=False.
Run: python -X utf8 validation/gse329052_k_sweep.py [comma-separated k list]
"""
import random
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

import signsat as ss
from signsat.graph import Graph
from signsat.sources import signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
SEED = 20261122
OUT = "validation/gse329052_ksweep_" + "_".join(sys.argv[1].split(",")) if len(sys.argv) > 1 else "validation/gse329052_ksweep"
KS = tuple(int(a) for a in sys.argv[1].split(",")) if len(sys.argv) > 1 else (5, 10, 15, 20, 30, 40)
QUERY_MAX_K = 10  # the MEK queries take up to 60 s each, so they are only run for small k
N_SHUF = 20
LIB = {"UPL25-007383": ("p1", "ctrl"), "UPL25-007384": ("p1", "MEKi"), "UPL25-007385": ("p1", "BCLxLi"),
       "UPL25-007386": ("p2", "ctrl"), "UPL25-007387": ("p2", "MEKi"), "UPL25-007388": ("p2", "BCLxLi"),
       "UPL25-007389": ("p3", "ctrl"), "UPL25-007390": ("p3", "MEKi"), "UPL25-007391": ("p3", "BCLxLi"),
       "UPL25-007392": ("p4", "ctrl"), "UPL25-007393": ("p4", "MEKi"), "UPL25-007394": ("p4", "BCLxLi")}

d = pd.read_csv("data/raw/gse329052/GSE329052_all.normalised_counts.tsv.gz", sep="\t")
d["gene_id"] = d.gene_id.str.split(".").str[0]
d = d[d.drop(columns="gene_id").sum(axis=1) >= 20]
m = pd.read_csv("data/raw/gse329052/ens2sym.csv")  # written by gse329052_organoids.py
x = d.merge(m, on="gene_id").drop(columns="gene_id").groupby("symbol").sum()
lg = np.log2(x + 1)
z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0).dropna()
ct = pd.read_csv("data/raw/collectri.tsv", sep="\t", dtype=str)
ct = ct[(ct.is_stimulation == "True") ^ (ct.is_inhibition == "True")]
ct["w"] = np.where(ct.is_stimulation == "True", 1.0, -1.0)
ct = ct.groupby(["source_genesymbol", "target_genesymbol"], as_index=False)["w"].mean()
sc = {}
for tf, g in ct.groupby("source_genesymbol"):
    g = g[g.target_genesymbol.isin(z.index)]
    if len(g) >= 8:
        sc[tf] = z.loc[g.target_genesymbol].mul(g.w.to_numpy(), axis=0).mean(axis=0)
S = pd.DataFrame(sc).T
pats = ("p1", "p2", "p3", "p4")
ctrl = {p: [c for c in LIB if LIB[c] == (p, "ctrl")][0] for p in pats}
mek = {p: [c for c in LIB if LIB[c] == (p, "MEKi")][0] for p in pats}
delta = pd.DataFrame({p: S[mek[p]] - S[ctrl[p]] for p in pats})
tt = stats.ttest_1samp(delta, 0, axis=1)
rank = pd.Series(np.abs(tt.statistic), index=delta.index).dropna().sort_values(ascending=False)

df, _ = signor_pkn("data/raw/signor_release.tsv")
P = ss.read_pkn(df)
G = Graph.from_pkn(P)
roots = set(G.structural_roots())
fwd = G.forward_reachable(sorted(roots))
usable = [g for g in rank.index if g in G.index and fwd[G.index[g]] and G.index[g] not in roots]
print("usable TFs (in SIGNOR, reachable, not roots):", len(usable), "; strongest:", usable[:10], flush=True)
rng = random.Random(SEED)
rows = []
for k in KS:
    tfs = usable[:k]
    inp = delta.loc[tfs].T  # patient x TF
    t0 = time.time()
    r = ss.infer(tf_scores=inp, pkn=df, solver="auto", classify=False, time_limit=120)
    st = r.status().to_dict()
    # shuffled-sign null: permute the signs of each patient across the k TFs
    shuf_sat = 0
    shuf_n = 0
    shuf_inc = 0
    sh = []
    for p in pats:
        for _ in range(N_SHUF):
            v = inp.loc[p].to_numpy().copy()
            rng.shuffle(v)
            sh.append(pd.Series(v, index=tfs, name=f"{p}_s{len(sh)}"))
    rs = ss.infer(tf_scores=pd.DataFrame(sh), pkn=df, solver="auto", classify=False, time_limit=120)
    sts = rs.status().value_counts().to_dict()
    row = {"k": k, "real": st, "shuffled_SAT": sts.get("SAT", 0), "shuffled_UNSAT": sts.get("UNSAT", 0),
           "shuffled_other": sum(v for a, v in sts.items() if a not in ("SAT", "UNSAT")), "shuffled_n": len(sh), "seconds": round(time.time() - t0)}
    for p in pats:
        n = r.network(p)
        if n.status == "UNSAT":
            row[f"{p}_conflict"] = n.diagnostics[n.diagnostics.code == "UNSAT_CONFLICT"].detail.iloc[0].replace("minimal conflicting TF observations: ", "")
        elif n.status == "SAT" and k <= QUERY_MAX_K:
            q = {lab: r.subnetwork(p, target=tgt, time_limit=60, edges=False) for lab, tgt in
                 (("minus", {"MAP2K1": -1, "MAP2K2": -1}), ("plus", {"MAP2K1": 1, "MAP2K2": 1}))}
            row[f"{p}_MEK-1_explains"] = f"{len(q['minus'].explained_tfs)}/{k}" + ("" if q["minus"].complete else " (incomplete)")
            row[f"{p}_MEK+1_explains"] = f"{len(q['plus'].explained_tfs)}/{k}" + ("" if q["plus"].complete else " (incomplete)")
    rows.append(row)
    print(row, flush=True)
    pd.DataFrame(rows).to_csv(OUT + ".csv", index=False)
