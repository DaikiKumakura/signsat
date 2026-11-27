"""Real-data check on paired patient-derived colorectal organoids: do SignSAT explanations of MEK-inhibitor-induced
TF changes route through MEK (MAP2K1/MAP2K2), and not for another drug?

Data: GEO GSE329052 (4 CRC patients; control, trametinib = MEK inhibitor, A1155463 = BCL-xL inhibitor; normalised
counts). The BCL-xL inhibitor is a different-mechanism drug used as a negative control for the MEK query.

Pipeline (simple on purpose):
 1. log2(normalised count + 1); per-gene z-score across the 12 samples.
 2. TF score per sample = mean over CollecTRI targets of (sign * z) (>= 8 targets); change of a drug sample = its score
    minus the control of the same patient.
 3. Per drug: observed TFs = paired-t FDR < 0.1 (Benjamini-Hochberg), in the SIGNOR PKN and reachable, at most TOP_TF by
    |t|; state per patient = sign of that patient's change. SignSAT on SIGNOR (human proteins).
 4. Query target={MAP2K1:-1, MAP2K2:-1}: number of observed TFs with a consistent explanation that passes through MEK;
    compared with the wrong direction (+1) and with random target pairs (out-degree 3-80 nodes that structurally
    reach an observed TF).
Run: python -X utf8 validation/gse329052_organoids.py
"""
import json
import random
import sys
import time
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

import signsat as ss
from signsat.graph import Graph
from signsat.sources import signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
SEED = 20261122
TOP_TF = 40
N_RANDOM = 60
Q_TIME = 20.0
OUT = "validation/gse329052_results"
# library name -> (patient, condition), from the GSM titles / descriptions
LIB = {"UPL25-007383": ("p1", "ctrl"), "UPL25-007384": ("p1", "MEKi"), "UPL25-007385": ("p1", "BCLxLi"),
       "UPL25-007386": ("p2", "ctrl"), "UPL25-007387": ("p2", "MEKi"), "UPL25-007388": ("p2", "BCLxLi"),
       "UPL25-007389": ("p3", "ctrl"), "UPL25-007390": ("p3", "MEKi"), "UPL25-007391": ("p3", "BCLxLi"),
       "UPL25-007392": ("p4", "ctrl"), "UPL25-007393": ("p4", "MEKi"), "UPL25-007394": ("p4", "BCLxLi")}

# ------------------------------------------------------------------------------------------ expression
d = pd.read_csv("data/raw/gse329052/GSE329052_all.normalised_counts.tsv.gz", sep="\t")
d["gene_id"] = d.gene_id.str.split(".").str[0]
d = d[d.drop(columns="gene_id").sum(axis=1) >= 20]
try:
    m = pd.read_csv("data/raw/gse329052/ens2sym.csv")
except FileNotFoundError:
    out = {}
    ids = list(d.gene_id)
    for i in range(0, len(ids), 1000):
        data = urllib.parse.urlencode({"q": ",".join(ids[i:i + 1000]), "scopes": "ensembl.gene", "fields": "symbol", "species": "human"}).encode()
        for h in json.load(urllib.request.urlopen(urllib.request.Request("https://mygene.info/v3/query", data=data), timeout=60)):
            if "symbol" in h:
                out[h["query"]] = h["symbol"]
        time.sleep(0.3)
    m = pd.DataFrame({"gene_id": list(out), "symbol": list(out.values())})
    m.to_csv("data/raw/gse329052/ens2sym.csv", index=False)
x = d.merge(m, on="gene_id").drop(columns="gene_id").groupby("symbol").sum()
lg = np.log2(x + 1)
ctrl = {p: [c for c in LIB if LIB[c] == (p, "ctrl")][0] for p in ("p1", "p2", "p3", "p4")}
chk = {g: round(float(np.mean([lg.loc[g, [c for c in LIB if LIB[c] == (p, "MEKi")][0]] - lg.loc[g, ctrl[p]] for p in ctrl])), 2)
       for g in ("DUSP6", "SPRY4", "ETV4", "ETV5", "DUSP4", "FOSL1") if g in lg.index}
print("MEKi - control paired change, log2 (MAPK output genes, expected negative):", chk, flush=True)
z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1).replace(0, np.nan), axis=0).dropna()

# ------------------------------------------------------------------------------------------ TF scores
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
print("TFs scored:", len(S), flush=True)

df, _ = signor_pkn("data/raw/signor_release.tsv")
P = ss.read_pkn(df)
G = Graph.from_pkn(P)
roots = set(G.structural_roots())
fwd = G.forward_reachable(sorted(roots))
rng = random.Random(SEED)


def bh(p):
    p = np.asarray(p)
    o = np.argsort(p)
    n = len(p)
    q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    return q


results = []
for drug in ("MEKi", "BCLxLi"):
    cols = {p: [c for c in LIB if LIB[c] == (p, drug)][0] for p in ctrl}
    delta = pd.DataFrame({p: S[cols[p]] - S[ctrl[p]] for p in ctrl})  # TF x patient
    tt = stats.ttest_1samp(delta, 0, axis=1)
    res = pd.DataFrame({"t": tt.statistic, "p": tt.pvalue, "mean_delta": delta.mean(axis=1)}, index=S.index).dropna()
    res["fdr"] = bh(res.p.values)
    cand = res[(res.fdr < 0.1) & res.index.isin(G.nodes)]
    cand = cand[[bool(fwd[G.index[g]]) and G.index[g] not in roots for g in cand.index]]
    cand = cand.reindex(cand.t.abs().sort_values(ascending=False).index).head(TOP_TF)
    print(f"\n== {drug}: observed TFs {len(cand)} (paired FDR<0.1); top {list(cand.index[:12])}", flush=True)
    if len(cand) < 3:
        results.append({"drug": drug, "observed_TFs": len(cand)})
        continue
    inp = delta.loc[cand.index].T  # patient x TF
    t0 = time.time()
    r = ss.infer(tf_scores=inp, pkn=df, solver="auto", classify=False, time_limit=300)
    print("infer", round(time.time() - t0, 1), "s; status:", r.status().to_dict(), flush=True)
    cand.to_csv(f"{OUT}_{drug}_observed_tfs.csv")
    for sample in inp.index:
        if r.network(sample).status != "SAT":
            continue
        nobs = len(r.network(sample).obs)
        for label, tgt in (("MEK -1", {"MAP2K1": -1, "MAP2K2": -1}), ("MEK +1", {"MAP2K1": 1, "MAP2K2": 1})):
            s = r.subnetwork(sample, target=tgt, time_limit=Q_TIME, edges=False)
            results.append({"drug": drug, "patient": sample, "query": label, "observed_TFs": nobs,
                            "explained_TFs": len(s.explained_tfs), "reason": s.empty_reason, "complete": s.complete})
            print(results[-1], flush=True)
    first = next((s for s in inp.index if r.network(s).status == "SAT"), None)
    if first is None:
        continue
    obs_idx = {G.index[t_] for t_ in r.network(first).obs}
    pool = [v for v in range(G.n) if 3 <= len(G.out_edges[v]) <= 80 and v not in obs_idx and G.nodes[v] not in ("MAP2K1", "MAP2K2")
            and any(G.forward_reachable([v])[f] for f in list(obs_idx)[:5])]
    rows = []
    t0 = time.time()
    for i in range(N_RANDOM):
        a, b = rng.sample(pool, 2)
        st = rng.choice((1, -1))
        s = r.subnetwork(first, target={G.nodes[a]: st, G.nodes[b]: st}, time_limit=Q_TIME, edges=False)
        rows.append({"targets": f"{G.nodes[a]},{G.nodes[b]}", "state": st, "explained_TFs": len(s.explained_tfs), "complete": s.complete})
    nd = pd.DataFrame(rows)
    nd.to_csv(f"{OUT}_{drug}_null.csv", index=False)
    done = nd[nd.complete]
    mek = next(x_["explained_TFs"] for x_ in results if x_.get("drug") == drug and x_.get("patient") == first and x_.get("query") == "MEK -1")
    print(f"{drug}/{first}: MEK -1 explains {mek} of {len(obs_idx)}; random pairs (complete n={len(done)}): median "
          f"{done.explained_TFs.median()}, max {done.explained_TFs.max()}, share >= MEK {(done.explained_TFs >= mek).mean():.3f}; "
          f"{round(time.time() - t0)} s", flush=True)
pd.DataFrame(results).to_csv(f"{OUT}_mek_queries.csv", index=False)
