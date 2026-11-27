"""How much does the sign of a path discriminate on a signed PKN? For random source nodes, the share of reachable nodes that
are reachable with BOTH signs (walk parity). Share 1.0 means every reachable target can be 'explained' with either state.
Run: python -X utf8 benchmarks/parity_ambiguity.py   (needs data/raw/signor_release.tsv)"""
import random
import sys
from collections import deque

import pandas as pd

import signsat as ss
from signsat.graph import Graph

sys.stdout.reconfigure(encoding="utf-8")
raw = pd.read_csv("data/raw/signor_release.tsv", sep="\t", dtype=str, quoting=3).apply(lambda c: c.str.strip('"'))
raw = raw[(raw.TAX_ID == "9606") & (raw.TYPEA == "protein") & (raw.TYPEB == "protein")]
raw["sign"] = raw.EFFECT.map(lambda e: -1 if "down-regulates" in str(e).lower() else (1 if "up-regulates" in str(e).lower() else 0))
raw = raw[raw.sign != 0]
variants = {
    "all signed protein-protein": raw,
    "DIRECT only": raw[raw.DIRECT == "t"],
    "direct phosphorylation / dephosphorylation": raw[(raw.DIRECT == "t") & raw.MECHANISM.isin(["phosphorylation", "dephosphorylation"])],
    "direct phospho + transcriptional regulation": raw[(raw.DIRECT == "t") & raw.MECHANISM.isin(["phosphorylation", "dephosphorylation", "transcriptional regulation"])],
    "transcriptional regulation only": raw[raw.MECHANISM == "transcriptional regulation"],
}
rng = random.Random(1)
rows = []
for name, d in variants.items():
    P = ss.read_pkn(pd.DataFrame({"source": d.ENTITYA, "target": d.ENTITYB, "sign": d.sign}))
    G = Graph.from_pkn(P)

    def reach(u):
        seen = {(u, 1)}
        dq = deque([(u, 1)])
        while dq:
            a, za = dq.popleft()
            for k in G.out_edges[a]:
                v = (int(G.dst[k]), int(G.sign[k]) * za)
                if v not in seen:
                    seen.add(v)
                    dq.append(v)
        return seen

    cand = [i for i in range(G.n) if len(G.out_edges[i]) >= 2]
    srcs = rng.sample(cand, min(400, len(cand)))
    one = both = 0
    for u in srcs:
        r = reach(u)
        for n in {n for n, _ in r if n != u}:
            c = ((n, 1) in r) + ((n, -1) in r)
            both += c == 2
            one += c == 1
    rows.append({"PKN variant": name, "edges": P.n_edges, "nodes": G.n, "reachable pairs, one sign": one, "both signs": both,
                 "share ambiguous": round(both / max(1, one + both), 3)})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv("benchmarks/parity_ambiguity.csv", index=False)
