"""Feasibility on real PKNs (OmniPath / SIGNOR) with SIMULATED but consistent TF states.

TF states: random root states propagated along a random BFS tree of the PKN (so every sample is SAT by
construction); k observed TFs are drawn from the CollecTRI regulators present in the PKN. This measures
runtime and result size only; it says nothing about biological validity. Each run has a wall-clock limit.
Run: python -X utf8 benchmarks/real_pkn.py [seconds_limit]
"""
import multiprocessing as mp
import random
import sys
import time
from collections import deque

import pandas as pd

import signsat as ss
from signsat.graph import Graph
from signsat.sources import collectri_tfs, omnipath_pkn, signor_pkn


def load(name):
    if name == "omnipath":
        df, _ = omnipath_pkn("data/raw/omnipath_interactions.tsv")
        df = df[~df.source.str.contains("_") & ~df.target.str.contains("_")]  # drop complexes (unsupported)
    else:
        df, _ = signor_pkn("data/raw/signor_release.tsv")
    return df


def simulate(df, tfs, k, seed):
    rng = random.Random(seed)
    P = ss.read_pkn(df); G = Graph.from_pkn(P)
    roots = G.structural_roots()
    state = {r: rng.choice((1, -1)) for r in roots}
    dq = deque(rng.sample(roots, len(roots)))
    while dq:
        u = dq.popleft()
        outs = list(G.out_edges[u]); rng.shuffle(outs)
        for e in outs:
            v = int(G.dst[e])
            if v not in state:
                state[v] = int(G.sign[e]) * state[u]
                dq.append(v)
    cand = [t for t in tfs if t in G.index and G.index[t] in state and G.index[t] not in set(roots)]
    pick = rng.sample(sorted(cand), k)
    return {t: float(state[G.index[t]]) for t in pick}


def worker(df, obs, solver, q, classify, tl):
    t0 = time.perf_counter()
    r = ss.infer(tf_scores=pd.DataFrame({t: [v] for t, v in obs.items()}, index=["s"]), pkn=df, solver=solver, classify=classify, time_limit=tl)
    net = r.network("s")
    m = net.model
    q.put({"status": net.status, "seconds": round(time.perf_counter() - t0, 2), "possible_edges": len(net.edges),
           "necessary": len(net.necessary_edges), "sat_vars": getattr(m, "nvar", None), "sat_clauses": len(m.clauses) if m else None})


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    limit = float(sys.argv[1]) if len(sys.argv) > 1 else 90
    solver = sys.argv[2] if len(sys.argv) > 2 else "glucose4"
    classify = (sys.argv[3] != "feasibility") if len(sys.argv) > 3 else True
    tfs = set(collectri_tfs("data/raw/collectri.tsv"))
    rows = []
    for name in ("signor", "omnipath"):
        df = load(name)
        for k in (2, 4, 8, 16, 32):
            for seed in range(3):
                obs = simulate(df, tfs, k, seed)
                q = mp.Queue()
                p = mp.Process(target=worker, args=(df, obs, solver, q, classify, limit))
                p.start(); p.join(limit + 60)
                if p.is_alive():
                    p.terminate(); p.join()
                    res = {"status": f"HARD-TIMEOUT>{limit + 60:.0f}s"}
                else:
                    res = q.get() if not q.empty() else {"status": "CRASH"}
                rows.append({"pkn": name, "k_TFs": k, "seed": seed, "solver": solver, "classify": classify, **res})
                print(rows[-1], flush=True)
                pd.DataFrame(rows).to_csv(f"benchmarks/real_pkn_{solver}_{'classify' if classify else 'feasibility'}.csv", index=False)
