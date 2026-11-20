"""Scale benchmark on synthetic PKNs (not biological data): runtime and memory of one-sample inference.
Run: python -X utf8 benchmarks/scale.py
The PKN is random (nodes with a few roots, random edges, some feedback edges); TF states come from a random
spanning tree propagated from the roots, so every sample is satisfiable by construction.
"""
import random
import sys
import time
import tracemalloc

import pandas as pd

import signsat as ss
from signsat.graph import Graph

sys.stdout.reconfigure(encoding="utf-8")


def make(n, m, n_tf, seed):
    rng = random.Random(seed)
    names = [f"G{i}" for i in range(n)]
    edges = set()
    roots = max(2, n // 100)
    tree = []
    for v in range(roots, n):  # random recursive tree guarantees reachability from the roots
        u = rng.randrange(0, v)
        edges.add((u, v, rng.choice((1, -1))))
    while len(edges) < m:
        u, v = rng.randrange(n), rng.randrange(n)
        if u != v and rng.random() < 0.9 and u < v or u != v and rng.random() < 0.1:
            edges.add((u, v, rng.choice((1, -1))))
    df = pd.DataFrame([(names[u], names[v], s) for u, v, s in sorted(edges)], columns=["source", "target", "sign"])
    # consistent states: propagate along a random in-edge tree
    state = {i: rng.choice((1, -1)) for i in range(roots)}
    by_target = {}
    for u, v, s in edges:
        by_target.setdefault(v, []).append((u, s))
    for v in range(roots, n):
        cands = [(u, s) for u, s in by_target.get(v, []) if u in state and u < v]
        u, s = rng.choice(cands)
        state[v] = s * state[u]
    tfs = rng.sample(range(roots, n), n_tf)
    return df, {names[t]: state[t] for t in tfs}


SOLVER = sys.argv[1] if len(sys.argv) > 1 else "glucose4"
ORDER = (sys.argv[2] != "lazy") if len(sys.argv) > 2 else True
SIZES = [(100, 400, 8), (200, 800, 10), (300, 1500, 15), (600, 3000, 20)]
rows = []
for n, m, k in SIZES:
    df, obs = make(n, m, k, seed=n)
    tf = pd.DataFrame({t: [float(s)] for t, s in obs.items()}, index=["s"])
    tracemalloc.start()
    t0 = time.perf_counter()
    r = ss.infer(tf_scores=tf, pkn=df, solver=SOLVER, conflict_budget=5_000_000, order_encoding=ORDER)
    dt = time.perf_counter() - t0
    peak = tracemalloc.get_traced_memory()[1] / 1e6
    tracemalloc.stop()
    net = r.network("s")
    mdl = net.model
    rows.append({"nodes": n, "edges": m, "observed_TFs": k, "solver": SOLVER, "order": ORDER, "status": net.status, "seconds": round(dt, 2), "peak_MB": round(peak, 1),
                 "possible_edges": len(net.edges), "necessary": len(net.necessary_edges),
                 "sat_vars": mdl.nvar if mdl else None, "sat_clauses": len(mdl.clauses) if mdl else None,
                 "cycle_clauses": mdl.n_cycle_clauses if mdl else None})
    print(rows[-1], flush=True)
