"""Does limiting the path length restore discrimination? Same hypothesis test as gse329052_hypothesis_test.py with
explanations restricted to paths of at most L edges (a departure from SignSAT's parameter-free design, tried here only to
see whether the saturation seen without a limit is the problem). Walk-parity BFS (validated against SAT on 120 cases).
Run: python -X utf8 validation/gse329052_depth_test.py [k]"""
import sys
from collections import deque

import numpy as np
import pandas as pd

sys.argv = [sys.argv[0], sys.argv[1] if len(sys.argv) > 1 else "15"]
src = open("validation/gse329052_hypothesis_test.py", encoding="utf-8").read().split("# ---- 2. hypotheses")[0]
exec(src)  # data, TF deltas, G, sets, all_tfs, rng, reach_size, ...


def reach_depth(root_names, sts, L):
    seen = {}
    dq = deque()
    for r, st in zip(root_names, sts):
        seen[(G.index[r], st)] = 0
        dq.append((G.index[r], st))
    while dq:
        u, zu = dq.popleft()
        d = seen[(u, zu)]
        if d >= L:
            continue
        for k in G.out_edges[u]:
            v = (int(G.dst[k]), int(G.sign[k]) * zu)
            if v not in seen:
                seen[v] = d + 1
                dq.append(v)
    return set(seen)


def cnt(seen, obs):
    return sum((G.index[tf], y) in seen for tf, y in obs.items())


pool = [g for g in G.nodes if g not in all_tfs and g not in ("MAP2K1", "MAP2K2") and len(G.out_edges[G.index[g]]) >= 2]
rows = []
N_NULL = 1000
sample_pairs = [tuple(rng.sample(pool, 2)) for _ in range(N_NULL)]
for L in (2, 3, 4, 5, 6, 8, 12):
    mek_m, mek_p = reach_depth(["MAP2K1", "MAP2K2"], (-1, -1), L), reach_depth(["MAP2K1", "MAP2K2"], (1, 1), L)
    null_seen = [reach_depth(list(pr), (-1, -1), L) for pr in sample_pairs]
    single = {g: reach_depth([g], (-1,), L) for g in G.nodes}
    for (variant, drug), (delta, top) in sets.items():
        pm = []
        for p in PATS:
            obs = {g: int(np.sign(delta.loc[g, p])) for g in top if delta.loc[g, p] != 0}
            n_m, n_p = cnt(mek_m, obs), cnt(mek_p, obs)
            null = np.array([cnt(t, obs) for t in null_seen])
            vals = np.array([cnt(s_, obs) for s_ in single.values()])
            r_m = float((vals > n_m).sum() + 1)
            rows.append({"L": L, "variant": variant, "drug": drug, "patient": p, "observed": len(obs), "MEK-1": n_m, "MEK+1": n_p,
                         "null_median": float(np.median(null)), "p_ge_MEK-1": float((null >= n_m).mean()),
                         "p_ge_MEK+1": float((null >= n_p).mean()), "rank_of_MEK_among_single_nodes(best=1)": r_m})
    d = pd.DataFrame(rows)
    sub = d[d.L == L].groupby(["variant", "drug"])[["MEK-1", "MEK+1", "null_median", "p_ge_MEK-1", "p_ge_MEK+1"]].mean().round(3)
    print(f"\nL = {L} (means over 4 patients)\n{sub.to_string()}", flush=True)
    d.to_csv("validation/gse329052_depth_k" + sys.argv[1] + ".csv", index=False)
