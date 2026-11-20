"""Query runtime on real PKNs with simulated consistent TF states (no classification at inference time).
Run: python -X utf8 benchmarks/real_query.py [limit_s] [solver]
"""
import multiprocessing as mp
import random
import sys
import time

import pandas as pd

sys.path.insert(0, "benchmarks")
import real_pkn as R
import signsat as ss
from signsat.graph import Graph
from signsat.sources import collectri_tfs


TL = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0


def work(df, obs, target, focus, solver, q):
    t0 = time.perf_counter()
    r = ss.infer(tf_scores=pd.DataFrame({t: [v] for t, v in obs.items()}, index=["s"]), pkn=df, solver=solver, classify=False, time_limit=TL)
    t1 = time.perf_counter()
    net = r.network("s")
    sub = r.subnetwork("s", target=target, focus=focus, time_limit=TL)
    q.put({"status": net.status, "infer_s": round(t1 - t0, 2), "query_s": round(time.perf_counter() - t1, 2),
           "query_edges": len(sub.edges), "explained": len(sub.explained_tfs), "complete": sub.complete, "reason": sub.empty_reason})


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    limit = float(sys.argv[1]) if len(sys.argv) > 1 else 60
    solver = sys.argv[2] if len(sys.argv) > 2 else "glucose4"
    tfs = set(collectri_tfs("data/raw/collectri.tsv"))
    rows = []
    for name in ("signor", "omnipath"):
        df = R.load(name)
        P = ss.read_pkn(df); G = Graph.from_pkn(P)
        for k in (4, 16):
            for seed in range(3):
                obs = R.simulate(df, tfs, k, seed)
                rng = random.Random(100 + seed)
                # target: a non-root node that reaches at least one observed TF; state arbitrary (+/-)
                reach = [n for n in G.nodes if n not in obs and G.out_edges[G.index[n]] and any(True for _ in [0])]
                tgt = rng.choice(reach)
                for mode in ("target", "target+focus", "focus"):
                    target = {tgt: rng.choice((1, -1))} if mode != "focus" else None
                    focus = [rng.choice(sorted(obs))] if mode != "target" else None
                    q = mp.Queue()
                    p = mp.Process(target=work, args=(df, obs, target, focus, solver, q))
                    p.start(); p.join(limit + 120)
                    if p.is_alive():
                        p.terminate(); p.join(); res = {"status": f"TIMEOUT>{limit:.0f}s"}
                    else:
                        res = q.get() if not q.empty() else {"status": "CRASH"}
                    rows.append({"pkn": name, "k_TFs": k, "seed": seed, "mode": mode, "solver": solver, **res})
                    print(rows[-1], flush=True)
                    pd.DataFrame(rows).to_csv(f"benchmarks/real_query_{solver}.csv", index=False)
