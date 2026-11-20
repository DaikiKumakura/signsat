"""How often are RANDOM TF states satisfiable on a real PKN? (a null for the 'UNSAT rate' question; the states are
random signs, not data). SIGNOR PKN, simulated k observed TFs, 20 samples per k, solver 'auto', 60 s limit."""
import random, sys, time
import pandas as pd
import signsat as ss
from signsat.graph import Graph
from signsat.sources import collectri_tfs, signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
df, _ = signor_pkn("data/raw/signor_release.tsv")
P = ss.read_pkn(df); G = Graph.from_pkn(P)
roots = set(G.structural_roots()); fwd = G.forward_reachable(sorted(roots))
tfs = sorted(t for t in collectri_tfs("data/raw/collectri.tsv") if t in G.index and fwd[G.index[t]] and G.index[t] not in roots)
rows = []
for k in (4, 8, 16, 32):
    rng = random.Random(k)
    states = pd.DataFrame([{t: float(rng.choice((1, -1))) for t in rng.sample(tfs, k)} for _ in range(20)], index=[f"r{i}" for i in range(20)])
    t0 = time.time()
    res = ss.infer(tf_states=(states > 0).astype(float).where(states.notna()), pkn=df, solver="auto", time_limit=60, classify=False)
    st = res.status().value_counts().to_dict()
    rows.append({"k": k, "n": 20, **st, "seconds": round(time.time() - t0, 1)})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv("benchmarks/random_states.csv", index=False)
print("TFs usable (reachable, non-root):", len(tfs))
