"""Structure of real PKNs: size, roots, SCCs, TF coverage. Inputs in data/raw (see docs/data.md)."""
import sys, time
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import signsat as ss
from signsat.graph import Graph
from signsat.sources import collectri_tfs, omnipath_pkn, signor_pkn

sys.stdout.reconfigure(encoding="utf-8")
tfs = set(collectri_tfs("data/raw/collectri.tsv"))
rows = []
for name, (df, bad) in {"omnipath": omnipath_pkn("data/raw/omnipath_interactions.tsv"),
                        "signor": signor_pkn("data/raw/signor_release.tsv")}.items():
    P = ss.read_pkn(df); G = Graph.from_pkn(P)
    roots = G.structural_roots(); fwd = G.forward_reachable(roots)
    k, lab = connected_components(csr_matrix((np.ones(len(G.src)), (G.src, G.dst)), shape=(G.n, G.n)), directed=True, connection="strong")
    sz = np.bincount(lab)
    tf_in = [t for t in tfs if t in G.index]
    rows.append({"pkn": name, "input_rows": len(df), "unusable_rows": len(bad), "edges": P.n_edges, "nodes": G.n,
                 "root_nodes": len(roots), "reachable_from_roots": int(fwd.sum()), "largest_SCC": int(sz.max()),
                 "SCCs_gt1": int((sz > 1).sum()), "CollecTRI_TFs": len(tfs), "TFs_in_PKN": len(tf_in),
                 "TFs_reachable": int(sum(fwd[G.index[t]] for t in tf_in))})
    print(rows[-1], flush=True)
    print("  unusable:", (bad["reason"].str.split(":").str[0].value_counts().head(4).to_dict() if len(bad) else {}))
    print("  root nodes (first 15):", [G.nodes[i] for i in roots[:15]])
pd.DataFrame(rows).to_csv("benchmarks/pkn_stats.csv", index=False)
