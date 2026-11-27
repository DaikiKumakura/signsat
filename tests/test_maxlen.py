"""max_path_length: SAT results equal brute force when explanation paths are limited to L edges."""
import random

import pandas as pd
import pytest

import signsat as ss
from signsat.graph import Graph

from bruteforce import classify, simple_paths
from test_bruteforce import random_case


@pytest.mark.parametrize("seed", range(300, 380))
def test_max_path_length_matches_bruteforce(seed):
    rng = random.Random(seed)
    df, obs = random_case(rng)
    L = rng.choice((1, 2, 3, 4))
    tf = pd.DataFrame({t: [float(s)] for t, s in obs.items()}, index=["s"])
    res = ss.infer(tf_scores=tf, pkn=df, max_path_length=L)
    net = res.network("s")
    G = Graph.from_pkn(res.pkn)
    roots = set(G.structural_roots())
    oidx = {G.index[t]: s for t, s in obs.items()}
    if any(len(simple_paths(G, roots, f, L)) == 0 for f in oidx):
        assert net.status == "UNREACHABLE", (seed, L)
        return
    ref = classify(G, roots, oidx, L)
    if ref is None:
        assert net.status == "UNSAT", (seed, L)
        return
    assert net.status == "SAT", (seed, L)
    eid = {e: i for i, e in enumerate(G.edge_ids)}
    assert {eid[e] for e in net.edges.edge_id} == ref["possible"], (seed, L)
    assert {eid[e] for e in net.necessary_edges.edge_id} == ref["necessary"], (seed, L)
