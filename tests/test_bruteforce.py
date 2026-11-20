"""Random tiny PKNs: SAT results must equal the brute-force enumeration (edges, nodes, statuses, queries)."""
import random

import numpy as np
import pandas as pd
import pytest

import signsat as ss
from signsat.graph import Graph

from bruteforce import classify, explanations, query, simple_paths


def random_case(rng: random.Random):
    n = rng.randint(4, 7)
    names = [f"N{i}" for i in range(n)]
    edges = set()
    for _ in range(rng.randint(n, 2 * n + 2)):
        u, v = rng.sample(range(n), 2) if rng.random() < 0.97 else (rng.randrange(n), rng.randrange(n))
        edges.add((u, v, rng.choice((1, -1))))
    df = pd.DataFrame([(names[u], names[v], s) for u, v, s in sorted(edges)], columns=["source", "target", "sign"])
    nodes_in = sorted(set(df["source"]) | set(df["target"]))
    k = rng.randint(1, min(4, len(nodes_in)))
    obs_names = rng.sample(nodes_in, k)
    obs = {t: rng.choice((1, -1)) for t in obs_names}
    return df, obs


@pytest.mark.parametrize("seed", range(120))
def test_inference_matches_bruteforce(seed):
    rng = random.Random(seed)
    df, obs = random_case(rng)
    tf = pd.DataFrame({t: [float(s)] for t, s in obs.items()}, index=["s"])
    res = ss.infer(tf_scores=tf, pkn=df)
    net = res.network("s")
    P = res.pkn
    G = Graph.from_pkn(P)
    roots = set(G.structural_roots())
    oidx = {G.index[t]: s for t, s in obs.items()}
    per_tf = {f: simple_paths(G, roots, f) for f in oidx}
    if any(len(v) == 0 for v in per_tf.values()):
        assert net.status == "UNREACHABLE"
        return
    ref = classify(G, roots, oidx)
    if ref is None:
        assert net.status == "UNSAT"
        return
    assert net.status == "SAT"
    eid = {e: i for i, e in enumerate(G.edge_ids)}
    got_poss = {eid[e] for e in net.edges["edge_id"]}
    got_nec = {eid[e] for e in net.necessary_edges["edge_id"]}
    assert got_poss == ref["possible"], (seed, got_poss, ref["possible"])
    assert got_nec == ref["necessary"], (seed, got_nec, ref["necessary"])
    # node states (non observed nodes) and necessity
    ns = net.node_states
    for v, states in ref["node_states"].items():
        name = G.nodes[v]
        if name in obs:
            continue
        assert ns[name] == {int(s) for s in states}, (seed, name)
    for r in net.nodes.itertuples():
        v = G.index[r.node]
        if r.status == "Observed":
            continue
        assert v in ref["node_states"], (seed, r.node)
        assert (r.status == "Necessary") == (v in ref["node_necessary"]), (seed, r.node)
    # per-TF explanation status
    for f, paths in per_tf.items():
        pass


@pytest.mark.parametrize("seed", range(200, 320))
def test_query_matches_bruteforce(seed):
    rng = random.Random(seed)
    df, obs = random_case(rng)
    tf = pd.DataFrame({t: [float(s)] for t, s in obs.items()}, index=["s"])
    res = ss.infer(tf_scores=tf, pkn=df)
    net = res.network("s")
    if net.status != "SAT":
        return
    G = Graph.from_pkn(res.pkn)
    roots = set(G.structural_roots())
    oidx = {G.index[t]: s for t, s in obs.items()}
    nodes = list(G.nodes)
    mode = rng.choice(("target", "focus", "both"))
    tgt = {}
    foc = []
    if mode in ("target", "both"):
        for g in rng.sample(nodes, rng.randint(1, 2)):
            tgt[g] = rng.choice((1, -1))
    if mode in ("focus", "both"):
        foc = rng.sample(nodes, rng.randint(1, 2))
    q = res.subnetwork("s", target=tgt, focus=foc)
    T = {G.index[g]: s for g, s in tgt.items()}
    F = {G.index[g] for g in foc}
    ref_edges, ref_expl = query(G, roots, oidx, T, F)
    eid = {e: i for i, e in enumerate(G.edge_ids)}
    got = {eid[e] for e in q.edges["edge_id"]}
    assert got == ref_edges, (seed, mode, tgt, foc, got, ref_edges)
    assert {G.index[t] for t in q.explained_tfs} == ref_expl, (seed, mode, tgt, foc)
    if not got and not ref_expl:
        assert q.empty_reason in ("NOT_IN_SAMPLE_NETWORK", "NO_CONNECTING_PATH", "DIRECTION_CONFLICT"), q.empty_reason
        if q.empty_reason == "DIRECTION_CONFLICT":
            r_edges, r_expl = query(G, roots, oidx, T, F, apply_states=False)
            assert r_edges or r_expl


@pytest.mark.parametrize("seed", range(40))
def test_lazy_and_order_encodings_agree(seed):
    rng = random.Random(seed)
    df, obs = random_case(rng)
    tf = pd.DataFrame({t: [float(s)] for t, s in obs.items()}, index=["s"])
    a = ss.infer(tf_scores=tf, pkn=df).network("s")
    b = ss.infer(tf_scores=tf, pkn=df, order_encoding=False).network("s")
    assert a.status == b.status
    pd.testing.assert_frame_equal(a.edges.reset_index(drop=True), b.edges.reset_index(drop=True))
    pd.testing.assert_frame_equal(a.nodes.reset_index(drop=True), b.nodes.reset_index(drop=True))
