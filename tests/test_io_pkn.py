import pandas as pd
import pytest

import signsat as ss
from signsat.pkn import read_pkn

PKN = pd.DataFrame(
    [("r", " a ", "+"), ("A", "T1", "activation"), ("R", "T2", "-1"), ("A", "T2", -1), ("A", "T2", -1),
     ("A", "T3", "?"), ("A", None, 1), ("A", "A", 1), ("A", "T2", 1)],
    columns=["source", "target", "sign"])
PKN["reference"] = ["p1", "p2", "p3", "p4", "p5", "p6", "p7", "p8", "p9"]


def test_pkn_normalisation_rules():
    p = read_pkn(PKN)
    e = p.edges
    assert list(e.columns)[:4] == ["edge_id", "source", "target", "sign"]
    # normalised names, duplicate merged with evidence kept, both signs of A->T2 kept as separate edges
    assert set(zip(e.source, e.target, e.sign)) == {("R", "A", 1), ("A", "T1", 1), ("R", "T2", -1), ("A", "T2", -1), ("A", "A", 1), ("A", "T2", 1)}
    dup = e[(e.source == "A") & (e.target == "T2") & (e.sign == -1)].iloc[0]
    assert dup["reference"] == "p4;p5"
    codes = set(p.diagnostics["reason"])
    assert {"UNKNOWN_SIGN", "MISSING_ENDPOINT", "DUPLICATE_MERGED"} <= codes
    with pytest.raises(ValueError):
        read_pkn(pd.DataFrame({"source": ["a"], "target": ["b"]}))


def test_tf_input_rules():
    pkn = pd.DataFrame([("R", "T1", 1)], columns=["source", "target", "sign"])
    with pytest.raises(ValueError):
        ss.infer(pkn=pkn)  # neither scores nor states
    with pytest.raises(ValueError):
        ss.infer(tf_scores=pd.DataFrame({"T1": [1]}), tf_states=pd.DataFrame({"T1": [1]}), pkn=pkn)
    with pytest.raises(ValueError):
        ss.infer(tf_states=pd.DataFrame({"T1": [2]}), pkn=pkn)  # not 0/1
    # score exactly 0 is "unknown", state 0 is "down"
    r = ss.infer(tf_scores=pd.DataFrame({"T1": [0.0]}, index=["s"]), pkn=pkn)
    assert r.network("s").status == "NO_OBSERVATIONS"
    r = ss.infer(tf_states=pd.DataFrame({"T1": [0]}, index=["s"]), pkn=pkn)
    assert r.network("s").status == "SAT"
    assert r.network("s").node_states["R"] == {-1}


def test_same_pattern_same_network_and_roundtrip(tmp_path):
    pkn = pd.DataFrame([("R", "A", 1), ("A", "T1", 1), ("A", "T2", -1), ("R", "B", 1), ("B", "T2", 1)], columns=["source", "target", "sign"])
    tf = pd.DataFrame({"T1": [1, 1, -1], "T2": [-1, -1, 1]}, index=["a", "b", "c"])
    r = ss.infer(tf_scores=tf, pkn=pkn)
    ea = r.network("a").edges.drop(columns="sample")
    eb = r.network("b").edges.drop(columns="sample")
    pd.testing.assert_frame_equal(ea, eb)
    r.save(tmp_path / "out")
    for name in ("edges.csv", "nodes.csv", "tf_explanations.csv", "diagnostics.csv", "metadata.json", "model/pkn_normalised.csv"):
        assert (tmp_path / "out" / name).exists(), name
    r2 = ss.Result.load(tmp_path / "out")
    pd.testing.assert_frame_equal(r.edges.reset_index(drop=True), r2.edges.reset_index(drop=True), check_dtype=False)
    q1 = r.subnetwork("a", target={"A": 1}, focus=["T1"])
    q2 = r2.subnetwork("a", target={"A": 1}, focus=["T1"])
    assert set(q1.edges["edge_id"]) == set(q2.edges["edge_id"]) and q1.explained_tfs == q2.explained_tfs
    q1.to_graphml(tmp_path / "q.graphml")
    assert "<graphml" in (tmp_path / "q.graphml").read_text(encoding="utf-8")
    meta = r.metadata()
    assert meta["roots"] == ["R"] and meta["pkn"]["sha256_normalised_edges"]
