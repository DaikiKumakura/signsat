"""SPEC section 6 worked example (7-edge artificial PKN): expected results fixed by hand and by brute force."""
import pandas as pd
import pytest

import signsat as ss

PKN = pd.DataFrame(
    [("R", "A", 1), ("R", "B", 1), ("R", "C", 1), ("A", "TF1", 1), ("A", "TF2", -1), ("B", "TF2", 1), ("C", "TF2", 1)],
    columns=["source", "target", "sign"])
EXPECT_IDS = {"e1": "R->A", "e2": "R->B", "e3": "R->C", "e4": "A->TF1", "e5": "A->TF2", "e6": "B->TF2", "e7": "C->TF2"}
TF = pd.DataFrame({"TF1": [1, 1, -1], "TF2": [1, -1, float("nan")]}, index=["S1", "S2", "S3"])


@pytest.fixture(scope="module")
def result():
    return ss.infer(tf_states=TF.replace({-1: 0}), pkn=PKN)


def _status(net):
    return dict(zip(net.edges["edge_id"], net.edges["status"]))


def test_roots_and_status(result):
    assert result.roots == ["R"]
    assert dict(result.status()) == {"S1": "SAT", "S2": "SAT", "S3": "SAT"}


def test_edge_classes(result):
    s1 = _status(result.network("S1"))
    assert {k for k, v in s1.items() if v == "Necessary"} == {"e1", "e4"}
    assert {k for k, v in s1.items() if v == "Alternative"} == {"e2", "e3", "e6", "e7"}
    assert len(s1) == 6
    s2 = _status(result.network("S2"))
    assert s2 == {"e1": "Necessary", "e4": "Necessary", "e5": "Necessary"}
    s3 = _status(result.network("S3"))
    assert s3 == {"e1": "Necessary", "e4": "Necessary"}


def test_node_states(result):
    n1 = result.network("S1").node_states
    assert n1["R"] == {1} and n1["A"] == {1} and n1["B"] == {1} and n1["C"] == {1}
    n3 = result.network("S3").node_states
    assert n3["R"] == {-1} and n3["A"] == {-1}
    assert "B" not in n3


def test_tf_explanations(result):
    t = result.tf_explanations.set_index(["sample", "tf"])
    assert t.loc[("S1", "TF2"), "explanation_status"] == "Multiple alternatives"
    assert t.loc[("S1", "TF1"), "explanation_status"] == "Unique"
    assert bool(t.loc[("S2", "TF2"), "explained"])
    assert ("S3", "TF2") not in t.index


def test_scores_vs_states_equivalent():
    r = ss.infer(tf_scores=TF * 2.5, pkn=PKN)
    r2 = ss.infer(tf_states=TF.replace({-1: 0}), pkn=PKN)
    pd.testing.assert_frame_equal(r.edges.reset_index(drop=True), r2.edges.reset_index(drop=True))


def test_query_table_of_spec():
    r = ss.infer(tf_states=TF.replace({-1: 0}), pkn=PKN)
    q = r.subnetwork("S1", focus=["B", "C"])
    assert set(q.edges["edge_id"]) == {"e2", "e3", "e6", "e7"}
    assert r.subnetwork("S2", focus=["B", "C"]).is_empty()
    q = r.subnetwork("S2", target={"A": 1}, focus=["TF2"])
    assert set(q.edges["edge_id"]) == {"e5"}
    assert r.subnetwork("S1", target={"A": 1}, focus=["TF2"]).is_empty()
    q = r.subnetwork("S3", target={"A": -1}, focus=["TF1"])
    assert set(q.edges["edge_id"]) == {"e4"}
    # A=-1 contradicts S1 (A must be +1 there): direction conflict, not "no path"
    e = r.subnetwork("S1", target={"A": -1}, focus=["TF1"])
    assert e.is_empty() and e.empty_reason == "DIRECTION_CONFLICT"
    assert r.subnetwork("S1", focus=["NOPE"]).empty_reason == "GENE_NOT_IN_PKN"
    # B is in the PKN but not part of the S3 explanation
    assert r.subnetwork("S3", focus=["B"]).empty_reason == "NOT_IN_SAMPLE_NETWORK"


def test_unsat_and_unreachable():
    pkn = pd.DataFrame([("R", "A", 1), ("A", "T1", 1), ("A", "T2", -1)], columns=["source", "target", "sign"])
    tf = pd.DataFrame({"T1": [1], "T2": [1], "X": [1]}, index=["s"])
    r = ss.infer(tf_scores=tf, pkn=pkn)
    n = r.network("s")
    assert n.status == "SAT" or n.status == "UNSAT"
    assert n.status == "UNSAT"  # T1=+1 needs A=+1, T2=+1 needs A=-1
    d = n.diagnostics
    assert (d["code"] == "UNSAT_CONFLICT").any() and (d["code"] == "TF_NOT_IN_PKN").any()
    assert r.subnetwork("s", focus=["A"]).empty_reason == "INFERENCE_UNAVAILABLE"
    pkn2 = pd.DataFrame([("R", "A", 1), ("Z", "Z2", 1)], columns=["source", "target", "sign"])
    r2 = ss.infer(tf_scores=pd.DataFrame({"A": [1], "Z2": [1]}, index=["s"]), pkn=pkn2, roots=["R"])
    assert r2.network("s").status == "UNREACHABLE"
    r3 = ss.infer(tf_scores=pd.DataFrame({"Q": [1]}, index=["s"]), pkn=pkn2)
    assert r3.network("s").status == "NO_OBSERVATIONS"


def test_cycle_not_part_of_explanation():
    # a 2-cycle between A and B hanging off the path must not appear
    pkn = pd.DataFrame([("R", "T", 1), ("A", "B", 1), ("B", "A", 1), ("R", "A", 1)], columns=["source", "target", "sign"])
    r = ss.infer(tf_scores=pd.DataFrame({"T": [1]}, index=["s"]), pkn=pkn)
    ids = set(r.network("s").edges["edge_id"])
    assert ids == {"e1"}


def test_root_states_fix_the_hypothesis():
    # R=+1 forced: TF2=+1 is explained via B or C but not A; R=-1 forced contradicts TF1=+1
    ok = ss.infer(tf_states=pd.DataFrame({"TF1": [1], "TF2": [1]}, index=["s"]), pkn=PKN, roots=["R"], root_states={"R": 1})
    assert ok.network("s").status == "SAT"
    bad = ss.infer(tf_states=pd.DataFrame({"TF1": [1]}, index=["s"]), pkn=PKN, roots=["R"], root_states={"R": -1})
    assert bad.network("s").status == "UNSAT"
    with pytest.raises(ValueError):
        ss.infer(tf_states=TF.replace({-1: 0}), pkn=PKN, root_states={"NOPE": 1})
