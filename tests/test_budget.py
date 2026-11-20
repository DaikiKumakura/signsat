import pandas as pd

import signsat as ss

PKN = pd.DataFrame([("R", "A", 1), ("R", "B", 1), ("A", "T1", 1), ("B", "T1", 1)], columns=["source", "target", "sign"])
TF = pd.DataFrame({"T1": [1.0]}, index=["s"])


def test_classify_false_still_answers_queries():
    r = ss.infer(tf_scores=TF, pkn=PKN, classify=False)
    n = r.network("s")
    assert n.status == "SAT" and not n.classified and len(n.edges) == 0
    q = r.subnetwork("s", focus=["A"])
    assert set(q.edges["source"] + ">" + q.edges["target"]) == {"R>A", "A>T1"}
    n.classify()
    assert n.classified and len(n.edges) == 4 and set(n.necessary_edges.edge_id) == set()


def test_zero_time_limit_is_incomplete_not_empty():
    r = ss.infer(tf_scores=TF, pkn=PKN, time_limit=0.0)
    n = r.network("s")
    assert n.status == "INCOMPLETE"
    assert (r.diagnostics["code"] == "INCOMPLETE").any()
    r2 = ss.infer(tf_scores=TF, pkn=PKN, classify=False)
    q = r2.subnetwork("s", focus=["A"], time_limit=0.0)
    assert q.empty_reason == "INCOMPLETE" and q.complete is False


def test_auto_solver_portfolio_matches_default():
    a = ss.infer(tf_scores=TF, pkn=PKN).network("s")
    b = ss.infer(tf_scores=TF, pkn=PKN, solver="auto", time_limit=30).network("s")
    assert b.status == "SAT"
    pd.testing.assert_frame_equal(a.edges.reset_index(drop=True), b.edges.reset_index(drop=True))
