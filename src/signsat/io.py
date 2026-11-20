"""Saving and loading a Result.

Tables are plain CSV. The SAT model itself is not stored: the normalised PKN and the TF states are stored, and the
constraints are regenerated (deterministically) when a query needs them.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .graph import Graph
from .pkn import PKN, read_pkn
from .results import DIAG_COLS, EDGE_COLS, NODE_COLS, TF_COLS, Result, SampleNetwork


def save_result(res: Result, directory) -> None:
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    res.edges.to_csv(d / "edges.csv", index=False)
    res.nodes.to_csv(d / "nodes.csv", index=False)
    res.tf_explanations.to_csv(d / "tf_explanations.csv", index=False)
    res.diagnostics.to_csv(d / "diagnostics.csv", index=False)
    (d / "model").mkdir(exist_ok=True)
    res.pkn.edges.to_csv(d / "model" / "pkn_normalised.csv", index=False)
    res.pkn.diagnostics.to_csv(d / "model" / "pkn_diagnostics.csv", index=False)
    res.tf_states.to_csv(d / "model" / "tf_states.csv")
    meta = res.metadata()
    meta["sample_status"] = {str(k): v for k, v in res.status().items()}
    meta["unmapped_tf_columns"] = res.unmapped_tfs
    meta["witness"] = {str(s): n.witness for s, n in res.networks.items()}
    (d / "metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def load_result(directory) -> Result:
    d = Path(directory)
    meta = json.loads((d / "metadata.json").read_text(encoding="utf-8"))
    pkn_edges = pd.read_csv(d / "model" / "pkn_normalised.csv", dtype={"edge_id": str})
    pkn = read_pkn(pkn_edges)
    pkn.edges = pkn_edges  # keep ids exactly as saved
    pkn.diagnostics = pd.read_csv(d / "model" / "pkn_diagnostics.csv")
    states = pd.read_csv(d / "model" / "tf_states.csv", index_col=0)
    G = Graph.from_pkn(pkn)
    roots = [G.index[r] for r in meta["roots"]]
    res = Result(pkn=pkn, tf_states=states, roots=meta["roots"], roots_source=meta["roots_definition"], solver=meta["solver"],
                 unmapped_tfs=meta.get("unmapped_tf_columns", []), input_form=meta.get("tf_input_form", "tf_scores"))
    tables = {n: pd.read_csv(d / f"{n}.csv") for n in ("edges", "nodes", "tf_explanations", "diagnostics")}
    status = meta["sample_status"]
    samples = list(states.index)
    for sample in samples:
        key = str(sample)
        obs = {t: int(v) for t, v in states.loc[sample].items() if pd.notna(v) and t in G.index}
        net = SampleNetwork(sample, status[key], G, roots, obs)
        for attr in ("edges", "nodes", "tf_explanations", "diagnostics"):
            df = tables[attr]
            sel = df[df["sample"].astype(str) == key].reset_index(drop=True)
            setattr(net, attr, sel)
        net.witness = meta.get("witness", {}).get(key, {})
        net.pkn_edges = pkn.edges
        res.add(net)
    return res
