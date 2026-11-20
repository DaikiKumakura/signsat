"""Result objects: per-sample network, whole-run result, and the tables they expose."""
from __future__ import annotations

import hashlib
import platform
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .graph import Graph
from .pkn import PKN

EDGE_COLS = ["sample", "edge_id", "source", "target", "sign", "status"]
NODE_COLS = ["sample", "node", "possible_states", "status"]
TF_COLS = ["sample", "tf", "observed_state", "explained", "explanation_status"]
DIAG_COLS = ["sample", "level", "code", "detail"]


def _fmt_states(states) -> str:
    return ",".join("+1" if s == 1 else "-1" for s in sorted(states, reverse=True))


class SampleNetwork:
    def __init__(self, sample, status, graph: Graph, roots: list[int], obs: dict[str, int]):
        self.sample = sample
        self.status = status
        self.graph = graph
        self.roots = roots
        self.obs = obs  # TF name -> +1/-1 (mapped TFs only)
        self.model = None
        self.diagnostics = pd.DataFrame(columns=DIAG_COLS)
        self.edges = pd.DataFrame(columns=EDGE_COLS)
        self.nodes = pd.DataFrame(columns=NODE_COLS)
        self.tf_explanations = pd.DataFrame(columns=TF_COLS)
        self.witness: dict[str, list[str]] = {}
        self.pkn_edges = None
        self.classified = True

    # ------------------------------------------------------------------ fill
    def _finish_empty(self):
        if not (self.status == "INCOMPLETE" and len(self.edges)):  # keep the proven partial edges of an unfinished run
            self.edges = pd.DataFrame(columns=EDGE_COLS)
        self.nodes = pd.DataFrame(columns=NODE_COLS)
        rows = [{"sample": self.sample, "tf": t, "observed_state": s, "explained": False,
                 "explanation_status": "Not explained (" + self.status + ")"} for t, s in self.obs.items()]
        self.tf_explanations = pd.DataFrame(rows, columns=TF_COLS)

    def _fill(self, model, out: dict):
        g = self.graph
        self.classified = not out.get("unclassified", False)
        if not self.classified:  # satisfiable, edges not classified: queries still work
            self.tf_explanations = pd.DataFrame([(self.sample, t, s_, True, "Not classified") for t, s_ in self.obs.items()], columns=TF_COLS)
            self.witness = {g.nodes[f]: [g.nodes[v] for v in path] for f, path in out["witness"].items()}
            return
        e = []
        for k, st in sorted(out["edge_status"].items()):
            if st == "Not inferred":
                continue
            e.append((self.sample, g.edge_ids[k], g.nodes[int(g.src[k])], g.nodes[int(g.dst[k])], int(g.sign[k]), st))
        self.edges = pd.DataFrame(e, columns=EDGE_COLS)
        obs_idx = {g.index[t]: s for t, s in self.obs.items()}
        nrows = {}
        for v, (states, st) in out["node_info"].items():
            if not states:  # node never lies on an explanation (all its edges are Not inferred)
                continue
            nrows[v] = (g.nodes[v], _fmt_states(states), st)
        for v, s in obs_idx.items():
            nrows[v] = (g.nodes[v], _fmt_states([s]), "Observed")
        self.nodes = pd.DataFrame([(self.sample,) + r for _, r in sorted(nrows.items(), key=lambda kv: g.nodes[kv[0]])],
                                  columns=NODE_COLS)
        tfrows = []
        for t, s_ in self.obs.items():
            f = g.index[t]
            tfrows.append((self.sample, t, s_, True, "Multiple alternatives" if out["tf_multi"].get(f, False) else "Unique"))
        self.tf_explanations = pd.DataFrame(tfrows, columns=TF_COLS)
        self.witness = {g.nodes[f]: [g.nodes[v] for v in path] for f, path in (out["witness"] or {}).items()}

    def _fill_partial(self, edge_idx):
        """Edges proven possible before the budget ran out (necessity and the rest are unresolved)."""
        g = self.graph
        rows = [(self.sample, g.edge_ids[k], g.nodes[int(g.src[k])], g.nodes[int(g.dst[k])], int(g.sign[k]), "Possible (partial)")
                for k in edge_idx]
        self.edges = pd.DataFrame(rows, columns=EDGE_COLS)

    def copy_for(self, sample, dropped: list[str]) -> "SampleNetwork":
        n = SampleNetwork(sample, self.status, self.graph, self.roots, dict(self.obs))
        n.model = self.model
        for attr in ("edges", "nodes", "tf_explanations"):
            df = getattr(self, attr).copy()
            if "sample" in df.columns:
                df["sample"] = sample
            setattr(n, attr, df)
        d = self.diagnostics[self.diagnostics["code"] != "TF_NOT_IN_PKN"].copy()
        d["sample"] = sample
        extra = pd.DataFrame([{"sample": sample, "level": "info", "code": "TF_NOT_IN_PKN", "detail": t} for t in dropped],
                             columns=DIAG_COLS)
        n.diagnostics = pd.concat([d, extra], ignore_index=True)
        n.witness = dict(self.witness)
        n.pkn_edges = self.pkn_edges
        n.classified = self.classified
        return n

    # --------------------------------------------------------------- tables
    @property
    def necessary_edges(self) -> pd.DataFrame:
        return self.edges[self.edges["status"] == "Necessary"].reset_index(drop=True)

    @property
    def alternative_edges(self) -> pd.DataFrame:
        return self.edges[self.edges["status"] == "Alternative"].reset_index(drop=True)

    @property
    def node_states(self) -> dict[str, set[int]]:
        return {r.node: {int(x) for x in r.possible_states.split(",")} for r in self.nodes.itertuples()}

    def classify(self) -> None:
        """Classify all edges now (for a network inferred with classify=False). Can be slow on large PKNs."""
        from .inference import analyze_sample
        if self.status != "SAT" or self.classified:
            return
        model = self.get_model()
        out = analyze_sample(model, self.graph, self.sample, {})
        if out["status"] != "SAT":
            self.status = out["status"]
            return
        self._fill(model, out)

    def get_model(self):
        if self.model is None:
            from .encoding import SampleModel
            obs = {self.graph.index[t]: s for t, s in self.obs.items()}
            self.model = SampleModel(self.graph, self.roots, obs)
        return self.model

    def to_csv(self, path) -> None:
        self.edges.to_csv(path, index=False)


@dataclass
class Result:
    pkn: PKN
    tf_states: pd.DataFrame
    roots: list[str]
    roots_source: str
    solver: str
    unmapped_tfs: list[str]
    input_form: str = "tf_scores"
    networks: dict = field(default_factory=dict)

    def add(self, net: SampleNetwork) -> None:
        self.networks[net.sample] = net

    @property
    def samples(self) -> list:
        return list(self.networks)

    def network(self, sample) -> SampleNetwork:
        if sample not in self.networks:
            raise KeyError(f"unknown sample: {sample!r}")
        return self.networks[sample]

    def status(self) -> pd.Series:
        return pd.Series({s: n.status for s, n in self.networks.items()}, name="status")

    def _cat(self, attr: str, cols) -> pd.DataFrame:
        dfs = [getattr(n, attr) for n in self.networks.values() if len(getattr(n, attr))]
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame(columns=cols)

    @property
    def edges(self) -> pd.DataFrame:
        return self._cat("edges", EDGE_COLS)

    @property
    def nodes(self) -> pd.DataFrame:
        return self._cat("nodes", NODE_COLS)

    @property
    def tf_explanations(self) -> pd.DataFrame:
        return self._cat("tf_explanations", TF_COLS)

    @property
    def diagnostics(self) -> pd.DataFrame:
        d = self._cat("diagnostics", DIAG_COLS)
        extra = []
        for t in self.unmapped_tfs:
            extra.append({"sample": None, "level": "info", "code": "TF_COLUMN_NOT_IN_PKN", "detail": t})
        for r in self.pkn.diagnostics.itertuples():
            extra.append({"sample": None, "level": "info", "code": "PKN_" + r.reason, "detail": r.detail})
        return pd.concat([d, pd.DataFrame(extra, columns=DIAG_COLS)], ignore_index=True) if extra else d

    def pkn_hash(self) -> str:
        return hashlib.sha256(self.pkn.edges.to_csv(index=False).encode("utf-8")).hexdigest()

    def metadata(self) -> dict:
        import pysat
        from . import __version__
        return {
            "signsat_version": __version__,
            "pkn": {"label": getattr(self.pkn, "source_label", None), "n_edges": int(self.pkn.n_edges),
                    "n_nodes": len(self.pkn.nodes), "sha256_normalised_edges": self.pkn_hash()},
            "tf_input_form": self.input_form,
            "roots": self.roots,
            "roots_definition": self.roots_source,
            "solver": self.solver,
            "pysat_version": getattr(pysat, "__version__", None),
            "python": platform.python_version(),
            "status_counts": self.status().value_counts().to_dict(),
            "n_samples": len(self.networks),
        }

    # ----------------------------------------------------------------- query
    def subnetwork(self, sample, target: dict | None = None, focus: list | None = None, time_limit: float | None = None):
        from .query import subnetwork
        return subnetwork(self.network(sample), target=target, focus=focus, time_limit=time_limit)

    def save(self, directory) -> None:
        from .io import save_result
        save_result(self, directory)

    @classmethod
    def load(cls, directory) -> "Result":
        from .io import load_result
        return load_result(directory)
