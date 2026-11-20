"""target / focus queries (SPEC.md section 5).

A query is answered on the explanations of ONE sample, never on the union graph, so that edges that cannot be
used together in a single consistent explanation are never joined:

* every explanation path P_f (root -> observed TF f) is a simple path of the sample's consistent explanation;
* ``target = {gene: +-1}``: the state z_gene = +-1 is imposed for every target gene (an intervention), and a path
  qualifies if it passes through at least one target gene; the returned part is the segment from that gene on;
* ``focus = [genes]``: a path qualifies if it passes through at least one focus gene; without target the whole
  path (root -> focus -> TF) is returned;
* both: the path must pass through a target gene and, at or after it, a focus gene (the same gene may be both);
  the segment from the target gene on is returned.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from xml.sax.saxutils import escape, quoteattr

import pandas as pd

from .encoding import cover
from .pkn import normalize_gene
from .results import SampleNetwork, _fmt_states

REASONS = ("GENE_NOT_IN_PKN", "NOT_IN_SAMPLE_NETWORK", "NO_CONNECTING_PATH", "DIRECTION_CONFLICT", "INFERENCE_UNAVAILABLE", "INCOMPLETE")


@dataclass
class Subnetwork:
    sample: object
    target: dict
    focus: list
    edges: pd.DataFrame
    nodes: list
    node_states: dict
    explained_tfs: list
    unexplained_tfs: list
    matched_targets: list
    matched_foci: list
    evidence: pd.DataFrame
    witness: list
    empty_reason: str | None = None
    missing_genes: list = field(default_factory=list)
    note: str | None = None
    complete: bool = True  # False: the time limit stopped the search; edges listed are proven, others unresolved

    def is_empty(self) -> bool:
        return len(self.edges) == 0

    def to_csv(self, path) -> None:
        self.edges.to_csv(path, index=False)

    def to_graphml(self, path) -> None:
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
                 '<key id="states" for="node" attr.name="possible_states" attr.type="string"/>',
                 '<key id="sign" for="edge" attr.name="sign" attr.type="int"/>',
                 '<key id="edge_id" for="edge" attr.name="edge_id" attr.type="string"/>',
                 f'<graph id={quoteattr(str(self.sample))} edgedefault="directed">']
        for n in self.nodes:
            lines.append(f'<node id={quoteattr(n)}><data key="states">{escape(_fmt_states(self.node_states.get(n, [])))}</data></node>')
        for r in self.edges.itertuples():
            lines.append(f'<edge source={quoteattr(r.source)} target={quoteattr(r.target)}>'
                         f'<data key="sign">{int(r.sign)}</data><data key="edge_id">{escape(str(r.edge_id))}</data></edge>')
        lines += ["</graph>", "</graphml>"]
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))


def _empty(net, target, focus, reason, missing=None, note=None) -> Subnetwork:
    return Subnetwork(net.sample, dict(target or {}), list(focus or []), pd.DataFrame(columns=["edge_id", "source", "target", "sign"]),
                      [], {}, [], list(net.obs), [], [], pd.DataFrame(), [], reason, list(missing or []), note)


def subnetwork(net: SampleNetwork, target: dict | None = None, focus: list | None = None,
               time_limit: float | None = None) -> Subnetwork:
    G = net.graph
    target = {normalize_gene(k): int(v) for k, v in (target or {}).items()}
    focus = [normalize_gene(g) for g in (focus or [])]
    if any(v not in (1, -1) for v in target.values()):
        raise ValueError("target values must be +1 or -1")
    if None in target or None in focus:
        raise ValueError("empty gene name in target/focus")
    if net.status != "SAT":
        return _empty(net, target, focus, "INFERENCE_UNAVAILABLE")
    if not target and not focus:
        return _whole(net)
    named = list(target) + [g for g in focus if g not in target]
    missing = [g for g in named if g not in G.index]
    if len(missing) == len(named):
        return _empty(net, target, focus, "GENE_NOT_IN_PKN", missing)
    if net.classified:
        in_net = set(net.nodes["node"]) if len(net.nodes) else set()
    else:  # not classified: any node on a possible root -> TF route counts as part of the sample's network
        mdl = net.get_model()
        in_net = {G.nodes[v] for reg in mdl.regions.values() for v in reg.nodes} | set(net.obs)
    present = [g for g in named if g in G.index]
    if not any(g in in_net for g in present):
        return _empty(net, target, focus, "NOT_IN_SAMPLE_NETWORK", missing)
    T = {G.index[g]: s for g, s in target.items() if g in G.index}
    F = {G.index[g] for g in focus if g in G.index}
    if (target and not T) or (focus and not F):
        return _empty(net, target, focus, "GENE_NOT_IN_PKN", missing)
    net.get_model().set_time_limit(time_limit)
    out = _run(net, T, F, apply_states=True)
    incomplete = net.get_model().incomplete
    net.get_model().set_time_limit(None)
    if incomplete and not out["edges"] and not out["explained"]:
        e = _empty(net, target, focus, "INCOMPLETE", missing, "time limit reached before any answer was found")
        e.complete = False
        return e
    if not out["edges"] and not out["explained"]:
        reason = "NO_CONNECTING_PATH"
        if T:
            relaxed = _run(net, T, F, apply_states=False)
            if relaxed["edges"] or relaxed["explained"]:
                reason = "DIRECTION_CONFLICT"
        return _empty(net, target, focus, reason, missing)
    sub = _assemble(net, target, focus, out, missing)
    sub.complete = not incomplete
    return sub


def _whole(net: SampleNetwork) -> Subnetwork:
    e = net.edges[["edge_id", "source", "target", "sign", "status"]].copy()
    nodes = list(net.nodes["node"]) if len(net.nodes) else []
    states = net.node_states
    P = _pkn_edges(net)
    ev = e.merge(P, on=["edge_id", "source", "target", "sign"], how="left") if len(e) else pd.DataFrame()
    return Subnetwork(net.sample, {}, [], e, nodes, states, list(net.obs), [], [], [], ev,
                      [{"tf": t, "path": p} for t, p in net.witness.items()])


def _pkn_edges(net: SampleNetwork) -> pd.DataFrame:
    pk = getattr(net, "pkn_edges", None)
    if pk is None:
        raise RuntimeError("network has no PKN attached")
    return pk


def _run(net: SampleNetwork, T: dict[int, int], F: set[int], apply_states: bool) -> dict:
    model = net.get_model()
    G = net.graph
    for t in T:
        model.zvar(t)
    for f in model.obs:
        model.zvar(f)
    s = model.fresh_solver()
    newv = model.new_var  # shared id counter; the definitional clauses below go to the private solver ``s`` only

    def add(cl):
        s.add_clause(cl)

    def or_def(var, lits):
        if not lits:
            add([-var])
            return
        add([-var] + lits)
        for l in lits:
            add([var, -l])

    def and_v(a, b):
        c = newv()
        add([-c, a])
        add([-c, b])
        add([c, -a, -b])
        return c

    desc_T = G_forward(G, list(T)) if T else None
    anc_F = G_backward(G, list(F)) if F else None
    desc_F = G_forward(G, list(F)) if F else None

    def structurally_possible(f):
        # exact necessary conditions on the PKN alone: a target must reach the TF; with focus genes the TF must lie
        # below a focus gene, which in turn must lie below a target gene
        if T and not desc_T[f]:
            return False
        if F and not desc_F[f]:
            return False
        if T and F and not any(desc_T[g] for g in F):
            return False
        return True

    active = {f: reg for f, reg in model.regions.items() if structurally_possible(f)}
    TRUE = newv()
    add([TRUE])
    obs = model.obs_lits()
    tl = [model.zvar(t) if st == 1 else -model.zvar(t) for t, st in T.items()] if apply_states else []
    onpath: dict[tuple[int, int], int] = {}
    A: dict[tuple[int, int], int] = {}
    Q: dict[int, int] = {}
    qe: dict[tuple[int, int], int] = {}
    for f, reg in active.items():
        eset = set(reg.edges)
        for v in reg.nodes:
            if v == f:
                onpath[(f, v)] = TRUE
            else:
                ov = newv()
                or_def(ov, [model.inv[(f, v)], model.outv[(f, v)]])
                onpath[(f, v)] = ov
        if T:
            for v in reg.nodes:
                A[(f, v)] = newv()
            for v in reg.nodes:
                terms = []
                if v in T:
                    terms.append(onpath[(f, v)])
                for k in G.in_edges[v]:
                    if k in eset:
                        terms.append(and_v(model.p[(f, k)], A[(f, int(G.src[k]))]))
                or_def(A[(f, v)], terms)
        qf = newv()
        Q[f] = qf
        if T and not F:
            or_def(qf, [A[(f, f)]])
        elif F and not T:
            or_def(qf, [onpath[(f, g)] for g in F if g in reg.nodes])
        else:
            or_def(qf, [A[(f, g)] for g in F if g in reg.nodes])
        for k in reg.edges:
            # exact structural filter: an edge of the answer lies after a target gene, or (focus only) on a path
            # through a focus gene, so it must be below a target / above-or-below a focus gene in the PKN
            if T and not desc_T[int(G.src[k])]:
                continue
            if F and not T and not (anc_F[int(G.dst[k])] or desc_F[int(G.src[k])]):
                continue
            if T:
                base = and_v(model.p[(f, k)], A[(f, int(G.src[k]))])
            else:
                base = model.p[(f, k)]
            qe[(f, k)] = and_v(base, qf)
    edge_var: dict[int, int] = {}
    for k in sorted({k for (_, k) in qe}):
        v = newv()
        or_def(v, [qe[(f, k)] for f in active if (f, k) in qe])
        edge_var[k] = v
    node_var: dict[int, int] = {}
    nodes_all = {int(G.src[k]) for k in edge_var} | {int(G.dst[k]) for k in edge_var}
    for n in nodes_all:
        v = newv()
        or_def(v, [edge_var[k] for k in G.out_edges[n] + G.in_edges[n] if k in edge_var])
        node_var[n] = v

    seen_edges: dict[int, dict] = {}

    def record(m):
        for k, ev in edge_var.items():
            if k not in seen_edges and ev in m:
                for f in active:
                    if (f, k) in qe and qe[(f, k)] in m:
                        seen_edges[k] = {"tf": f, "path": model.path_nodes(m, f)}
                        break

    def sat(lits):
        r = model.solve(obs + tl + lits, solver=s)
        if r:
            record(model.model)
        return r

    explained, matched_t, matched_f = [], set(), set()
    for f in active:
        r = sat([Q[f]])
        if r:
            explained.append(f)
            if T:
                for t in T:
                    if t in model.regions[f].nodes and t not in matched_t and sat([Q[f], onpath[(f, t)]]):
                        matched_t.add(t)
            for gnode in F:
                if gnode in model.regions[f].nodes and gnode not in matched_f:
                    key = A[(f, gnode)] if T else onpath[(f, gnode)]
                    if sat([Q[f], key]):
                        matched_f.add(gnode)
    base = obs + tl
    cover(model, base, dict(edge_var), record, solver=s)  # every edge that can be in the answer (final UNSAT proves the rest cannot)
    wnode = {}
    for n, nv_ in node_var.items():
        for st in (1, -1):
            w = newv()
            zl = model.z[n] if st == 1 else -model.z[n]
            add([-w, nv_]); add([-w, zl]); add([w, -nv_, -zl])
            wnode[(n, st)] = w
    covn = cover(model, base, wnode, None, solver=s) or set()
    seen_nodes: dict[int, set[int]] = {}
    for (n, st) in covn:
        seen_nodes.setdefault(n, set()).add(st)
    return {"edges": sorted(seen_edges), "witness": seen_edges, "nodes": seen_nodes, "explained": explained,
            "matched_t": matched_t, "matched_f": matched_f}


def G_forward(G, starts):
    return G.forward_reachable(list(starts))


def G_backward(G, targets):
    import numpy as np
    seen = np.zeros(G.n, dtype=bool)
    for t in targets:
        seen |= G.backward_reachable(t)
    return seen


def _assemble(net: SampleNetwork, target, focus, out, missing) -> Subnetwork:
    G = net.graph
    P = _pkn_edges(net)
    rows = [{"edge_id": G.edge_ids[k], "source": G.nodes[int(G.src[k])], "target": G.nodes[int(G.dst[k])], "sign": int(G.sign[k])}
            for k in out["edges"]]
    edges = pd.DataFrame(rows, columns=["edge_id", "source", "target", "sign"])
    nodes = sorted(G.nodes[n] for n in out["nodes"])
    states = {G.nodes[n]: set(s) for n, s in out["nodes"].items()}
    expl = sorted(G.nodes[f] for f in out["explained"])
    unexpl = sorted(set(net.obs) - set(expl))
    ev = edges.merge(P, on=["edge_id", "source", "target", "sign"], how="left") if len(edges) else pd.DataFrame()
    wit = [{"edge_id": G.edge_ids[k], "tf": G.nodes[w["tf"]], "path": [G.nodes[v] for v in w["path"]]}
           for k, w in sorted(out["witness"].items())]
    note = None
    if not len(edges):
        note = "the qualifying paths contain no edge after the target gene (the target is an observed TF itself)"
    return Subnetwork(net.sample, target, focus, edges, nodes, states, expl, unexpl,
                      sorted(G.nodes[t] for t in out["matched_t"]), sorted(G.nodes[g] for g in out["matched_f"]),
                      ev, wit, None, missing, note)
