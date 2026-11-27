"""Inference: TF states + signed PKN -> per-sample possible / necessary causal networks."""
from __future__ import annotations

import time

import numpy as np
import pandas as pd

from .encoding import SampleModel, cover, shrink
from .graph import Graph
from .pkn import PKN, normalize_gene, read_pkn
from .results import Result, SampleNetwork

STATUS_SAT = "SAT"
STATUS_UNSAT = "UNSAT"
STATUS_UNREACHABLE = "UNREACHABLE"
STATUS_NO_OBS = "NO_OBSERVATIONS"
STATUS_INCOMPLETE = "INCOMPLETE"


def encode_tf_input(tf_scores=None, tf_states=None) -> pd.DataFrame:
    """Samples x TFs table of +1 / -1 / NaN.

    ``tf_scores``: continuous scores; > 0 -> +1, < 0 -> -1, 0 or missing -> NaN (no threshold).
    ``tf_states``: already binarised; 1 = up (+1), 0 = down (-1), NA = unknown.
    """
    if (tf_scores is None) == (tf_states is None):
        raise ValueError("give exactly one of tf_scores or tf_states")
    df = pd.DataFrame(tf_scores if tf_scores is not None else tf_states).copy()
    df.columns = [normalize_gene(c) for c in df.columns]
    if df.columns.isna().any() or df.columns.duplicated().any():
        raise ValueError("TF names must be non-empty and unique after normalisation (strip, upper case)")
    if df.index.duplicated().any():
        raise ValueError("sample names must be unique")
    vals = df.apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    if tf_scores is not None:
        out = np.where(vals > 0, 1.0, np.where(vals < 0, -1.0, np.nan))
    else:
        bad = ~np.isnan(vals) & ~np.isin(vals, (0.0, 1.0))
        if bad.any():
            raise ValueError("tf_states must be 1, 0 or NA")
        out = np.where(vals == 1, 1.0, np.where(vals == 0, -1.0, np.nan))
    return pd.DataFrame(out, index=df.index, columns=df.columns)


def _sgn(s: int) -> str:
    return "+1" if s == 1 else "-1"


def analyze_sample(model: SampleModel, graph: Graph, name: str, sample_obs_names: dict[int, str]) -> dict:
    """Full possible / necessary classification for one built model."""
    obs = model.obs_lits()
    path_sigs: dict[int, set] = {f: set() for f in model.regions}
    witness = None

    def record(m):
        nonlocal witness
        for f in model.regions:
            path_sigs[f].add(tuple(sorted(model.used_edges(m, f))))
        if witness is None:
            witness = {f: model.path_nodes(m, f) for f in model.regions}

    r = model.solve(obs)
    if r is None:
        return {"status": STATUS_INCOMPLETE}
    if r is False:
        return {"status": STATUS_UNSAT, "conflict": model.conflicting_observations()}
    first = model.model
    record(first)
    xs = dict(model.x)
    # ---- edges: Possible = true in some model, Necessary = true in all models
    poss = cover(model, obs, xs, record)
    if model.incomplete:
        return {"status": STATUS_INCOMPLETE, "partial_edges": sorted(poss)}
    nec = shrink(model, obs, {k: xs[k] for k in poss}, {k for k in poss if xs[k] in first}, record)
    if model.incomplete:
        return {"status": STATUS_INCOMPLETE, "partial_edges": sorted(poss)}
    edge_status = {k: ("Necessary" if k in nec else "Alternative") if k in poss else "Not inferred" for k in xs}
    # ---- is the explanation of each TF unique?
    tf_multi: dict[int, bool] = {}
    for f in model.regions:
        if len(path_sigs[f]) > 1:
            tf_multi[f] = True
            continue
        sig = next(iter(path_sigs[f]))
        act = model.new_var()
        model.add([-act] + [-model.p[(f, k)] for k in sig] + [var for k, var in model.p_by_f[f] if k not in sig])
        res = model.solve(obs + [act])
        if res is None:
            return {"status": STATUS_INCOMPLETE, "partial_edges": sorted(poss)}
        if res:
            record(model.model)
        tf_multi[f] = bool(res)
        model.add([-act])
    # ---- nodes: possible states and necessity
    nodes = list(model.u)
    wgoals = {(v, st): model.node_state_var(v, st) for v in nodes for st in (1, -1)}
    cov = cover(model, obs, wgoals, record)
    if model.incomplete:
        return {"status": STATUS_INCOMPLETE, "partial_edges": sorted(poss)}
    ugoals = {v: model.u[v] for v in nodes if any((v, st) in cov for st in (1, -1))}
    nec_nodes = shrink(model, obs, ugoals, {v for v in ugoals if ugoals[v] in first}, record)
    if model.incomplete:
        return {"status": STATUS_INCOMPLETE, "partial_edges": sorted(poss)}
    node_info = {}
    for v in nodes:
        states = {st for st in (1, -1) if (v, st) in cov}
        if states:
            node_info[v] = (states, "Necessary" if v in nec_nodes else "Alternative")
    return {"status": STATUS_SAT, "edge_status": edge_status, "tf_multi": tf_multi, "node_info": node_info,
            "witness": witness, "cycle_clauses": model.n_cycle_clauses}


def _portfolio_model(G, root_idx, obs, conflict_budget, order_encoding, time_limit, fixed=None, max_path_length=None):
    """Satisfiability check with alternating solvers (maplechrono, glucose4, then cadical153 without interruption). On the
    real PKNs tried no solver was best everywhere: each one failed within 60 s on some formula another solved in seconds. Returns the model whose solver finished, or None."""
    total = float("inf") if time_limit is None else float(time_limit)
    start = time.monotonic()
    end = start + total
    schedule = [("maplechrono", 10), ("glucose4", 10), ("maplechrono", 30), ("cadical153", 1e9)]  # cadical cannot be interrupted: its time limit is soft
    for name, sl in schedule:
        rem = end - time.monotonic()
        if rem <= 0:
            break
        model = SampleModel(G, root_idx, obs, solver=name, conflict_budget=conflict_budget, order_encoding=order_encoding,
                            fixed_states=fixed, max_path_length=max_path_length)
        model.set_time_limit(min(sl, end - time.monotonic()))
        r = model.solve(model.obs_lits())
        if r is not None:
            model._auto_deadline = end
            model.chosen_solver = name
            model.set_time_limit(None)
            return model
    return None


def check_only(model: SampleModel) -> dict:
    """Satisfiability only (no edge classification): status, conflict diagnosis, one witness."""
    obs = model.obs_lits()
    r = model.solve(obs)
    if r is None:
        return {"status": STATUS_INCOMPLETE}
    if r is False:
        return {"status": STATUS_UNSAT, "conflict": model.conflicting_observations()}
    return {"status": STATUS_SAT, "unclassified": True, "witness": {f: model.path_nodes(model.model, f) for f in model.regions}}


def infer(tf_scores=None, pkn=None, *, tf_states=None, roots: list[str] | None = None, solver: str = "auto",
          conflict_budget: int | None = None, order_encoding: bool = True, classify: bool = True,
          time_limit: float | None = None, root_states: dict[str, int] | None = None,
          max_path_length: int | None = None) -> Result:
    """Reconstruct, for every sample, the set of causal networks that explain all observed TFs at once."""
    if pkn is None:
        raise ValueError("pkn is required")
    P: PKN = read_pkn(pkn)
    states = encode_tf_input(tf_scores, tf_states)
    G = Graph.from_pkn(P)
    if roots is None:
        root_idx = G.structural_roots()
        roots_source = "structural (source strongly connected components of the PKN)"
    else:
        names = [normalize_gene(r) for r in roots]
        unknown = [r for r in names if r not in G.index]
        if unknown:
            raise ValueError(f"roots not in the PKN: {unknown}")
        root_idx = sorted(G.index[r] for r in names)
        roots_source = "user-supplied"
    tf_names = list(states.columns)
    fixed = None
    if root_states:
        bad = [g for g in root_states if normalize_gene(g) not in G.index]
        if bad:
            raise ValueError(f"root_states genes not in the PKN: {bad}")
        fixed = {G.index[normalize_gene(g)]: int(v) for g, v in root_states.items()}
        if any(v not in (1, -1) for v in fixed.values()):
            raise ValueError("root_states values must be +1 or -1")
    unmapped = [t for t in tf_names if t not in G.index]
    result = Result(pkn=P, tf_states=states, roots=[G.nodes[i] for i in root_idx], roots_source=roots_source,
                    solver=solver, unmapped_tfs=unmapped,
                    input_form="tf_states" if tf_states is not None else "tf_scores")
    cache: dict[tuple, SampleNetwork] = {}
    for sample, row in states.iterrows():
        obs_named = {t: int(v) for t, v in row.items() if not np.isnan(v) and t in G.index}
        dropped = [t for t, v in row.items() if not np.isnan(v) and t not in G.index]
        key = tuple(sorted(obs_named.items()))
        if key in cache:
            base = cache[key]
            net = base.copy_for(sample, dropped)
            net.pkn_edges = P.edges
            result.add(net)
            continue
        net = _infer_one(sample, obs_named, dropped, G, root_idx, solver, conflict_budget, order_encoding, classify, time_limit, fixed, max_path_length)
        net.pkn_edges = P.edges
        cache[key] = net
        result.add(net)
    return result


def _infer_one(sample, obs_named: dict[str, int], dropped: list[str], G: Graph, root_idx: list[int], solver: str,
               conflict_budget, order_encoding: bool = True, classify: bool = True,
               time_limit: float | None = None, fixed=None, max_path_length: int | None = None) -> SampleNetwork:
    diag = [{"sample": sample, "level": "info", "code": "TF_NOT_IN_PKN", "detail": t} for t in dropped]
    obs = {G.index[t]: s for t, s in obs_named.items()}
    net = SampleNetwork(sample=sample, status=STATUS_NO_OBS, graph=G, roots=root_idx, obs=obs_named)
    if not obs:
        net.diagnostics = pd.DataFrame(diag, columns=["sample", "level", "code", "detail"])
        net._finish_empty()
        return net
    if solver == "auto":
        model = _portfolio_model(G, root_idx, obs, conflict_budget, order_encoding, time_limit, fixed, max_path_length)
        if model is not None and classify and model.chosen_solver != "cadical153":
            # classification issues thousands of incremental calls; on the real PKNs cadical153 proved about 3x more
            # edges per minute than maplechrono. The satisfiability result is already known, so rebuild with it.
            end = model._auto_deadline
            model = SampleModel(G, root_idx, obs, solver="cadical153", conflict_budget=conflict_budget, order_encoding=order_encoding,
                                fixed_states=fixed, max_path_length=max_path_length)
            model._auto_deadline = end
            model.chosen_solver = "cadical153"
        time_limit = None if model is None or time_limit is None else max(0.0, model._auto_deadline - time.monotonic())
        if model is None:
            net.status = STATUS_INCOMPLETE
            diag.append({"sample": sample, "level": "error", "code": "INCOMPLETE", "detail": "no solver finished the satisfiability check in the time limit"})
            net.diagnostics = pd.DataFrame(diag, columns=["sample", "level", "code", "detail"])
            net._finish_empty()
            return net
    else:
        model = SampleModel(G, root_idx, obs, solver=solver, conflict_budget=conflict_budget, order_encoding=order_encoding,
                            fixed_states=fixed, max_path_length=max_path_length)
    if model.unreachable:
        net.status = STATUS_UNREACHABLE
        for f in model.unreachable:
            diag.append({"sample": sample, "level": "error", "code": "UNREACHABLE_TF",
                         "detail": f"{G.nodes[f]} is not reachable from any candidate root"})
    else:
        model.set_time_limit(time_limit)
        out = analyze_sample(model, G, sample, {}) if classify else check_only(model)
        model.set_time_limit(None)
        net.status = out["status"]
        if net.status == STATUS_UNSAT:
            diag.append({"sample": sample, "level": "error", "code": "UNSAT_CONFLICT",
                         "detail": "minimal conflicting TF observations: " + ",".join(G.nodes[f] + _sgn(obs[f]) for f in out["conflict"])})
        elif net.status == STATUS_INCOMPLETE:
            part = out.get("partial_edges")
            diag.append({"sample": sample, "level": "error", "code": "INCOMPLETE",
                         "detail": "solver budget/time exhausted" + ("" if part is None else f"; {len(part)} edges proven possible, classification unfinished")})
            if part:
                net._fill_partial(part)
        else:
            net._fill(model, out)
    net.model = model
    net.diagnostics = pd.DataFrame(diag, columns=["sample", "level", "code", "detail"])
    if net.status != STATUS_SAT:
        net._finish_empty()
    return net
