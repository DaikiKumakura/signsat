"""SAT encoding of one sample (SPEC.md section 3).

Variables
  z[v]      node v has state +1 (True) or -1 (False)
  p[f,k]    the explanation path of observed TF f uses edge k
  x[k]      edge k is used by at least one explanation path   (x_k <-> OR_f p[f,k])
  u[v]      node v lies on at least one explanation path (through an edge)

Per observed TF f (a single simple path from a root node to f; trivial if f is itself a root)
  in-degree <= 1 and out-degree <= 1 on every node, f has no used out-edge,
  non-root nodes other than f: in <-> out, root nodes: in -> out (a path may start there),
  f has an in-edge unless f is a root (empty path).
These force "one path plus possibly disjoint cycles". Disjoint cycles are excluded by an order encoding: inside every
nontrivial strongly connected component of the TF's region each node gets a binary position and a used edge must
increase it strictly (positions 0..|SCC|-1 always suffice for a simple path). As a safety net, if a model still
contains a cycle a blocking clause is added; it is implied by the definition of a simple path, so no explanation is lost.

Edge signs: p[f,(u,v)] -> z[v] = sign * z[u]; observations: z[f] = y_f (passed as assumptions so that an
UNSAT core names the conflicting observations).
"""
from __future__ import annotations

import threading
import time
from dataclasses import dataclass, field

import numpy as np
from pysat.solvers import Solver
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .graph import Graph


class ModelBits:
    """Truth values of a model as a numpy bool array indexed by variable id (variables created later are False)."""

    __slots__ = ("arr",)

    def __init__(self, lits, nvar: int):
        m = np.fromiter(lits, dtype=np.int64, count=len(lits))
        arr = np.zeros(max(nvar, int(np.abs(m).max()) if len(m) else 0) + 1, dtype=bool)
        arr[m[m > 0]] = True
        self.arr = arr

    def __contains__(self, var: int) -> bool:
        return 0 < var < len(self.arr) and bool(self.arr[var])

    def get(self, var, default=False):
        return var in self

    def many(self, vars_) -> np.ndarray:
        idx = np.asarray(vars_, dtype=np.int64)
        out = np.zeros(len(idx), dtype=bool)
        ok = idx < len(self.arr)
        out[ok] = self.arr[idx[ok]]
        return out


@dataclass
class Region:
    nodes: set[int]
    edges: list[int]  # edge indices usable by this TF's path (no self loops)


class SampleModel:
    def __init__(self, graph: Graph, roots: list[int], obs: dict[int, int], *, solver: str = "cadical153",
                 conflict_budget: int | None = None, order_encoding: bool = True):
        self.g = graph
        self.roots = set(roots)
        self.obs = dict(obs)  # node index -> +1/-1
        self.solver_name = solver
        self.conflict_budget = conflict_budget
        self.order_encoding = order_encoding
        self.deadline: float | None = None  # time.monotonic() value after which solving is interrupted
        self.nvar = 0
        self.clauses: list[list[int]] = []
        self.z: dict[int, int] = {}
        self.p: dict[tuple[int, int], int] = {}
        self.x: dict[int, int] = {}
        self.u: dict[int, int] = {}
        self.inv: dict[tuple[int, int], int] = {}
        self.outv: dict[tuple[int, int], int] = {}
        self.regions: dict[int, Region] = {}
        self.unreachable: list[int] = []
        self.p_by_f: dict[int, list[tuple[int, int]]] = {}
        self._wvars: dict[tuple[int, int], int] = {}
        self.n_cycle_clauses = 0
        self.incomplete = False
        self._solver = None
        self._build()

    # ------------------------------------------------------------------ helpers
    def new_var(self) -> int:
        self.nvar += 1
        return self.nvar

    def add(self, clause: list[int]) -> None:
        self.clauses.append(clause)
        if self._solver is not None:
            self._solver.add_clause(clause)

    def zvar(self, v: int) -> int:
        if v not in self.z:
            self.z[v] = self.new_var()
        return self.z[v]

    def obs_lit(self, f: int) -> int:
        return self.zvar(f) if self.obs[f] == 1 else -self.zvar(f)

    def obs_lits(self, subset=None) -> list[int]:
        return [self.obs_lit(f) for f in (self.obs if subset is None else subset)]

    def _amo(self, lits: list[int]) -> None:
        n = len(lits)
        if n <= 1:
            return
        if n <= 5:
            for i in range(n):
                for j in range(i + 1, n):
                    self.add([-lits[i], -lits[j]])
            return
        s = [self.new_var() for _ in range(n - 1)]
        self.add([-lits[0], s[0]])
        for i in range(1, n - 1):
            self.add([-lits[i], s[i]])
            self.add([-s[i - 1], s[i]])
            self.add([-lits[i], -s[i - 1]])
        self.add([-lits[n - 1], -s[n - 2]])

    def _or_def(self, var: int, lits: list[int]) -> None:
        """var <-> OR(lits)  (var False if lits is empty)."""
        if not lits:
            self.add([-var])
            return
        self.add([-var] + lits)
        for l in lits:
            self.add([var, -l])

    def and_var(self, a: int, b: int) -> int:
        c = self.new_var()
        self.add([-c, a])
        self.add([-c, b])
        self.add([c, -a, -b])
        return c

    # -------------------------------------------------------------------- build
    def _build(self) -> None:
        g = self.g
        fwd = g.forward_reachable(sorted(self.roots))
        for f in self.obs:
            self.zvar(f)
            if not fwd[f]:
                self.unreachable.append(f)
                continue
            back = g.backward_reachable(f)
            nodes = {v for v in range(g.n) if fwd[v] and back[v]}
            edges = [k for k in range(len(g.src)) if g.src[k] != g.dst[k] and int(g.src[k]) in nodes and int(g.dst[k]) in nodes]
            self.regions[f] = Region(nodes, edges)
        if self.unreachable:
            return
        for f, reg in self.regions.items():
            eset = set(reg.edges)
            for k in reg.edges:
                self.p[(f, k)] = self.new_var()
            self.p_by_f[f] = [(k, self.p[(f, k)]) for k in reg.edges]
            for v in reg.nodes:
                ins = [self.p[(f, k)] for k in g.in_edges[v] if k in eset]
                outs = [self.p[(f, k)] for k in g.out_edges[v] if k in eset]
                self._amo(ins)
                self._amo(outs)
                iv, ov = self.new_var(), self.new_var()
                self.inv[(f, v)], self.outv[(f, v)] = iv, ov
                self._or_def(iv, ins)
                self._or_def(ov, outs)
                if v == f:
                    self.add([-ov])
                    if f not in self.roots:
                        self.add([iv])
                elif v in self.roots:
                    self.add([-iv, ov])
                else:
                    self.add([-iv, ov])
                    self.add([iv, -ov])
            for k in reg.edges:
                a, b = int(g.src[k]), int(g.dst[k])
                za, zb, pk = self.zvar(a), self.zvar(b), self.p[(f, k)]
                if g.sign[k] == 1:
                    self.add([-pk, -za, zb])
                    self.add([-pk, za, -zb])
                else:
                    self.add([-pk, za, zb])
                    self.add([-pk, -za, -zb])
            if self.order_encoding:
                self._order_encoding(f, reg)
            self._two_cycles(f, reg)
        all_edges = sorted({k for reg in self.regions.values() for k in reg.edges})
        for k in all_edges:
            xk = self.new_var()
            self.x[k] = xk
            self._or_def(xk, [self.p[(f, k)] for f in self.regions if (f, k) in self.p])
        nodes_used = {int(g.src[k]) for k in all_edges} | {int(g.dst[k]) for k in all_edges}
        for v in nodes_used:
            uv = self.new_var()
            self.u[v] = uv
            inc = [self.x[k] for k in g.out_edges[v] + g.in_edges[v] if k in self.x]
            self._or_def(uv, inc)


    def _two_cycles(self, f: int, reg: Region) -> None:
        """Forbid u->v together with v->u for one TF (a simple path never uses both)."""
        g = self.g
        by_pair: dict[tuple[int, int], list[int]] = {}
        for k in reg.edges:
            by_pair.setdefault((int(g.src[k]), int(g.dst[k])), []).append(k)
        for (a, b), ks in by_pair.items():
            if a < b and (b, a) in by_pair:
                for k1 in ks:
                    for k2 in by_pair[(b, a)]:
                        self.add([-self.p[(f, k1)], -self.p[(f, k2)]])

    def _order_encoding(self, f: int, reg: Region) -> None:
        """Forbid cycles among the edges used by TF f: strictly increasing binary positions inside each nontrivial SCC."""
        g = self.g
        nodes = sorted(reg.nodes)
        pos = {v: i for i, v in enumerate(nodes)}
        if not reg.edges:
            return
        r = np.array([pos[int(g.src[k])] for k in reg.edges])
        c = np.array([pos[int(g.dst[k])] for k in reg.edges])
        m = csr_matrix((np.ones(len(r)), (r, c)), shape=(len(nodes), len(nodes)))
        _, label = connected_components(m, directed=True, connection="strong")
        sizes = np.bincount(label)
        big = [v for v in nodes if sizes[label[pos[v]]] > 1]
        if not big:
            return
        nbits = max(1, (int(max(sizes)) - 1).bit_length())
        bits = {v: [self.new_var() for _ in range(nbits)] for v in big}
        for k in reg.edges:
            a, b = int(g.src[k]), int(g.dst[k])
            if label[pos[a]] != label[pos[b]] or sizes[label[pos[a]]] <= 1:
                continue
            pk = self.p[(f, k)]
            lt = [self.new_var() for _ in range(nbits)] + [None]
            self.add([-pk, lt[0]])
            for i in range(nbits):
                ai, bi = bits[a][i], bits[b][i]
                self.add([-lt[i], -ai, bi])  # lt_i -> not (a_i and not b_i)
                if lt[i + 1] is None:
                    self.add([-lt[i], -ai])
                    self.add([-lt[i], bi])
                else:
                    self.add([-lt[i], -ai, lt[i + 1]])
                    self.add([-lt[i], bi, lt[i + 1]])

    def node_state_var(self, v: int, st: int) -> int:
        """Variable w <-> (node v lies on an explanation path AND has state st). Created on first use."""
        key = (v, st)
        if key not in self._wvars:
            w = self.new_var()
            zl = self.zvar(v) if st == 1 else -self.zvar(v)
            uv = self.u[v]
            self.add([-w, uv])
            self.add([-w, zl])
            self.add([w, -uv, -zl])
            self._wvars[key] = w
        return self._wvars[key]

    # -------------------------------------------------------------------- solve
    @property
    def solver(self):
        if self._solver is None:
            self._solver = self.fresh_solver()
        return self._solver

    def fresh_solver(self, name: str | None = None):
        name = name or self.solver_name
        s = Solver(name=name, bootstrap_with=self.clauses)
        s.signsat_name = name
        try:  # prefer 'edge not used': keeps spurious cycles rare (the result does not depend on this)
            s.set_phases([-v for v in self.p.values()] + [-v for v in self.x.values()])
        except (NotImplementedError, AttributeError):
            pass
        return s

    def set_time_limit(self, seconds: float | None) -> None:
        self.deadline = None if seconds is None else time.monotonic() + float(seconds)
        self.incomplete = False

    def solve(self, assumptions: list[int], solver=None) -> bool | None:
        """SAT with lazy removal of disjoint cycles. True / False / None (budget exhausted)."""
        s = solver or self.solver
        sname = getattr(s, "signsat_name", self.solver_name)
        while True:
            if self.deadline is not None and sname.startswith("cadical"):
                # CaDiCaL in PySAT cannot be interrupted: the time limit is checked between calls only (soft limit)
                if self.deadline - time.monotonic() <= 0:
                    self.incomplete = True
                    return None
                r = s.solve(assumptions=assumptions)
            elif self.deadline is not None:
                rem = self.deadline - time.monotonic()
                if rem <= 0:
                    self.incomplete = True
                    return None
                timer = threading.Timer(rem, s.interrupt)
                timer.start()
                try:
                    r = s.solve_limited(assumptions=assumptions, expect_interrupt=True)
                finally:
                    timer.cancel()
                    s.clear_interrupt()
            elif self.conflict_budget is not None and not sname.startswith("cadical"):
                s.conf_budget(self.conflict_budget)
                r = s.solve_limited(assumptions=assumptions)
            else:
                r = s.solve(assumptions=assumptions)
            if r is None:
                self.incomplete = True
                return None
            if not r:
                return False
            model = self._model_dict(s.get_model())
            cycles = [] if self.order_encoding else self._find_cycles(model)
            if not cycles:
                self.model = model
                return True
            for cl in cycles:
                self.clauses.append(cl)
                if self._solver is not None and self._solver is not s:
                    self._solver.add_clause(cl)
                s.add_clause(cl)
                self.n_cycle_clauses += 1

    def _model_dict(self, model: list[int]) -> ModelBits:
        return ModelBits(model, self.nvar)

    def used_edges(self, model, f: int) -> list[int]:
        return [k for k, var in self.p_by_f[f] if var in model]

    def _find_cycles(self, model) -> list[list[int]]:
        g = self.g
        out = []
        for f in self.regions:
            used = self.used_edges(model, f)
            if not used:
                continue
            succ = {int(g.src[k]): k for k in used}
            has_in = {int(g.dst[k]) for k in used}
            starts = [int(g.src[k]) for k in used if int(g.src[k]) not in has_in]
            seen: set[int] = set()
            for s0 in starts:
                v = s0
                while v in succ:
                    k = succ[v]
                    seen.add(k)
                    v = int(g.dst[k])
            rest = [k for k in used if k not in seen]
            done: set[int] = set()
            for k in rest:
                if k in done:
                    continue
                cyc = [k]
                done.add(k)
                v = int(g.dst[k])
                while v in succ and succ[v] not in done:
                    kk = succ[v]
                    cyc.append(kk)
                    done.add(kk)
                    v = int(g.dst[kk])
                out.append([-self.p[(f, kk)] for kk in cyc])
                for f2, reg2 in self.regions.items():
                    if f2 != f and all((f2, kk) in self.p for kk in cyc):
                        out.append([-self.p[(f2, kk)] for kk in cyc])
        return out

    def path_nodes(self, model, f: int) -> list[int]:
        """Node sequence (root ... f) of the explanation of f in ``model``."""
        g = self.g
        used = self.used_edges(model, f)
        if not used:
            return [f]
        succ = {int(g.src[k]): k for k in used}
        has_in = {int(g.dst[k]) for k in used}
        start = next(int(g.src[k]) for k in used if int(g.src[k]) not in has_in)
        seq, v = [start], start
        while v in succ:
            v = int(g.dst[succ[v]])
            seq.append(v)
        return seq

    # ----------------------------------------------------------- UNSAT diagnosis
    def conflicting_observations(self) -> list[int]:
        """A minimal (deletion-based) subset of observed TFs that is already unsatisfiable."""
        s = self.solver
        lits = self.obs_lits()
        lit2f = {l: f for l, f in zip(lits, self.obs)}
        if self.solve(lits) is not False:
            return []
        core = [l for l in s.get_core() if l in lit2f]
        keep = list(core)
        for l in list(core):
            trial = [m for m in keep if m != l]
            if self.solve(trial) is False:
                keep = trial
        return sorted(lit2f[l] for l in keep)


# ---------------------------------------------------------------------------------------------- set-cover loops
def _add(model: SampleModel, solver, clause: list[int]) -> None:
    """Add to the private query solver if given, otherwise to the model (and its incremental solver)."""
    if solver is not None:
        solver.add_clause(clause)
    else:
        model.add(clause)


def _set_phase(solver, lits):
    try:
        solver.set_phases(lits)
    except (NotImplementedError, AttributeError):
        pass


def cover(model: SampleModel, base: list[int], goals: dict, on_model=None, solver=None):
    """Keys whose goal variable is True in at least one model (under ``base`` assumptions).

    Repeatedly asks for a model that makes at least one still-uncovered goal true (clause guarded by a fresh
    activation variable), so the number of SAT calls is the number of models needed to cover the goals, not the
    number of goals. The final UNSAT call proves the remaining goals can never be true. Returns None if the
    solver budget is exhausted.
    """
    s = solver or model.solver
    remaining = dict(goals)
    covered: set = set()
    while remaining:
        act = model.new_var()
        cl = [-act] + list(remaining.values())
        _add(model, solver, cl)
        r = model.solve(base + [act], solver=s)
        _add(model, solver, [-act])
        if r is None:  # budget/time exhausted: what was covered so far is proven, the rest is unresolved
            return covered
        if not r:
            break
        m = model.model
        if on_model:
            on_model(m)
        found = [k for k, v in remaining.items() if v in m]
        covered.update(found)
        for k in found:
            del remaining[k]
    return covered


def shrink(model: SampleModel, base: list[int], goals: dict, start: set, on_model=None, solver=None):
    """Keys whose goal variable is True in EVERY model: start from ``start`` (true in one model) and keep
    asking for a model that falsifies at least one remaining key. Returns None if the budget is exhausted."""
    s = solver or model.solver
    cand = set(start)
    while cand:
        act = model.new_var()
        cl = [-act] + [-goals[k] for k in cand]
        _add(model, solver, cl)
        r = model.solve(base + [act], solver=s)
        _add(model, solver, [-act])
        if r is None:
            return cand  # unresolved: model.incomplete is True and the caller must not treat cand as necessary
        if not r:
            break
        m = model.model
        if on_model:
            on_model(m)
        cand = {k for k in cand if goals[k] in m}
    return cand
