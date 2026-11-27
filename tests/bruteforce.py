"""Independent brute-force reference (no SAT): enumerate every simple root->TF path per observed TF, every
combination of paths, and every node-state assignment. Used only by tests on tiny graphs.

Definitions (SPEC.md): an explanation H is a union of one simple path per observed TF (trivial if the TF is a
root) such that some assignment z of all nodes gives z_f = y_f for observed f and z_v = sign * z_u on every
edge of H. Target states are additional constraints on z.
"""
from __future__ import annotations

import itertools

from signsat.graph import Graph


def simple_paths(G: Graph, roots: set[int], f: int, max_len: int | None = None):
    """All simple edge-paths from a root node to f (list of edge-index tuples; () if f is a root)."""
    out = []
    if f in roots:
        out.append(())

    def dfs(v, visited, edges):
        if v == f and edges:
            out.append(tuple(edges))
            return
        if max_len is not None and len(edges) >= max_len:
            return
        for k in G.out_edges[v]:
            w = int(G.dst[k])
            if w in visited:
                continue
            dfs(w, visited | {w}, edges + [k])

    for r in sorted(roots):
        if r == f:
            # a path that starts at a root and ends at f != trivial: must leave and return -> not simple; skip
            continue
        dfs(r, {r}, [])
    # paths ending at f that start at a root other than f; paths that start at f are impossible (simple)
    return out


def explanations(G: Graph, roots: set[int], obs: dict[int, int], targets: dict[int, int] | None = None, max_len: int | None = None):
    """Yield (paths_per_tf, edge_set, z_assignments) for consistent combinations."""
    per_tf = {f: simple_paths(G, roots, f, max_len) for f in obs}
    if any(len(v) == 0 for v in per_tf.values()):
        return []
    res = []
    fs = list(obs)
    for combo in itertools.product(*[per_tf[f] for f in fs]):
        edges = set(k for p in combo for k in p)
        zs = []
        for z in itertools.product((1, -1), repeat=G.n):
            if any(z[f] != s for f, s in obs.items()):
                continue
            if targets and any(z[t] != s for t, s in targets.items()):
                continue
            if all(z[int(G.dst[k])] == G.sign[k] * z[int(G.src[k])] for k in edges):
                zs.append(z)
        if zs:
            res.append((dict(zip(fs, combo)), edges, zs))
    return res


def classify(G, roots, obs, max_len=None):
    ex = explanations(G, roots, obs, None, max_len)
    if not ex:
        return None
    sets = [e for _, e, _ in ex]
    poss = set().union(*sets)
    nec = set.intersection(*sets)
    node_states = {}
    node_used_all = None
    for paths, edges, zs in ex:
        used = {int(G.src[k]) for k in edges} | {int(G.dst[k]) for k in edges}
        node_used_all = used if node_used_all is None else node_used_all & used
        for v in used:
            node_states.setdefault(v, set()).update(z[v] for z in zs)
    return {"possible": poss, "necessary": nec, "node_states": node_states, "node_necessary": node_used_all,
            "n_explanations": len(ex)}


def query(G, roots, obs, targets: dict[int, int], focus: set[int], apply_states=True):
    """Reference answer for edges of target/focus queries (semantics of query.py docstring)."""
    ex = explanations(G, roots, obs, targets if apply_states else None)
    edges = set()
    explained = set()
    for paths, eset, zs in ex:
        for f, pe in paths.items():
            nodes = _path_nodes(G, f, pe)
            if targets:
                idx = [i for i, v in enumerate(nodes) if v in targets]
                if not idx:
                    continue
                i0 = idx[0]
                # segment from the target gene; qualifies if (focus empty or a focus gene at/after it)
                for i in idx:
                    seg_nodes = nodes[i:]
                    if focus and not any(v in focus for v in seg_nodes):
                        continue
                    explained.add(f)
                    seg_edges = [k for k in pe if nodes.index(int(G.src[k])) >= i]
                    edges.update(seg_edges)
            else:
                if focus and not any(v in focus for v in nodes):
                    continue
                explained.add(f)
                edges.update(pe)
    return edges, explained


def _path_nodes(G, f, pe):
    if not pe:
        return [f]
    nodes = [int(G.src[pe[0]])]
    for k in pe:
        nodes.append(int(G.dst[k]))
    return nodes
