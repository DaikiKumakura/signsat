"""Directed-graph helpers: indexing, strongly connected components, candidate roots, reachability."""
from __future__ import annotations

from collections import deque
from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .pkn import PKN


@dataclass
class Graph:
    nodes: list[str]
    index: dict[str, int]
    src: np.ndarray  # edge source node index
    dst: np.ndarray  # edge target node index
    sign: np.ndarray  # +1/-1
    edge_ids: list[str]
    out_edges: list[list[int]]  # node -> edge indices
    in_edges: list[list[int]]

    @classmethod
    def from_pkn(cls, pkn: PKN) -> "Graph":
        nodes = list(pkn.nodes)
        index = {n: i for i, n in enumerate(nodes)}
        e = pkn.edges
        src = np.array([index[s] for s in e["source"]], dtype=np.int64)
        dst = np.array([index[t] for t in e["target"]], dtype=np.int64)
        out_edges = [[] for _ in nodes]
        in_edges = [[] for _ in nodes]
        for k, (u, v) in enumerate(zip(src, dst)):
            out_edges[u].append(k)
            in_edges[v].append(k)
        return cls(nodes, index, src, dst, e["sign"].to_numpy(dtype=np.int64), e["edge_id"].tolist(), out_edges, in_edges)

    @property
    def n(self) -> int:
        return len(self.nodes)

    # ------------------------------------------------------------------ roots
    def structural_roots(self) -> list[int]:
        """Nodes of the strongly connected components that receive no edge from another component."""
        n = self.n
        if n == 0:
            return []
        m = csr_matrix((np.ones(len(self.src)), (self.src, self.dst)), shape=(n, n))
        k, label = connected_components(m, directed=True, connection="strong")
        has_in = np.zeros(k, dtype=bool)
        for u, v in zip(self.src, self.dst):
            if label[u] != label[v]:
                has_in[label[v]] = True
        return [i for i in range(n) if not has_in[label[i]]]

    # ------------------------------------------------------------ reachability
    def forward_reachable(self, starts: list[int]) -> np.ndarray:
        seen = np.zeros(self.n, dtype=bool)
        dq = deque(starts)
        for s in starts:
            seen[s] = True
        while dq:
            u = dq.popleft()
            for k in self.out_edges[u]:
                v = int(self.dst[k])
                if not seen[v]:
                    seen[v] = True
                    dq.append(v)
        return seen

    def distances(self, starts: list[int], forward: bool = True) -> np.ndarray:
        """Unweighted BFS distance from the start nodes (along edges if forward, against them otherwise); large if unreachable."""
        INF = 10**9
        dist = np.full(self.n, INF, dtype=np.int64)
        dq = deque()
        for s in starts:
            dist[s] = 0
            dq.append(s)
        adj = self.out_edges if forward else self.in_edges
        nxt = self.dst if forward else self.src
        while dq:
            u = dq.popleft()
            for k in adj[u]:
                v = int(nxt[k])
                if dist[v] == INF:
                    dist[v] = dist[u] + 1
                    dq.append(v)
        return dist

    def backward_reachable(self, target: int) -> np.ndarray:
        seen = np.zeros(self.n, dtype=bool)
        seen[target] = True
        dq = deque([target])
        while dq:
            v = dq.popleft()
            for k in self.in_edges[v]:
                u = int(self.src[k])
                if not seen[u]:
                    seen[u] = True
                    dq.append(u)
        return seen
