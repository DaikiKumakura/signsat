"""Signed prior-knowledge network (PKN): reading and fixed normalisation rules.

Rules (SPEC.md section 2):
* gene symbols are stripped and upper-cased;
* identical (source, target, sign) rows are merged and their evidence is kept;
* a (source, target) pair with both signs is kept as two candidate edges;
* rows with an unknown sign or a missing endpoint are not used and are reported;
* self loops are kept in the table but never used in an explanation (a simple path has no self loop).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

_POS = {"1", "+1", "+", "positive", "activation", "stimulation", "up", "increase", "true"}
_NEG = {"-1", "-", "negative", "inhibition", "repression", "down", "decrease", "false"}
REQUIRED = ("source", "target", "sign")


def normalize_gene(name) -> str | None:
    if name is None or (isinstance(name, float) and np.isnan(name)):
        return None
    s = str(name).strip().upper()
    return s or None


def parse_sign(value) -> int | None:
    """Return +1 / -1, or None when the sign is unknown."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, (int, np.integer, float, np.floating)):
        if value == 1:
            return 1
        if value == -1:
            return -1
        return None
    s = str(value).strip().lower()
    if s in _POS:
        return 1
    if s in _NEG:
        return -1
    return None


@dataclass
class PKN:
    """Normalised network. ``edges`` has one row per candidate edge (unique ``edge_id``)."""

    edges: pd.DataFrame
    diagnostics: pd.DataFrame
    nodes: list[str] = field(init=False)

    def __post_init__(self):
        self.nodes = sorted(set(self.edges["source"]) | set(self.edges["target"]))

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def to_csv(self, path) -> None:
        self.edges.to_csv(path, index=False)


def read_pkn(pkn, *, source_label: str | None = None) -> PKN:
    """Build a :class:`PKN` from a DataFrame or a CSV path with columns source, target, sign."""
    if isinstance(pkn, PKN):
        return pkn
    if isinstance(pkn, (str, Path)):
        df = pd.read_csv(pkn, dtype=str, keep_default_na=False, na_values=[""])
        source_label = source_label or str(pkn)
    else:
        df = pd.DataFrame(pkn).copy()
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"PKN is missing required column(s): {missing}")
    extra = [c for c in df.columns if c not in REQUIRED]
    diag: list[dict] = []
    rows = []
    for i, r in df.iterrows():
        s, t, sg = normalize_gene(r["source"]), normalize_gene(r["target"]), parse_sign(r["sign"])
        if s is None or t is None:
            diag.append({"row": i, "reason": "MISSING_ENDPOINT", "detail": f"{r['source']}->{r['target']}"})
            continue
        if sg is None:
            diag.append({"row": i, "reason": "UNKNOWN_SIGN", "detail": f"{s}->{t} sign={r['sign']!r}"})
            continue
        rows.append((i, s, t, sg, {c: r[c] for c in extra}))
    # merge identical (source, target, sign); keep all evidence
    merged: dict[tuple, dict] = {}
    for i, s, t, sg, ex in rows:
        key = (s, t, sg)
        if key not in merged:
            merged[key] = {"source": s, "target": t, "sign": sg, "_ex": [ex], "_row": i}
        else:
            merged[key]["_ex"].append(ex)
            diag.append({"row": i, "reason": "DUPLICATE_MERGED", "detail": f"{s}->{t} sign={sg:+d}"})
    out = []
    used_ids: set[str] = set()
    for k, (key, m) in enumerate(sorted(merged.items(), key=lambda kv: kv[1]["_row"]), start=1):
        rec = {"source": m["source"], "target": m["target"], "sign": m["sign"]}
        for c in extra:
            vals = [str(e[c]) for e in m["_ex"] if e[c] is not None and str(e[c]) != "nan"]
            uniq = list(dict.fromkeys(vals))
            rec[c] = ";".join(uniq) if uniq else None
        out.append(rec)
    edges = pd.DataFrame(out, columns=["source", "target", "sign"] + extra)
    if "edge_id" in extra:
        ids = edges["edge_id"].tolist()
        if any(x is None for x in ids) or len(set(ids)) != len(ids) or any(";" in str(x) for x in ids):
            raise ValueError("edge_id must be present, unique, and not merged across duplicate rows")
    else:
        edges.insert(0, "edge_id", [f"e{k}" for k in range(1, len(edges) + 1)])
    for c in ("edge_id", "source", "target", "sign"):
        pass
    cols = ["edge_id", "source", "target", "sign"] + [c for c in edges.columns if c not in ("edge_id", "source", "target", "sign")]
    edges = edges[cols].reset_index(drop=True)
    edges["sign"] = edges["sign"].astype(int)
    d = pd.DataFrame(diag, columns=["row", "reason", "detail"])
    p = PKN(edges=edges, diagnostics=d)
    p.source_label = source_label  # type: ignore[attr-defined]
    return p
