"""Readers that turn public resource dumps into the (source, target, sign) table SignSAT expects.

The package does not download anything. Fetch the files yourself (commands are in docs/data.md) and pass the paths.
Rows that cannot be used are returned in a second table with the reason, never dropped silently.
"""
from __future__ import annotations

import pandas as pd


def omnipath_pkn(path_or_df, *, directed_only: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """OmniPath ``interactions`` table (genesymbols=1, tsv). A row that is both stimulation and inhibition yields
    two candidate edges (kept as alternatives); a directed row with neither sign is unusable (reported)."""
    d = pd.read_csv(path_or_df, sep="\t", dtype=str) if not isinstance(path_or_df, pd.DataFrame) else path_or_df.astype(str)
    for c in ("is_directed", "is_stimulation", "is_inhibition"):
        d[c] = d[c].str.lower().eq("true")
    keep = d[d["is_directed"]] if directed_only else d
    rows, bad = [], []
    ev_cols = [c for c in ("sources", "references", "curation_effort") if c in d.columns]
    for r in keep.itertuples(index=False):
        base = {"source": r.source_genesymbol, "target": r.target_genesymbol}
        for c in ev_cols:
            base[c] = getattr(r, c)
        signs = [s for s, on in ((1, r.is_stimulation), (-1, r.is_inhibition)) if on]
        if not signs:
            bad.append({**base, "reason": "UNKNOWN_SIGN"})
        for s in signs:
            rows.append({**base, "sign": s})
    out = pd.DataFrame(rows, columns=["source", "target", "sign"] + ev_cols)
    return out, pd.DataFrame(bad)


def signor_pkn(path_or_df, *, entity_types=("protein",)) -> tuple[pd.DataFrame, pd.DataFrame]:
    """SIGNOR ``getLatestRelease`` TSV. Only edges between the listed entity types are used (complexes, protein
    families, chemicals and phenotypes are reported as unsupported); the sign comes from EFFECT
    (``up-regulates*`` -> +1, ``down-regulates*`` -> -1, anything else unknown)."""
    if not isinstance(path_or_df, pd.DataFrame):
        # the dump contains stray quote characters inside fields: read without quote handling, then strip quotes
        d = pd.read_csv(path_or_df, sep="\t", dtype=str, quoting=3)
        d = d.apply(lambda c: c.str.strip('"') if c.dtype == object else c)
    else:
        d = path_or_df
    rows, bad = [], []
    for r in d.itertuples(index=False):
        a, b = getattr(r, "ENTITYA"), getattr(r, "ENTITYB")
        eff = str(getattr(r, "EFFECT")).lower()
        base = {"source": a, "target": b, "mechanism": getattr(r, "MECHANISM"), "reference": getattr(r, "PMID"),
                "signor_id": getattr(r, "SIGNOR_ID")}
        if getattr(r, "TYPEA") not in entity_types or getattr(r, "TYPEB") not in entity_types:
            bad.append({**base, "reason": f"UNSUPPORTED_ENTITY_TYPE:{getattr(r, 'TYPEA')}/{getattr(r, 'TYPEB')}"})
        elif "down-regulates" in eff:
            rows.append({**base, "sign": -1})
        elif "up-regulates" in eff:
            rows.append({**base, "sign": 1})
        else:
            bad.append({**base, "reason": f"UNKNOWN_SIGN:{getattr(r, 'EFFECT')}"})
    return pd.DataFrame(rows, columns=["source", "target", "sign", "mechanism", "reference", "signor_id"]), pd.DataFrame(bad)


def collectri_tfs(path_or_df) -> list[str]:
    """TF symbols (regulators) of a CollecTRI table obtained from OmniPath (datasets=collectri)."""
    d = pd.read_csv(path_or_df, sep="\t", dtype=str) if not isinstance(path_or_df, pd.DataFrame) else path_or_df
    return sorted(set(d["source_genesymbol"].str.upper()))
