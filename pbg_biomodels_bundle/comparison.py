"""Shared simulator-comparison math.

Used by both the post-hoc HTML report builder (``analysis.py``) and the
runtime :class:`SimulatorComparisonStep` (``steps/simulator_comparison.py``).

Inputs are per-engine ``numeric_result`` payloads — dicts of shape
``{time: list[float], columns: list[str], values: list[list[float]]}``,
emitted by the COPASI / Tellurium UTC steps.
"""
from __future__ import annotations

import math
from typing import Any, Dict, Iterable, List, Optional


# Normalized-RMSE thresholds (mean across shared species).
BUCKET_THRESHOLDS = [
    ("good",       0.01,     "Good (≤1%)"),
    ("borderline", 0.10,     "Borderline (1–10%)"),
    ("large",      math.inf, "Large diff (>10%)"),
]
NO_COMPARISON_BUCKET = ("none", "No comparison")


def union_species(engines: Dict[str, Any]) -> List[str]:
    """Return the ordered union of species columns across all engines."""
    species: List[str] = []
    seen: set = set()
    for payload in engines.values():
        if not payload:
            continue
        for c in payload.get("columns", []):
            if c not in seen:
                seen.add(c)
                species.append(c)
    return species


def compare_two_engines(
    engine_a: Optional[Dict[str, Any]],
    engine_b: Optional[Dict[str, Any]],
    name_a: str = "a",
    name_b: str = "b",
) -> Dict[str, Any]:
    """Compare two numeric_result payloads species-by-species.

    Computes RMSE and normalized RMSE (RMSE / max|y|) for every species that
    appears in BOTH engines' columns, then a mean nRMSE across shared species
    and a coarse bucket label.

    Returns the shape consumed by the HTML report's metadata card and the
    dashboard Step. Returns the ``none`` bucket if either engine is missing
    or no species are shared.
    """
    engines = {name_a: engine_a, name_b: engine_b}
    species = union_species(engines)

    rmse_by_species: Dict[str, float] = {}
    nrmse_by_species: Dict[str, float] = {}
    n_shared = 0

    for sp in species:
        ys_per_engine: Dict[str, List[float]] = {}
        for eng_name, payload in engines.items():
            if not payload:
                continue
            cols = payload.get("columns", [])
            if sp not in cols:
                continue
            j = cols.index(sp)
            ys_per_engine[eng_name] = [row[j] for row in payload.get("values", [])]
        if len(ys_per_engine) < 2:
            continue
        n_shared += 1
        keys = sorted(ys_per_engine.keys())
        y1, y2 = ys_per_engine[keys[0]], ys_per_engine[keys[1]]
        n = min(len(y1), len(y2))
        if n == 0:
            continue
        rmse = math.sqrt(sum((y1[k] - y2[k]) ** 2 for k in range(n)) / n)
        rmse_by_species[sp] = rmse
        denom = max((abs(v) for v in y1[:n]), default=0.0)
        denom = max(denom, max((abs(v) for v in y2[:n]), default=0.0))
        if denom > 0:
            nrmse_by_species[sp] = rmse / denom

    if nrmse_by_species:
        mean_nrmse: Optional[float] = sum(nrmse_by_species.values()) / len(nrmse_by_species)
    else:
        mean_nrmse = None

    if mean_nrmse is None:
        bucket_id, bucket_label = NO_COMPARISON_BUCKET
    else:
        for bid, threshold, label in BUCKET_THRESHOLDS:
            if mean_nrmse <= threshold:
                bucket_id, bucket_label = bid, label
                break

    return {
        "n_shared": n_shared,
        "rmse_by_species": rmse_by_species,
        "nrmse_by_species": nrmse_by_species,
        "mean_nrmse": mean_nrmse,
        "bucket": bucket_id,
        "bucket_label": bucket_label,
    }
