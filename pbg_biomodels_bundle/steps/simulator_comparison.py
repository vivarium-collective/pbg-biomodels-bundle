"""``SimulatorComparisonStep`` — compare two simulator trajectories at runtime.

Inputs are two ``numeric_result`` payloads (the shape emitted by the COPASI
and Tellurium UTC steps in ``pbsim_common``):

    {"time": list[float], "columns": list[str], "values": list[list[float]]}

The Step delegates to ``pbg_biomodels_bundle.comparison.compare_two_engines``
so the math is shared with the post-hoc HTML report builder in
``analysis.py``. Output ``comparison`` is a free-form dict containing per-
species RMSE, normalized RMSE, mean nRMSE, and a coarse quality bucket.
"""
from __future__ import annotations

from typing import Any, Dict

from process_bigraph import Step

from pbg_biomodels_bundle.comparison import compare_two_engines


class SimulatorComparisonStep(Step):
    """Compare two simulator trajectories species-by-species.

    Configure with the engine names you want labelled in the output dict.
    Defaults match the canonical pair used by the bundle's regression flow.
    """

    config_schema = {
        "engine_a_name": {"_type": "string", "_default": "copasi"},
        "engine_b_name": {"_type": "string", "_default": "tellurium"},
    }

    def inputs(self) -> Dict[str, str]:
        return {
            "engine_a_result": "numeric_result",
            "engine_b_result": "numeric_result",
        }

    def outputs(self) -> Dict[str, Any]:
        # Concrete struct shape matches the dict returned by
        # `pbg_biomodels_bundle.comparison.compare_two_engines`. Declared
        # inline (rather than registered as a named type) so the bundle has
        # no module-init order requirement on workspace `register_types`.
        return {
            "comparison": {
                "n_shared":         "integer",
                "rmse_by_species":  "map[float]",
                "nrmse_by_species": "map[float]",
                "mean_nrmse":       "maybe[float]",
                "bucket":           "string",
                "bucket_label":     "string",
            }
        }

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        summary = compare_two_engines(
            engine_a=state.get("engine_a_result"),
            engine_b=state.get("engine_b_result"),
            name_a=self.config["engine_a_name"],
            name_b=self.config["engine_b_name"],
        )
        return {"comparison": summary}
