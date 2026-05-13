"""``LoadBiomodelStep`` — fetch a BioModel by ID and emit SBML path + UTC spec.

Reads ``biomodel_id`` from its single input port, calls
:func:`pbg_biomodels_bundle.run_biomodels.load_biomodel` (which queries the
``biomodels`` REST API, caches the SBML + SED-ML files locally under
``models/<id>/``, and parses the first UniformTimeCourse simulation out of
the SED-ML), and writes the resolved SBML path, simulation duration, and
output-point count out on three output ports.

Designed to be the entry point of a single-composite pipeline:

    biomodel_id (store)
        ↓ LoadBiomodelStep
    sbml_path / time / n_points (stores)
        ↓ LocalCopasiUTCStep + LocalTelluriumUTCStep (in parallel)
    results[<engine>] (stores)
        ↓ SimulatorComparisonStep / CompareOverlay
    comparison + viz_html
"""
from __future__ import annotations

from typing import Any, ClassVar, Dict

from process_bigraph import Step


class LoadBiomodelStep(Step):
    """Resolve a BioModels identifier to a local SBML file and UTC spec.

    Inputs:
        biomodel_id: BioModels identifier, e.g. ``"BIOMD0000000001"``.

    Outputs:
        sbml_path: Absolute path to the cached SBML XML file.
        time: SED-ML UniformTimeCourse duration (``output_end_time
              - output_start_time``).
        n_points: SED-ML UniformTimeCourse number_of_points.

    Side effects:
        Caches the SBML + SED-ML files under ``models/<biomodel_id>/`` in
        the workspace (or current working directory) — same convention
        the bundle's ``run_biomodels`` uses, so the cache is shared.
    """

    config_schema: ClassVar[Dict[str, str]] = {}

    def inputs(self) -> Dict[str, str]:
        return {"biomodel_id": "string"}

    def outputs(self) -> Dict[str, str]:
        return {
            "sbml_path": "string",
            "time": "float",
            "n_points": "integer",
        }

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        # Lazy import so loading this module doesn't pull in biomodels'
        # heavy chain (pooch, pydantic-xml, …) when only the type is being
        # inspected (e.g. dashboard registry introspection).
        import os

        import biomodels

        from pbg_biomodels_bundle.run_biomodels import load_biomodel

        biomodel_id = state.get("biomodel_id") or ""
        if not biomodel_id:
            raise ValueError(
                "LoadBiomodelStep: input port `biomodel_id` is empty; "
                "set it in the composite state before running."
            )

        meta = biomodels.get_metadata(biomodel_id)
        result = load_biomodel(biomodel_id, meta)
        return {
            "sbml_path": os.path.abspath(result.sbml_path),
            "time":     float(result.utc.duration),
            "n_points": int(result.utc.number_of_points),
        }
