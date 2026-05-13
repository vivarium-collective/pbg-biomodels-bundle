"""Local one-shot UTC simulator Steps — COPASI and Tellurium.

Wraps ``copasi-basico`` and ``tellurium`` directly (no ``pbsim_common``
indirection). Both Steps consume three runtime inputs:

- ``model_source`` — absolute path to an SBML file.
- ``time``         — UniformTimeCourse duration (model-time seconds).
- ``n_points``     — number of output samples (incl. endpoints).

They emit a single ``numeric_result`` (``{time, columns, values}``) on
the ``result`` output port. Because the inputs are runtime ports
(not config), upstream Steps such as
:class:`pbg_biomodels_bundle.steps.load_biomodel.LoadBiomodelStep` can
populate them dynamically — supporting a true single-composite pipeline
where ``biomodel_id`` is the only external parameter.
"""
from __future__ import annotations

from typing import Any, ClassVar, Dict, List

from process_bigraph import Step


_UTC_INPUTS: Dict[str, str] = {
    "model_source": "string",
    "time":         "float",
    "n_points":     "integer",
}


def _resolve_n_points(n: Any, step: str) -> int:
    n = int(n)
    if n < 2:
        raise ValueError(f"{step}: n_points must be >= 2, got {n}")
    return n


class LocalCopasiUTCStep(Step):
    """Run an SBML model in COPASI for a fixed UniformTimeCourse window."""

    config_schema: ClassVar[Dict[str, Any]] = {}

    def inputs(self) -> Dict[str, str]:
        return dict(_UTC_INPUTS)

    def outputs(self) -> Dict[str, str]:
        return {"result": "numeric_result"}

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        from basico import load_model, run_time_course  # type: ignore

        model_source = state["model_source"]
        n_points = _resolve_n_points(state["n_points"], "LocalCopasiUTCStep")

        dm = load_model(model_source)
        tc = run_time_course(
            start_time=0.0,
            duration=float(state["time"]),
            intervals=n_points - 1,  # COPASI parameterizes by intervals
            update_model=True,
            use_sbml_id=True,
            model=dm,
        )
        return {
            "result": {
                "time":    tc.index.to_list(),
                "columns": list(tc.columns),
                "values":  tc.values.tolist(),
            }
        }


class LocalTelluriumUTCStep(Step):
    """Run an SBML model in Tellurium (libroadrunner) for a UTC window."""

    config_schema: ClassVar[Dict[str, Any]] = {}

    def inputs(self) -> Dict[str, str]:
        return dict(_UTC_INPUTS)

    def outputs(self) -> Dict[str, str]:
        return {"result": "numeric_result"}

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        import tellurium as te  # type: ignore

        model_source = state["model_source"]
        n_points = _resolve_n_points(state["n_points"], "LocalTelluriumUTCStep")

        rr = te.loadSBMLModel(model_source)
        result = rr.simulate(0.0, float(state["time"]), n_points)

        col_names: List[str] = list(result.colnames)
        # Strip the [...] wrapper so species names match SBML IDs (and COPASI).
        normalized = [c.strip("[]") for c in col_names]
        try:
            t_idx = normalized.index("time")
        except ValueError as exc:
            raise RuntimeError(
                f"LocalTelluriumUTCStep: no 'time' column in simulate() result; "
                f"got columns {col_names}"
            ) from exc

        times = [float(x) for x in result[:, t_idx]]
        species_cols = [(i, name) for i, name in enumerate(normalized) if i != t_idx]
        species_names = [name for _, name in species_cols]
        values = [
            [float(result[r, i]) for i, _ in species_cols]
            for r in range(result.shape[0])
        ]
        return {
            "result": {
                "time":    times,
                "columns": species_names,
                "values":  values,
            }
        }
