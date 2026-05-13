"""Local one-shot UTC simulator Steps — COPASI and Tellurium.

Wraps ``copasi-basico`` (the COPASI Python bindings) and ``tellurium``
directly, without going through ``pbsim_common``. The key behavioral
difference: these Steps declare ``inputs() == {}``, so process-bigraph's
runtime fires them on composite startup. ``pbsim_common``'s UTC Steps
wire their inputs to ``species_concentrations`` / ``species_counts``
stores that no upstream process writes to, and never trigger.

Both Steps load an SBML file from ``config['model_source']``, simulate
from time 0 for ``config['time']`` seconds at ``config['n_points']``
samples, and emit a single ``numeric_result``
(``{time, columns, values}``) on the ``result`` output port.

Designed for the ``build_compare_biomodel`` composite — one COPASI Step
and one Tellurium Step in parallel, feeding ``SimulatorComparisonStep``.
"""
from __future__ import annotations

from typing import Any, ClassVar, Dict, List

from process_bigraph import Step


class LocalCopasiUTCStep(Step):
    """Run an SBML model in COPASI for a fixed UniformTimeCourse window."""

    config_schema: ClassVar[Dict[str, Any]] = {
        "model_source": "string",
        "time": "float",
        "n_points": "integer",
    }

    def inputs(self) -> Dict[str, str]:
        # Empty — fires once on composite startup, no input dependencies.
        return {}

    def outputs(self) -> Dict[str, str]:
        return {"result": "numeric_result"}

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        # Lazy import so a broken COPASI install doesn't sink module import.
        from basico import load_model, run_time_course  # type: ignore

        n_points = int(self.config["n_points"])
        if n_points < 2:
            raise ValueError(f"LocalCopasiUTCStep: n_points must be >= 2, got {n_points}")

        dm = load_model(self.config["model_source"])
        tc = run_time_course(
            start_time=0.0,
            duration=float(self.config["time"]),
            intervals=n_points - 1,  # COPASI parameterizes by intervals
            update_model=True,
            use_sbml_id=True,
            model=dm,
        )
        return {
            "result": {
                "time": tc.index.to_list(),
                "columns": list(tc.columns),
                "values": tc.values.tolist(),
            }
        }


class LocalTelluriumUTCStep(Step):
    """Run an SBML model in Tellurium (libroadrunner) for a UTC window."""

    config_schema: ClassVar[Dict[str, Any]] = {
        "model_source": "string",
        "time": "float",
        "n_points": "integer",
    }

    def inputs(self) -> Dict[str, str]:
        return {}

    def outputs(self) -> Dict[str, str]:
        return {"result": "numeric_result"}

    def update(self, state: Dict[str, Any]) -> Dict[str, Any]:
        import tellurium as te  # type: ignore

        n_points = int(self.config["n_points"])
        if n_points < 2:
            raise ValueError(f"LocalTelluriumUTCStep: n_points must be >= 2, got {n_points}")

        rr = te.loadSBMLModel(self.config["model_source"])
        result = rr.simulate(0.0, float(self.config["time"]), n_points)

        col_names: List[str] = list(result.colnames)
        # Tellurium column names look like 'time', '[S1]', '[S2]'...
        # Strip the brackets so species names match SBML IDs (and COPASI's labels).
        normalized = [c.strip("[]") for c in col_names]
        # Pull out the time column; remaining columns are species values.
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
        # Rows of [s1_val, s2_val, ...] in species_names order.
        values = [
            [float(result[r, i]) for i, _ in species_cols]
            for r in range(result.shape[0])
        ]
        return {
            "result": {
                "time": times,
                "columns": species_names,
                "values": values,
            }
        }
