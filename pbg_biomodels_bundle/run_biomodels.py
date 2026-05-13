"""
Load BioModels entries, extract UniformTimeCourse settings from SED-ML,
resolve SBML source, and emit/run process-bigraph documents for UTC steps.
"""

from __future__ import annotations

import asyncio
import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

import biomodels
import libsbml
import libsedml


# Workaround for biomodels<=0.x: ``get_metadata`` opens the cached JSON without
# specifying an encoding, so non-ASCII metadata (e.g. BIOMD0000000192) crashes
# under any C locale. Re-read with utf-8 when the platform default fails.
_original_get_metadata = biomodels.get_metadata


def _utf8_safe_get_metadata(model_id, **kwargs):
    try:
        return _original_get_metadata(model_id, **kwargs)
    except UnicodeDecodeError:
        from biomodels import metadata as _bm_meta
        from biomodels.common import base_url, cache_path, fix_id, pooch
        mid = fix_id(model_id)
        path = pooch.retrieve(
            f"{base_url}/model/files/{mid}?format=json",
            known_hash=kwargs.get("known_hash"),
            fname=mid,
            path=Path(cache_path, "model", "files"),
            progressbar=kwargs.get("progress_bar", False),
        )
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        meta = _bm_meta.Metadata(
            model_id=mid, files=[*data["main"], *data["additional"]]
        )
        meta._file = Path(path)
        meta.files = [meta.files[0], *sorted(meta.files[1:], key=lambda x: x.name)]
        return meta


biomodels.get_metadata = _utf8_safe_get_metadata

from pbest.globals import get_loaded_core
from process_bigraph import allocate_core, Composite
from process_bigraph.emitter import add_emitter_to_composite, gather_emitter_results


# ----------------------------
# Data structures
# ----------------------------

@dataclass(frozen=True)
class UniformTimeCourseSpec:
    initial_time: float
    output_start_time: float
    output_end_time: float
    number_of_points: int

    @property
    def duration(self) -> float:
        return float(self.output_end_time - self.output_start_time)


@dataclass(frozen=True)
class BiomodelLoadResult:
    biomodel_id: str
    sbml_path: str
    sedml_path: str
    utc: UniformTimeCourseSpec
    units: Optional[Dict[str, Any]] = None  # {"time_unit": str|None, "species_units": dict[str, str|None]}


# ----------------------------
# Helpers: picking files
# ----------------------------

_SBML_RE = re.compile(r"\.(xml|sbml)$", re.IGNORECASE)
_SEDML_RE = re.compile(r"\.sedml$", re.IGNORECASE)


def _iter_entry_files(entry: Any) -> Iterable[Any]:
    if entry is None:
        return []
    if isinstance(entry, (list, tuple)):
        return entry
    if isinstance(entry, dict):
        for key in ("files", "main_files", "model_files"):
            v = entry.get(key)
            if isinstance(v, (list, tuple)):
                return v
        return []
    try:
        return list(entry)
    except TypeError:
        return []


def _file_name(obj: Any) -> str:
    return getattr(obj, "name", str(obj))


def find_first_sedml(entry_files: Iterable[Any]) -> Optional[Any]:
    for f in entry_files:
        if _SEDML_RE.search(_file_name(f)):
            return f
    return None


def find_first_sbml(entry_files: Iterable[Any]) -> Optional[Any]:
    candidates = []
    for f in entry_files:
        name = _file_name(f)
        if _SEDML_RE.search(name):
            continue
        if _SBML_RE.search(name):
            candidates.append(f)

    # Prefer SBML-ish names
    for key in ("sbml", "model"):
        for c in candidates:
            if key in _file_name(c).lower():
                return c

    return candidates[0] if candidates else None


# ----------------------------
# SED-ML parsing
# ----------------------------

def read_sedml_doc(sedml_path: str) -> libsedml.SedDocument:
    doc = libsedml.readSedMLFromFile(str(sedml_path))
    if doc is None:
        raise RuntimeError(f"libsedml returned None reading: {sedml_path}")
    if doc.getNumErrors() > 0:
        msg = doc.getErrorLog().toString()
        raise RuntimeError(f"SED-ML parse errors in {sedml_path}:\n{msg}")
    return doc


def extract_first_uniform_time_course(sed_doc: libsedml.SedDocument) -> UniformTimeCourseSpec:
    n_sims = int(sed_doc.getNumSimulations())
    for i in range(n_sims):
        sim = sed_doc.getSimulation(i)
        if sim is None:
            continue

        is_utc = False
        if hasattr(sim, "isSedUniformTimeCourse"):
            try:
                is_utc = bool(sim.isSedUniformTimeCourse())
            except Exception:
                is_utc = False

        if not is_utc:
            needed = ("getInitialTime", "getOutputStartTime", "getOutputEndTime", "getNumberOfPoints")
            is_utc = all(hasattr(sim, m) for m in needed)

        if not is_utc:
            continue

        return UniformTimeCourseSpec(
            initial_time=float(sim.getInitialTime()),
            output_start_time=float(sim.getOutputStartTime()),
            output_end_time=float(sim.getOutputEndTime()),
            number_of_points=int(sim.getNumberOfPoints()),
        )

    raise ValueError("No UniformTimeCourse simulation found in SED-ML.")


def resolve_sbml_source_from_sedml(
    sed_doc: libsedml.SedDocument,
    sedml_dir: str,
    fallback_sbml_path: str,
) -> str:
    if sed_doc.getNumModels() == 0:
        return fallback_sbml_path

    model = sed_doc.getModel(0)
    if model is None:
        return fallback_sbml_path

    src = model.getSource()
    if not src:
        return fallback_sbml_path

    if src.startswith(("http://", "https://", "urn:", "biomodels:", "BIOMD")):
        return fallback_sbml_path

    candidate = os.path.abspath(os.path.join(sedml_dir, src))
    return candidate if os.path.exists(candidate) else fallback_sbml_path


# ----------------------------
# SBML unit extraction
# ----------------------------

def _format_unit_definition(ud: "libsbml.UnitDefinition") -> str:
    """Render a UnitDefinition as a short string like 'mole * second^-1'."""
    parts = []
    for j in range(ud.getNumUnits()):
        u = ud.getUnit(j)
        kind = libsbml.UnitKind_toString(u.getKind())
        e = u.getExponent()
        parts.append(kind if e == 1 else f"{kind}^{e}")
    return " * ".join(parts)


def _resolve_unit(model: "libsbml.Model", unit_id: str) -> Optional[str]:
    """Resolve an SBML unit reference to a printable string, or None if unknown."""
    if not unit_id:
        return None
    # SBML base unit kind names ('second', 'mole', ...) are valid as-is
    if libsbml.UnitKind_forName(unit_id) != libsbml.UNIT_KIND_INVALID:
        return unit_id
    ud = model.getUnitDefinition(unit_id)
    if ud is None:
        return unit_id  # unknown reference; pass through verbatim
    rendered = _format_unit_definition(ud)
    return rendered or unit_id


def extract_sbml_metadata(sbml_path: str) -> Dict[str, Any]:
    """Pull display-friendly metadata from an SBML model.

    Returns ``{"name", "n_species", "n_reactions", "n_parameters",
    "compartments"}`` with empty/0 values when the document is missing.
    """
    doc = libsbml.SBMLReader().readSBMLFromFile(str(sbml_path))
    model = doc.getModel() if doc is not None else None
    if model is None:
        return {
            "name": None,
            "n_species": 0,
            "n_reactions": 0,
            "n_parameters": 0,
            "compartments": [],
        }
    compartments = []
    for i in range(model.getNumCompartments()):
        c = model.getCompartment(i)
        compartments.append(c.getName() or c.getId())
    return {
        "name": model.getName() or model.getId() or None,
        "n_species": int(model.getNumSpecies()),
        "n_reactions": int(model.getNumReactions()),
        "n_parameters": int(model.getNumParameters()),
        "compartments": compartments,
    }


def extract_sbml_units(sbml_path: str) -> Dict[str, Any]:
    """Pull time + per-species units from an SBML model.

    Returns ``{"time_unit": str|None, "species_units": {species_id: str|None}}``.
    Values are ``None`` when the model leaves them unset (SBML defaults apply).
    """
    doc = libsbml.SBMLReader().readSBMLFromFile(str(sbml_path))
    model = doc.getModel() if doc is not None else None
    if model is None:
        return {"time_unit": None, "species_units": {}}

    time_unit = _resolve_unit(model, model.getTimeUnits())
    default_substance = model.getSubstanceUnits() if hasattr(model, "getSubstanceUnits") else ""

    species_units: Dict[str, Optional[str]] = {}
    for i in range(model.getNumSpecies()):
        sp = model.getSpecies(i)
        unit_id = sp.getUnits() or sp.getSubstanceUnits() or default_substance
        species_units[sp.getId()] = _resolve_unit(model, unit_id)

    return {"time_unit": time_unit, "species_units": species_units}


# ----------------------------
# BioModels fetching
# ----------------------------

def fetch_biomodel_files_to_dir(biomodel_file_entry: Any, out_dir: str) -> str:
    f = biomodels.get_file(biomodel_file_entry)

    if isinstance(f, (str, os.PathLike)) and os.path.exists(str(f)):
        return str(f)

    name = _file_name(biomodel_file_entry)
    out_path = os.path.join(out_dir, name)

    if isinstance(f, bytes):
        Path(out_path).write_bytes(f)
        return out_path

    Path(out_path).write_text(str(f), encoding="utf-8")
    return out_path


def load_biomodel(biomodel_id: str, metadata_or_entry: Any) -> BiomodelLoadResult:
    entry_files = list(_iter_entry_files(metadata_or_entry))

    sedml_entry = find_first_sedml(entry_files)
    sbml_entry = find_first_sbml(entry_files)

    if sedml_entry is None:
        raise ValueError(f"{biomodel_id}: could not find a .sedml file in entry.")
    if sbml_entry is None:
        raise ValueError(f"{biomodel_id}: could not find an SBML (.xml/.sbml) file in entry.")

    with tempfile.TemporaryDirectory(prefix=f"biomodel_{biomodel_id}_") as tmp:
        sedml_path = fetch_biomodel_files_to_dir(sedml_entry, tmp)
        sbml_path = fetch_biomodel_files_to_dir(sbml_entry, tmp)

        sed_doc = read_sedml_doc(sedml_path)
        utc = extract_first_uniform_time_course(sed_doc)

        sedml_dir = os.path.dirname(sedml_path)
        resolved_sbml = resolve_sbml_source_from_sedml(sed_doc, sedml_dir, sbml_path)

        stable_dir = os.path.join("models", biomodel_id)
        os.makedirs(stable_dir, exist_ok=True)

        stable_sedml = os.path.join(stable_dir, os.path.basename(sedml_path))
        stable_sbml = os.path.join(stable_dir, os.path.basename(resolved_sbml))

        Path(stable_sedml).write_bytes(Path(sedml_path).read_bytes())
        Path(stable_sbml).write_bytes(Path(resolved_sbml).read_bytes())

    units = extract_sbml_units(stable_sbml)

    return BiomodelLoadResult(
        biomodel_id=biomodel_id,
        sbml_path=stable_sbml,
        sedml_path=stable_sedml,
        utc=utc,   # TODO also support steady state
        units=units,
    )


# ----------------------------
# Document creation (matches your UTC Step demos)
# ----------------------------

def make_utc_step_state(
    step_name: str,
    step_address: str,
    sbml_path: str,
    utc: UniformTimeCourseSpec,
) -> Dict[str, Any]:
    return {
        f"{step_name}_step": {
            "_type": "step",
            "address": step_address,
            "config": {
                "model_source": sbml_path,
                "time": float(utc.duration),
                "n_points": int(utc.number_of_points),
            },
            # ✅ paths are a single list, not list-of-lists
            "inputs": {
                "species_concentrations": ["species_concentrations"],
                "species_counts": ["species_counts"],
            },
            "outputs": {
                "result": ["results", step_name],
            },
        },
    }

def make_multi_biomodel_document(
        # biomodel_id: list[str],
        # sbml_path: list[str],
        # utc: UniformTimeCourseSpec,
        # steps: Dict[str, str],
        biomodel_info = list[Dict]
):
    doc = {"state": {}, "schema": {}}
    for biomodel in biomodel_info:
        biomodel_id = biomodel["biomodel_id"]
        sbml_path = biomodel["sbml_path"]
        utc = biomodel["utc"]
        steps = biomodel["steps"]
        single_model_doc = make_biomodel_document(biomodel_id, sbml_path, utc, steps)
        single_model_state = single_model_doc["state"]
        single_model_schema = single_model_doc["schema"]
        doc["state"][biomodel_id] = single_model_state
        doc["schema"][biomodel_id] = single_model_schema
    # TODO add comparison step
    # TODO add emitter
    return doc

def make_biomodel_document(
    biomodel_id: str,
    sbml_path: str,
    utc: UniformTimeCourseSpec,
    steps: Dict[str, str],
    with_emitter: bool = False,
    emitter_address: str = "local:RAMEmitter",
) -> Dict[str, Any]:
    """
    Store schemas align with step contracts; numeric_result is assumed to exist already.
    """
    state: Dict[str, Any] = {
        "species_concentrations": {},
        "species_counts": {},
        "results": {},
    }

    schema: Dict[str, Any] = {
        "species_concentrations": "map[float]",
        "species_counts": "map[float]",
        # results[step_name] is numeric_result
        "results": "map[numeric_result]",
    }

    engine_inputs: Dict[str, Any] = {}
    for engine_name, engine_address in steps.items():
        step_key = f"{biomodel_id}_{engine_name}"
        state.update(
            make_utc_step_state(
                step_name=step_key,
                step_address=engine_address,
                sbml_path=sbml_path,
                utc=utc,
            )
        )
        engine_inputs[step_key] = ["results", step_key]

    if with_emitter:
        wires = {**engine_inputs, "global_time": ["global_time"]}
        state["emitter"] = {
            "_type": "step",
            "address": emitter_address,
            "config": {"emit": {k: "node" for k in wires}},
            "inputs": wires,
        }

    return {"schema": schema, "state": state}


# ----------------------------
# Runner
# ----------------------------
import zipfile


async def submit_composite_document(
    document: Dict[str, Any],
    core,
    name: Optional[str] = None,
    outdir: str = "out",
    time: Optional[float] = None,
    save: bool = True,
    max_retries: int = 3,
    retry_delay: float = 5.0,
    sbml_path: Optional[str] = None,
):
    # Lazy-imported because the rest of this module (load_biomodel,
    # make_biomodel_document, run_composite_document) doesn't need them, and
    # ExecutionProgramArguments has historically moved around in pbest.
    import pbest as pb
    from pbest.utils.input_types import ExecutionProgramArguments

    outdir = Path(outdir)
    # Create Omex that gets sent to the server
    ts = datetime.now().strftime('%Y%m%d_%H%M%S')
    ts_name = f"{name}_{ts}"
    biomodel_pbg = outdir / "models" / ts_name / f"{ts_name}_url"
    omex_file = outdir / f"submit_{ts_name}.omex"
    biomodel_pbg.parent.mkdir(parents=True, exist_ok=True)
    omex_file = str(omex_file)
    # sbml_path = Path(sbml_path).name

    with open(biomodel_pbg, "w") as f:
        json.dump(document, f)

    with zipfile.ZipFile(omex_file, "w") as f:
        f.write(filename=biomodel_pbg, arcname=f"{ts_name}.pbg")
        if sbml_path:
            # sbml_arcname = Path(sbml_path).name
            f.write(filename=sbml_path, arcname=sbml_path)

    # get the runtime
    if time is None:
        times = []
        for node in document["state"].values():
            if isinstance(node, dict) and node.get("_type") == "step":
                t = node.get("config", {}).get("time")
                if isinstance(t, (int, float)):
                    times.append(float(t))
        time = max(times) if times else 10.0

    # -----------------------------------------------------#
    # Specify the input file, time length of experiment,  #
    # and where it gets saved locally                     #
    # -----------------------------------------------------#
    args: ExecutionProgramArguments = ExecutionProgramArguments(
        input_file_path=omex_file,
        interval=time,
        output_directory=outdir
    )

    last_error = None
    for attempt in range(1, max_retries + 1):
        try:
            # todo -- call http client from here? in multiprocessing

            await pb.run_remote_experiment(prog_args=args)
            print(f"All done executing {name}.")
            return
        except Exception as e:
            last_error = e
            if attempt < max_retries:
                print(f"Attempt {attempt}/{max_retries} failed for {name}: {e}")
                print(f"Retrying in {retry_delay}s ...")
                await asyncio.sleep(retry_delay)
            else:
                print(f"All {max_retries} attempts failed for {name}: {e}")
    raise last_error


def standardize_emitter_results(
    history: List[Dict[str, Any]],
    biomodel_id: str,
    engines: Iterable[str],
    duration: Optional[float] = None,
    units: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Reduce raw emitter history to one standardized JSON payload per biomodel.

    The UTC step writes its `numeric_result` once at the end of the composite
    step, so the canonical per-engine series lives in the final history entry
    under `results.{biomodel_id}_{engine}` with shape {time, columns, values}.
    """
    final = history[-1] if history else {}
    if not isinstance(final, dict):
        final = {}

    engine_payloads: Dict[str, Any] = {}
    for engine in engines:
        key = f"{biomodel_id}_{engine}"
        engine_payloads[engine] = final.get(key)

    units = units or {"time_unit": None, "species_units": {}}
    return {
        "biomodel_id": biomodel_id,
        "duration": duration,
        "global_time": final.get("global_time") if isinstance(final, dict) else None,
        "time_unit": units.get("time_unit"),
        "species_units": units.get("species_units") or {},
        "engines": engine_payloads,
    }


def run_composite_document(
    document: Dict[str, Any],
    core,
    name: Optional[str] = None,
    outdir: str = "out_biomodels",
    time: Optional[float] = None,
    save: bool = True,
    biomodel_id: Optional[str] = None,
    engines: Optional[Iterable[str]] = None,
    units: Optional[Dict[str, Any]] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> Composite:
    os.makedirs(outdir, exist_ok=True)

    if "state" not in document or not isinstance(document.get("state"), dict):
        document = {"state": document}

    if name is None:
        name = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"

    if time is None:
        times = []
        for node in document["state"].values():
            if isinstance(node, dict) and node.get("_type") == "step":
                t = node.get("config", {}).get("time")
                if isinstance(t, (int, float)):
                    times.append(float(t))
        time = max(times) if times else 10.0

    sim = Composite(document, core=core)

    if save:
        Path(os.path.join(outdir, f"{name}.json")).write_text(
            json.dumps(document, indent=2), encoding="utf-8"
        )
        Path(os.path.join(outdir, f"{name}_schema.json")).write_text(
            json.dumps(core.render(sim.schema), indent=2), encoding="utf-8"
        )

    print(f"⏱ Running {name} for {time}s ...")
    sim.run(time)
    print(f"✅ Done: {name}")

    if save and biomodel_id and engines:
        emitter_results = gather_emitter_results(sim)
        history = next(iter(emitter_results.values()), [])
        payload = standardize_emitter_results(
            history=history,
            biomodel_id=biomodel_id,
            engines=engines,
            duration=time,
            units=units,
        )
        if metadata:
            payload["metadata"] = metadata
        Path(os.path.join(outdir, f"{name}_results.json")).write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )

    return sim


async def run_biomodels(
        core,
        number_of_models: int = 2,
        mode: str = "local",
) -> List[BiomodelLoadResult]:
    if mode not in ("local", "remote"):
        raise ValueError(f"mode must be 'local' or 'remote', got {mode!r}")

    biomodel_ids = biomodels.get_all_identifiers()[:number_of_models]

    # Addresses match your discovered step classes
    steps = {
        "copasi": "local:pbsim_common.simulators.copasi_process.CopasiUTCStep",
        "tellurium": "local:pbsim_common.simulators.tellurium_process.TelluriumUTCStep",
        # "copasi": "local:CopasiUTCStep",
        # "tellurium": "local:TelluriumUTCStep",
    }

    os.makedirs("documents", exist_ok=True)
    loaded: List[BiomodelLoadResult] = []
    failed: List[str] = []

    for biomodel_id in biomodel_ids:
        try:
            meta = biomodels.get_metadata(biomodel_id)
            result = load_biomodel(biomodel_id, meta)
            loaded.append(result)

            # Local COPASI/Tellurium steps resolve relative paths against the
            # package install dir, so absolutise here. For remote we keep the
            # relative path because it doubles as the OMEX arcname.
            doc_sbml_path = (
                os.path.abspath(result.sbml_path) if mode == "local" else result.sbml_path
            )
            doc = make_biomodel_document(
                biomodel_id=biomodel_id,
                sbml_path=doc_sbml_path,
                utc=result.utc,
                steps=steps,
                with_emitter=(mode == "local"),
            )

            # Save doc for inspection
            Path(os.path.join("documents", f"{biomodel_id}.json")).write_text(
                json.dumps(doc, indent=2), encoding="utf-8"
            )


            if mode == "remote":
                await submit_composite_document(
                    doc,
                    core=core,
                    name=f"{biomodel_id}_utc",
                    outdir="out_biomodels",
                    time=None,
                    save=True,
                    sbml_path=result.sbml_path,
                )
            else:
                metadata = extract_sbml_metadata(result.sbml_path)
                metadata["n_points"] = result.utc.number_of_points
                metadata["biomodels_url"] = (
                    f"https://www.ebi.ac.uk/biomodels/{biomodel_id}"
                )
                run_composite_document(
                    doc,
                    core=core,
                    name=f"{biomodel_id}_utc",
                    outdir="out_biomodels",
                    time=None,
                    save=True,
                    biomodel_id=biomodel_id,
                    engines=list(steps.keys()),
                    units=result.units,
                    metadata=metadata,
                )
        except Exception as e:
            print(f"FAILED {biomodel_id}: {e}")
            failed.append(biomodel_id)

    if failed:
        print(f"\n{len(failed)}/{len(biomodel_ids)} model(s) failed: {failed}")

    return loaded

async def run_one_biomodel(core, biomodel_id: str, mode: str = "local") -> None:
    """Run a single biomodel by ID, end-to-end, in this process."""
    meta = biomodels.get_metadata(biomodel_id)
    result = load_biomodel(biomodel_id, meta)

    steps = {
        "copasi": "local:pbsim_common.simulators.copasi_process.CopasiUTCStep",
        "tellurium": "local:pbsim_common.simulators.tellurium_process.TelluriumUTCStep",
    }
    os.makedirs("documents", exist_ok=True)

    doc_sbml_path = (
        os.path.abspath(result.sbml_path) if mode == "local" else result.sbml_path
    )
    doc = make_biomodel_document(
        biomodel_id=biomodel_id,
        sbml_path=doc_sbml_path,
        utc=result.utc,
        steps=steps,
        with_emitter=(mode == "local"),
    )
    Path(os.path.join("documents", f"{biomodel_id}.json")).write_text(
        json.dumps(doc, indent=2), encoding="utf-8"
    )

    if mode == "remote":
        await submit_composite_document(
            doc,
            core=core,
            name=f"{biomodel_id}_utc",
            outdir="out_biomodels",
            time=None,
            save=True,
            sbml_path=result.sbml_path,
        )
    else:
        metadata = extract_sbml_metadata(result.sbml_path)
        metadata["n_points"] = result.utc.number_of_points
        metadata["biomodels_url"] = f"https://www.ebi.ac.uk/biomodels/{biomodel_id}"
        run_composite_document(
            doc,
            core=core,
            name=f"{biomodel_id}_utc",
            outdir="out_biomodels",
            time=None,
            save=True,
            biomodel_id=biomodel_id,
            engines=list(steps.keys()),
            units=result.units,
            metadata=metadata,
        )


def supervise_biomodels(
    number_of_models: int,
    mode: str = "local",
    timeout: float = 300.0,
) -> None:
    """Run each biomodel in its own subprocess so segfaults / hangs are isolated."""
    import subprocess
    import sys

    biomodel_ids = biomodels.get_all_identifiers()[:number_of_models]
    failed: List[str] = []
    n = len(biomodel_ids)
    for i, bid in enumerate(biomodel_ids, start=1):
        print(f"[{i}/{n}] {bid}", flush=True)
        try:
            r = subprocess.run(
                [sys.executable, __file__, "--single-model", bid, "--mode", mode],
                timeout=timeout,
            )
            if r.returncode != 0:
                print(f"FAILED {bid}: subprocess exit {r.returncode}", flush=True)
                failed.append(bid)
        except subprocess.TimeoutExpired:
            print(f"FAILED {bid}: timed out after {timeout}s", flush=True)
            failed.append(bid)

    if failed:
        print(
            f"\n{len(failed)}/{n} model(s) failed: {failed[:20]}"
            f"{'...' if len(failed) > 20 else ''}"
        )
    print(f"Loaded {n - len(failed)} biomodel(s).")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run biomodels regression.")
    parser.add_argument("-n", "--number-of-models", type=int, default=5)
    parser.add_argument("--mode", choices=("local", "remote"), default="local")
    parser.add_argument(
        "--single-model",
        metavar="BIOMD_ID",
        help="Run exactly one model in this process (used by the supervisor).",
    )
    parser.add_argument(
        "--no-isolate",
        action="store_true",
        help="Run all models in a single process (faster, no segfault isolation).",
    )
    parser.add_argument(
        "--per-model-timeout",
        type=float,
        default=300.0,
        help="Seconds per model when supervising subprocesses (default 300).",
    )
    args = parser.parse_args()

    if args.single_model:
        core = get_loaded_core()
        asyncio.run(run_one_biomodel(core, args.single_model, mode=args.mode))
    elif args.no_isolate:
        core = get_loaded_core()
        loaded = asyncio.run(
            run_biomodels(core, number_of_models=args.number_of_models, mode=args.mode)
        )
        print(f"Loaded {len(loaded)} biomodel(s).")
    else:
        supervise_biomodels(
            number_of_models=args.number_of_models,
            mode=args.mode,
            timeout=args.per_model_timeout,
        )
