# pbg-biomodels-bundle (superseded)

This package has been **superseded** by [pbg-biomodels](https://github.com/vivarium-collective/pbg-biomodels), which now contains all functionality that used to live here:

- `LoadBiomodelStep` (BioModels fetcher) → `pbg_biomodels.steps.load_biomodel.LoadBiomodelStep`
- `SimulatorComparisonStep` (per-species nRMSE) → `pbg_biomodels.steps.simulator_comparison.SimulatorComparisonStep`
- `comparison.py` (shared comparison math) → `pbg_biomodels.comparison`
- `run_biomodels.py` (BioModels loader + batch UTC runner) → `pbg_biomodels.run_biomodels`
- `analysis.py` (offline HTML report builder) → `pbg_biomodels.analysis`
- `LocalCopasiUTCStep` / `LocalTelluriumUTCStep` → REPLACED by [`BiomodelsCopasiStep` / `BiomodelsTelluriumStep`](https://github.com/vivarium-collective/pbg-biomodels/blob/main/pbg_biomodels/steps/simulators.py), thin wrappers that delegate to [pbg-copasi](https://github.com/vivarium-collective/pbg-copasi) and [pbg-tellurium](https://github.com/vivarium-collective/pbg-tellurium) — eliminating duplicated simulator code.

**No further development happens in this repository.** New work continues in [pbg-biomodels](https://github.com/vivarium-collective/pbg-biomodels).
