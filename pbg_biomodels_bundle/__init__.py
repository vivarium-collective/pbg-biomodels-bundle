sed_types = {
    'result': {
        'time': 'list[float]',
        'species_concentrations': 'map[list[float]]',
    },
    'results': 'map[result]'
}

standard_types = {
    'numeric_result': {
        'time': 'list[float]',
        'columns': 'list[string]',
        'values': 'list[list[float]]',
        # 'n_spacial_dimensions': 'tuple[int, int]'
    },
    'numeric_results': 'map[numeric_result]',
    'columns_of_interest': 'list[string]'
}


TYPES_DICT = {
    **standard_types,
    **sed_types
}


# -----------------------------------------------------------------------------
# Public API
# -----------------------------------------------------------------------------

def register_types(core):
    core.register_types(TYPES_DICT)
    return core


# Re-export the public Step so callers can do
# `from pbg_biomodels_bundle import SimulatorComparisonStep`.
# Imported here (not at the top) so that `register_types` keeps working even
# when downstream dependencies of the steps subpackage are missing.
from pbg_biomodels_bundle.steps import SimulatorComparisonStep  # noqa: E402

__all__ = ["TYPES_DICT", "register_types", "SimulatorComparisonStep"]
