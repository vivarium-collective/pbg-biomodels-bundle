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
