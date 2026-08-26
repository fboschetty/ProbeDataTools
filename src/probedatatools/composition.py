from .probe import ProbeData


def convert_units(probe_data: ProbeData, from_unit: str, to_unit: str) -> ProbeData:
    """Convert analytical concentrations between ppm and wt.%."""

    valid_units = {"ppm", "wt%"}

    if from_unit not in valid_units or to_unit not in valid_units:
        raise ValueError(f"Unsupported units: {from_unit} -> {to_unit}")

    if from_unit == to_unit:
        return ProbeData(probe_data.data.copy(), list(probe_data.species))

    factor = 1e-4 if from_unit == "ppm" else 1e4
    data = probe_data.data.copy()
    data[list(probe_data.species)] *= factor

    return ProbeData(data, list(probe_data.species))


def element_to_oxide(probe_data: ProbeData) -> ProbeData:
    """Convert elemental concentrations to equivalent oxide concentrations."""

    elements = list(probe_data.species)
    if "Fe" in elements:
        raise ValueError("Iron must be specified as Fe2 or Fe3, not Fe.")

    mapping = {}
    for element in elements:
        oxides = ProbeData.cat_str[ProbeData.cat_str == element].index.tolist()

        if not oxides:
            raise ValueError(f"No oxide found for element: {element}")

        if len(oxides) > 1:
            raise ValueError(f"Multiple oxides found for element: {element}")

        mapping[element] = oxides[0]

    oxides = list(mapping.values())
    factors = ProbeData.MR.loc[oxides] / (
        ProbeData.cat_num.loc[oxides] * ProbeData.MR_El.loc[oxides]
    )

    data = probe_data.data.copy()
    data[oxides] = data[elements].mul(factors.to_numpy(), axis=1)
    data.drop(columns=elements, inplace=True)

    return ProbeData(data, oxides)


def oxide_to_element(probe_data: ProbeData) -> ProbeData:
    """Convert oxide concentrations to valence-specific cation concentrations."""

    invalid = probe_data.cat_num_use[probe_data.cat_num_use <= 0].index.tolist()
    if invalid:
        raise ValueError(f"Species cannot be converted to elements: {invalid}")

    oxides = list(probe_data.species)
    elements = probe_data.cat_str_use.tolist()
    factors = probe_data.cat_num_use * probe_data.MR_El_use / probe_data.MR_use

    data = probe_data.data.copy()
    data[elements] = data[oxides].mul(factors, axis=1).set_axis(elements, axis=1)
    data.drop(columns=oxides, inplace=True)

    return ProbeData(data, elements)