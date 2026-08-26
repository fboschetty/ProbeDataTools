"""
Calculate composition variables and atoms per formula unit (APFU).

The functions in this module operate on ``ProbeData`` objects and support
calculation of molar proportions, oxygen proportions, cation proportions,
oxygen-basis renormalisation, and atoms per formula unit.
"""

from __future__ import annotations

import pandas as pd

from .probe import ProbeData


def _cation_columns(probe_data: ProbeData) -> pd.Index:
    """Return APFU column names corresponding to analysed cations."""

    cation_species = probe_data.cat_num_use[probe_data.cat_num_use > 0].index

    return probe_data.cat_str.loc[cation_species]


def calc_mol_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate molar proportions from analytical concentrations.

    Analytical concentrations are divided by the corresponding molar
    masses defined by ``ProbeData.MR_use``.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and species
        definitions.

    Returns
    -------
    pandas.DataFrame
        Molar proportions for each analysed species.
    """

    return probe_data.data[list(probe_data.species)].div(probe_data.MR_use, axis=1)


def calc_mol_frac(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate molar fractions from analytical concentrations.

    Molar proportions are normalised to sum to one for each analysis.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations.

    Returns
    -------
    pandas.DataFrame
        Molar fractions for each analysed species.
    """

    mol_prop = calc_mol_prop(probe_data)
    return mol_prop.div(mol_prop.sum(axis=1, skipna=True), axis=0)


def calc_oxygen_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate oxygen proportions from analytical concentrations.

    Each species' molar proportion is multiplied by its number of oxygen
    atoms defined by ``ProbeData.ox_num_use``.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and oxygen
        stoichiometries.

    Returns
    -------
    pandas.DataFrame
        Oxygen proportions contributed by each analysed species.
    """

    return calc_mol_prop(probe_data).mul(probe_data.ox_num_use, axis=1)


def calc_cat_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate unnormalised cation proportions.

    Each species' molar proportion is multiplied by its number of cations
    defined by ``ProbeData.cat_num_use``.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and cation
        stoichiometries.

    Returns
    -------
    pandas.DataFrame
        Unnormalised cation proportions for each analysed species.
    """

    return calc_mol_prop(probe_data).mul(probe_data.cat_num_use, axis=1)


def calc_cat_frac(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate cation fractions.

    Unnormalised cation proportions are normalised so that the cation
    fractions sum to one for each analysis.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations.

    Returns
    -------
    pandas.DataFrame
        Cation fractions for each analysed species.
    """

    cat_prop = calc_cat_prop(probe_data)
    return cat_prop.div(cat_prop.sum(axis=1, skipna=True), axis=0)


def calc_ORF(probe_data: ProbeData, afu: float) -> pd.Series:
    """Calculate the oxygen renormalisation factor.

    The oxygen renormalisation factor (ORF) converts the analysed oxygen
    total to a specified number of oxygens per formula unit.

    The ORF is calculated as:

        ORF = afu / sum(O)

    where ``afu`` is the target number of oxygens per formula unit.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing the analytical concentrations.
    afu : float
        Target number of oxygens per formula unit.

    Returns
    -------
    pandas.Series
        Oxygen renormalisation factor for each analysis.
    """

    oxygen_tot = calc_oxygen_prop(probe_data).sum(axis=1, skipna=True)

    return afu / oxygen_tot


def calc_nonO_anions(probe_data: ProbeData, afu: float) -> pd.DataFrame:
    """Calculate analysed non-oxygen anions per formula unit.

    Non-oxygen species identified as anions by ``ProbeData.an_num_use``
    are converted to atoms per formula unit using the oxygen
    renormalisation factor. Oxygen itself is excluded because it defines
    the normalisation basis.

    For example, F and Cl are calculated on the same oxygen basis as the
    cations but do not contribute to the oxygen renormalisation factor.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and anion
        stoichiometries.
    afu : float
        Target number of oxygens per formula unit.

    Returns
    -------
    pandas.DataFrame
        Analysed non-oxygen anions in atoms per formula unit. An empty
        dataframe is returned when no non-oxygen anions are present.
    """

    anion_species = probe_data.an_num_use[probe_data.an_num_use > 0].index
    nonO_anions = anion_species.intersection(probe_data.data.columns).difference(["O"])

    if not len(nonO_anions):
        return pd.DataFrame(index=probe_data.data.index)

    anions = calc_mol_prop(probe_data)[nonO_anions].mul(
        probe_data.an_num_use[nonO_anions],
        axis=1,
    )

    return anions.mul(calc_ORF(probe_data, afu), axis=0)


def _change_headers_cfu(cations: pd.DataFrame, probe_data: ProbeData) -> pd.DataFrame:
    """Rename analytical species to their cation labels.

    Parameters
    ----------
    cations : pandas.DataFrame
        Dataframe containing cation proportions with analytical species
        as column names.
    probe_data : ProbeData
        Probe-data object providing the mapping from analytical species
        to cation names.

    Returns
    -------
    pandas.DataFrame
        Cation dataframe with cation names as column headers.
    """

    headers = probe_data.cat_str.loc[cations.columns]
    return cations.rename(columns=headers)


def calc_apfu(
    probe_data: ProbeData,
    afu: float,
    change_head: bool = True,
) -> pd.DataFrame:
    """Calculate atoms per formula unit (APFU).

    Cations are normalised to a specified number of oxygens per formula
    unit. Analysed non-oxygen anions, such as F and Cl, are then
    calculated using the same oxygen renormalisation factor and appended
    to the result.

    Thus, non-oxygen anions are reported in the APFU dataframe but do not
    contribute to the oxygen normalisation.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and
        stoichiometric information.
    afu : float
        Target number of oxygens per formula unit.
    change_head : bool, default=True
        If True, rename cation columns from analytical species names to
        cation names. Non-oxygen anion names are unchanged.

    Returns
    -------
    pandas.DataFrame
        Atoms per formula unit for all analysed cations and non-oxygen
        anions.
    """

    orf = calc_ORF(probe_data, afu)

    cations = calc_cat_prop(probe_data)
    cations = cations.loc[:, probe_data.cat_num_use > 0]
    cations = cations.mul(orf, axis=0)

    if change_head:
        cations = _change_headers_cfu(cations, probe_data)

    anions = calc_nonO_anions(probe_data, afu)

    return pd.concat([cations, anions], axis=1)


def calc_cat_tot(probe_data: ProbeData, afu: float) -> pd.DataFrame:
    """Calculate APFU and total cations.

    The returned dataframe contains all values from :func:`calc_apfu` and
    an additional ``cat_tot`` column. The cation total includes only
    species identified as cations by ``ProbeData.cat_num_use``; non-oxygen
    anions such as F and Cl are excluded.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and
        stoichiometric information.
    afu : float
        Target number of oxygens per formula unit.

    Returns
    -------
    pandas.DataFrame
        Atoms per formula unit with an additional total-cation column.
    """

    apfu = calc_apfu(probe_data, afu)
    cation_columns = _cation_columns(probe_data)

    apfu["cat_tot"] = apfu[cation_columns].sum(axis=1, skipna=True)

    return apfu


def check_cat_tot(probe_data: ProbeData, cfu: float, afu: float, tolerance: float = 0.005) -> list[bool]:
    """Check whether calculated cation totals match a target value.

    The calculated cation total is compared with ``cfu``. An analysis
    passes when its relative deviation from the target cation total is
    less than or equal to ``tolerance``.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object containing analytical concentrations and
        stoichiometric information.
    cfu : float
        Target number of cations per formula unit.
    afu : float
        Target number of oxygens per formula unit used for APFU
        normalisation.
    tolerance : float, default=0.005
        Maximum allowed relative deviation from ``cfu``.

    Returns
    -------
    list[bool]
        Boolean result for each analysis, where ``True`` indicates that
        the calculated cation total is within the specified tolerance.
    """

    cat_tot = calc_cat_tot(probe_data, afu)["cat_tot"]
    deviation = (cat_tot - cfu) / cfu

    return (deviation.abs() <= tolerance).tolist()