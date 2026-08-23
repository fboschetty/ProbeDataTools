"""
Calculate composition variables and cations per formula unit from oxide analyses.

The functions in this module operate on ``ProbeData`` objects.
"""

from __future__ import annotations

import pandas as pd

from .probe import ProbeData


def calc_mol_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate molar proportions from wt% oxides."""

    return probe_data.data[list(probe_data.species)].div(probe_data.MR_use, axis=1)


def calc_mol_frac(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate molar fractions from wt% oxides."""

    mol_prop = calc_mol_prop(probe_data)
    return mol_prop.div(mol_prop.sum(axis=1, skipna=True), axis=0)


def calc_oxygen_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate oxygen proportions from oxide molar proportions."""

    return calc_mol_prop(probe_data).mul(probe_data.ox_num_use, axis=1)


def calc_cat_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate unnormalized cation proportions from wt% oxides."""

    return calc_mol_prop(probe_data).mul(probe_data.cat_num_use, axis=1)


def calc_cat_frac(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate cation fractions from wt% oxides."""

    cat_prop = calc_cat_prop(probe_data)
    return cat_prop.div(cat_prop.sum(axis=1, skipna=True), axis=0)


def calc_ARF(
    probe_data: ProbeData,
    afu: float,
    anion_species: tuple[str, ...] = ("O",),
) -> pd.Series:
    """Calculate the anion renormalisation factor.

    Parameters
    ----------
    afu : float
        Target number of anions per formula unit.

    anion_species : tuple[str, ...]
        Anions included in the normalisation basis.
    """

    anion_tot = pd.Series(0.0, index=probe_data.data.index)

    if "O" in anion_species:
        anion_tot += calc_oxygen_prop(probe_data).sum(axis=1, skipna=True)

    other_anions = pd.Index(anion_species).difference(["O"])
    other_anions = other_anions.intersection(probe_data.data.columns)

    if len(other_anions):
        anion_tot += (
            calc_mol_prop(probe_data)[other_anions]
            .mul(probe_data.an_num_use[other_anions], axis=1)
            .sum(axis=1, skipna=True)
        )

    return afu / anion_tot


def calc_nonO_anions(
    probe_data: ProbeData,
    afu: float,
    anion_species: tuple[str, ...] = ("O",),
) -> pd.DataFrame:
    """Calculate non-oxygen anions per formula unit."""

    other_anions = pd.Index(anion_species).difference(["O"])
    other_anions = other_anions.intersection(probe_data.data.columns)

    if not len(other_anions):
        return pd.DataFrame(index=probe_data.data.index)

    anions = calc_mol_prop(probe_data)[other_anions].mul(
        probe_data.an_num_use[other_anions],
        axis=1,
    )

    return anions.mul(
        calc_ARF(probe_data, afu, anion_species),
        axis=0,
    )


def _change_headers_cfu(cations: pd.DataFrame, probe_data: ProbeData,
) -> pd.DataFrame:
    """Change headers from analytical oxides to cation names."""

    headers = probe_data.cat_str.loc[cations.columns]
    return cations.rename(columns=headers)


def calc_cations(
    probe_data: ProbeData,
    afu: float,
    anion_species: tuple[str, ...] = ("O",),
    change_head: bool = True,
) -> pd.DataFrame:
    """Calculate cations per formula unit on a specified anion basis.

    Parameters
    ----------
    afu : float
        Target number of anions per formula unit.

    anion_species : tuple[str, ...]
        Analytical anions included in the normalisation basis.

    change_head : bool
        If True, change analytical headers to cation names.
    """

    cations = calc_cat_prop(probe_data)
    cations = cations.loc[:, probe_data.cat_num_use > 0]
    cations = cations.mul(
        calc_ARF(probe_data, afu, anion_species),
        axis=0,
    )

    return _change_headers_cfu(cations, probe_data) if change_head else cations


def calc_cat_tot(
    probe_data: ProbeData,
    afu: float,
    anion_species: tuple[str, ...] = ("O",),
) -> pd.DataFrame:
    """Calculate cations per formula unit and their total."""

    cations = calc_cations(probe_data, afu, anion_species)
    cations["cat_tot"] = cations.sum(axis=1, skipna=True)
    return cations


def check_cat_tot(
    probe_data: ProbeData,
    cfu: float,
    afu: float,
    anion_species: tuple[str, ...] = ("O",),
    tolerance: float = 0.005,
) -> list[bool]:
    """Check whether cation totals fall within the specified tolerance."""

    cat_tot = calc_cat_tot(probe_data, afu, anion_species)["cat_tot"]
    deviation = (cat_tot - cfu) / cfu
    return (deviation.abs() <= tolerance).tolist()