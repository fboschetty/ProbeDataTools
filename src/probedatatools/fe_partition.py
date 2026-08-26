"""
Functions for calculating Fe2+/Fe3+ in iron-bearing mineral analyses.

Includes:
    (1) Droop (1987) stoichiometric Fe2+/Fe3+ calculations.
    (2) Stormer (1983) recalculation for spinel-group minerals.
    (3) Generic charge-balance Fe2+/Fe3+ calculation.

All Fe-partitioning functions expect total iron to be reported as FeO.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .cations import _cation_columns, calc_apfu
from .probe import ProbeData


def _calc_cat_tot(probe_data: ProbeData, apfu: pd.DataFrame) -> pd.Series:
    """Calculate total cations from an APFU dataframe.

    Parameters
    ----------
    probe_data : ProbeData
        Probe-data object defining cation stoichiometries and labels.
    apfu : pandas.DataFrame
        APFU dataframe containing cations and any analysed non-O anions.

    Returns
    -------
    pandas.Series
        Total cations for each analysis.
    """

    columns = _cation_columns(probe_data)

    return apfu[columns].sum(axis=1, skipna=True)


def partition_Fe(probe_data: ProbeData, cations: pd.DataFrame, Fe3ideal: pd.Series) -> ProbeData:
    """Partition total Fe, reported as FeO, into FeO and Fe2O3.

    Parameters
    ----------
    probe_data:
        Probe data with total Fe reported as FeO.
    cations : pandas.DataFrame
        Cation APFU dataframe in which the ``Fe2`` column represents
        total analysed Fe temporarily treated as Fe2+.
    Fe3ideal:
        Calculated Fe3+ content in atoms per formula unit.

    Returns
    -------
    ProbeData
        A new ProbeData object with FeO and Fe2O3.

    Notes
    -----
    The calculated Fe3+ content is clipped between zero and the
    available total Fe.
    """

    Fe3 = np.clip(Fe3ideal, 0, cations["Fe2"])
    Fe2 = cations["Fe2"] - Fe3

    Fe2FeT = (Fe2 / (Fe2 + Fe3)).fillna(0)

    new_data = probe_data.data.copy()
    new_data["FeO"] = probe_data.data["FeO"] * Fe2FeT
    Fe2O3_factor = probe_data.MR["Fe2O3"] / (2 * probe_data.MR["FeO"])

    new_data["Fe2O3"] = probe_data.data["FeO"] * (1 - Fe2FeT) * Fe2O3_factor

    new_oxides = list(probe_data.species)
    if "Fe2O3" not in new_oxides:
        new_oxides.append("Fe2O3")

    return ProbeData(new_data, new_oxides)


def calc_Fe2O3_Droop(probe_data: ProbeData, cfu: float, afu: float, norm: bool = False) -> ProbeData:
    """Calculate Fe2+/Fe3+ stoichiometrically using Droop (1987).

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.
    cfu : float
        Ideal number of cations per formula unit.
    afu : float
        Target number of oxygens per formula unit.
    norm : bool, default=False
        If True, normalise the APFU cations to ``cfu`` before partitioning
        Fe. If False, retain the original analysed concentrations when
        partitioning Fe.

    Returns
    -------
    ProbeData
        New probe data with total Fe partitioned between FeO and Fe2O3.
    """

    apfu = calc_apfu(probe_data, afu=afu)
    S = _calc_cat_tot(probe_data, apfu)

    factor = cfu / S
    Fe3ideal = 2 * afu * (1 - factor)

    if norm:
        cation_columns = _cation_columns(probe_data)
        apfu[cation_columns] = apfu[cation_columns].mul(factor, axis=0)

    return partition_Fe(probe_data, apfu, Fe3ideal)


def calc_Fe2O3_Droop_Eq4(probe_data: ProbeData, norm: bool = False) -> ProbeData:
    """Calculate Fe2+/Fe3+ using Droop (1987), Equation 4.

    Equation 4 applies to amphiboles with full A-sites:

        Fe3+ = 46 * (1 - 16 / S)

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.

    norm : bool, default=False
        If True, normalize the oxide analysis so that the calculated
        cation total equals ``cfu``. If False, retain the original
        analysed oxide concentrations.
    """

    return calc_Fe2O3_Droop(probe_data, cfu=16.0, afu=23.0, norm=norm)


def calc_Fe2O3_Droop_Eq5(probe_data: ProbeData, norm: bool = False) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
    Equation 5 tailored for Fe-Mg amphiboles. Requires calculation of cations on the basis of 23.0 oxygens, anhydrous.
    Using 15.0 cations exclusive of Na and K.

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.

    norm : bool, default=False
        If True, normalize the oxide analysis so that the calculated
        cation total equals ``cfu``. If False, retain the original
        analysed oxide concentrations.
    """

    apfu = calc_apfu(probe_data, afu=23)

    rel_cat = ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca"]

    S = apfu[rel_cat].sum(axis=1, skipna=True)
    factor = 15.0 / S
    Fe3ideal = 46.0 * (1 - factor)

    if norm:
        cation_columns = _cation_columns(probe_data)
        apfu[cation_columns] = apfu[cation_columns].mul(factor, axis=0)

    return partition_Fe(probe_data, apfu, Fe3ideal)


def calc_Fe2O3_Droop_Eq6(probe_data: ProbeData, norm: bool = False) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
    Equation 6 tailored for many calcic amphiboles. Requires cations on the basis of 23.0 oxygens, anhydrous.
    Using 13.0 cations, exclusive of Ca, Na and K.

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.
    """

    apfu = calc_apfu(probe_data, afu=23)

    rel_cat = ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg"]

    S = apfu[rel_cat].sum(axis=1, skipna=True)
    factor = 13.0 / S
    Fe3ideal = 46.0 * (1.0 - factor)

    if norm:
        cation_columns = _cation_columns(probe_data)
        apfu[cation_columns] = apfu[cation_columns].mul(factor, axis=0)

    return partition_Fe(probe_data, apfu, Fe3ideal)


def _normalize_cations(cations: pd.DataFrame, cfu: float) -> pd.DataFrame:
    """Normalise cation proportions to a specified cation total.

    Parameters
    ----------
    cations : pandas.DataFrame
        Dataframe containing cations only.
    cfu : float
        Target number of cations per formula unit.

    Returns
    -------
    pandas.DataFrame
        Cation proportions normalised to ``cfu``.
    """

    return cations.mul(cfu / cations.sum(axis=1, skipna=True), axis=0)


def calc_Fe2O3_Stormer(probe_data: ProbeData) -> ProbeData:
    """Calculate Fe2+/Fe3+ for spinel using Stormer (1983) [1].

    The cations are normalized to 3 cation sites, and Fe2+/Fe3+
    is adjusted to balance the total cation charge to 8.

    References
    ----------
    [1] https://pubs.geoscienceworld.org/msa/ammin/article/68/5-6/586/104818/
    """

    apfu = calc_apfu(probe_data, afu=4.0, change_head=False)
    cation_species = probe_data.cat_num_use[probe_data.cat_num_use > 0].index
    cations = apfu[cation_species]
    cations = _normalize_cations(cations, cfu=3.0)

    charge = cations.mul(probe_data.cat_chrg_use, axis=1).sum(axis=1, skipna=True)
    Fe3ideal = 8.0 - charge
    cations = cations.rename(columns={"FeO": "Fe2"})

    return partition_Fe(probe_data, cations, Fe3ideal)


def _calc_Fe3_charge_balance(cations: pd.DataFrame, cat_chrg: pd.Series, afu: float) -> pd.Series:
    """Calculate Fe3+ required for charge balance."""

    Fe_total = cations["FeO"]

    charge_other = (
        cations.drop(columns="FeO")
        .mul(cat_chrg.drop("FeO"), axis=1)
        .sum(axis=1, skipna=True)
    )

    return np.clip(2 * afu - charge_other - 2 * Fe_total, 0, Fe_total)


def calc_Fe2O3_charge_balance(probe_data: ProbeData, afu: float) -> ProbeData:
    """Calculate Fe2+/Fe3+ by charge balance.

    All measured Fe is initially treated as Fe2+. Fe3+ is then
    calculated from charge neutrality on the specified oxygen basis.

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.

    afu : float
        Number of oxygens per formula unit.

    Returns
    -------
    ProbeData
        Probe data with Fe partitioned between FeO and Fe2O3.
    """

    apfu = calc_apfu(probe_data, afu, change_head=False)

    cation_species = probe_data.cat_num_use[probe_data.cat_num_use > 0].index
    cations = apfu[cation_species]

    Fe3ideal = _calc_Fe3_charge_balance(cations, probe_data.cat_chrg_use, afu)

    cations = cations.rename(columns={"FeO": "Fe2"})

    return partition_Fe(probe_data, cations, Fe3ideal)