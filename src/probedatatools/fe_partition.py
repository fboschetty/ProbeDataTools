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

from .cations import calc_cations
from .probe import ProbeData


def partition_Fe(probe_data: ProbeData, cations: pd.DataFrame, Fe3ideal: pd.Series) -> ProbeData:
    """Partition total Fe, reported as FeO, into FeO and Fe2O3.

    Parameters
    ----------
    probe_data:
        Probe data with total Fe reported as FeO.

    cations:
        Cations per formula unit, with total Fe represented by ``Fe2``.

    Fe3ideal:
        Ideal Fe3+ content per formula unit.

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

    new_oxides = list(probe_data.oxides)
    if "Fe2O3" not in new_oxides:
        new_oxides.append("Fe2O3")

    return ProbeData(new_data, new_oxides)


def calc_Fe2O3_Droop(probe_data: ProbeData, cfu: float, afu: float, norm: bool = False) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop (1987).

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.

    cfu : float
        Ideal number of cations per formula unit (T in Droop, 1987).

    afu : float
        Number of anions per formula unit (X in Droop, 1987).

    norm : bool, default=False
        If True, normalize the oxide analysis so that the calculated
        cation total equals ``cfu``. If False, retain the original
        analysed oxide concentrations.
    """

    cations = calc_cations(probe_data, afu=afu)
    S = cations.sum(axis=1, skipna=True)

    factor = cfu / S
    Fe3ideal = 2 * afu * (1 - (cfu / S))

    if norm:
        cations = cations.mul(factor, axis=0)

    return partition_Fe(probe_data, cations, Fe3ideal)


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

    cations = calc_cations(probe_data, afu=23)

    rel_cat = ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca"]

    S = cations[rel_cat].sum(axis=1, skipna=True)
    factor = 15.0 / S
    Fe3ideal = 46.0 * (1 - factor)

    if norm:
        cations = cations.mul(factor, axis=0)

    return partition_Fe(probe_data, cations, Fe3ideal)


def calc_Fe2O3_Droop_Eq6(probe_data: ProbeData, norm: bool = False) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
    Equation 6 tailored for many calcic amphiboles. Requires cations on the basis of 23.0 oxygens, anhydrous.
    Using 13.0 cations, exclusive of Ca, Na and K.

    Parameters
    ----------
    probe_data : ProbeData
        Probe data with total Fe reported as FeO.
    """

    cations = calc_cations(probe_data, afu=23)

    rel_cat = ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg"]

    S = cations[rel_cat].sum(axis=1, skipna=True)
    factor = 13.0 / S
    Fe3ideal = 46.0 * (1.0 - factor)

    if norm:
        cations = cations.mul(factor, axis=0)

    return partition_Fe(probe_data, cations, Fe3ideal)


def _normalize_cations(cations: pd.DataFrame, cfu: float) -> pd.DataFrame:
    """Normalize cations to a specified total cation number."""
    return cations.mul(cfu / cations.sum(axis=1, skipna=True), axis=0)


def calc_Fe2O3_Stormer(probe_data: ProbeData) -> ProbeData:
    """Calculate Fe2+/Fe3+ for spinel using Stormer (1983) [1].

    The cations are normalized to 3 cation sites, and Fe2+/Fe3+
    is adjusted to balance the total cation charge to 8.

    References
    ----------
    [1] https://pubs.geoscienceworld.org/msa/ammin/article/68/5-6/586/104818/
    """

    cations = calc_cations(probe_data, afu=4.0, change_head=False)
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

    cations = calc_cations(probe_data, afu, change_head=False)

    Fe3ideal = _calc_Fe3_charge_balance(cations, probe_data.cat_chrg_use, afu)

    cations = cations.rename(columns={"FeO": "Fe2"})

    return partition_Fe(probe_data, cations, Fe3ideal)