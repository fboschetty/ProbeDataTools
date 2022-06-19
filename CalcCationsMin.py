#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains a group of functions used to calculate the cations per formula unit of any mineral,
given an ideal number of cations and anions per formula unit.
It also contains functions for filtering out non-ideal analyses.

@author: FelixBoschetty
"""

import pandas as pd
from ProbeData import probedata


def calc_mol_prop(probe_data: probedata) -> pd.DataFrame:
    """Calculate molar proportions from wt% oxides."""
    return probe_data.data.div(probe_data.MR, axis=1)


def calc_ox_prop(probe_data: probedata) -> pd.DataFrame:
    """Calculate anion proportions from wt% oxides."""
    return calc_mol_prop(probe_data).mul(probe_data.ox_num, axis=1)


def calc_ORF(probe_data: probedata, afu: float) -> pd.Series:
    """Calculate Oxygen Renormalisation Factor.

    Parameters
    ----------
    afu : float
        Anions per formula unit e.g. for olivine, afu = 4.
    """
    ox_tot = calc_ox_prop(probe_data).sum(axis=1, skipna=True)
    return afu / ox_tot


def calc_anions(probe_data: probedata, afu: float) -> pd.DataFrame:
    """Calculate the number of anions per formula unit.

    Parameters
    ----------
    afu : float
        Anions per formula unit e.g. for olivine, afu = 4.
    """
    return calc_ox_prop(probe_data).mul(calc_ORF(probe_data, afu), axis=0)


def change_headers_cfu(df: pd.DataFrame, probe_data: probedata) -> pd.DataFrame:
    """Change headers on pd.DataFrame from wt% oxide to cfu."""
    cat_str = probe_data.cat_str[df.columns].to_dict()
    return df.rename(cat_str, axis=1)


def calc_cations(
    probe_data: probedata, afu: float, change_head: bool = True
) -> pd.DataFrame:
    """Calculate cations from wt% oxides.

    Parameters
    ----------
    afu : float
        ideal nions per formula unit e.g. Ol = 4.

    change_head : bool
        True: Change headers to cations (default).
        False: Retain oxide headers (requirement for Pyrolite log transforms).
    """
    cations = calc_anions(probe_data, afu).mul(probe_data.cat_num, axis=1)
    cations = cations.div(probe_data.ox_num, axis=1)

    # remove non-relavent oxide headers
    cations = cations[probe_data.oxides]

    if change_head:
        return change_headers_cfu(cations, probe_data)
    return cations


def calc_cat_tot(probe_data: probedata, afu: float) -> pd.Series:
    """Calculate the cation total per analysis. Add to column cat_tot."""
    cations = calc_cations(probe_data, afu)
    cations["cat_tot"] = cations.sum(axis=1, skipna=True)
    return cations


def check_cat_tot(probe_data: probedata,
                  cfu: float,
                  afu: float,
                  wiggle: float = 0.005) -> list:
    """Check whether the cation total of each analysis lies within a range.

    Returns a list of booleans for boolean indexing.

    Parameters
    ----------
    cfu : float
        ideal cation per formula unit total e.g. for olivine = 3.

    afu : float
        ideal anions per formula unit e.g. for olivine = 4.

    wiggle : float
        fraction of ideal cfu either side of which is acceptable (default = 0.005)
    """
    cat_tot = calc_cat_tot(probe_data, afu)
    upper, lower = cfu + cfu * wiggle, cfu - cfu * wiggle

    cat_tot = [
        True if (x <= upper) & (x >= lower) else False for x in cat_tot["cat_tot"]
    ]

    return cat_tot
