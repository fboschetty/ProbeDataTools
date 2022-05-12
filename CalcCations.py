#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 15:46:58 2022

@author: FelixBoschetty
"""

import pandas as pd

from ProbeData import ProbeData


def calc_mol_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate molar proportions from wt% oxides."""

    return probe_data.data.div(probe_data.MR, axis=1)

def calc_ox_prop(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate anion proportions from wt% oxides."""

    return calc_mol_prop(probe_data).mul(probe_data.ox_num, axis=1)

def calc_ORF(probe_data: ProbeData, afu: float) -> pd.Series:
    """Calculate Oxygen Renormalisation Factor.

    Parameters
    ----------
    afu : float
        Anions per formula unit e.g. for olivine, afu = 4.
    """

    ox_tot = calc_ox_prop(probe_data).sum(axis=1, skipna=True)
    return afu / ox_tot

def calc_anions(probe_data: ProbeData, afu: float) -> pd.DataFrame:
    """Calculate the number of anions per formula unit

    Parameters
    ----------
    afu : float
        Anions per formula unit e.g. for olivine, afu = 4.
    """

    return calc_ox_prop(probe_data).mul(calc_ORF(probe_data, afu), axis=0)

def change_headers_cfu(df: pd.DataFrame, probe_data: ProbeData) -> pd.DataFrame:
    """Change headers on pd.DataFrame from wt% oxide to cfu"""

    cat_str = probe_data.cat_str[df.columns].to_dict()
    return df.rename(cat_str, axis=1)

def calc_cations(probe_data: ProbeData, afu: float, change_head: bool = True) -> pd.DataFrame:
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

def calc_cat_tot(probe_data: ProbeData, afu: float) -> pd.Series:
    """ Calculate the cation total per analysis. Add to column cat_tot"""
    
    cations = calc_cations(probe_data, afu)
    cations["cat_tot"] = cations.sum(axis=1, skipna=True)
    return cations

def check_cat_tot(probe_data: ProbeData,
                  cfu: float,
                  afu: float,
                  wiggle: float = 0.005) -> pd.DataFrame:
    """ Check whether the cation total of each analysis lies within a range. Add column to reflect this.
    
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
    upper, lower = cfu + cfu*wiggle, cfu - cfu*wiggle
    
    cat_tot["cat_good"] = [True if (x <= upper) & (x >= lower) else False
                           for x in cat_tot["cat_tot"]]
    return cat_tot

# %% Test
# Path to Master xls
PathXLS = "/Users/FelixBoschetty/Desktop/Vill_DB_SM_1.xlsx"

# Read in sheets replacing strings as nans where needed
Ol_E   = pd.read_excel(PathXLS, "Olivine",
                       engine="openpyxl", na_values=["<", "-"])
Feld_E = pd.read_excel(PathXLS, "Feldspar",
                       engine="openpyxl", na_values=["<", "-"])
Cpx_E  = pd.read_excel(PathXLS, "Clinopyroxene",
                       engine="openpyxl", na_values=["<", "-"])

# Oxides for each dataset
Ol_ox = ["SiO2", "FeO", "Cr2O3", "MgO", "MnO", "NiO", "CaO"]
Feld_ox = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O"]
Cpx_ox = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MgO", "MnO", "CaO", "Na2O"]

# Construct ProbeData objects for each dataset
Ol = ProbeData(Ol_E, Ol_ox)
Feld = ProbeData(Feld_E, Feld_ox)
Cpx = ProbeData(Cpx_E, Cpx_ox)

# Calc Cations using Generic Class
Ol_Cat = calc_cations(probe_data=Ol, afu=4., change_head=False)   # Keep Oxide Headers
Feld_Cat = calc_cations(Feld, 32.)  # Change Headers
Cpx_Cat_Fe2 = calc_cations(Cpx, afu=6.) 

# Calc Cations using Droop Fe3+ calc
# Cpx_test = CalcCationsFe3(Cpx_E[Cpx_ox], Cpx_ox)

# Calc Cations using Papike Fe3+ calc
