#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 16:42:11 2022

@author: FelixBoschetty
"""

import pandas as pd

from ProbeData import ProbeData
import CalcCationsMin


def calc_cations_Droop(self, cfu, afu) -> pd.DataFrame:
    """Calculate Fe2/3 ratio using the method of Droop, 1987
      Parameters
    ----------
    cfu : int or float
        ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

    afu : int or float
        ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

    """

    if "FeO" not in oxides:
        raise Exception(
            "Dataset must include a column called FeO. Try calc_cations instead."
        )

    if "Fe2O3" in oxides:
        raise Exception(
            "Dataset already contains Fe2O3. Try calc_cations instead.")
    
    # Calculate cation and oxygen proportions 
    mol_prop = cat_num.mul(calc_mol_prop())
    mol_prop_tot = mol_prop.sum(axis=1, skipna=True)
    mol_fact = cfu / mol_prop_tot
    
    ox_no = mol_prop.mul(ox_num, axis=1)
    ox_tot = ox_no.sum(axis=1, skipna=True)
    ox_fact = afu / ox_tot

    F1 = calc_mol_prop().mul(mol_fact, axis=0)
    F2 = calc_mol_prop().mul(ox_fact, axis=0)
    
    # Calculate Droop Components S, T and X THESE ARE NOT GIVING RIGHT VALUES
    S = F1.sum(axis=1, skipna=True)
    T = F2.sum(axis=1, skipna=True)
    OxNum6 = F2.mul(ox_num/cat_num, axis=1)
    X = OxNum6.sum(axis=1, skipna=True)
    # Calculate N
    
    # check S/T and X/N are equal
    
    # Calculate new FeO and Fe2O3 contents
    
    
    
    # Fe3, Fe2 = np.zeros(len(S)), np.zeros(len(S))
    # for idx in range(len(S)):
    #     if 2*X[idx]*(1-(T[idx]/S[idx])) > 0.:
    #         Fe3[idx] = 2*X[idx]*(1-T[idx]/S[idx])
    #     else:
    #         Fe3[idx] = 0.
    # Fe2[idx] = F1[4][idx] - Fe3[idx]

    # Fe2FeT = Fe2/(Fe2+Fe3)

    # FeO = CpxOrg[4] * Fe2FeT
    # Fe2O3 = CpxOrg[4] * (1-Fe2FeT) * 1.1113
    # Recalculate Cations

    return F1, F2, S, T, X

def calc_cations_Papike(self, divisor):
    """Calculate Fe2/3 ratio using the method of Papike (1947)"""
    return calc_mols().mul(calc_cat_frac(divisor), axis=0)

def calc_cations_Lindsey(self, divisor):
    """Calculate Fe2/3 ratio using the method of Lindsley"""
    return calc_mols().mul(calc_cat_frac(divisor), axis=0)