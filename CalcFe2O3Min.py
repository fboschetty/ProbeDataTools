#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 16:42:11 2022

@author: FelixBoschetty
"""
from copy import deepcopy

from ProbeData import probedata
from CalcCationsMin import calc_mol_prop, calc_cations


def calc_New_FeO_Fe2O3(probe_data: probedata, Fe3: list) -> probedata:
    """Makes a deepcopy of the probedata object and appends new FeO and Fe2O3 values"""
    Fe2 = probe_data.data['FeO'] - Fe3

    Fe2FeT = Fe2/(Fe2+Fe3)

    FeO = probe_data.data['FeO'] * Fe2FeT
    Fe2O3 = probe_data.data['FeO'] * (1-Fe2FeT) * 1.1113

    probe_data_new = deepcopy(probe_data)
    probe_data_new.data['FeO'] = FeO
    probe_data_new.data['Fe2O3'] = Fe2O3
    probe_data_new.oxides.append('Fe2O3')

    return probe_data_new


def calc_Fe2O3_Droop(probe_data: probedata, cfu: float, afu: float) -> probedata:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
      Parameters
    ----------
    cfu : float
        ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

    afu : float
        ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

    """

    # Calculate cation and oxygen proportions
    mol_prop = probe_data.cat_num.mul(calc_mol_prop(probe_data))
    mol_prop_tot = mol_prop.sum(axis=1, skipna=True)
    mol_fact = cfu / mol_prop_tot

    rel_ox = probe_data.ox_num / probe_data.cat_num  # Relative number of oxygens
    ox_no = mol_prop.mul(rel_ox, axis=1)
    ox_tot = ox_no.sum(axis=1, skipna=True)
    ox_fact = afu / ox_tot

    # Calculate Factor 1 and 2
    F1 = mol_prop.mul(mol_fact, axis=0)
    F2 = mol_prop.mul(ox_fact, axis=0)

    # Calculate Droop Components S, T, X and N
    S = F2.sum(axis=1, skipna=True)
    T = F1.sum(axis=1, skipna=True)

    OxNum6 = F2.mul(rel_ox, axis=1)
    X = OxNum6.sum(axis=1, skipna=True)

    OxNum4 = F1.mul(rel_ox, axis=1)
    N = OxNum4.sum(axis=1, skipna=True)

    # check S/T and X/N are equal
    if all(S/T - X/N) < 0.00001:
        pass
    else:
        raise Exception("S/T not equal to X/N, something has gone wrong!")

    # Calculate new FeO and Fe2O3 contents
    test = 2*X*(1-(T/S))
    Fe3 = [x if x > 0. else 0. for x in test]

    probe_data_new = calc_New_FeO_Fe2O3(probe_data, Fe3)

    return probe_data_new


def calc_Fe2O3_Papike(probe_data: probedata, afu: float) -> probedata:
    """Calculate Fe2/3 ratio of Pyroxene using the method of Papike et al., (1947)"""

    cations = calc_cations(probe_data, afu=afu)

    AlIV = 2. - cations.Si
    AlVI = [x if x > 0. else 0. for x in cations.Al - AlIV]
    Fe3 = [x if x > 0 else 0. for x in cations.Na + AlIV - AlVI - cations.Ti - cations.Cr]

    probe_data_new = calc_New_FeO_Fe2O3(probe_data, Fe3)

    return probe_data_new


def calc_sp_Fe3_Stormer(probe_data: probedata, afu: float) -> probedata:
    """Calculate Fe3+ for Spinel cations using the method of Stormer, 1983"""

    cations = calc_cations(probe_data, afu=afu)

    charge = cations * probedata.cat_chrg
    charge_tot = charge.sum(axis=1, skipna=True)

    charge_diff = 8. - (charge_tot - charge.Fe2)
    Fe2 = 3*cations.Fe2 - charge_diff
    Fe3 = cations.Fe2 - Fe2

    probe_data_new = calc_New_FeO_Fe2O3(probe_data, Fe3)

    return probe_data_new
