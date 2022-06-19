#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 19 11:31:13 2022

@author: FelixBoschetty
"""

import pandas as pd
from ProbeData import probedata
from CalcCationsMin import calc_mol_prop


def calc_mol_frac(probe_data: probedata) -> pd.DataFrame:
    """Calculate molar fractions from wt% oxides."""
    mol_prop = calc_mol_prop(probe_data)[probe_data.oxides]
    return mol_prop.div(mol_prop.sum(axis=1, skipna=True), axis=0)


def calc_cat_frac(probe_data: probedata) -> pd.DataFrame:
    """Calculate cation fractions from wt% oxides."""
    cat_prop = calc_mol_prop(probe_data) * probe_data.cat_num
    return cat_prop.div(cat_prop.sum(axis=1, skipna=True), axis=0)[probe_data.oxides]
