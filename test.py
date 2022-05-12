#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 16:50:53 2022

@author: FelixBoschetty
"""

import pandas as pd

from ProbeData import ProbeData
from CalcCationsMin import calc_cations, check_cat_tot
from CalcCationsMinFe3 import calc_Fe2O3_Droop, calc_Fe2O3_Papike

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
Ol_Cat = calc_cations(probe_data=Ol, afu=4.)                          # Change Headers (default)
Feld_Cat = calc_cations(probe_data=Feld, afu=32., change_head=False)  # Keep Oxide Headers           
Cpx_Cat_Fe2 = calc_cations(Cpx, 6.) 

# Check whether Ol Cfu are good
Ol_Cfu = check_cat_tot(probe_data=Ol, afu=4., cfu=3.)

# Calc Cpx Cations using Droop Fe3+ calc
Cpx_Droop = calc_Fe2O3_Droop(Cpx, cfu=4., afu=6.)
Cpx_Cat_Droop = calc_cations(Cpx_Droop, afu=6.)

# Calc CpxCations using Papike Fe3+ calc
Cpx_Papike = calc_Fe2O3_Papike(Cpx, afu=6.)
Cpx_Cat_Papike = calc_cations(Cpx_Papike, afu=6.)