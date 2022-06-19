#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This pulls mineral data from Vill_DB_SM_1.xlsx.
Then calculates the cations per formula unit of olivine, plagioclase and clinopyroxene data.

@author: FelixBoschetty
"""

import pandas as pd

from ProbeData import probedata
from CalcCationsMin import calc_cations, check_cat_tot
from CalcFe2O3Min import calc_Fe2O3_Droop, calc_Fe2O3_Papike
import CalcEMMin


# IMPORT AND PREPARE DATA #

# Path to Master xls
PathXLS = "Vill_DB_SM_1.xlsx"

# Read in sheets replacing strings as nans where needed
Ol_E = pd.read_excel(PathXLS, "Olivine", engine="openpyxl", na_values=["<", "-"])
Feld_E = pd.read_excel(PathXLS, "Feldspar", engine="openpyxl", na_values=["<", "-"])
Cpx_E = pd.read_excel(PathXLS, "Clinopyroxene", engine="openpyxl", na_values=["<", "-"])

# Oxides for each dataset
Ol_ox = ["SiO2", "FeO", "Cr2O3", "MgO", "MnO", "NiO", "CaO"]
Feld_ox = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O"]
Cpx_ox = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MgO", "MnO", "CaO", "Na2O", "K2O"]

# Construct ProbeData objects for each dataset
Ol = probedata(Ol_E, Ol_ox)
Feld = probedata(Feld_E, Feld_ox)
Cpx = probedata(Cpx_E, Cpx_ox)


# CALCULATE MINERAL CATIONS #

# Calc Cations using Generic Class
Ol_Cat = calc_cations(probe_data=Ol, afu=4.)
Feld_Cat = calc_cations(probe_data=Feld, afu=32.)
Cpx_Cat_Fe2 = calc_cations(Cpx, 6.)

# Calc Cpx Cations using Droop Fe3+ calc
Cpx_Droop = calc_Fe2O3_Droop(Cpx, cfu=4., afu=6.)
Cpx_Cat_Droop = calc_cations(Cpx_Droop, afu=6.)

# Calc CpxCations using Papike Fe3+ calc
Cpx_Papike = calc_Fe2O3_Papike(Cpx, afu=6.)
Cpx_Cat_Papike = calc_cations(Cpx_Papike, afu=6.)

# Check which Cpx Cfu Totals are good
Cpx_Cfu_Good = check_cat_tot(probe_data=Cpx_Droop, cfu=4., afu=6., wiggle=0.005)


# ASSIGN CATIONS TO METAL SITES #

Cpx_Sites = CalcEMMin.assign_cpx_sites(Cpx_Cat_Droop)


# CALCULATE MINERAL ENDMEMBERS #

# Calculate Ol EM
Ol_EM = CalcEMMin.calc_ol_EM(Ol_Cat)

# Calculate Feld EM
Feld_EM = CalcEMMin.calc_feld_EM(Feld_Cat)

# Calculate Cpx Quad EM
Cpx_Em_Quad = CalcEMMin.calc_cpx_EM_quad(Cpx_Cat_Droop)

# Calculagte Cpx EM using Putirka Method
Cpx_EM_Putirka = CalcEMMin.calc_cpx_EM_Putirka(Cpx_Cat_Fe2)

# Calculate Cpx EM using Dietrech method
Cpx_EM_Ditrich = CalcEMMin.calc_cpx_EM_Dietrich(Cpx_Cat_Droop)
