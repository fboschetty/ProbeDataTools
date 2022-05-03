#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 15:46:58 2022

@author: FelixBoschetty
"""

from typing import List
import pandas as pd
# from dataclasses import dataclass


class ProbeData(object):
    """
        Class that relates input probe data, oxides analysed and useful data about all oxides.

        Parameters
        ----------
        data: pd.DataFrame
            The probe data that you want to process. Ensure headers are wt% oxides.

        oxides: List[str] 
            A list of oxides that you want to use for analysis.

        Returns
        -------
        None.

        """

    # Do I want this to be unchangeable?
    oxide_info = pd.read_csv('oxides.csv', index_col=0)

    def __init__(self, data: pd.DataFrame, oxides: List[str]):

        self.data = data
        self.oxides = oxides

        # Extract useful info from oxides_info using oxides [THIS IS PROBABLY BAD PYTHON. HOW DO I MAKE IT BETTER?]
        self.MR = self.oxide_info.loc['MR'][self.oxides].astype('float')
        self.cat_num = self.oxide_info.loc['cations'][self.oxides].astype(
            'float')
        self.ox_num = self.oxide_info.loc['oxygens'][self.oxides].astype(
            'float')
        self.cat_str = self.oxide_info.loc['cat_str'][self.oxides].astype(
            'str')


class CalcCationsMin(ProbeData):
    """Class containing the functions needed to calculate and check goodnes of cations per formula unit of a mineral"""

    def calc_mol_prop(self) -> pd.DataFrame:
        """Calculate molar proportions from wt% oxides."""

        return self.data.div(self.MR, axis=1)

    def calc_ox_prop(self) -> pd.DataFrame:
        """Calculate anion proportions from wt% oxides."""

        return self.calc_mol_prop().mul(self.ox_num, axis=1)

    def calc_ORF(self, afu) -> pd.Series:
        """Calculate Oxygen Renormalisation Factor.

        Parameters
        ----------
        afu : int or float
            Anions per formula unit e.g. for olivine, afu = 4.
        """

        ox_tot = self.calc_ox_prop().sum(axis=1, skipna=True)
        return afu / ox_tot

    def calc_anions(self, afu) -> pd.DataFrame:
        """Calculate the number of anions per formula unit

        Parameters
        ----------
        afu : int or float
            Anions per formula unit e.g. for olivine, afu = 4.
        """

        return self.calc_ox_prop().mul(self.calc_ORF(afu), axis=0)

    def change_headers_cfu(self, df) -> pd.DataFrame:
        """Change headers on pd.DataFrame from wt% oxide to cfu"""

        cat_str = self.cat_str[self.data.columns].to_dict()
        return df.rename(cat_str, axis=1)

    def calc_cations(self, afu, change_head=True) -> pd.DataFrame:
        """Calculate cations from wt% oxides.

        Parameters
        ----------
        afu : int or float
            ideal nions per formula unit e.g. Ol = 4.

        change_head : bool
            True: Change headers to cations (default).
            False: Retain oxide headers (requirement for Pyrolite log transforms).
        """

        cations = self.calc_anions(afu).mul(self.cat_num, axis=1)
        cations = cations.div(self.ox_num, axis=1)

        if change_head:
            return self.change_headers_cfu(cations)
        return cations

    def calc_cat_tot(self, afu) -> pd.Series:
        """ Calculate the cation total per analysis. Add to column cat_tot"""
        
        cations = self.calc_cations(afu)
        cations["cat_tot"] = cations.sum(axis=1, skipna=True)
        return cations

    def check_cat_tot(self, cfu, afu, wiggle=0.005) -> pd.DataFrame:
        """ Check whether the cation total of each analysis lies within a range. Add column to reflect this.
        
        Parameters
        ----------
        
        cfu : int or float
            ideal cation per formula unit total e.g. for olivine = 3.
        
        afu : int or float
            ideal anions per formula unit e.g. for olivine = 4.
        
        wiggle : float
            fraction of ideal cfu either side of which is acceptable (default = 0.005)
            
        """
        
        cat_tot = self.calc_cat_tot(afu)
        upper, lower = cfu + cfu*wiggle, cfu - cfu*wiggle
        
        cat_tot["cat_good"] = [True if (x <= upper) & (x >= lower) else False for x in cat_tot["cat_tot"]]
        return cat_tot


class CalcCationsFe3(CalcCationsMin):
    """Class containing functions to calculate cations per formula unit of a mineral with Fe3+ estimation methods"""

    def calc_cations_Droop(self, cfu, afu) -> pd.DataFrame:
        """Calculate Fe2/3 ratio using the method of Droop, 1987
         Parameters
        ----------
        cfu : int or float
            ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

        afu : int or float
            ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

        """

        if "FeO" not in self.oxides:
            raise Exception(
                "Dataset must include a column called FeO. Try calc_cations instead."
            )

        if "Fe2O3" in self.oxides:
            raise Exception(
                "Dataset already contains Fe2O3. Try calc_cations instead.")

        mol_fact = self.calc_cat_frac(afu)
        ox_fact = self.calc_cat_frac(cfu)

        F1 = mol_fact.mul(self.calc_cat_prop())
        F2 = self.calc_ox_prop()

        self.data["S"] = F2.sum(axis=1, skipna=True)
        # self.data['T'] = F1.sum(axis=1, skipna=True)
        # self.data['X'] = np.nansum(F2*self.ox_num, axis=1)  # remove nansum

        return self.data

    def calc_cations_Papike(self, divisor):
        """Calculate Fe2/3 ratio using the method of Papike (1947)"""
        return self.calc_mols().mul(self.calc_cat_frac(divisor), axis=0)

    def calc_cations_Lindsey(self, divisor):
        """Calculate Fe2/3 ratio using the method of Lindsley"""
        return self.calc_mols().mul(self.calc_cat_frac(divisor), axis=0)


# %% Test
# Path to Master xls
PathXLS = "/Users/FelixBoschetty/Desktop/Vill_DB_SM_1.xlsx"

# Read in oxides csv
oxide_info = pd.read_csv('oxides.csv', index_col=0)
MR = oxide_info.loc['MR']
cat_num = oxide_info.loc['cations']
rel_ox = oxide_info.loc['oxygens']
cat_str = oxide_info.loc['cat_str']

# Read in sheets replacing strings as nans where needed
Ol_E = pd.read_excel(PathXLS, "Olivine",
                     engine="openpyxl", na_values=["<", "-"])
Feld_E = pd.read_excel(PathXLS, "Feldspar",
                       engine="openpyxl", na_values=["<", "-"])

Ol_ox = ["SiO2", "FeO", "MgO", "MnO", "CaO"]
Feld_ox = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O"]

Feld = ProbeData(Feld_E[Feld_ox], Feld_ox)
Ol = ProbeData(Ol_E[Ol_ox], Ol_ox)

Ol_Cat = CalcCationsMin(Ol_E[Ol_ox], Ol_ox).calc_cations(4, change_head=False)
Feld_Cat = CalcCationsMin(Feld_E[Feld_ox], Feld_ox).calc_cations(32)

# Ol_test = CalcCationsFe3(Ol_E[Ol_ox], Ol_ox).calc_cations_Droop(4, 6)
