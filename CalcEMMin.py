#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
(In Progress) contains relevant functions to calculate the endmembers of common minerals
e.g. olivine, feldspar, clinopyroxene, amphibole and spinels.
@author: FelixBoschetty
"""

import pandas as pd
import numpy as np
from copy import deepcopy
from typing import List


def calc_ol_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate olivine endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Mg, Fe, Ca, and Mn headers

    Returns
    -------
    EM : pd.DataFrame
        contains Forsterite, Fayalite, Tephroite and Monticellite olivine endmembers

    """
    # Ensure No Nans
    cations[cations.isnull()] = 0.

    Fo = 100 * cations.Mg/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)    # Forsterite [Mg2 SiO4]
    Fay = 100 * cations.Fe2/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)  # Fayalite [Fe2 SiO4]
    Teph = 100 * cations.Mn/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)  # Tephroite [Mn2 SiO4]
    Mont = 100 * cations.Ca/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)  # Monticellite [Ca2 SiO4]

    EM = pd.DataFrame({'Fo': Fo, 'Fay': Fay, 'Teph': Teph, 'Mont': Mont})

    return EM


def calc_feld_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate feldspar endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Ca, Na and K headers

    Returns
    -------
    EM : pd.DataFrame
        contains Anorthite, Albite and Orthoclase feldspar endmembers

    """
    # Ensure No Nans
    cations[cations.isnull()] = 0.

    An = 100 * cations.Ca/(cations.Ca+cations.Na+cations.K)  # Anorthite [CaAl2 Si2O8]
    Ab = 100 * cations.Na/(cations.Ca+cations.Na+cations.K)  # Albite [NaAl Si2O8]
    Or = 100 * cations.K/(cations.Ca+cations.Na+cations.K)   # Orthoclase [KAl Si2O8]

    EM = pd.DataFrame({'An': An, 'Ab': Ab, 'Or': Or})

    return EM


def calc_cpx_EM_quad(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinopyroxene QUAD endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Mg, Fe, Ca headers

    Returns
    -------
    EM : pd.DataFrame
        contains Ferrosilite, Enstatite and Wollastonite clinopyroxene endmembers.

    """
    # Ensure No Nans
    cations[cations.isnull()] = 0.

    Fs = 100 * cations.Fe2/(cations.Fe2+cations.Mg+cations.Ca)  # Ferrosilite [Fe2 Si2O6]
    En = 100 * cations.Mg/(cations.Fe2+cations.Mg+cations.Ca)   # Enstatite [Mg2 Si2O6]
    Wo = 100 * cations.Ca/(cations.Fe2+cations.Mg+cations.Ca)   # Wollastonite [Ca2 Si2O6]

    EM = pd.DataFrame({'Fs': Fs, 'En': En, 'Wo': Wo})

    return EM


def calc_cpx_EM_Putirka(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinopyroxene endmembers using the method of Putirka (2008).

    Also calculates Fe3+ so only use on datasets that haven't calculated stoichiometrically.

    """
    EM = deepcopy(cations)

    # Ensure No Nans
    EM[EM.isnull()] = 0.

    # Calculate T-site Cations after Papike et al., 1974
    EM["Al_IV"] = 2 - EM.Si
    EM["Al_VI"] = EM.Al - EM.Al_IV
    EM["Fe3"] = EM.Na + EM.Al_IV - EM.Al_VI - 2*EM.Ti - EM.Cr

    # Calculate Endmembers
    EM["Jd"] = [EM.loc[row].Na if EM.loc[row].Na <= EM.loc[row].Al_VI
                else EM.loc[row].Al_VI
                for row in range(len(EM))]
    EM["CaTs"] = EM.Al_VI - EM.Jd
    EM["CaTi"] = (EM.Al_IV - EM.CaTs)/2
    EM["CrCaTs"] = EM.Cr/2
    EM["DiHd"] = EM.Ca - EM.CaTi - EM.CaTs - EM.CrCaTs
    EM["EnFs"] = (EM.Fe2 + EM.Mg - EM.DiHd)/2

    return EM


def sites(cations: pd.DataFrame,
          total: pd.Series,
          cat_sites: List[str],
          con: str,
          rem: str) -> pd.DataFrame:
    """Assign cations to metal sites. Generic function.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe containing calculated cations w/ cation headers

    total :pd.Series
        The total number of cations in the site, as a series of length Cations

    cat_sites : list[str]
        list of cations that may be contained in site in order of occupancy

    con : str
        string to form cations in site

    Returns
    -------
    df : pd.DataFrame
        dataframe with new columns for sites.

    """
    # Create lists of strings for column headers
    constituents = [cat + con for cat in cat_sites]
    remainder = [cat + rem for cat in cat_sites]
    sites = deepcopy(cations)

    for idx, cat in enumerate(cat_sites):  # for each relevant oxide

        # add columns to sites w/ default 0.
        sites[constituents[idx]] = 0.
        sites[remainder[idx]] = 0.

        for row in range(len(cations)):  # for each row

            if total.loc[row] > 0.:  # if cations still need to be assigned
                sites[constituents[idx]].loc[row] = total.loc[row]

                if sites[constituents[idx]].loc[row] > sites[cat].loc[row]:
                    sites[constituents[idx]].loc[row] = sites[cat].loc[row]

            # update remaining total
            total.loc[row] = total.loc[row] - sites[constituents[idx]].loc[row]

        sites[remainder[idx]] = sites[cat] - sites[constituents[idx]]

    # Prevent long column names
    for col in sites.columns:
        if len(col.split("_")) > 2.:
            sites = sites.drop(col, axis=1)

    return sites


def assign_cpx_sites(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign clinopyroxene cations to metal sites in unit cell after Morimoto (1998) [1].

    Parameters
    ----------
    cations: pd.DataFrame
        dataframe containing calculated cations w/ cation headers.

    Returns
    -------
    M2: pd.DataFrame
        dataframe with cations and filled metal cation sites.

    References
    ----------
    [1] https://doi.org/10.1007/BF01226262

    """
    # Ensure No Nans
    cations[cations.isnull()] = 0.

    # check Si is good
    for row in range(len(cations)):
        if cations.Si[row] < 1.0:
            print(r'Analysis no. %i may be bad, has Si < 1.0.' % row)

        if cations.Si[row] > 2.0:
            print(r'Analysis no. %i may be bad, has Si > 2.0.' % row)

    T = sites(cations=cations,
              total=2.0-cations.Si,
              cat_sites=['Al', 'Fe3', 'Cr'],
              con="_VI",
              rem="_IV")

    M1_tot = pd.Series(np.ones(len(T)))

    M1 = sites(cations=T,
               total=M1_tot,
               cat_sites=['Al_IV', 'Fe3_IV', 'Ti', 'Cr_IV', 'Mg', 'Fe2', 'Mn'],
               con="_M1",
               rem="_M2")

    M2_tot = pd.Series(np.ones(len(T)))

    M2 = sites(cations=M1,
               total=M2_tot,
               cat_sites=['Mg_M2', 'Fe2_M2', 'Mn_M2', 'Ca', 'Na'],
               con='_M2',
               rem='_Ex')

    return M2


def calc_cpx_EM_Dietrich(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinpyroxene endmembers using the method of Dietrich & Petrakasis (1996) [1].

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers.

    Returns
    -------
    EM : pd.DataFrame
        dataframe containing the following 11 linearly independent components.

        Jd: Jadeite [NaAl Si2O6].

        Ae: Aegirine [NaFe3+ Si2O6].

        Ur: Ureyite [NaCr Si2O6].

        Ti-Ts: Ti-Tschermak's [CaTiAlAlO6].

        Ca-Ts: Ca-Tschermak's [CaAlAlSiO6].

        Fe-Ts: Fe-Tschermak's [CaFe3+Fe3+SiO6].

        Cr-Ts: Cr-Tschermak's [CaCrCrSiO6].

        Pm: Pyroxmangite [Mn2 Si2O6].

        Fs: Ferrosilite [Fe2 Si2O6].

        En: Enstatite [Mg2 Si2O6].

        Wo: Wollastonite [Ca2 Si2O6].

    References
    ----------
    [1] https://doi.org/10.1007/BF01191990

    """
    # Ensure No Nans
    cations[cations.isnull()] = 0.

    # Calculate Sites
    SitesT = sites(cations=cations,
                   total=2.0-cations.Si,
                   cat_sites=['Al', 'Fe3', 'Cr'],
                   con="_VI",
                   rem="_IV")

    # Extract Relevant Columns
    SitesT["Na+K"] = SitesT["Na"] + SitesT["K"]
    Cols = ["Al_IV", "Cr_IV", "Al_VI", "Ti", "Cr_VI",
            "Fe3", "Mn", "Fe2", "Mg", "Ca", "Na+K"]
    Sites = SitesT[Cols]

    # Al(IV), Cr(IV), Al(VI), Ti, Cr(VI), Fe3+, Mn, Fe2+, Mg, Ca, Na+k
    Coeff = np.array([[0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0],   # Jd
                      [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],   # Ae
                      [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],   # Ur
                      [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],   # Ti-Ts
                      [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],   # Fe-Ts
                      [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],   # Cr-Ts
                      [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],   # Ca-Ts
                      [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],   # Pm
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],   # Fs
                      [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2],   # En
                      [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0]])  # Wo

    ICoeff = np.linalg.inv(Coeff)

    # Vectorise linear algebra using np.einsum rather than looping through array
    LinComp = np.einsum('ij,kj->ik', Sites, ICoeff)

    EM_head = ['Jd', 'Ae', 'Ur', 'Ti-Ts', 'Fe-Ts',
               'Cr-Ts', 'Ca-Ts', 'Pm', 'Fs', 'En', 'Wo']

    EM = pd.DataFrame(LinComp, columns=EM_head)

    return EM


def assign_amph_sites_Leake1978(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign Amphibole Cations to Metal Sites according to Leake et al., 1978 [1].

    Use when Fe3 has already been calculated.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers.

    Returns
    -------
    M2: pd.DataFrame
         dataframe with cations and filled metal cation sites.

    References
    ----------
    [1] https://pubs.geoscienceworld.org/msa/ammin/article/63/11-12/1023/40814/

    """
    T = sites(cations,
              total=8.0-cations.Si,
              cat_sites=["Al", "Cr", "Fe3", "Ti"],
              con="_T",
              rem="_C")

    C_tot = pd.Series(5*np.ones(len(T)))

    C = sites(T,
              total=C_tot,
              cat_sites=["Al_C", "Cr_C", "Fe3_C", "Ti_C", "Mg", "Fe2", "Mn"],
              con="_C",
              rem="_B")

    B_tot = pd.Series(2*np.ones(len(T)))

    B = sites(C,
              total=B_tot,
              cat_sites=["Mg_B", "Fe2_B", "Mn_B", "Ca", "Na"],
              con="_B",
              rem="_A")

    A_tot = pd.Series(np.ones(len(T)))

    A = sites(B,
              total=A_tot,
              cat_sites=["Na_A", "K", "Vac"],
              con="_A",
              rem="_Ex")

    return A


def assign_amph_sites_Leake1997(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign Amphibole Cations to Metal Sites according to Leake et al., 1997 [1].

    Use when Fe3 is known.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers.

    Returns
    -------
    M2: pd.DataFrame
         dataframe with cations and filled metal cation sites.

    References
    ----------
    [1] https://doi.org/10.1180/minmag.1997.061.405.13

    """
    T = sites(cations,
              total=8.0-cations.Si,
              cat_sites=["Al", "Ti"],
              con="_T",
              rem="_C")

    C_tot = pd.Series(5*np.ones(len(T)))

    C = sites(T,
              total=C_tot,
              cat_sites=["Al_C", "Ti_C", "Cr", "Fe3", "Mn", "Mg", 'Fe2'],
              con="_C",
              rem="_B")

    B_tot = pd.Series(2*np.ones(len(T)))

    B = sites(C,
              total=B_tot,
              cat_sites=["Mg_B", "Fe2_B", "Mn_B", "Ca", "Na"],
              con="_B",
              rem="_A")

    A_tot = pd.Series(np.ones(len(T)))

    A = sites(B,
              total=A_tot,
              cat_sites=["Na_A", "K", "Vac"],
              con="_A",
              rem="_Ex")

    return A


def assign_amph_sites_Leake1997_Fe3(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign Amphibole Cations to Metal Sites according to Leake et al., 1997.

    Will calculate Fe3 using stoichiometry at expense of Fe2. See Rock and Leake et al., 1984.

    """
    return NotImplementedError


def assign_amph_name(sites: pd.DataFrame) -> List[str]:
    """Assign Amphibole Name to Amphibole cations."""
    CANA = sites["Ca_B"] + sites["Na_B"]
    NAB = sites["Na_B"]

    names = []
    for row in range(len(sites)):
        if CANA.loc[row] < 1. and ...:
            names[row] = "Mg-Fe_Mn"
        elif CANA.loc[row] >= 1. and NAB.loc[row] < 0.5:
            names[row] = "Calcic"
        elif CANA.loc[row] >= 1. and 0.5 <= NAB.loc[row] < 1.5:
            names[row] = "Sodic-Calcic"
        elif NAB.loc[row] >= 1.5:
            names[row] = "Sodic"

    return NotImplementedError


def calc_sp_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate spinel endmembers as in Ferracutti et al., (2014).

    Requirements
        1. Cations have been calculated with 3 cations per formula unit.
        2. Fe3+ has been calculated with the method of Droop 1987 unless already known.

    Calculates proportion of the following 19 endmembers:

    Al subgroup [X Al2O4]:
        Sp: Spinel s.s. [Mg Al2O4]
        He: Hercynite [Fe Al2O4]
        Gal: Galaxite [Mn Al2O4]
        Gah: Gahnite [Zn Al2O4]

    Iron subgroup [X Fe3+O2O4]
        MFe: Magnesioferrite [Mg Fe2O4]
        Mag: Magnetite [Fe2+ Fe2O4]
        Jac: Jacobsite [Mn Fe2O4]
        Fra: Franklinite [Zn Fe2O4]
        Tev: Tevolite [Ni Fe2O4]

    Chromium subgroup [X Cr2O4]
        MCr: Magnesiochromite [Mg Cr2O4]
        Cr: Chromite [Fe2+ Cr2O4]
        MnCr: Manganochromite [Mn Cr2O4]
        ZnCr: Zincochromite [Zn Cr2O4]
        NiCr: Nichromite [Ni Cr2O4]

    Vanadium subgroup [X Ni2O4]
        MgCo: Magnesiocoulsonite [Mg V2O4]
        Cou: Coulsonite [Fe2+ V2O4]
        MnCo: Vuorelainenite [Mn V2O4]

    Titanium subgroup [X2 TiO4]
        Qa: Qandilite [Mg2TiO4]
        Ulv: Ulvöspinel [Fe2TiO4]

    """
    raise NotImplementedError
