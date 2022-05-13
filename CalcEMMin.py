#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 12 13:12:55 2022

@author: FelixBoschetty
"""

import pandas as pd
import numpy as np


def calc_ol_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """ Calculate olivine endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Mg, Fe, Ca, and Mn headers

    Returns
    -------
    EM : pd.DataFrame
        contains Forsterite, Fayalite, Tephroite and Monticellite olivine endmembers
  
    """
   
    Fo   = 100 * cations.Mg/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)   # Forsterite [Mg2 SiO4]
    Fay  = 100 * cations.Fe2/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)  # Fayalite [Fe2 SiO4]
    Teph = 100 * cations.Mn/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)   # Tephroite [Mn2 SiO4]
    Mont = 100 * cations.Ca/(cations.Mg+cations.Fe2+cations.Mn+cations.Ca)   # Monticellite [Ca2 SiO4]
    
    EM = pd.DataFrame({'Fo': Fo, 'Fay': Fay, 'Teph': Teph, 'Mont': Mont})
    
    return EM


def calc_feld_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """ Calculate feldspar endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Ca, Na and K headers

    Returns
    -------
    EM : pd.DataFrame
        contains Anorthite, Albite and Orthoclase feldspar endmembers
  
    """
    
    An = 100 * cations.Ca/(cations.Ca+cations.Na+cations.K)  # Anorthite [CaAl2 Si2O8]
    Ab = 100 * cations.Ca/(cations.Ca+cations.Na+cations.K)  # Albite [NaAl Si2O8]
    Or = 100 * cations.Ca/(cations.Ca+cations.Na+cations.K)  # Orthoclase [KAl Si2O8]
    
    EM = pd.DataFrame({'An': An, 'Ab': Ab, 'Or': Or})
    
    return EM 


def calc_cpx_EM_quad(cations: pd.DataFrame) -> pd.DataFrame:
    """ Calculate clinopyroxene QUAD endmembers.

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers. requires Mg, Fe, Ca headers

    Returns
    -------
    EM : pd.DataFrame
        contains Ferrosilite, Enstatite, Wollastonite, Diopside and Hedenbergite clinopyroxene endmembers
  
    """
    
    Fs = 100 * cations.Fe/(cations.Fe+cations.Mg+cations.Ca)  # Ferrosilite [Fe2 Si2O6]
    En = 100 * cations.Mg/(cations.Fe+cations.Mg+cations.Ca)  # Enstatite [Mg2 Si2O6]
    Wo = 100 * cations.Ca/(cations.Fe+cations.Mg+cations.Ca)  # Wollastonite [Ca2 Si2O6]
    
    Di = 100 * 0.5*(cations.Ca+cations.Mg)/(cations.Fe+cations.Mg+cations.Ca)  # Diopside [CaMg Si2O6]
    Hd = 100 * 0.5*(cations.Ca+cations.Fe)/(cations.Fe+cations.Mg+cations.Ca)  # Hedenbergite [CaFe Si2O6]
    
    EM = pd.DataFrame({'Fs': Fs, 'En': En, 'Wo': Wo, 'Di': Di, 'Hd': Hd})
    
    return EM

def calc_cpx_EM_Putirka(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinopyroxene endmembers using the method of Putirka (2008)"""
    raise NotImplementedError

def assign_cpx_sites(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign clinopyroxene cations to metal sites in unit cell"""
    raise NotImplementedError

def calc_cpx_EM_Dietrich(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinpyroxene endmembers using the method of Dietrich & Petrakasis (1996).
    
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
    
    """
    
    # Calculate Sites
    Sites = assign_cpx_sites().to_numpy()
    
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
    # return EM        
    raise NotImplementedError


def calc_sp_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate spinel endmembers as in Ferracutti et al., (2014)
    
    Calculates proportion of 19 endmembers.
    
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