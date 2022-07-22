#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains relevant functions to calculate the Fe2/3 ratio of analyses of iron-containing
minerals via the method of Droop (1987).

It also contains functions to calculate the Fe2/3 of specific minerals:
    (1) pyroxene using the method of Papike 1947.
    (2) spinel using the method of Stormer 1983.

@author: FelixBoschetty
"""
from copy import deepcopy

from ProbeData import probedata
from CalcCationsMin import calc_mol_prop, calc_cations


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 16:42:11 2022

@author: FelixBoschetty
"""
from copy import deepcopy

from CalcCationsMin import calc_cations, calc_mol_prop
from ProbeData import ProbeData


def calc_Fe2O3_Droop(probe_data: ProbeData, cfu: float, afu: float) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
      Parameters
    ----------
    cfu : float
        ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

    afu : float
        ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

    """

    # Calculate for 6 oxygens
    cpxCations = calc_cations(probe_data, afu=afu)
    S = cpxCations.sum(axis=1)

    # Calculate for 4 Cations
    MolProp = calc_mol_prop(probe_data)
    MolPropTot = MolProp.sum(axis=1, skipna=True)
    Factor = cfu/MolPropTot
    normMolProp = MolProp.multiply(Factor, axis=0)
    T = normMolProp.sum(axis=1)

    # Calculate new FeO and Fe2O3 contents
    Fe3ideal = 2*afu*(1-(T/S))
    Fe3 = [x if x > 0. else 0. for x in Fe3ideal]  # Check amount of Fe3+ isn't negative

    Fe2 = cpxCations['Fe2'] - Fe3  # Calculate new Fe2
    Fe2FeT = Fe2/(Fe2+Fe3)
    FeO = probe_data.data['FeO'] * Fe2FeT
    Fe2O3 = probe_data.data['FeO'] * (1-Fe2FeT) * 1.1113

    probe_data_new = deepcopy(probe_data)
    probe_data_new.data['FeO'] = FeO
    probe_data_new.data['Fe2O3'] = Fe2O3
    probe_data_new.oxides.append('Fe2O3')

    return probe_data_new



def calc_Fe2O3_Droop_Eq5(probe_data: ProbeData, cfu: float = 15.0, afu: float = 23.0) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
    Equation 5 tailored for Fe-Mg amphiboles. Requires calculation of cations on the basis of 23.0 oxygens, anhydrous.
    Using 15.0 cations exclusive of Na and K.

      Parameters
    ----------
    cfu : float
        ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

    afu : float
        ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

    """

    raise NotImplementedError()

    # Calculate for 6 oxygens
    cpxCations = calc_cations(probe_data, afu=afu)
    S = cpxCations.sum(axis=1)

    # Calculate for 4 Cations
    catnorm = cfu/S
    cpx4Cations = cpxCations.mul(catnorm, axis=0)
    T = cpx4Cations.sum(axis=1)

    # Calculate new FeO and Fe2O3 contents
    Fe3ideal = 2*afu*(1-(T/S))
    Fe3 = [x if x > 0. else 0. for x in Fe3ideal]  # Check amount of Fe3+ isn't negative

    Fe2 = cpxCations['Fe2'] - Fe3  # Calculate new Fe2
    Fe2FeT = Fe2/(Fe2+Fe3)
    FeO = probe_data.data['FeO'] * Fe2FeT
    Fe2O3 = probe_data.data['FeO'] * (1-Fe2FeT) * 1.1113

    probe_data_new = deepcopy(probe_data)
    probe_data_new.data['FeO'] = FeO
    probe_data_new.data['Fe2O3'] = Fe2O3
    probe_data_new.oxides.append('Fe2O3')

    return probe_data_new


def calc_Fe2O3_Droop_Eq6(probe_data: ProbeData, cfu: float = 13.0, afu: float = 23.0) -> ProbeData:
    """Calculate Fe2/3 ratio stoichiometrically using the method of Droop, 1987.
    Equation 6 tailored for many calcic amphiboles. Requires cations on the basis of 23.0 oxygens, anhydrous.
    Using 13.0 cations, exclusive of Ca, Na and K.

      Parameters
    ----------
    cfu : float
        ideal cations per formula unit e.g. for clinopyroxene, cfu = 4.

    afu : float
        ideal anions per formula unit e.g. for clinopyroxene, afu = 6.

    """

    raise NotImplementedError()

    # Calculate for 6 oxygens
    cpxCations = calc_cations(probe_data, afu=afu)
    S = cpxCations.sum(axis=1)

    # Calculate for 4 Cations
    catnorm = cfu/S
    cpx4Cations = cpxCations.mul(catnorm, axis=0)
    T = cpx4Cations.sum(axis=1)

    # Calculate new FeO and Fe2O3 contents
    Fe3ideal = 2*afu*(1-(T/S))
    Fe3 = [x if x > 0. else 0. for x in Fe3ideal]  # Check amount of Fe3+ isn't negative

    Fe2 = cpxCations['Fe2'] - Fe3  # Calculate new Fe2
    Fe2FeT = Fe2/(Fe2+Fe3)
    FeO = probe_data.data['FeO'] * Fe2FeT
    Fe2O3 = probe_data.data['FeO'] * (1-Fe2FeT) * 1.1113

    probe_data_new = deepcopy(probe_data)
    probe_data_new.data['FeO'] = FeO
    probe_data_new.data['Fe2O3'] = Fe2O3
    probe_data_new.oxides.append('Fe2O3')

    return probe_data_new

# def calc_Fe2O3_Papike(probe_data: ProbeData, afu: float) -> ProbeData:
#     """Calculate Fe2/3 ratio of Pyroxene using the method of Papike et al., (1947)"""

#     cations = calc_cations(probe_data, afu=afu)

#     AlIV = 2. - cations.Si
#     AlVI = [x if x > 0. else 0. for x in cations.Al - AlIV]
#     Fe3 = [x if x > 0 else 0. for x in cations.Na + AlIV - AlVI - cations.Ti - cations.Cr]

#     probe_data_new = calc_New_FeO_Fe2O3(probe_data, Fe3)

#     return probe_data_new

