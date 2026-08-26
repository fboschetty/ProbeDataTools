"""
Contains functions to calculate the endmembers of common minerals
e.g. olivine, feldspar, clinopyroxene and spinel.
"""
import numpy as np
import pandas as pd

from .cations import calc_apfu
from .fe_partition import calc_Fe2O3_Droop
from .probe import ProbeData


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
    cations = cations.fillna(0.0)
    total = cations[["Mg", "Fe2", "Mn", "Ca"]].sum(axis=1)

    return pd.DataFrame({
        "Fo": 100 * cations["Mg"] / total,    # Forsterite [Mg2 SiO4]
        "Fay": 100 * cations["Fe2"] / total,  # Fayalite [Fe2 SiO4]
        "Teph": 100 * cations["Mn"] / total,  # Tephroite [Mn2 SiO4]
        "Mont": 100 * cations["Ca"] / total,  # Monticellite [Ca2 SiO4]
    })


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
    cations = cations.fillna(0.0)
    total = cations[["Ca", "Na", "K"]].sum(axis=1)

    return pd.DataFrame({
        "An": 100 * cations["Ca"] / total,  # Anorthite [CaAl2 Si2O8]
        "Ab": 100 * cations["Na"] / total,  # Albite [NaAl Si2O8]
        "Or": 100 * cations["K"] / total    # Orthoclase [KAl Si2O8]
    })


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
    cations = cations.fillna(0.0)
    total = cations[["Fe2", "Mg", "Ca"]].sum(axis=1)

    return pd.DataFrame({
        "Fs": 100 * cations.Fe2 / total,   # Ferrosilite [Fe2 Si2O6]
        "En": 100 * cations.Mg / total,    # Enstatite [Mg2 Si2O6]
        "Wo": 100 * cations.Ca / total     # Wollastonite [Ca2 Si2O6]
    })


def calc_cpx_EM_Putirka(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinopyroxene components using the method of Putirka (2008).

    Fe3+ is calculated following Putirka (2008), Table 3, rather than
    requiring a prior stoichiometric Fe2+/Fe3+ partition.
    """
    components = cations.fillna(0.0).copy()

    # Al site allocation after Papike et al. (1974)
    components["Al_IV"] = 2.0 - components["Si"]
    components["Al_VI"] = components["Al"] - components["Al_IV"]

    components["Fe3"] = components["Na"] + components["Al_IV"] - components["Al_VI"] - 2.0 * components["Ti"] - components["Cr"]

    # Calculate Endmembers
    components["Jd"] = np.minimum(components["Na"], components["Al_VI"])
    components["CaTs"] = components["Al_VI"] - components["Jd"]
    components["CaTi"] = (components["Al_IV"] - components["CaTs"]) / 2.0
    components["CrCaTs"] = components["Cr"] / 2.0
    components["DiHd"] = components["Ca"] - components["CaTi"] - components["CaTs"] - components["CrCaTs"]
    components["EnFs"] = (components["Fe2"] + components["Mg"] - components["DiHd"]) / 2.0

    return components


def _assign_site_sequentially(
    cations: pd.DataFrame,
    site_total: pd.Series,
    cat_sites: list[str],
    con: str,
    rem: str,
    ) -> pd.DataFrame:
    """Assign cations sequentially to a site."""

    components = cations.copy()
    remaining = site_total.copy()

    for cat in cat_sites:
        assigned = np.minimum(remaining, components[cat])
        components[f"{cat}{con}"] = assigned
        components[f"{cat}{rem}"] = components[cat] - assigned
        remaining -= assigned

    return components


def calc_cpx_EM_Neave(probe_data: ProbeData) -> pd.DataFrame:
    """Calculate clinopyroxene endmembers using the method of Neave (2024).
        Uses Droop 1987 to estimate Fe2O3 first.

    Args:
       probe_data (ProbeData): clinopyroxene ProbeData object.

    Returns:
        pd.DataFrame: dataframe including endmembers.
    """

    # Calc Cpx components using Droop Fe3+ calc
    droop_data = calc_Fe2O3_Droop(probe_data, cfu=4.0, afu=6.0, norm=True)
    components = calc_apfu(droop_data, afu=6.0)

    components["Al_IV"] = np.clip(2 - components["Si"], 0., None)
    components["Al_VI"] = np.clip(components["Al"] - components["Al_IV"], 0., None)
    components["Jd"] = np.minimum(components["Na"], components["Al_VI"])

    Na_Remaining = components["Na"] - components["Jd"]
    components["Ae"] = np.minimum(components["Fe3"], Na_Remaining)
    components["Np"] = Na_Remaining - components["Ae"]
    Fe2_remaining = components["Fe2"] - components["Np"] / 2.0

    components["Es"] = components["Fe3"] - components["Ae"]
    components["CaTs"] = components["Al_VI"] - components["Jd"]
    components["CaTi"] = np.clip(components["Ti"] - components["Np"] / 2.0, 0.0, None)
    components["CrAlTs"] = components["Cr"]

    components["DiHd"] = components["Ca"] - components["Es"] - components["CaTs"] - components["CaTi"] - components["CrAlTs"]
    components["EnFs"] = (components["Mg"] + Fe2_remaining + components["Mn"] - components["DiHd"]) / 2.0

    return components


def assign_cpx_sites(cations: pd.DataFrame) -> pd.DataFrame:
    """Assign clinopyroxene cations to crystallographic sites.

    Site allocation follows Morimoto (1998).

    Parameters
    ----------
    cations : pd.DataFrame
        Cations per formula unit with cation headers.

    Returns
    -------
    pd.DataFrame
        Copy of the input cations with T, M1, M2 and excess
        site assignments added.

    References
    ----------
    Morimoto (1998).
    https://doi.org/10.1007/BF01226262
    """

    components = cations.fillna(0.0).copy()

    t_sites = _assign_site_sequentially(
        components,
        2.0 - components["Si"],
        ["Al", "Fe3", "Cr"],
        "_VI",
        "_IV",
    )

    m1_sites = _assign_site_sequentially(
        t_sites,
        pd.Series(1.0, index=t_sites.index),
        ["Al_IV", "Fe3_IV", "Ti", "Cr_IV", "Mg", "Fe2", "Mn"],
        "_M1",
        "_M2",
    )

    return _assign_site_sequentially(
        m1_sites,
        pd.Series(1.0, index=m1_sites.index),
        ["Mg_M2", "Fe2_M2", "Mn_M2", "Ca", "Na"],
        "_M2",
        "_Ex",
    )


_DIETRICH_COMPONENTS = [
    "Al_IV", "Cr_IV", "Al_VI", "Ti", "Cr_VI",
    "Fe3", "Mn", "Fe2", "Mg", "Ca", "Na+K",
]

_DIETRICH_ENDMEMBERS = [
    "Jd", "Ac", "Ur", "Ti-Ts", "Fe-Ts",
    "Cr-Ts", "Ca-Ts", "Pm", "Fs", "En", "Wo",
]

_DIETRICH_COEFFICIENTS = np.array([
    # Jd Ac Ur TiTs FeTs CrTs CaTs Pm Fs En Wo
    [0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0],  # Al_IV
    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],  # Cr_IV
    [1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],  # Al_VI
    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],  # Ti
    [0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0],  # Cr_VI
    [0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0],  # Fe3
    [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0],  # Mn
    [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0],  # Fe2
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0],  # Mg
    [0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2],  # Ca
    [1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],  # Na+K
], dtype=float)


def calc_cpx_EM_Dietrich(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate clinpyroxene endmembers using the method of Dietrich & Petrakakis (1986) [1].

    Parameters
    ----------
    cations : pd.DataFrame
        dataframe of cations with cfu headers.

    Returns
    -------
    EM : pd.DataFrame
        dataframe containing the following 11 linearly independent components.

        Jd: Jadeite [NaAl Si2O6].
        Ac: Acmite [NaFe3+ Si2O6].
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

    components = cations.fillna(0.0).copy()
    components = _assign_site_sequentially(components, 2.0 - cations["Si"], ["Al", "Fe3", "Cr"], "_VI", "_IV")
    components["Na+K"] = components["Na"] + components["K"]

    endmembers = np.linalg.solve(_DIETRICH_COEFFICIENTS, components[_DIETRICH_COMPONENTS].T).T
    return pd.DataFrame(endmembers, columns=_DIETRICH_ENDMEMBERS, index=cations.index)

# Spinel endmembers

_SPINEL_R2_R3 = {
    "MgAl2O4": ("Mg", "Al"),
    "FeAl2O4": ("Fe2", "Al"),
    "MnAl2O4": ("Mn2", "Al"),
    "ZnAl2O4": ("Zn", "Al"),
    "NiAl2O4": ("Ni", "Al"),
    "CuAl2O4": ("Cu", "Al"),
    "MgFe2O4": ("Mg", "Fe3"),
    "FeFe2O4": ("Fe2", "Fe3"),
    "MnFe2O4": ("Mn2", "Fe3"),
    "ZnFe2O4": ("Zn", "Fe3"),
    "NiFe2O4": ("Ni", "Fe3"),
    "CuFe2O4": ("Cu", "Fe3"),
    "FeMn2O4": ("Fe2", "Mn3"),
    "MgMn2O4": ("Mg", "Mn3"),
    "MnMn2O4": ("Mn2", "Mn3"),
    "ZnMn2O4": ("Zn", "Mn3"),
    "MgCr2O4": ("Mg", "Cr"),
    "FeCr2O4": ("Fe2", "Cr"),
    "MnCr2O4": ("Mn2", "Cr"),
    "ZnCr2O4": ("Zn", "Cr"),
    "NiCr2O4": ("Ni", "Cr"),
    "CoCr2O4": ("Co", "Cr"),
    "MgV2O4": ("Mg", "V"),
    "FeV2O4": ("Fe2", "V"),
    "MnV2O4": ("Mn2", "V"),
    "CoCo2O3": ("Co", "Co3"),
}

_SPINEL_R2_R4 = {
    "TiMg2O4": ("Mg", "Ti"),
    "TiFe2O4": ("Fe2", "Ti"),
    "SiMg2O4": ("Mg", "Si"),
    "SiFe2O4": ("Fe2", "Si"),
    "GeFe2O4": ("Fe2", "Ge"),
}

_SPINEL_ENDMEMBERS = {
    **_SPINEL_R2_R3,
    **_SPINEL_R2_R4,
}


def _spinel_cations(cations: pd.DataFrame, species: list[str]) -> pd.DataFrame:
    """Return requested cations with missing species set to zero."""

    return cations.reindex(columns=species, fill_value=0.0).fillna(0.0)


def _calc_spinel_raw_components(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate unnormalised OxyEMG spinel components."""

    components = cations.fillna(0.0).copy()

    return pd.DataFrame({
        name: 4.0 * components.get(r2, 0.0) * components.get(r3, 0.0)
        for name, (r2, r3) in _SPINEL_ENDMEMBERS.items()
        },
        index=components.index,
    )


def calc_sp_EM(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate spinel endmembers using OxyEMG (2024).

    Parameters
    ----------
    cations : pd.DataFrame
        Cations on a 3-cation-per-formula-unit basis with variable-
        valence cations already discriminated.

    Returns
    -------
    pd.DataFrame
        Normalised OxyEMG endmember proportions.
    """

    components = _calc_spinel_raw_components(cations)
    total = components.sum(axis=1)

    if total.eq(0).any():
        raise ValueError("Cannot calculate spinel endmembers from zero cations.")

    return components.mul(8.0 / total, axis=0)


def _calc_spinel_prism_basis(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate R3+, R2+ and R3+/R2+ for prism classification.

    R3+ = Al3+ + Fe3+ + Cr3+
    R2+ = Mg2+ + Fe2+
    """

    R3 = _spinel_cations(cations, ["Al", "Fe3", "Cr"]).sum(axis=1)
    R2 = _spinel_cations(cations, ["Mg", "Fe2"]).sum(axis=1)

    return pd.DataFrame({"R3": R3, "R2": R2, "R3_R2": R3.div(R2.replace(0.0, np.nan))}, index=cations.index)


def _classify_spinel_prisms(cations: pd.DataFrame, ratio_tolerance: float = 0.0) -> tuple[pd.Series, pd.Series]:
    """Classify analyses as magnetite- or ulvospinel-type prisms."""

    basis = _calc_spinel_prism_basis(cations)
    ratio = basis["R3_R2"]

    magnetite = basis["R3"].gt(1.0) & basis["R2"].lt(1.5) & ratio.ge(2.0 / 3.0) & ratio.le(2.0 + ratio_tolerance)
    ulvospinel = basis["R3"].lt(1.0) & basis["R2"].gt(1.5) & ratio.ge(0.0) & ratio.le(2.0 / 3.0)

    return magnetite, ulvospinel


def _normalised_spinel_pairs(cations: pd.DataFrame, divalent: list[str], trivalent: list[str]) -> pd.DataFrame:
    """Calculate normalized divalent/trivalent cation pair proportions."""

    div = _spinel_cations(cations, divalent)
    tri = _spinel_cations(cations, trivalent)

    div = div.div(div.sum(axis=1).replace(0.0, np.nan), axis=0)
    tri = tri.div(tri.sum(axis=1).replace(0.0, np.nan), axis=0)

    return pd.concat([div.add_suffix("_R2"), tri.add_suffix("_R3")], axis=1)


def _spinel_magnetite_components(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate magnetite-prism redistribution components."""

    pairs = _normalised_spinel_pairs(cations, ["Mg", "Fe2"], ["Al", "Fe3", "Cr"])

    return pd.DataFrame({
        "MgAl2O4": pairs["Mg_R2"] * pairs["Al_R3"],
        "FeAl2O4": pairs["Fe2_R2"] * pairs["Al_R3"],
        "MgFe2O4": pairs["Mg_R2"] * pairs["Fe3_R3"],
        "FeFe2O4": pairs["Fe2_R2"] * pairs["Fe3_R3"],
        "MgCr2O4": pairs["Mg_R2"] * pairs["Cr_R3"],
        "FeCr2O4": pairs["Fe2_R2"] * pairs["Cr_R3"],
        },
        index=cations.index,
    ).fillna(0.0)


def _spinel_ulvospinel_components(cations: pd.DataFrame) -> pd.DataFrame:
    """Calculate ulvospinel-prism redistribution components."""

    pairs = _normalised_spinel_pairs(cations, ["Mg", "Fe2"], ["Al", "Ti", "Cr"])

    return pd.DataFrame({
        "MgAl2O4": pairs["Mg_R2"] * pairs["Al_R3"],
        "FeAl2O4": pairs["Fe2_R2"] * pairs["Al_R3"],
        "TiMg2O4": pairs["Mg_R2"] * pairs["Ti_R3"],
        "TiFe2O4": pairs["Fe2_R2"] * pairs["Ti_R3"],
        "MgCr2O4": pairs["Mg_R2"] * pairs["Cr_R3"],
        "FeCr2O4": pairs["Fe2_R2"] * pairs["Cr_R3"],
        },
        index=cations.index,
    ).fillna(0.0)


def calc_sp_prism_EM(cations: pd.DataFrame, ratio_tolerance: float = 0.0) -> pd.DataFrame:
    """Calculate OxyEMG magnetite/ulvospinel prism components."""

    magnetite, ulvospinel = _classify_spinel_prisms(cations, ratio_tolerance=ratio_tolerance)
    magnetite_components = _spinel_magnetite_components(cations)
    ulvospinel_components = _spinel_ulvospinel_components(cations)
    assigned = magnetite | ulvospinel

    return pd.DataFrame({
        "MgAl2O4": np.where(magnetite, magnetite_components["MgAl2O4"], np.where(ulvospinel, ulvospinel_components["MgAl2O4"], 0.0)),
        "FeAl2O4": np.where(magnetite, magnetite_components["FeAl2O4"], np.where(ulvospinel, ulvospinel_components["FeAl2O4"], 0.0)),
        "MgFe2O4": np.where(magnetite, magnetite_components["MgFe2O4"], 0.0),
        "FeFe2O4": np.where(magnetite, magnetite_components["FeFe2O4"], 0.0),
        "TiMg2O4": np.where(ulvospinel, ulvospinel_components["TiMg2O4"], 0.0),
        "TiFe2O4": np.where(ulvospinel, ulvospinel_components["TiFe2O4"], 0.0),
        "MgCr2O4": np.where(assigned, np.where(magnetite, magnetite_components["MgCr2O4"], ulvospinel_components["MgCr2O4"]), 0.0),
        "FeCr2O4": np.where(assigned, np.where(magnetite, magnetite_components["FeCr2O4"], ulvospinel_components["FeCr2O4"]), 0.0),
        },
        index=cations.index
    )
