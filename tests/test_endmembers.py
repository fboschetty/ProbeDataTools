from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from probedatatools import ProbeData
from probedatatools.cations import calc_apfu
from probedatatools.endmembers import (
    _DIETRICH_COEFFICIENTS,
    _DIETRICH_COMPONENTS,
    _DIETRICH_ENDMEMBERS,
    _assign_site_sequentially,
    _calc_spinel_prism_basis,
    _classify_spinel_prisms,
    assign_cpx_sites,
    calc_cpx_EM_Dietrich,
    calc_cpx_EM_Neave,
    calc_cpx_EM_Putirka,
    calc_cpx_EM_quad,
    calc_feld_EM,
    calc_ol_EM,
    calc_sp_EM,
    calc_sp_prism_EM,
)
from probedatatools.fe_partition import calc_Fe2O3_Droop

DATA_DIR = Path(__file__).parent / "test_data"

CATION_COLS = [
    "Si", "Ti", "Ge", "Al", "Cr", "V",
    "Fe3+", "Mn3+", "Co3+",
    "Fe2+", "Mn2+", "Mg", "Ca", "Zn",
    "Ni", "Cu", "Co",
]

CATION_RENAME = {
    "Fe2+": "Fe2",
    "Fe3+": "Fe3",
    "Mn2+": "Mn2",
    "Mn3+": "Mn3",
    "Co3+": "Co3",
}

PUTIRKA_COMPONENTS = [
    "Al_IV", "Al_VI", "Fe3", "Jd", "CaTs",
    "CaTi", "CrCaTs", "DiHd", "EnFs",
]

NEAVE_OXIDES = [
    "SiO2", "TiO2", "Al2O3", "Cr2O3",
    "FeO", "MnO", "MgO", "CaO",
    "Na2O", "K2O", "P2O5", "NiO",
]

NEAVE_CATIONS = ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca", "Na"]
NEAVE_DROOP = ["Si", "Ti", "Al", "Cr", "Fe3", "Fe2", "Mn", "Mg", "Ca", "Na"]
NEAVE_EM = ["Jd", "Ae", "Np", "Es", "CaTs", "CaTi", "CrAlTs", "DiHd", "EnFs"]


# Simple mineral endmembers

@pytest.mark.parametrize(
    "function, cations, expected",
    [
        (
            calc_ol_EM,
            {"Mg": 0.60, "Fe2": 0.25, "Mn": 0.05, "Ca": 0.10},
            {"Fo": 60, "Fay": 25, "Teph": 5, "Mont": 10},
        ),
        (
            calc_feld_EM,
            {"Ca": 0.50, "Na": 0.30, "K": 0.20},
            {"An": 50, "Ab": 30, "Or": 20},
        ),
        (
            calc_cpx_EM_quad,
            {"Fe2": 0.20, "Mg": 0.50, "Ca": 0.30},
            {"Fs": 20, "En": 50, "Wo": 30},
        ),
    ],
)
def test_simple_endmembers(function, cations, expected):
    result = function(pd.DataFrame([cations]))

    assert np.allclose(result.loc[0, list(expected)], list(expected.values()))
    assert result.loc[0].sum() == pytest.approx(100.0)


@pytest.mark.parametrize(
    "function, cations",
    [
        (calc_ol_EM, {"Mg": 0.6, "Fe2": 0.3, "Mn": np.nan, "Ca": 0.1}),
        (calc_feld_EM, {"Ca": 0.5, "Na": np.nan, "K": 0.5}),
        (calc_cpx_EM_quad, {"Fe2": 0.2, "Mg": np.nan, "Ca": 0.8}),
    ],
)
def test_simple_endmembers_handle_nan(function, cations):
    result = function(pd.DataFrame([cations]))

    assert np.isfinite(result.to_numpy()).all()
    assert result.loc[0].sum() == pytest.approx(100.0)


# Clinopyroxene site assignment

def test_assign_cpx_sites_T_site():
    cations = pd.DataFrame({
        "Si": [1.80], "Ti": [0.10], "Al": [0.20], "Cr": [0.05], "Fe3": [0.10],
        "Fe2": [0.20], "Mn": [0.05], "Mg": [0.70], "Ca": [0.60], "Na": [0.10],
    })

    result = assign_cpx_sites(cations)

    assert np.allclose(result.loc[0, ["Al_VI", "Fe3_VI", "Cr_VI"]], [0.20, 0.0, 0.0])


def test_assign_cpx_sites_conserves_cations():
    cations = pd.DataFrame({
        "Si": [1.95], "Ti": [0.02], "Al": [0.20], "Cr": [0.01], "Fe3": [0.03],
        "Fe2": [0.20], "Mn": [0.02], "Mg": [0.70], "Ca": [0.80], "Na": [0.06],
    })

    result = assign_cpx_sites(cations)

    assert result.loc[0, "Mg_M1"] + result.loc[0, "Mg_M2"] == pytest.approx(cations.loc[0, "Mg"])
    assert result.loc[0, "Fe2_M1"] + result.loc[0, "Fe2_M2"] == pytest.approx(cations.loc[0, "Fe2"])
    assert result.loc[0, "Mn_M1"] + result.loc[0, "Mn_M2"] == pytest.approx(cations.loc[0, "Mn"])
    assert result.loc[0, "Ca_M2"] + result.loc[0, "Ca_Ex"] == pytest.approx(cations.loc[0, "Ca"])
    assert result.loc[0, "Na_M2"] + result.loc[0, "Na_Ex"] == pytest.approx(cations.loc[0, "Na"])


# OxyEMG spinel endmembers

@pytest.fixture
def oxyemg_tables():
    s2 = pd.read_excel(DATA_DIR / "Ferracutti2024_Table S2.xlsx", header=2)
    s3 = pd.read_excel(DATA_DIR / "Ferracutti2024_Table S3.xlsx", header=2)
    s5 = pd.read_excel(DATA_DIR / "Ferracutti2024_Table S5.xlsx", header=3).rename(columns={"Unnamed: 0": "Sample"})
    return s2, s3, s5


def make_oxyemg_cations(s2):
    return s2[CATION_COLS].rename(columns=CATION_RENAME)


def test_calc_sp_EM_against_oxyemg(oxyemg_tables):
    s2, s3, _ = oxyemg_tables
    result = calc_sp_EM(make_oxyemg_cations(s2))

    assert np.allclose(
        result.to_numpy(),
        s3[result.columns].to_numpy(),
        atol=1e-6,
    )


def test_calc_sp_EM_sums_to_8(oxyemg_tables):
    s2, _, _ = oxyemg_tables
    result = calc_sp_EM(make_oxyemg_cations(s2))

    assert np.allclose(result.sum(axis=1), 8.0, atol=1e-10)


def test_calc_sp_EM_zero_cations_raises():
    cations = pd.DataFrame({"Mg": [0.0], "Al": [0.0], "Fe2": [0.0], "Fe3": [0.0]})

    with pytest.raises(ValueError, match="zero cations"):
        calc_sp_EM(cations)


# Spinel prism classification

def test_spinel_prism_basis():
    cations = pd.DataFrame({
        "Al": [1.0], "Fe3": [0.5], "Cr": [0.5],
        "Mg": [0.4], "Fe2": [0.6],
        "Mn2": [1.0], "Mn3": [1.0], "V": [1.0],
    })

    result = _calc_spinel_prism_basis(cations)

    assert result.loc[0, "R3"] == pytest.approx(2.0)
    assert result.loc[0, "R2"] == pytest.approx(1.0)
    assert result.loc[0, "R3_R2"] == pytest.approx(2.0)


@pytest.mark.parametrize(
    "R3, R2, expected_mag, expected_ulv",
    [
        (1.5, 1.0, True, False),
        (2.0, 1.0, True, False),
        (1.0, 1.0, False, False),
        (1.5, 1.5, False, False),
        (0.5, 2.0, False, True),
        (0.9, 1.6, False, True),
        (1.0, 1.6, False, False),
        (0.9, 1.5, False, False),
    ],
)
def test_spinel_prism_masks(R3, R2, expected_mag, expected_ulv):
    cations = pd.DataFrame({"Al": [R3], "Mg": [R2]})

    magnetite, ulvospinel = _classify_spinel_prisms(cations)

    assert magnetite.iloc[0] == expected_mag
    assert ulvospinel.iloc[0] == expected_ulv


def test_spinel_prism_excludes_extra_divalent_cations():
    cations = pd.DataFrame({
        "Al": [0.0], "Fe3": [2.0], "Cr": [0.0],
        "Mg": [0.0], "Fe2": [0.57], "Mn2": [0.43],
    })

    basis = _calc_spinel_prism_basis(cations)
    magnetite, ulvospinel = _classify_spinel_prisms(cations)

    assert basis.loc[0, "R2"] == pytest.approx(0.57)
    assert basis.loc[0, "R3_R2"] == pytest.approx(2.0 / 0.57)
    assert not magnetite.iloc[0]
    assert not ulvospinel.iloc[0]


def test_spinel_prism_magnetite():
    cations = pd.DataFrame({
        "Mg": [0.5], "Fe2": [0.5], "Al": [1.0],
        "Fe3": [1.0], "Cr": [0.0], "Ti": [0.0],
    })

    result = calc_sp_prism_EM(cations)

    assert np.allclose(
        result.loc[0, ["MgAl2O4", "FeAl2O4", "MgFe2O4", "FeFe2O4"]],
        0.25,
    )
    assert np.allclose(result.loc[0, ["TiMg2O4", "TiFe2O4"]], 0.0)
    assert result.loc[0].sum() == pytest.approx(1.0)


def test_spinel_prism_ulvospinel():
    cations = pd.DataFrame({
        "Mg": [1.0], "Fe2": [1.0], "Al": [0.5],
        "Fe3": [0.0], "Cr": [0.0], "Ti": [0.5],
    })

    result = calc_sp_prism_EM(cations)

    assert np.allclose(
        result.loc[0, ["MgAl2O4", "FeAl2O4", "TiMg2O4", "TiFe2O4"]],
        0.25,
    )
    assert np.allclose(result.loc[0, ["MgFe2O4", "FeFe2O4"]], 0.0)
    assert result.loc[0].sum() == pytest.approx(1.0)


def test_spinel_prism_unassigned():
    cations = pd.DataFrame({
        "Mg": [0.0], "Fe2": [0.5], "Mn2": [0.5],
        "Al": [0.0], "Fe3": [2.0], "Cr": [0.0], "Ti": [0.0],
    })

    result = calc_sp_prism_EM(cations)

    assert np.allclose(result.loc[0], 0.0)


def test_spinel_prism_against_oxyemg(oxyemg_tables):
    s2, _, s5 = oxyemg_tables
    cations = make_oxyemg_cations(s2)

    magnetite, ulvospinel = _classify_spinel_prisms(cations, ratio_tolerance=0.01)
    result = calc_sp_prism_EM(cations, ratio_tolerance=0.01)

    expected = pd.DataFrame(0.0, index=s2.index, columns=result.columns)

    mag_columns = {
        "MgAl2O4": "MgAl2O4",
        "FeAl2O4": "FeAl2O4",
        "MgFe2O4": "MgFe2O4",
        "FeFe2O4": "FeFe2O4",
        "MgCr2O4": "MgCr2O4",
        "FeCr2O4": "FeCr2O4",
    }
    ulv_columns = {
        "MgAl2O4": "MgAl2O4.1",
        "FeAl2O4": "FeAl2O4.1",
        "TiMg2O4": "TiMg2O4",
        "TiFe2O4": "TiFe2O4",
        "MgCr2O4": "MgCr2O4.1",
        "FeCr2O4": "FeCr2O4.1",
    }

    for output, source in mag_columns.items():
        expected.loc[magnetite, output] = s5.loc[magnetite, source]

    for output, source in ulv_columns.items():
        expected.loc[ulvospinel, output] = s5.loc[ulvospinel, source]

    eligible = magnetite | ulvospinel
    assert np.allclose(result.loc[eligible], expected.loc[eligible], atol=1e-6)


def test_spinel_ulvospinel_oktyabrsky(oxyemg_tables):
    s2, _, s5 = oxyemg_tables
    sample = s2["Sample"].eq("#5 Oktyabrsky et al. (1992)")
    result = calc_sp_prism_EM(
        make_oxyemg_cations(s2.loc[sample]),
        ratio_tolerance=0.01,
    )

    expected = s5.loc[sample, ["MgAl2O4.1", "FeAl2O4.1", "TiMg2O4", "TiFe2O4", "MgCr2O4.1", "FeCr2O4.1"]].to_numpy()
    actual = result[["MgAl2O4", "FeAl2O4", "TiMg2O4", "TiFe2O4", "MgCr2O4", "FeCr2O4"]].to_numpy()

    assert np.allclose(actual, expected, atol=1e-6)


# Neave (2024) clinopyroxene components

@pytest.fixture
def neave_data():
    data = pd.read_excel(DATA_DIR / "Neave2024_Supplement.xlsx", header=[0, 1])
    data.columns = [
        "_".join(str(part).strip() for part in column if str(part) != "nan")
        for column in data.columns
    ]
    return data


def make_neave_probe(row: pd.Series) -> ProbeData:
    data = pd.DataFrame([[row[f"Oxides_{oxide}"] for oxide in NEAVE_OXIDES]], columns=NEAVE_OXIDES)
    return ProbeData(data, NEAVE_OXIDES)


def test_neave_cations_against_supplement(neave_data):
    results = []

    for _, row in neave_data.iterrows():
        cations = calc_apfu(make_neave_probe(row), afu=6.0)
        results.append(cations.loc[0, NEAVE_CATIONS])

    result = pd.DataFrame(results, columns=NEAVE_CATIONS)
    expected = neave_data[[f"Cations_{column}" for column in NEAVE_CATIONS]].copy()
    expected.columns = NEAVE_CATIONS

    assert np.allclose(result.to_numpy(), expected.to_numpy(), atol=0.001)


def test_neave_droop_against_supplement(neave_data):
    results = []

    for _, row in neave_data.iterrows():
        corrected = calc_Fe2O3_Droop(make_neave_probe(row), cfu=4.0, afu=6.0, norm=True)
        cations = calc_apfu(corrected, afu=6.0)
        results.append(cations.loc[0, NEAVE_DROOP])

    result = pd.DataFrame(results, columns=NEAVE_DROOP)
    expected = neave_data[[f"Droop_{column}" for column in NEAVE_DROOP]].copy()
    expected.columns = NEAVE_DROOP

    assert np.allclose(result.to_numpy(), expected.to_numpy(), atol=0.002)


def test_neave_droop_cations_sum_to_four(neave_data):
    for _, row in neave_data.iterrows():
        corrected = calc_Fe2O3_Droop(make_neave_probe(row), cfu=4.0, afu=6.0, norm=True)
        cations = calc_apfu(corrected, afu=6.0)
        assert cations.loc[0, NEAVE_DROOP].sum() == pytest.approx(4.0, abs=0.005)


def test_neave_endmembers_against_supplement(neave_data):
    results = []

    for _, row in neave_data.iterrows():
        components = calc_cpx_EM_Neave(make_neave_probe(row))
        results.append(components.loc[0, NEAVE_EM])

    result = pd.DataFrame(results, columns=NEAVE_EM)
    expected = neave_data[[f"Endmembers_{column}" for column in NEAVE_EM]].copy()
    expected.columns = NEAVE_EM

    assert np.allclose(result.to_numpy(), expected.to_numpy(), atol=0.002)


# Dietrich & Petrakakis (1996) clinopyroxene components

@pytest.fixture
def dietrich_data():
    return pd.read_excel(DATA_DIR / "DietrichPetrakakis1986_T3.xlsx", header=[0, 1])


def test_dietrich_coefficients_against_table(dietrich_data):
    components = pd.DataFrame({
        "Al_IV": dietrich_data[("T", "Al_IV")],
        "Cr_IV": dietrich_data[("T", "CR_IV")],
        "Al_VI": dietrich_data[("M1", "Al_VI")],
        "Ti": dietrich_data[("M1", "Ti")],
        "Cr_VI": dietrich_data[("M1", "Cr_VI")],
        "Fe3": dietrich_data[("M1", "Fe3")],
        "Mn": dietrich_data[("M1", "Mn")] + dietrich_data[("M2", "Mn")],
        "Fe2": dietrich_data[("M1", "Fe2")] + dietrich_data[("M2", "Fe2")],
        "Mg": dietrich_data[("M1", "Mg")] + dietrich_data[("M2", "Mg")],
        "Ca": dietrich_data[("M2", "Ca")],
        "Na+K": dietrich_data[("M2", "Na")] + dietrich_data[("M2", "K")],
    })

    expected = dietrich_data[
        [
            ("Components", "Jd"),
            ("Components", "Ac"),
            ("Components", "Ur"),
            ("Components", "Ti Ts"),
            ("Components", "Fe Ts"),
            ("Components", "Cr Ts"),
            ("Components", "CaTs"),
            ("Components", "Pm"),
            ("Components", "Fs"),
            ("Components", "En"),
            ("Components", "Wo"),
        ]
    ].to_numpy()

    endmembers = np.linalg.solve(
        _DIETRICH_COEFFICIENTS,
        components[_DIETRICH_COMPONENTS].T,
    ).T

    assert np.allclose(endmembers, expected, atol=0.005, rtol=0.0)


def test_calc_cpx_EM_Dietrich_reconstructs_components():
    cations = pd.DataFrame({
        "Si": [2.0],
        "Ti": [0.1],
        "Al": [0.5],
        "Cr": [0.05],
        "Fe3": [0.4],
        "Fe2": [0.2],
        "Mn": [0.05],
        "Mg": [0.7],
        "Ca": [0.9],
        "Na": [0.3],
        "K": [0.0],
    })

    components = _assign_site_sequentially(
        cations,
        2.0 - cations["Si"],
        ["Al", "Fe3", "Cr"],
        "_VI",
        "_IV",
    )
    components["Na+K"] = components["Na"] + components["K"]

    endmembers = calc_cpx_EM_Dietrich(cations)

    assert np.allclose(
        endmembers.to_numpy() @ _DIETRICH_COEFFICIENTS.T,
        components[_DIETRICH_COMPONENTS].to_numpy(),
        atol=1e-10,
    )


def test_calc_cpx_EM_Dietrich_handles_nan():
    cations = pd.DataFrame({
        "Si": [2.017], "Ti": [0.096], "Al": [0.010], "Cr": [0.000],
        "Fe3": [0.689], "Fe2": [0.167], "Mn": [np.nan], "Mg": [0.011],
        "Ca": [0.048], "Na": [0.923], "K": [np.nan],
    })

    result = calc_cpx_EM_Dietrich(cations)

    assert np.isfinite(result[_DIETRICH_ENDMEMBERS].to_numpy()).all()


# Putirka (2008) clinopyroxene components

def make_putirka_cations() -> pd.DataFrame:
    return pd.DataFrame({
        "Si": [1.9553], "Ti": [0.0092], "Al": [0.0694], "Fe2": [0.3080],
        "Mn": [0.0076], "Mg": [0.7379], "Ca": [0.8978], "Na": [0.0267],
        "K": [0.0000], "Cr": [0.0015],
    })


def test_calc_cpx_EM_Putirka_sites():
    result = calc_cpx_EM_Putirka(make_putirka_cations())

    assert result.loc[0, "Al_IV"] == pytest.approx(0.0447, abs=0.0001)
    assert result.loc[0, "Al_VI"] == pytest.approx(0.0247, abs=0.0001)
    assert result.loc[0, "Fe3"] == pytest.approx(0.0268, abs=0.0001)


def test_calc_cpx_EM_Putirka_against_table3():
    result = calc_cpx_EM_Putirka(make_putirka_cations())
    expected = [0.0247, 0.0000, 0.0223, 0.0007, 0.8747, 0.0856]

    assert np.allclose(
        result.loc[0, ["Jd", "CaTs", "CaTi", "CrCaTs", "DiHd", "EnFs"]].to_numpy(),
        expected,
        atol=0.0001,
        rtol=0.0,
    )


def test_calc_cpx_EM_Putirka_returns_components():
    result = calc_cpx_EM_Putirka(make_putirka_cations())
    assert set(PUTIRKA_COMPONENTS) <= set(result.columns)


def test_calc_cpx_EM_Putirka_does_not_modify_input():
    cations = make_putirka_cations()
    original = cations.copy()

    calc_cpx_EM_Putirka(cations)

    pd.testing.assert_frame_equal(cations, original)


def test_calc_cpx_EM_Putirka_handles_nan():
    cations = make_putirka_cations()
    cations.loc[0, "Cr"] = np.nan

    result = calc_cpx_EM_Putirka(cations)

    assert np.isfinite(result.to_numpy()).all()