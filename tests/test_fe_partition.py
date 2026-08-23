import numpy as np
import pandas as pd
import pytest

from probedatatools import ProbeData
from probedatatools.cations import calc_cations
from probedatatools.fe_partition import (
    _calc_Fe3_charge_balance,
    calc_Fe2O3_charge_balance,
    calc_Fe2O3_Droop,
    calc_Fe2O3_Droop_Eq4,
    calc_Fe2O3_Droop_Eq5,
    calc_Fe2O3_Droop_Eq6,
    calc_Fe2O3_Stormer,
    partition_Fe,
)

OXIDES = ["SiO2", "TiO2", "Al2O3", "Cr2O3", "FeO", "MnO", "MgO", "CaO", "ZnO", "Na2O", "K2O"]

DROOP = {
    "GARNET": {
        "afu": 12, "cfu": 8,
        "A": [35.94, 0.12, 0.28, 18.11, 11.37, 0.03, 0, 33.33, 0.19, 0, 0],
        "B": [0.18, 12.43],
    },
    "PYROXENE": {
        "afu": 6, "cfu": 4,
        "A": [56.21, 0, 10.56, 0, 6.10, 0.01, 7.33, 12.13, 0, 7.64, 0.04],
        "B": [3.05, 3.39],
    },
    "SPINEL": {
        "afu": 32, "cfu": 24,
        "A": [0, 0, 61.99, 0, 24.68, 0.11, 13.70, 0, 0, 0, 0],
        "B": [21.00, 4.09],
    },
    "SAPPHIRINE": {
        "afu": 20, "cfu": 14,
        "A": [12.93, 0.03, 60.99, 0, 8.75, 0.04, 17.05, 0, 0.04, 0, 0],
        "B": [5.13, 4.02],
    },
    "AMPHIBOLE 1": {
        "afu": 23, "cfu": 13,
        "A": [42.89, 1.80, 11.55, 0, 14.96, 0, 11.12, 11.49, 0, 1.87, 1.06],
        "B": [12.50, 2.74],
    },
    "AMPHIBOLE 2": {
        "afu": 23, "cfu": 15,
        "A": [43.43, 0.20, 18.26, 0, 12.85, 0.21, 19.67, 0.64, 0, 1.97, 0],
        "B": [12.13, 0.80],
    },
}


def make_probe(name: str) -> ProbeData:
    entry = DROOP[name]
    return ProbeData(pd.DataFrame([entry["A"]], columns=OXIDES), OXIDES)


def make_fe_probe() -> ProbeData:
    return ProbeData(pd.DataFrame({"FeO": [100.0]}), ["FeO"])


def assert_fe_result(result: ProbeData, expected: list[float]) -> None:
    assert result.data.loc[0, "FeO"] == pytest.approx(expected[0], abs=0.03)
    assert result.data.loc[0, "Fe2O3"] == pytest.approx(expected[1], abs=0.03)


@pytest.mark.parametrize(
    "mineral, elements, values",
    [
        ("GARNET", ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca", "Zn"], [3.093, 0.008, 0.028, 1.233, 0.818, 0.002, 0, 3.074, 0.012]),
        ("PYROXENE", ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca", "Na", "K"], [2.013, 0, 0.446, 0, 0.183, 0, 0.391, 0.465, 0.531, 0.002]),
        ("SPINEL", ["Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Zn"], [0, 15.509, 0, 4.381, 0.020, 4.335, 0]),
        ("SAPPHIRINE", ["Si", "Ti", "Al", "Cr", "Fe2", "Mn", "Mg", "Ca", "Zn"], [1.553, 0.003, 8.633, 0, 0.879, 0.004, 3.052, 0, 0.004]),
    ],
)
def test_calc_cations_against_droop_table3(mineral, elements, values):
    result = calc_cations(make_probe(mineral), afu=DROOP[mineral]["afu"])
    assert np.allclose(result.loc[0, elements], values, atol=0.002)


@pytest.mark.parametrize("mineral", ["GARNET", "PYROXENE", "SPINEL", "SAPPHIRINE"])
def test_droop_eq3_against_table3(mineral):
    entry = DROOP[mineral]
    result = calc_Fe2O3_Droop(make_probe(mineral), afu=entry["afu"], cfu=entry["cfu"], norm=True)
    assert_fe_result(result, entry["B"])


def test_droop_eq4_matches_eq3():
    probe_data = make_probe("AMPHIBOLE 1")
    eq3 = calc_Fe2O3_Droop(probe_data, afu=23, cfu=16)
    eq4 = calc_Fe2O3_Droop_Eq4(probe_data)

    assert np.allclose(eq3.data[["FeO", "Fe2O3"]], eq4.data[["FeO", "Fe2O3"]])


@pytest.mark.parametrize("function, mineral", [(calc_Fe2O3_Droop_Eq5, "AMPHIBOLE 2"), (calc_Fe2O3_Droop_Eq6, "AMPHIBOLE 1")])
def test_droop_amphibole_equations(function, mineral):
    assert_fe_result(function(make_probe(mineral), norm=True), DROOP[mineral]["B"])


@pytest.mark.parametrize("mineral", ["GARNET", "PYROXENE", "SPINEL", "SAPPHIRINE"])
def test_droop_preserves_total_fe(mineral):
    probe_data = make_probe(mineral)
    entry = DROOP[mineral]
    result = calc_Fe2O3_Droop(probe_data, afu=entry["afu"], cfu=entry["cfu"])

    original = probe_data.data.loc[0, "FeO"] / probe_data.MR["FeO"]
    corrected = result.data.loc[0, "FeO"] / result.MR["FeO"] + 2 * result.data.loc[0, "Fe2O3"] / result.MR["Fe2O3"]

    assert corrected == pytest.approx(original, rel=1e-4)


def test_droop_norm_false_preserves_oxide_scale():
    probe_data = make_probe("GARNET")
    result = calc_Fe2O3_Droop(probe_data, afu=12, cfu=8)

    for oxide in probe_data.species:
        if oxide != "FeO":
            assert result.data.loc[0, oxide] == pytest.approx(probe_data.data.loc[0, oxide])


def test_droop_norm_true_normalizes_cations():
    result = calc_Fe2O3_Droop(make_probe("GARNET"), afu=12, cfu=8, norm=True)
    assert calc_cations(result, afu=12).sum(axis=1).iloc[0] == pytest.approx(8.0)


def test_droop_norm_changes_fe_partition():
    probe_data = make_probe("GARNET")
    raw = calc_Fe2O3_Droop(probe_data, afu=12, cfu=8)
    norm = calc_Fe2O3_Droop(probe_data, afu=12, cfu=8, norm=True)

    raw_ratio = raw.data.loc[0, "Fe2O3"] / raw.data.loc[0, "FeO"]
    norm_ratio = norm.data.loc[0, "Fe2O3"] / norm.data.loc[0, "FeO"]

    assert norm_ratio != pytest.approx(raw_ratio)


@pytest.mark.parametrize("Fe3ideal, FeO, Fe3_fraction", [(0.25, 75.0, 0.25), (-0.25, 100.0, 0.0), (1.5, 0.0, 1.0)])
def test_partition_Fe(Fe3ideal, FeO, Fe3_fraction):
    probe_data = make_fe_probe()
    result = partition_Fe(probe_data, pd.DataFrame({"Fe2": [1.0]}), pd.Series([Fe3ideal]))
    factor = probe_data.MR["Fe2O3"] / (2 * probe_data.MR["FeO"])

    assert result.data.loc[0, "FeO"] == pytest.approx(FeO)
    assert result.data.loc[0, "Fe2O3"] == pytest.approx(100 * Fe3_fraction * factor)
    assert "Fe2O3" in result.species


@pytest.mark.parametrize(
    "cations, charges, afu, expected",
    [
        ({"SiO2": 1.0, "MgO": 0.5, "FeO": 1.0}, {"SiO2": 4.0, "MgO": 2.0, "FeO": 2.0}, 4.0, 1.0),
        ({"MgO": 2.0, "FeO": 1.0}, {"MgO": 2.0, "FeO": 2.0}, 2.0, 0.0),
        ({"SiO2": 0.0, "FeO": 1.0}, {"SiO2": 4.0, "FeO": 2.0}, 2.0, 1.0),
    ],
)
def test_calc_Fe3_charge_balance(cations, charges, afu, expected):
    result = _calc_Fe3_charge_balance(pd.DataFrame([cations]), pd.Series(charges), afu)
    assert result.loc[0] == pytest.approx(expected)


def test_calc_Fe2O3_charge_balance_pure_FeO():
    result = calc_Fe2O3_charge_balance(make_fe_probe(), afu=1.0)
    assert result.data.loc[0, "FeO"] == pytest.approx(100.0)
    assert result.data.loc[0, "Fe2O3"] == pytest.approx(0.0)


def test_calc_Fe2O3_charge_balance_returns_probe_data():
    result = calc_Fe2O3_charge_balance(make_fe_probe(), afu=1.0)
    assert isinstance(result, ProbeData)
    assert {"FeO", "Fe2O3"} <= set(result.data.columns)
    assert "Fe2O3" in result.species


def test_stormer_pure_FeO():
    probe_data = make_fe_probe()
    result = calc_Fe2O3_Stormer(probe_data)
    factor = probe_data.MR["Fe2O3"] / (2 * probe_data.MR["FeO"])

    assert result.data.loc[0, "FeO"] == pytest.approx(100 / 3)
    assert result.data.loc[0, "Fe2O3"] == pytest.approx(200 / 3 * factor)


def test_stormer_ulvospinel():
    FeO, TiO2 = ProbeData.MR["FeO"], ProbeData.MR["TiO2"]
    total = 2 * FeO + TiO2

    probe_data = ProbeData(
        pd.DataFrame({"TiO2": [100 * TiO2 / total], "FeO": [100 * 2 * FeO / total]}),
        ["TiO2", "FeO"],
    )

    result = calc_Fe2O3_Stormer(probe_data)
    assert result.data.loc[0, "Fe2O3"] == pytest.approx(0.0)


def test_stormer_returns_probe_data():
    result = calc_Fe2O3_Stormer(make_fe_probe())
    assert isinstance(result, ProbeData)
    assert "Fe2O3" in result.data.columns
    assert "Fe2O3" in result.species