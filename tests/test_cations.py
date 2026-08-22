import numpy as np
import pandas as pd
import pytest

from probedatatools import ProbeData
from probedatatools.cations import (
    calc_ARF,
    calc_cat_frac,
    calc_cat_prop,
    calc_cat_tot,
    calc_cations,
    calc_mol_frac,
    calc_mol_prop,
    calc_nonO_anions,
    calc_oxygen_prop,
    check_cat_tot,
)


@pytest.fixture
def olivine_probe():
    """Exact Mg-endmember olivine: Mg2SiO4."""

    mgO, siO2 = ProbeData.MR["MgO"], ProbeData.MR["SiO2"]
    formula_mass = 2 * mgO + siO2

    data = pd.DataFrame({
        "SiO2": [100 * siO2 / formula_mass],
        "MgO": [100 * 2 * mgO / formula_mass],
    })

    return ProbeData(data, ["SiO2", "MgO"])


@pytest.fixture
def halogen_probe():
    """Mg2SiO4 composition with F and Cl added."""

    mgO, siO2 = ProbeData.MR["MgO"], ProbeData.MR["SiO2"]
    formula_mass = 2 * mgO + siO2

    data = pd.DataFrame({
        "SiO2": [98 * siO2 / formula_mass],
        "MgO": [98 * 2 * mgO / formula_mass],
        "F": [1.0],
        "Cl": [1.0],
    })

    return ProbeData(data, ["SiO2", "MgO", "F", "Cl"])


def test_calc_mol_prop(olivine_probe):
    result = calc_mol_prop(olivine_probe)
    assert np.allclose(result["SiO2"], olivine_probe.data["SiO2"] / olivine_probe.MR["SiO2"])
    assert np.allclose(result["MgO"], olivine_probe.data["MgO"] / olivine_probe.MR["MgO"])


def test_calc_mol_frac(olivine_probe):
    result = calc_mol_frac(olivine_probe)
    assert np.allclose(result.sum(axis=1), 1.0)


def test_calc_oxygen_prop(olivine_probe):
    result = calc_oxygen_prop(olivine_probe)
    mol = calc_mol_prop(olivine_probe)

    assert np.allclose(result["SiO2"], mol["SiO2"] * olivine_probe.ox_num["SiO2"])
    assert np.allclose(result["MgO"], mol["MgO"] * olivine_probe.ox_num["MgO"])


def test_calc_cat_prop(olivine_probe):
    result = calc_cat_prop(olivine_probe)
    mol = calc_mol_prop(olivine_probe)

    assert np.allclose(result["SiO2"], mol["SiO2"] * olivine_probe.cat_num["SiO2"])
    assert np.allclose(result["MgO"], mol["MgO"] * olivine_probe.cat_num["MgO"])


def test_calc_cat_frac(olivine_probe):
    result = calc_cat_frac(olivine_probe)
    assert np.allclose(result.sum(axis=1), 1.0)


def test_calc_ARF_oxygen_basis(olivine_probe):
    arf = calc_ARF(olivine_probe, afu=4)
    oxygen_total = calc_oxygen_prop(olivine_probe).sum(axis=1)

    assert np.allclose(arf, 4 / oxygen_total)


def test_calc_ARF_includes_F_and_Cl(halogen_probe):
    O = calc_ARF(halogen_probe, afu=4)
    OFCl = calc_ARF(halogen_probe, afu=4, anion_species=("O", "F", "Cl"))

    assert np.all(OFCl < O)


def test_calc_ARF_default_is_oxygen(olivine_probe):
    assert np.allclose(
        calc_ARF(olivine_probe, 4),
        calc_ARF(olivine_probe, 4, ("O",)),
    )


def test_calc_nonO_anions(halogen_probe):
    result = calc_nonO_anions(halogen_probe, 4, ("O", "F", "Cl"))
    arf = calc_ARF(halogen_probe, 4, ("O", "F", "Cl"))

    assert np.allclose(result["F"], halogen_probe.data["F"] / halogen_probe.MR["F"] * arf)
    assert np.allclose(result["Cl"], halogen_probe.data["Cl"] / halogen_probe.MR["Cl"] * arf)


def test_calc_nonO_anions_oxygen_only_is_empty(olivine_probe):
    assert calc_nonO_anions(olivine_probe, 4).empty


def test_calc_cations_exact_olivine(olivine_probe):
    result = calc_cations(olivine_probe, 4)

    assert result.loc[0, "Mg"] == pytest.approx(2.0)
    assert result.loc[0, "Si"] == pytest.approx(1.0)


def test_calc_cations_total_exact_olivine(olivine_probe):
    result = calc_cat_tot(olivine_probe, 4)
    assert result.loc[0, "cat_tot"] == pytest.approx(3.0)


def test_F_and_Cl_do_not_affect_O_only_cations(olivine_probe, halogen_probe):
    olivine = calc_cations(olivine_probe, 4)
    halogen = calc_cations(halogen_probe, 4)

    assert halogen.loc[0, "Mg"] == pytest.approx(olivine.loc[0, "Mg"])
    assert halogen.loc[0, "Si"] == pytest.approx(olivine.loc[0, "Si"])


def test_F_and_Cl_change_cation_normalization(halogen_probe):
    O = calc_cations(halogen_probe, 4)
    OFCl = calc_cations(halogen_probe, 4, ("O", "F", "Cl"))

    assert OFCl.loc[0, "Mg"] < O.loc[0, "Mg"]
    assert OFCl.loc[0, "Si"] < O.loc[0, "Si"]


def test_calc_cations_can_keep_oxide_headers(olivine_probe):
    result = calc_cations(olivine_probe, 4, change_head=False)
    assert list(result.columns) == ["SiO2", "MgO"]


def test_check_cat_tot_exact_analysis_passes(olivine_probe):
    assert check_cat_tot(olivine_probe, cfu=3, afu=4) == [True]


def test_check_cat_tot_wrong_cfu_fails(olivine_probe):
    assert check_cat_tot(olivine_probe, cfu=3.1, afu=4) == [False]