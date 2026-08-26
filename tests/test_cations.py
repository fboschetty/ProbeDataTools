import numpy as np
import pandas as pd
import pytest

from probedatatools import ProbeData
from probedatatools.cations import (
    calc_apfu,
    calc_cat_frac,
    calc_cat_prop,
    calc_cat_tot,
    calc_mol_frac,
    calc_mol_prop,
    calc_nonO_anions,
    calc_ORF,
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

    assert np.allclose(
        result["SiO2"],
        olivine_probe.data["SiO2"] / olivine_probe.MR["SiO2"],
    )
    assert np.allclose(
        result["MgO"],
        olivine_probe.data["MgO"] / olivine_probe.MR["MgO"],
    )


def test_calc_mol_frac(olivine_probe):
    result = calc_mol_frac(olivine_probe)

    assert np.allclose(result.sum(axis=1), 1.0)


def test_calc_oxygen_prop(olivine_probe):
    result = calc_oxygen_prop(olivine_probe)
    mol = calc_mol_prop(olivine_probe)

    assert np.allclose(
        result["SiO2"],
        mol["SiO2"] * olivine_probe.ox_num["SiO2"],
    )
    assert np.allclose(
        result["MgO"],
        mol["MgO"] * olivine_probe.ox_num["MgO"],
    )


def test_calc_cat_prop(olivine_probe):
    result = calc_cat_prop(olivine_probe)
    mol = calc_mol_prop(olivine_probe)

    assert np.allclose(
        result["SiO2"],
        mol["SiO2"] * olivine_probe.cat_num["SiO2"],
    )
    assert np.allclose(
        result["MgO"],
        mol["MgO"] * olivine_probe.cat_num["MgO"],
    )


def test_calc_cat_frac(olivine_probe):
    result = calc_cat_frac(olivine_probe)

    assert np.allclose(result.sum(axis=1), 1.0)


def test_calc_ORF_oxygen_basis(olivine_probe):
    orf = calc_ORF(olivine_probe, afu=4)
    oxygen_total = calc_oxygen_prop(olivine_probe).sum(axis=1)

    assert np.allclose(orf, 4 / oxygen_total)


def test_calc_ORF_ignores_F_and_Cl(halogen_probe):
    """F and Cl do not affect oxygen-basis normalisation."""

    oxygen_total = calc_oxygen_prop(halogen_probe).sum(axis=1)
    expected = 4 / oxygen_total

    assert np.allclose(
        calc_ORF(halogen_probe, 4),
        expected,
    )


def test_calc_nonO_anions(halogen_probe):
    result = calc_nonO_anions(halogen_probe, 4)

    orf = calc_ORF(halogen_probe, 4)

    expected_F = (
        halogen_probe.data["F"]
        / halogen_probe.MR["F"]
        * halogen_probe.an_num_use["F"]
        * orf
    )
    expected_Cl = (
        halogen_probe.data["Cl"]
        / halogen_probe.MR["Cl"]
        * halogen_probe.an_num_use["Cl"]
        * orf
    )

    assert np.allclose(result["F"], expected_F)
    assert np.allclose(result["Cl"], expected_Cl)


def test_calc_nonO_anions_oxygen_only_is_empty(olivine_probe):
    assert calc_nonO_anions(olivine_probe, 4).empty


def test_calc_nonO_anions_only_returns_anions_present(
    olivine_probe,
    halogen_probe,
):
    olivine_result = calc_nonO_anions(olivine_probe, 4)
    halogen_result = calc_nonO_anions(halogen_probe, 4)

    assert olivine_result.empty
    assert list(halogen_result.columns) == ["Cl", "F"]


def test_calc_nonO_anions_use_oxygen_ORF(halogen_probe):
    """F and Cl are normalised using the oxygen-only ORF."""

    result = calc_nonO_anions(halogen_probe, 4)

    orf = calc_ORF(halogen_probe, 4)

    expected = (
        calc_mol_prop(halogen_probe)[["F", "Cl"]]
        .mul(
            halogen_probe.an_num_use[["F", "Cl"]],
            axis=1,
        )
        .mul(orf, axis=0)
    )

    pd.testing.assert_frame_equal(
        result.sort_index(axis=1),
        expected.sort_index(axis=1),
    )


def test_calc_nonO_anions_expected_values(halogen_probe):
    result = calc_nonO_anions(halogen_probe, 4)
    orf = calc_ORF(halogen_probe, 4).iloc[0]

    assert result.loc[0, "F"] == pytest.approx(
        1.0 / ProbeData.MR["F"] * orf
    )
    assert result.loc[0, "Cl"] == pytest.approx(
        1.0 / ProbeData.MR["Cl"] * orf
    )


def test_calc_apfu_exact_olivine(olivine_probe):
    result = calc_apfu(olivine_probe, 4)

    assert result.loc[0, "Mg"] == pytest.approx(2.0)
    assert result.loc[0, "Si"] == pytest.approx(1.0)


def test_calc_apfu_automatically_reports_anions(halogen_probe):
    result = calc_apfu(halogen_probe, 4)

    assert "F" in result.columns
    assert "Cl" in result.columns

    assert result.loc[0, "F"] > 0
    assert result.loc[0, "Cl"] > 0


def test_calc_apfu_anions_do_not_affect_cation_normalization(halogen_probe):
    result = calc_apfu(halogen_probe, 4)

    assert result.loc[0, "Mg"] == pytest.approx(2.0)
    assert result.loc[0, "Si"] == pytest.approx(1.0)


def test_calc_apfu_without_anions_has_only_cations(olivine_probe):
    result = calc_apfu(olivine_probe, 4)

    assert list(result.columns) == ["Si", "Mg"]


def test_calc_cat_tot_exact_olivine(olivine_probe):
    result = calc_cat_tot(olivine_probe, 4)

    assert result.loc[0, "cat_tot"] == pytest.approx(3.0)


def test_calc_cat_tot_includes_anions_but_excludes_them_from_total(
    halogen_probe,
):
    result = calc_cat_tot(halogen_probe, 4)

    assert "F" in result.columns
    assert "Cl" in result.columns

    cation_total = result.loc[0, ["Si", "Mg"]].sum()

    assert result.loc[0, "cat_tot"] == pytest.approx(cation_total)

    assert result.loc[0, "cat_tot"] != pytest.approx(
        cation_total + result.loc[0, "F"] + result.loc[0, "Cl"]
    )


def test_calc_apfu_can_keep_oxide_headers(olivine_probe):
    result = calc_apfu(
        olivine_probe,
        4,
        change_head=False,
    )

    assert list(result.columns) == ["SiO2", "MgO"]


def test_calc_apfu_with_anions_can_keep_oxide_headers(halogen_probe):
    result = calc_apfu(
        halogen_probe,
        4,
        change_head=False,
    )

    assert "SiO2" in result.columns
    assert "MgO" in result.columns
    assert "F" in result.columns
    assert "Cl" in result.columns


def test_check_cat_tot_exact_analysis_passes(olivine_probe):
    assert check_cat_tot(
        olivine_probe,
        cfu=3,
        afu=4,
    ) == [True]


def test_check_cat_tot_wrong_cfu_fails(olivine_probe):
    assert check_cat_tot(
        olivine_probe,
        cfu=3.1,
        afu=4,
    ) == [False]