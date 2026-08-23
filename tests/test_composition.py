import pandas as pd
import pytest

from probedatatools import ProbeData
from probedatatools.composition import convert_units, element_to_oxide, oxide_to_element

# Unit conversion

@pytest.mark.parametrize(
    "from_unit, to_unit, factor",
    [("ppm", "wt%", 1e-4), ("wt%", "ppm", 1e4)],
)
def test_convert_units(from_unit, to_unit, factor):
    data = pd.DataFrame({"SiO2": [50.0], "MgO": [10.0], "Sample": ["A"]})
    probe = ProbeData(data, ["SiO2", "MgO"])

    result = convert_units(probe, from_unit, to_unit)

    assert result.data.loc[0, "SiO2"] == pytest.approx(50.0 * factor)
    assert result.data.loc[0, "MgO"] == pytest.approx(10.0 * factor)
    assert result.data.loc[0, "Sample"] == "A"
    assert result.species == probe.species


def test_convert_units_same_units_returns_copy():
    data = pd.DataFrame({"SiO2": [50.0], "Sample": ["A"]})
    probe = ProbeData(data, ["SiO2"])

    result = convert_units(probe, "wt%", "wt%")

    assert result is not probe
    pd.testing.assert_frame_equal(result.data, probe.data)
    assert result.species == probe.species


@pytest.mark.parametrize(
    "from_unit, to_unit",
    [("ppb", "wt%"), ("ppm", "ppb"), ("wt", "ppm")],
)
def test_convert_units_rejects_unsupported_units(from_unit, to_unit):
    probe = ProbeData(pd.DataFrame({"SiO2": [50.0]}), ["SiO2"])

    with pytest.raises(ValueError, match="Unsupported units"):
        convert_units(probe, from_unit, to_unit)


def test_convert_units_does_not_modify_input():
    data = pd.DataFrame({"SiO2": [500000.0], "MgO": [100000.0], "Sample": ["A"]})
    probe = ProbeData(data, ["SiO2", "MgO"])
    original = probe.data.copy()

    convert_units(probe, "ppm", "wt%")

    pd.testing.assert_frame_equal(probe.data, original)


# Element -> oxide

def test_element_to_oxide():
    data = pd.DataFrame({"Si": [50.0], "Mg": [10.0]})
    probe = ProbeData(data, ["Si", "Mg"])

    result = element_to_oxide(probe)

    assert list(result.data.columns) == ["SiO2", "MgO"]
    assert result.species == ["SiO2", "MgO"]
    assert result.data.loc[0, "SiO2"] == pytest.approx(
        50.0 * ProbeData.MR["SiO2"] / ProbeData.MR_El["SiO2"]
    )
    assert result.data.loc[0, "MgO"] == pytest.approx(
        10.0 * ProbeData.MR["MgO"] / ProbeData.MR_El["MgO"]
    )


@pytest.mark.parametrize(
    "element, oxide",
    [("Fe2", "FeO"), ("Fe3", "Fe2O3")],
)
def test_element_to_oxide_iron(element, oxide):
    probe = ProbeData(pd.DataFrame({element: [10.0]}), [element])

    result = element_to_oxide(probe)

    assert result.species == [oxide]
    assert result.data.loc[0, oxide] == pytest.approx(
        10.0 * ProbeData.MR[oxide]
        / (ProbeData.cat_num[oxide] * ProbeData.MR_El[oxide])
    )


def test_total_iron_species_raises():
    with pytest.raises(ValueError, match="Species not found in species.csv"):
        ProbeData(pd.DataFrame({"Fe": [10.0]}), ["Fe"])


def test_element_to_oxide_preserves_metadata():
    data = pd.DataFrame({"Si": [50.0], "Sample": ["A"]})
    probe = ProbeData(data, ["Si"])

    result = element_to_oxide(probe)

    assert result.data.loc[0, "Sample"] == "A"


def test_element_to_oxide_does_not_modify_input():
    data = pd.DataFrame({"Si": [50.0], "Sample": ["A"]})
    probe = ProbeData(data, ["Si"])
    original = probe.data.copy()

    element_to_oxide(probe)

    pd.testing.assert_frame_equal(probe.data, original)


# Oxide -> element

def test_oxide_to_element():
    data = pd.DataFrame({"SiO2": [50.0], "MgO": [10.0]})
    probe = ProbeData(data, ["SiO2", "MgO"])

    result = oxide_to_element(probe)

    assert list(result.data.columns) == ["Si", "Mg"]
    assert result.species == ["Si", "Mg"]
    assert result.data.loc[0, "Si"] == pytest.approx(
        50.0 * ProbeData.MR_El["SiO2"] / ProbeData.MR["SiO2"]
    )
    assert result.data.loc[0, "Mg"] == pytest.approx(
        10.0 * ProbeData.MR_El["MgO"] / ProbeData.MR["MgO"]
    )


@pytest.mark.parametrize(
    "oxide, element",
    [("FeO", "Fe2"), ("Fe2O3", "Fe3")],
)
def test_oxide_to_element_iron(oxide, element):
    probe = ProbeData(pd.DataFrame({oxide: [10.0]}), [oxide])

    result = oxide_to_element(probe)

    assert result.species == [element]
    assert result.data.loc[0, element] == pytest.approx(
        10.0 * ProbeData.cat_num[oxide] * ProbeData.MR_El[oxide] / ProbeData.MR[oxide]
    )


def test_oxide_to_element_rejects_non_cation_species():
    probe = ProbeData(pd.DataFrame({"F": [1.0]}), ["F"])

    with pytest.raises(ValueError, match="Species cannot be converted"):
        oxide_to_element(probe)


def test_oxide_to_element_does_not_modify_input():
    data = pd.DataFrame({"SiO2": [50.0]})
    probe = ProbeData(data, ["SiO2"])
    original = probe.data.copy()

    oxide_to_element(probe)

    pd.testing.assert_frame_equal(probe.data, original)


# Round trips

@pytest.mark.parametrize(
    "oxide, element",
    [("SiO2", "Si"), ("Al2O3", "Al"), ("FeO", "Fe2"), ("Fe2O3", "Fe3")],
)
def test_oxide_element_round_trip(oxide, element):
    probe = ProbeData(pd.DataFrame({oxide: [10.0]}), [oxide])

    elemental = oxide_to_element(probe)
    recovered = element_to_oxide(elemental)

    assert recovered.species == [oxide]
    assert recovered.data.loc[0, oxide] == pytest.approx(10.0)