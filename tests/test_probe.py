import pandas as pd
import pytest

from probedatatools import ProbeData


@pytest.fixture
def probe():
    """Create a simple ProbeData object."""

    return ProbeData(
        pd.DataFrame({"SiO2": [50.0], "MgO": [30.0], "FeO": [20.0]}),
        ["SiO2", "MgO", "FeO"],
    )


def test_probe_data_stores_input(probe):
    assert list(probe.data.columns) == ["SiO2", "MgO", "FeO"]
    assert probe.species == ["SiO2", "MgO", "FeO"]


@pytest.mark.parametrize("attribute", ["MR_use", "cat_num_use", "ox_num_use", "an_num_use", "cat_str_use", "cat_chrg_use"])
def test_metadata_filtered_to_species(probe, attribute):
    assert getattr(probe, attribute).index.tolist() == probe.species


@pytest.mark.parametrize("oxide, expected", [("SiO2", 4.0), ("MgO", 2.0), ("Al2O3", 3.0), ("FeO", 2.0), ("Fe2O3", 3.0), ("TiO2", 4.0)])
def test_cation_charge(oxide, expected):
    assert ProbeData.cat_chrg[oxide] == pytest.approx(expected)


def test_filtered_metadata_values(probe):
    assert probe.cat_num_use["SiO2"] == pytest.approx(1.0)
    assert probe.ox_num_use["SiO2"] == pytest.approx(2.0)
    assert probe.cat_str_use["SiO2"] == "Si"
    assert probe.cat_num_use["MgO"] == pytest.approx(1.0)
    assert probe.ox_num_use["MgO"] == pytest.approx(1.0)
    assert probe.cat_str_use["MgO"] == "Mg"


def test_non_cation_species():
    probe = ProbeData(pd.DataFrame({"F": [1.0], "Cl": [0.5]}), ["F", "Cl"])
    assert (probe.cat_num_use == 0).all()
    assert (probe.cat_chrg_use == 0).all()


def test_halogen_anion_metadata():
    probe = ProbeData(pd.DataFrame({"F": [1.0], "Cl": [0.5]}), ["F", "Cl"])
    assert probe.an_num_use["F"] == pytest.approx(1.0)
    assert probe.an_num_use["Cl"] == pytest.approx(1.0)
    assert probe.ox_num_use["F"] == pytest.approx(0.0)
    assert probe.ox_num_use["Cl"] == pytest.approx(0.0)


def test_metadata_preserves_oxide_order():
    probe = ProbeData(
        pd.DataFrame({"MgO": [30.0], "SiO2": [50.0], "FeO": [20.0]}),
        ["FeO", "SiO2", "MgO"],
    )
    assert probe.MR_use.index.tolist() == ["FeO", "SiO2", "MgO"]


def test_unknown_species_raises_error():
    with pytest.raises(ValueError, match="Species not found in species.csv"):
        ProbeData(pd.DataFrame({"SiO2": [50.0], "UnobtainiumO2": [1.0]}), ["SiO2", "UnobtainiumO2"])

@pytest.mark.parametrize(
    "species, source",
    [("Si", "SiO2"), ("Mg", "MgO"), ("Fe2", "FeO"), ("Fe3", "Fe2O3")],
)
def test_elemental_species_use_linked_metadata(species, source):
    probe = ProbeData(pd.DataFrame({species: [1.0]}), [species])

    assert probe.MR_use[species] == pytest.approx(probe.MR_El[source])
    assert probe.cat_num_use[species] == pytest.approx(1.0)
    assert probe.ox_num_use[species] == pytest.approx(0.0)
    assert probe.cat_chrg_use[species] == pytest.approx(probe.cat_chrg[source])