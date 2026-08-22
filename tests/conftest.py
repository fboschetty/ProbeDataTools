from pathlib import Path

import pandas as pd
import pytest

from probedatatools import ProbeData


DATA_DIR = Path(__file__).parent / "test_data"


@pytest.fixture
def olivine_probe():
    mgO = ProbeData.MR["MgO"]
    siO2 = ProbeData.MR["SiO2"]
    formula_mass = 2 * mgO + siO2

    data = pd.DataFrame({
        "SiO2": [100 * siO2 / formula_mass],
        "MgO": [100 * (2 * mgO) / formula_mass],
    })

    return ProbeData(data, ["SiO2", "MgO"])


@pytest.fixture
def halogen_probe():
    mgO = ProbeData.MR["MgO"]
    siO2 = ProbeData.MR["SiO2"]
    formula_mass = 2 * mgO + siO2

    return ProbeData(
        pd.DataFrame({
            "SiO2": [98 * siO2 / formula_mass],
            "MgO": [98 * (2 * mgO) / formula_mass],
            "F": [1.0],
            "Cl": [1.0],
        }),
        ["SiO2", "MgO", "F", "Cl"],
    )


@pytest.fixture
def fe_probe():
    return ProbeData(
        pd.DataFrame({"FeO": [100.0]}),
        ["FeO"],
    )