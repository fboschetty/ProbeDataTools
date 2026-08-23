import numpy as np
import pandas as pd
import pytest

from probedatatools.utils import (
    aggregate_repeats,
    convert_headers_thermobar,
    gen_kd_line,
    revert_headers_thermobar,
)


@pytest.mark.parametrize("numeric_agg", ["mean", "median"])
def test_aggregate_repeats_numeric_aggregation(numeric_agg):
    data = pd.DataFrame({"Sample": ["A", "A", "A"], "Value": [1.0, 3.0, 5.0]})
    result = aggregate_repeats(data, ["Sample"], numeric_agg)
    assert result.loc[0, "Value"] == pytest.approx(3.0)


def test_aggregate_repeats_preserves_constant_metadata():
    data = pd.DataFrame({
        "Sample": ["A", "A", "B"],
        "Value": [1.0, 3.0, 10.0],
        "Comment": ["Good", "Good", "Other"],
    })

    result = aggregate_repeats(data, ["Sample"])

    assert result.loc[result["Sample"] == "A", "Comment"].iloc[0] == "Good"
    assert result.loc[result["Sample"] == "B", "Comment"].iloc[0] == "Other"


def test_aggregate_repeats_marks_conflicting_metadata():
    data = pd.DataFrame({
        "Sample": ["A", "A"],
        "Value": [1.0, 3.0],
        "Comment": ["Good", "Bad"],
    })

    result = aggregate_repeats(data, ["Sample"])

    assert pd.isna(result.loc[0, "Comment"])


def test_aggregate_repeats_multiple_grouping_columns():
    data = pd.DataFrame({
        "Sample": ["A", "A", "A", "B"],
        "Type": ["x", "x", "y", "x"],
        "Value": [1.0, 3.0, 5.0, 7.0],
    })

    result = aggregate_repeats(data, ["Sample", "Type"])

    ax = result.set_index(["Sample", "Type"])
    assert ax.loc[("A", "x"), "Value"] == pytest.approx(2.0)
    assert ax.loc[("A", "y"), "Value"] == pytest.approx(5.0)


def test_aggregate_repeats_ignores_missing_metadata():
    data = pd.DataFrame({
        "Sample": ["A", "A"],
        "Position": ["core", pd.NA],
        "Value": [1.0, 3.0],
    })

    result = aggregate_repeats(data, ["Sample"])

    assert result.loc[0, "Position"] == "core"


def test_aggregate_repeats_validates_grouping_columns():
    data = pd.DataFrame({"Sample": ["A"], "Value": [1.0]})

    with pytest.raises(KeyError):
        aggregate_repeats(data, ["Missing"])


def test_convert_headers_thermobar():
    data = pd.DataFrame({
        "Sample": ["A"],
        "SiO2": [50.0],
        "FeO": [10.0],
        "MgO": [5.0],
        "Extra": [99.0],
    })

    result = convert_headers_thermobar(
        data, "Liq", ["SiO2", "FeO", "MgO"], sample_id="Sample"
    )

    assert list(result.columns) == ["Sample_ID_Liq", "SiO2_Liq", "FeOt_Liq", "MgO_Liq"]
    assert "Extra" not in result


def test_convert_headers_thermobar_without_sample_id():
    data = pd.DataFrame({"SiO2": [50.0], "MgO": [5.0]})
    result = convert_headers_thermobar(data, "Liq", ["SiO2", "MgO"])

    assert list(result.columns) == ["SiO2_Liq", "MgO_Liq"]


def test_convert_headers_thermobar_fe_total_false():
    data = pd.DataFrame({"FeO": [10.0]})
    result = convert_headers_thermobar(data, "Liq", ["FeO"], fe_total=False)

    assert list(result.columns) == ["FeO_Liq"]


def test_convert_headers_thermobar_does_not_modify_input():
    data = pd.DataFrame({"Sample": ["A"], "SiO2": [50.0]})
    original = data.copy()

    convert_headers_thermobar(data, "Liq", ["SiO2"], sample_id="Sample")

    pd.testing.assert_frame_equal(data, original)


def test_revert_headers_thermobar():
    data = pd.DataFrame({
        "Sample_ID_Liq": ["A"],
        "SiO2_Liq": [50.0],
        "FeOt_Liq": [10.0],
        "MgO_Liq": [5.0],
    })

    result = revert_headers_thermobar(data)

    assert list(result.columns) == ["Sample", "SiO2", "FeO", "MgO"]


def test_revert_headers_thermobar_does_not_modify_input():
    data = pd.DataFrame({"SiO2_Liq": [50.0], "FeOt_Liq": [10.0]})
    original = data.copy()

    revert_headers_thermobar(data)

    pd.testing.assert_frame_equal(data, original)


def test_convert_revert_thermobar_round_trip():
    data = pd.DataFrame({
        "Sample": ["A"],
        "SiO2": [50.0],
        "FeO": [10.0],
        "MgO": [5.0],
    })

    converted = convert_headers_thermobar(data, "Liq", ["SiO2", "FeO", "MgO"], sample_id="Sample")
    reverted = revert_headers_thermobar(converted)

    pd.testing.assert_frame_equal(reverted, data)


@pytest.fixture
def kd_lines():
    return gen_kd_line(0.30, 0.20, 0.40)


def test_gen_kd_line_shape(kd_lines):
    assert kd_lines.shape == (4, 100)


def test_gen_kd_line_liquid_axis(kd_lines):
    assert kd_lines[0, 0] == pytest.approx(0.01)
    assert kd_lines[0, -1] == pytest.approx(1.0)
    assert np.all(np.diff(kd_lines[0]) > 0)


def test_gen_kd_line_limits(kd_lines):
    assert np.allclose(kd_lines[1:, -1], 1.0)
    assert np.all(kd_lines[1:, 0] < 0.1)


def test_gen_kd_line_ordering(kd_lines):
    assert np.all(kd_lines[3] <= kd_lines[2])
    assert np.all(kd_lines[2] <= kd_lines[1])


def test_gen_kd_line_matches_formula(kd_lines):
    liq = kd_lines[0]
    expected = [
        1 / ((kd / liq) + (1 - kd))
        for kd in (0.20, 0.30, 0.40)
    ]

    assert np.allclose(kd_lines[1], expected[0])
    assert np.allclose(kd_lines[2], expected[1])
    assert np.allclose(kd_lines[3], expected[2])