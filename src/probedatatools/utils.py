"""Small general-purpose utilities used by ProbeDataTools."""

from __future__ import annotations

import numpy as np
import pandas as pd


def groupby_str(data: pd.DataFrame, cols: list[str], num: str) -> pd.DataFrame:
    """Group a dataframe while preserving constant string columns.

    Numeric columns are aggregated using ``num``. Object columns are
    retained when a group contains a single unique value; otherwise
    they are set to NaN.
    """

    aggregations = {}

    for column in data.select_dtypes(include="number"):
        if column not in cols:
            aggregations[column] = (column, num)

    for column in data.select_dtypes(include="object"):
        if column not in cols:
            aggregations[column] = (
                column,
                lambda s: s.iloc[0] if s.nunique(dropna=False) == 1 else np.nan,
            )

    return data.groupby(cols).agg(**aggregations)


def convert_headers_thermobar(
    data: pd.DataFrame,
    data_type: str,
    oxides: list[str],
    sample_id: str | None = None,
    fe_total: bool = True,
) -> pd.DataFrame:
    """Add Thermobar suffixes to selected oxide columns."""

    columns = list(oxides)
    if sample_id is not None and sample_id not in columns:
        columns.insert(0, sample_id)

    rename = {col: f"{col}_{data_type}" for col in columns}

    if fe_total and "FeO" in columns:
        rename["FeO"] = f"FeOt_{data_type}"

    if sample_id is not None:
        rename[sample_id] = f"Sample_ID_{data_type}"

    return data.loc[:, columns].rename(columns=rename).copy()


def revert_headers_thermobar(data: pd.DataFrame) -> pd.DataFrame:
    """Remove Thermobar suffixes from column headers."""

    output = data.copy()
    output.columns = [column.split("_", 1)[0] for column in output.columns]
    return output.rename(columns={"FeOt": "FeO"})


def gen_kd_line(Kd_mid: float, Kd_min: float, Kd_max: float) -> np.ndarray:
    """
    Produce Kd Lines and calculate no. of equilibrium points

    Parameters
    ----------
    Kd_mid : Partition constant
    Kd_min : minimum Partition constant
    Kd_max : maximum Partition constant

    Returns
    -------
    Kd_Total : contains arrays to plot Kd Lines

    """

    liq = np.linspace(0.01, 1.0, 100)
    mid = 1 / ((Kd_mid / liq) + (1 - Kd_mid))
    min_ = 1 / ((Kd_min / liq) + (1 - Kd_min))
    max_ = 1 / ((Kd_max / liq) + (1 - Kd_max))
    return np.vstack([liq, min_, mid, max_])
