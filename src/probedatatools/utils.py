"""Small general-purpose utilities used by ProbeDataTools."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd


def aggregate_repeats(
    data: pd.DataFrame,
    group_by: Sequence[str],
    numeric_agg: str | Callable = "mean",
) -> pd.DataFrame:
    """Aggregate repeated probe analyses while retaining constant metadata.

    Numeric columns are aggregated with ``numeric_agg``. Non-numeric metadata
    columns are retained only when a group contains one unique non-null value;
    otherwise they are returned as ``pd.NA``.

    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing repeated analyses.

    group_by : Sequence[str]
        Columns defining which analyses belong to the same repeat group.

    numeric_agg : str or Callable, default="mean"
        Aggregation used for numeric columns, e.g. ``"mean"``, ``"median"``,
        or a callable.

    Returns
    -------
    pd.DataFrame
        One row per repeat group, with aggregated numerical analyses and
        constant metadata retained.
    """

    if not group_by:
        raise ValueError("`group_by` must contain at least one column")

    missing = [col for col in group_by if col not in data.columns]
    if missing:
        raise KeyError(f"Unknown grouping columns: {missing}")

    group_by = list(group_by)
    value_cols = [col for col in data.columns if col not in group_by]

    def keep_constant(series: pd.Series):
        values = series.dropna()
        return values.iloc[0] if not values.empty and values.nunique() == 1 else pd.NA

    aggregations: dict[str, str | Callable] = {}

    numeric_cols = data[value_cols].select_dtypes(include="number").columns
    metadata_cols = data[value_cols].select_dtypes(
        include=["object", "string", "category"]
    ).columns

    for col in numeric_cols:
        aggregations[col] = numeric_agg
    for col in metadata_cols:
        aggregations[col] = keep_constant

    return data.groupby(
        group_by, sort=False, dropna=False
    ).agg(aggregations).reset_index()


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
