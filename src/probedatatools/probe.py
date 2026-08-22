"""
Core data model for electron-microprobe data.

`ProbeData` associates a probe-analysis DataFrame with the oxide
metadata required by the calculation modules.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from importlib.resources import files

import pandas as pd


def load_oxide_info() -> pd.DataFrame:
    """Load the bundled oxide metadata table."""

    path = files("probedatatools").joinpath("oxides.csv")
    return pd.read_csv(path, index_col=0)


@dataclass
class ProbeData:
    """Probe analyses and metadata for the selected species."""

    data: pd.DataFrame
    oxides: Sequence[str]

    oxide_info = load_oxide_info()

    MR = oxide_info.loc["MR"].astype(float)
    cat_num = oxide_info.loc["cations"].astype(float)
    ox_num = oxide_info.loc["oxygens"].astype(float)
    an_num = oxide_info.loc["anions"].astype(float)
    cat_str = oxide_info.loc["cat_str"].astype(str)
    MR_El = oxide_info.loc["MR_El"].astype(float)
    cat_chrg = (2.0 * ox_num / cat_num).where(cat_num > 0, 0.0)

    def __post_init__(self) -> None:
        """Validate species and select their metadata."""

        self.oxides = list(self.oxides)

        unknown = sorted(set(self.oxides) - set(self.oxide_info.columns))
        if unknown:
            raise ValueError(
                f"Species not found in oxides.csv: {unknown}"
            )

        self._refresh_metadata()

    def _refresh_metadata(self) -> None:
        """Refresh metadata for the current species list."""

        self.MR_use = self.MR.loc[self.oxides]
        self.cat_num_use = self.cat_num.loc[self.oxides]
        self.ox_num_use = self.ox_num.loc[self.oxides]
        self.an_num_use = self.an_num.loc[self.oxides]
        self.cat_str_use = self.cat_str.loc[self.oxides]
        self.cat_chrg_use = self.cat_chrg.loc[self.oxides]
        self.MR_El_use = self.MR_El.loc[self.oxides]