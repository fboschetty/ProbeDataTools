"""
Core data model for electron-microprobe data.

`ProbeData` associates a probe-analysis DataFrame with metadata required
to process its analytical species.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from importlib.resources import files

import pandas as pd


def load_species_info() -> pd.DataFrame:
    """Load the bundled species metadata table."""

    path = files("probedatatools").joinpath("species.csv")
    return pd.read_csv(path, index_col=0)


@dataclass
class ProbeData:
    """Probe analyses and metadata for the selected species."""

    data: pd.DataFrame
    species: Sequence[str]

    species_info = load_species_info()

    MR = species_info.loc["MR"].astype(float)
    cat_num = species_info.loc["cations"].astype(float)
    ox_num = species_info.loc["oxygens"].astype(float)
    an_num = species_info.loc["anions"].astype(float)
    cat_str = species_info.loc["cat_str"].astype(str)
    MR_El = species_info.loc["MR_El"].astype(float)
    cat_chrg = (2.0 * ox_num / cat_num).where(cat_num > 0, 0.0)

    def __post_init__(self) -> None:
        """Validate species and select their metadata."""

        self.species = list(self.species)

        valid = set(self.species_info.columns) | set(self.cat_str)  # Allow headers to be elements as well as oxides
        unknown = sorted(set(self.species) - valid)

        if unknown:
            raise ValueError(f"Species not found in species.csv: {unknown}")

        self._refresh_metadata()

    def _metadata_source(self, species: str) -> str:
        """Return the metadata column associated with a species."""

        if species in self.species_info.columns:
            return species

        matches = self.cat_str[self.cat_str == species].index
        if len(matches) != 1:
            raise ValueError(f"Species has no unique metadata source: {species}")

        return matches[0]

    def _refresh_metadata(self) -> None:
        """Refresh metadata for the current species."""

        source = [self._metadata_source(species) for species in self.species]

        self.MR_use = self.MR.loc[source].set_axis(self.species)
        self.cat_num_use = self.cat_num.loc[source].set_axis(self.species)
        self.ox_num_use = self.ox_num.loc[source].set_axis(self.species)
        self.an_num_use = self.an_num.loc[source].set_axis(self.species)
        self.cat_str_use = self.cat_str.loc[source].set_axis(self.species)
        self.cat_chrg_use = self.cat_chrg.loc[source].set_axis(self.species)
        self.MR_El_use = self.MR_El.loc[source].set_axis(self.species)

        elemental = ~pd.Index(self.species).isin(self.species_info.columns)

        self.MR_use.loc[elemental] = self.MR_El_use.loc[elemental]
        self.cat_num_use.loc[elemental] = 1.0
        self.ox_num_use.loc[elemental] = 0.0
        self.an_num_use.loc[elemental] = 0.0