# ProbeDataTools

ProbeDataTools is a Python package for processing mineral and liquid compositions acquired using an electron microprobe.

It provides an alternative to specialised spreadsheet-based workflows, where calculations are often hidden in cell formulae and need to be modified when different oxides are analysed. ProbeDataTools instead provides transparent, reusable calculations that are largely independent of mineral type and analytical setup.

The package includes data for 46 commonly analysed oxides in `oxides.csv`. Additional analytical species can be added without modifying the calculation code.

## Citation

ProbeDataTools is open-source software released under the MIT License. You are welcome to use, modify, and adapt the code.

If you use ProbeDataTools in research or published work, please cite the software using the citation information provided in CITATION.cff.

When using calculations that implement published methods, please also cite the original methodological source.

## Functionality

- Convert analytical wt% to atoms per formula unit (APFU).
- Calculate molar and cation proportions and fractions.
- Estimate Fe2+/Fe3+ using published stoichiometric methods.
- Calculate endmember components for olivine, feldspar, clinopyroxene and spinel.
- Identify potentially anomalous analyses from apfu totals.

## Basic Usage

ProbeDataTools uses pandas DataFrames for analytical data.

```python
import pandas as pd

from probedatatools import ProbeData
from probedatatools.cations import calc_apfu

data = pd.read_excel("olivine_data.xlsx")
species = ["SiO2", "FeO", "MgO", "MnO", "CaO"]

probe = ProbeData(data, species)
apfu = calc_apfu(probe, afu=4)
```

`calc_apfu()` is not specific to olivine: the same function can be used for any mineral given the appropriate number of oxygens per formula unit.

For example, an amphibole analysis can be calculated on a 23-oxygen basis:

```python
species = [
    "SiO2", "TiO2", "Al2O3", "FeO", "MgO",
    "CaO", "Na2O", "K2O", "F", "Cl",
]

probe = ProbeData(data, species)
apfu = calc_apfu(probe, afu=23)
```

Analysed non-oxygen anions such as F and Cl are identified from the species definitions and reported automatically in the resulting APFU dataframe. They do not contribute to the oxygen renormalisation factor.
