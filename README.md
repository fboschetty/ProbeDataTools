# ProbeDataTools

ProbeDataTools is a Python package for processing mineral and liquid compositions acquired using an electron microprobe.

It provides an alternative to specialised spreadsheet-based workflows, where calculations are often hidden in cell formulae and need to be modified when different oxides are analysed. ProbeDataTools instead provides transparent, reusable calculations that are largely independent of mineral type and analytical setup.

The package includes data for 46 commonly analysed oxides in `oxides.csv`. Additional analytical species can be added without modifying the calculation code.

## Functionality

- Convert oxide wt% to cations per formula unit.
- Calculate molar and cation proportions and fractions.
- Estimate Fe2+/Fe3+ using published stoichiometric methods.
- Calculate endmember components for olivine, feldspar, clinopyroxene and spinel.
- Identify potentially anomalous analyses from cation totals.
- Convert data to and from Thermobar column conventions.

## Basic Usage

ProbeDataTools uses pandas DataFrames for analytical data.

```python
import pandas as pd

from probedatatools import ProbeData
from probedatatools.cations import calc_cations

data = pd.read_excel("olivine_data.xlsx")
oxides = ["SiO2", "FeO", "MgO", "MnO", "CaO"]

probe = ProbeData(data, oxides)
cations = calc_cations(probe, afu=4)
```

`calc_cations()` is not specific to olivine: the same function can be used for any mineral given the appropriate anion basis.