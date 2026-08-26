# ProbeDataTools

**ProbeDataTools** is a Python package for processing mineral and liquid compositions acquired using an electron microprobe.

The package provides reusable Python tools for common mineral-chemical calculations without relying on mineral-specific spreadsheet templates.

## Features

ProbeDataTools can:

* calculate atoms per formula unit (APFU) on a specified oxygen basis;
* calculate molar and cation proportions and fractions;
* handle analysed non-oxygen anions such as F and Cl;
* estimate Fe²⁺/Fe³⁺ using published stoichiometric methods;
* calculate endmember components for olivine, feldspar, clinopyroxene, and spinel;
* identify potentially anomalous analyses from cation totals; and
* convert data to and from Thermobar column conventions.

## Installation

ProbeDataTools can be installed using pip:

```bash
pip install probedatatools
```

## Updating

To update an existing installation to the latest release:

```bash
pip install --upgrade probedatatools
```

The installed version can be checked with:

```python
import probedatatools

probedatatools.__version__
```

## Quick Example

ProbeDataTools operates on analytical data stored in pandas DataFrames.

```python
import pandas as pd

from probedatatools import ProbeData
from probedatatools.cations import calc_apfu

data = pd.read_excel("olivine_data.xlsx")

species = ["SiO2", "FeO", "MgO", "MnO", "CaO"]

probe = ProbeData(data, species)

apfu = calc_apfu(probe, afu=4)
```

The same calculation framework can be applied to different minerals by specifying the appropriate number of oxygens per formula unit.

Analysed non-oxygen anions such as F and Cl are recognised automatically from the analytical species and carried through to the APFU calculation.

## Documentation

See [Getting Started](getting_started.ipynb) for an introduction to the basic ProbeDataTools workflow.

Further documentation covers:

* atoms per formula unit and compositional calculations;
* Fe²⁺/Fe³⁺ partitioning;
* mineral endmember calculations; and
* conversion to and from Thermobar column conventions.

## Citation

If you use ProbeDataTools in published work, please cite the software.

**Software citation:** *[citation and DOI to be added]*

Calculations that implement published methods should also be accompanied by a citation to the original source. Relevant references are provided in the documentation for each method.
