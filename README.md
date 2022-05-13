# ProbeDataTools
Python code for dealing with probe data.

Contains: 

**(1) ProbeData.py** which contains the probe data class. This is used to relate mineral data to useful information about their consitutent oxides.
This pulls information from *oxides.csv*. This can be edited to add more obscure oxides.

**(2) CalcCationsMin.py** which contains a group of functions used to calculate the cations per formula unit of any mineral, given an ideal number of cations and anions per formula unit.
It also contains functions for filtering out non-ideal analyses.

**(3) CalcFe2O3Min.py** which contains relevant functions to calculate the Fe2/3 ratio of a mineral via two methods. (1) the stoichiometric method of Droop (1987). (2) The method of Papike et al., (1947).

**(4) CalcEmMin.py** (In Progress) contains relevant functions to calculate the endmembers of common minerals e.g. olivine, feldspar, clinopyroxene, amphibole and spinels.

**(5) Test.py** This pulls mineral data from *Vill_DB_SM_1.xlsx* and calculates the cations per formula unit of olivine, plagioclase and clinopyroxene data.
