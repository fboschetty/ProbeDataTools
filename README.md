# ProbeDataTools
This package contains a series of python scripts which process mineral or liquid data acquired using an Electron MicroProbe. It presents an alternative to traditional complex spreadsheets which are prone to user error as they often need to be adapted to cater for the oxides measured.  

This package contains a comma-seperated-variable file (*oxides.csv*) which has information about 46 commonly measured oxides. This can easily be added to, to cater for other oxides. No changes to the code itself are required, minimising user error.  

This package has functionality for the following:  

- Convert mineral data from wt% oxide to cations per formula unit for any mineral given an ideal number of cations and anions per formula unit.
- Estimate Fe2O3 contents from FeOT using different stoichiometric methods.
- Calculate endmember proportions for common minerals. 
- Identify potentially anomalous analyses based on cation per formula unit totals.
- Convert liquid data from wt% oxide to molar or cation fractions.  

Detailed examples of the above can be found in the *test.py* script using a compilation of analyses of volcanic products from Villarrica volcano, Chile.

# Basic Usage

This package utilises the commonly used data manipulation package pandas (<https://pandas.pydata.org/>). Therefore you should load your data into python using one of their *reader* functions e.g., pandas.read\_csv() or pandas.read\_excel(). The data must be formated so that each row represents an analysis, and columns have headers for each oxide e.g. SiO2 or Al2O3.

```python
# read in raw data from an excel spreadsheet
Ol_raw = pd.read_excel("olivine_data.xlsx")
```  

The package contains a data class *probedata* that contains  the raw data you want to process and all the relevant data for processing it. This requires a pandas DataFrame containing the raw data and a list of the oxides you want to use to process it.  

```python
# create a probedata object for the data and specify oxides to analyse
oxides = ['SiO2', 'FeO', 'MgO', 'MnO', 'CaO']
Ol_probe = probedata(probe_data=Ol_raw, oxides=oxides)
```  

Now we can process the data. For example we could calculate the cations per formula unit of each data point, given an ideal anions-per-formula-unit of four:

```python
Ol_cations = calc_cations(probe_data=Ol_probe, afu=4.)
```  

The above line returns a pandas DataFrame containing the cations per formula unit from the raw data in wt% oxides. This could be saved as another excel spreadsheet or used in further calculations.

The code is written in such a way that it doesn't know what mineral is being input, so we could calculate the cations per formula unit for any given number of anions per formula unit. This makes it applicable to analyses of any mineral, given we have a constraint on the structure of the mineral's  unit lattice.

# Suggestions

If you have ideas for functionality that you would like to see, or bugs and errors you have found, please add them by opening an issue. Or me an email at eefob@leeds.ac.uk
