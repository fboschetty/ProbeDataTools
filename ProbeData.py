#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains the probedata class.
This is used to relate mineral data to useful information about their consitutent oxides.
This pulls information from oxides.csv. This can be edited to add more obscure oxides.

@author: FelixBoschetty
"""

import os
import pandas as pd
from dataclasses import dataclass

@dataclass
class probedata:
    """Class that relates input probe data, oxides analysed and useful data about oxides.

    Parameters
    ----------
    data: pd.DataFrame
        The probe data that you want to process. Ensure headers are wt% oxides.

    oxides: List[str]
        A list of oxides that you want to use for analysis.

    Returns
    -------
    None.

    """

    data: pd.DataFrame
    oxides: list[str]
    here = os.path.dirname(os.path.abspath(__file__))

    filename = os.path.join(here, 'oxides.csv')
    oxide_info = pd.read_csv(filename, index_col=0)

    # Extract useful information from oxides_info
    MR = oxide_info.loc['MR'].astype('float')
    cat_num = oxide_info.loc['cations'].astype('float')
    ox_num = oxide_info.loc['oxygens'].astype('float')
    cat_str = oxide_info.loc['cat_str'].astype('str')
    cat_chrg = 2.*ox_num/cat_num
