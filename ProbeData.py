#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 14:38:13 2022

@author: FelixBoschetty
"""

import pandas as pd
from dataclasses import dataclass


@dataclass
class probedata:
    """
        Class that relates input probe data, oxides analysed and useful data about all oxides.

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

    # Do I want this to be unchangeable?
    oxide_info = pd.read_csv(r'oxides.csv', index_col=0)

    # Extract useful information from dataframe
    MR = oxide_info.loc['MR'].astype('float')
    cat_num = oxide_info.loc['cations'].astype('float')
    ox_num = oxide_info.loc['oxygens'].astype('float')
    cat_str = oxide_info.loc['cat_str'].astype('str')
    cat_chrg = 2.*ox_num/cat_num
