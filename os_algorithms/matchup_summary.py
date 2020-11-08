#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:07:22 2020

@author: verwirrt
"""

import osoda_global_settings;
from netCDF4 import Dataset;
import numpy as np;

settings = osoda_global_settings.get_default_settings();


startYear = 1990;
stopYear = 2019+1;


totalRows = 0;
totalDIC = 0;
totalAT = 0;
for year in range(startYear, stopYear):
    matchupPath = settings["matchupDatasetTemplate"].safe_substitute(YYYY=year);
    
    matchup = Dataset(matchupPath, 'r');
    totalRows += len(matchup.variables["lat"][:]);
    totalDIC += np.sum(np.isfinite(matchup.variables["DIC_mean"][:]).data);
    totalAT += np.sum(np.isfinite(matchup.variables["AT_mean"][:]).data);


print("Total rows:", totalRows);
print("Total DIC: ", totalDIC);
print("Total AT:  ", totalAT);