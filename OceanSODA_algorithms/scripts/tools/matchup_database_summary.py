#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 22:40:21 2020

@author: verwirrt
"""

#import pandas as pd;
import numpy as np;
from netCDF4 import Dataset;
import os;
os.chdir(os.path.join(os.getcwd(), ".."));

import osoda_global_settings;
import utilities;
from string import Template;
from os import path;


settings = osoda_global_settings.get_default_settings();

regionNC = Dataset(settings["regionMasksPath"], 'r');
regionNames = ["global", "oceansoda_amazon_plume", "oceansoda_st_lawrence", "oceansoda_mississippi", "oceansoda_congo", "oceansoda_mediterranean"];
#regionNames = ["oceansoda_amazon_plume"];

#settings["matchupDatasetTemplate"] = Template(path.join("../../matchup_datasets/v1_incomplete", "${YYYY}_soda_mdb.nc"));
matchupData = utilities.load_matchup_to_dataframe(settings, commonNames=["date", "AT", "DIC", "lat", "lon"])#, "SSS", "SST"]);

print("Number of data points by region:");
for regionName in regionNames:
    subset = utilities.subset_from_mask(matchupData, regionNC, regionName);
    print(regionName+":\t", len(subset));
    matrix = np.array(subset.values[:,1:], dtype=float); #all but date
    w = np.where(np.isfinite(matrix)==False);
    print("sst missing:", np.sum(w[1]==2), str(np.sum(w[1]==2)/len(subset)*100)+"%");
    print("sss missing:", np.sum(w[1]==3), str(np.sum(w[1]==3)/len(subset)*100)+"%");
    print("oc missing:", np.sum(w[1]==4), str(np.sum(w[1]==4)/len(subset)*100)+"%");
    print("dic missing:", np.sum(w[1]==5), str(np.sum(w[1]==5)/len(subset)*100)+"%");
    print("at missing:", np.sum(w[1]==6), str(np.sum(w[1]==6)/len(subset)*100)+"%");
    print("");

print("\n\nNumber of data points with DIC by region:");
for regionName in regionNames:
    subset = utilities.subset_from_mask(matchupData, regionNC, regionName);
    subset = subset[np.isfinite(subset["DIC"])]
    print(regionName+":\t", len(subset));

print("\n\nNumber of data points with AT by region:");
for regionName in regionNames:
    subset = utilities.subset_from_mask(matchupData, regionNC, regionName);
    subset = subset[np.isfinite(subset["AT"])]
    print(regionName+":\t", len(subset));



