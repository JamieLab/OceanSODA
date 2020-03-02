#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:41:50 2020

Uses reef base 'reef location' data to produce summaries of number of reefs in each OSODA region.

@author: tom holding
"""

from os import path;
import sys;
osodaBaseDir = path.abspath(path.dirname(path.dirname(__file__)));
if osodaBaseDir not in sys.path:
    sys.path.append(osodaBaseDir); #tmp hack. TODO: reorganise folders

import pandas as pd;
import numpy as np;
from netCDF4 import Dataset, stringtochar;
from scripts.utilities import subset_from_mask;


regionMaskPath = path.join("../region_masks/osoda_region_masks_v2.nc");
inputPath = "ReefLocations.csv";
outputPath = "reef_locations.nc";

df = pd.read_csv(inputPath, sep=",", encoding="ISO-8859-1");
df.columns = [key.lower() for key in df.keys()];


regionMaskNC = Dataset(regionMaskPath, 'r');
#dfAll = subset_from_mask(df, regionMaskNC, "oceansoda_all");
dfAmazon = subset_from_mask(df, regionMaskNC, "oceansoda_amazon_plume");
dfCongo = subset_from_mask(df, regionMaskNC, "oceansoda_congo");
dfMississippi = subset_from_mask(df, regionMaskNC, "oceansoda_mississippi");
dfStLawrence = subset_from_mask(df, regionMaskNC, "oceansoda_st_lawrence");
dfMediterranean = subset_from_mask(df, regionMaskNC, "oceansoda_mediterranean");

dfBarents = subset_from_mask(df, regionMaskNC, "oceansoda_barents");
dfBenguelaSystem = subset_from_mask(df, regionMaskNC, "oceansoda_benguela_system");
dfCalifornianSystem = subset_from_mask(df, regionMaskNC, "oceansoda_californian_system");
dfCanarySystem = subset_from_mask(df, regionMaskNC, "oceansoda_canary_system");


### Print summaries
print("Number of reefs (ReefBase) by OceanSODA region:");
print("\tGlobal:", len(df));
print("\tAmazon:", len(dfAmazon));
print("\tCongo:", len(dfCongo));
print("\tBarents:", len(dfBarents));
print("\tBenguela System:", len(dfBenguelaSystem));
print("\tCalifornian System:", len(dfCalifornianSystem));
print("\tCanary System:", len(dfCanarySystem));
print("\tCMeditteranean:", len(dfMediterranean));
print("\tMississippi:", len(dfMississippi));
print("\tSt Lawrence:", len(dfStLawrence));


