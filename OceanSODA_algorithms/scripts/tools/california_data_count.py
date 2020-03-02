#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 15:03:05 2020

@author: tom holding
"""

import os;
from os import path;
if path.basename(os.getcwd()) != "20190816_OceanSODA":
    os.chdir(path.join(os.getcwd(), "../.."));


import pandas as pd;
from string import Template;
import numpy as np;
from netCDF4 import Dataset;
import geopy.distance;


def make_corners(north, south, east, west):
    return [[north, east],
            [north, west],
            [south, east],
            [south, west]];


regionTemplate = Template("matchup_datasets/matchup_dataset_v2/${YYYY}regions");

#years = np.arange(1957, 2019+1);
years = np.arange(2010, 2019+1);

#combine all dataframes into one
allRegionDFs = [];
for year in years:
    try:
        df = pd.read_csv(regionTemplate.safe_substitute(YYYY=year));
        allRegionDFs.append(df);
    except FileNotFoundError:
        print(year, "not found");
regionData = pd.concat(allRegionDFs, ignore_index=True);


if True:
    bounds = pd.DataFrame();
    bounds["region"] = regionData["REGION"];
    norths = [];
    souths = [];
    easts = [];
    wests = [];
    for r in range(len(regionData)):
        print(r, "of", len(regionData));
        row = regionData.loc[r];
        centre = geopy.Point(row["LAT"], row["LON"]);
        radius = row["RADIUS(KM)"];
        distObj = geopy.distance.geodesic(kilometers=100);
        
        norths.append(distObj.destination(point=centre, bearing=0.0).latitude);
        easts.append(distObj.destination(point=centre, bearing=90.0).longitude);
        souths.append(distObj.destination(point=centre, bearing=180.0).latitude);
        wests.append(distObj.destination(point=centre, bearing=270.0).longitude);
    
    bounds["north"] = norths;
    bounds["south"] = souths;
    bounds["east"] = easts;
    bounds["west"] = wests;


masksToUse = ["oceansoda_amazon_plume", "oceansoda_congo", "oceansoda_st_lawrence", "oceansoda_mississippi", "oceansoda_mediterranean", "oceansoda_californian_system"];
masknc = Dataset("region_masks/region_masks.nc", 'r');

maskCounts = {};
for maskName in masksToUse:
    maskCounts[maskName] = 0;
    mask = masknc.variables[maskName][:];
#    lons = np.floor(masknc.variables["lon"][:]);
#    lats = np.floor(masknc.variables["lat"][:]);
    
    for r in range(len(bounds)):
        print(maskName, r, "of", len(bounds), "n =", maskCounts[maskName]);
        row = bounds.loc[r];
        corners = make_corners(row["north"], row["south"], row["east"], row["west"]);
        corners = np.floor(corners);
        corners[:,0] += 90;
        corners[:,1] += 180;
        corners = np.array(corners, dtype=int);
        
        keep = True;
        for i in range(0, 4):
            if mask[corners[i,0], corners[i,1]] == 0:
                keep=False;
                break;
        if keep:
            maskCounts[maskName] += 1;
    

for maskName in maskCounts:
    print(maskName+":\t", maskCounts[maskName]);
    
    
    
    
    
