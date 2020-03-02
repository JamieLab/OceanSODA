#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 09:29:06 2020

Uses the NASA distance-from-coast dataset (https://oceancolor.gsfc.nasa.gov/docs/distfromcoast/)
to create a mask and saves the mask to a netCDF file

@author: Tom Holding
"""

import numpy as np;
import pandas as pd;
from netCDF4 import Dataset;
from os import path;
import matplotlib.pyplot as plt;
from scipy.stats import binned_statistic_2d;

#Calculate the minimum of subgrids. Subgrids defined by reshape to outputShape
def min_of_submatrices(matrix, outputShape):
    sh = outputShape[0], matrix.shape[0]//outputShape[0], outputShape[1], matrix.shape[1]//outputShape[1];
    return matrix.reshape(sh).min(3).min(1);


generateFromScratch = False;


distanceToCoastPath = path.expanduser("~/data/NASA_distance_to_coast/dist2coast.txt"); #0.04 degree distance to land grid
inputRes = 0.04; #degrees
outputResolution = 1.0; #degrees
#targetRes = (180, 360); #output resolution
minDistance = 250; #km
outputPath = path.join("../../region_masks", "distance_to_land_mask_"+str(int(minDistance))+"km.nc");
binnedMinPath = path.join(path.dirname(outputPath),"rebinned_dist_to_land_matrix_"+str(outputResolution)+"x"+str(outputResolution)+"deg.csv")



if generateFromScratch == True:
    #Read input
    inputDF = pd.read_csv(distanceToCoastPath, header=None, names=["lon", "lat", "distance"], sep="\t");
    subset = inputDF.iloc[0:10000];
    
    binnedMinDistances = np.zeros((int(180), int(360)), dtype=float);
    for y in range(-90,90):
        for x in range(-180, 180):
            print(x, ",", y);
            subset = inputDF[(inputDF["lon"] > x*1.0) & (inputDF["lon"] <= (x+1)*1.0) &
                             (inputDF["lat"] > y*1.0) & (inputDF["lat"] <= (y+1)*1.0)];
            binnedMinDistances[y+90,x+180] = subset["distance"].min();

    np.savetxt(binnedMinPath, binnedMinDistances, delimiter=",");
    
if generateFromScratch == False: #Use the bined min distance values
    binnedMinDistances = np.genfromtxt(binnedMinPath, delimiter=",");

#Apply threshold to create mask
mask = np.zeros(binnedMinDistances.shape, dtype=int);
mask[binnedMinDistances >= minDistance] = 1;
#plt.figure(); plt.imshow(mask);


#Write output to netCDF
nc = Dataset(outputPath, 'w');
nc.createDimension("lat", 180);
nc.createDimension("lon", 360);

#dimension variables
var = nc.createVariable("lat", float, ("lat",));
var.units = "lat (degrees North)";
var[:] = np.arange(-90, 90)+0.5;
var = nc.createVariable("lon", float, ("lon",));
var.units = "lon (degrees East)";
var[:] = np.arange(-180, 180)+0.5;


#full distance to coast
var = nc.createVariable("dist_to_coast", int, ("lat", "lon"));
var.long_name = "Minimum distance to the nearest coast for a 1x1 deg grid.";
var.units = "km";
var[:] = binnedMinDistances;

#mask
var = nc.createVariable("dist_to_coast_mask", int, ("lat", "lon"));
var.long_name = "Distance-to-coast mask indicating regions where the distance to the nearest coast is >= "+str(minDistance)+" km";
var.units = "integer mask";
var[:] = mask;

nc.close();
