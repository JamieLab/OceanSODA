#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 08:39:06 2020

Uses the gebco elevation dataset to create a minimum-depth mask netCDF file

@author: Tom Holding
"""

import numpy as np;
from netCDF4 import Dataset;
from os import path;
import matplotlib.pyplot as plt;

#Calculate the minimum of subgrids. Subgrids defined by reshape to outputShape
def min_of_submatrices(matrix, outputShape):
    sh = outputShape[0], matrix.shape[0]//outputShape[0], outputShape[1], matrix.shape[1]//outputShape[1];
    return matrix.reshape(sh).min(3).min(1);


gebcoElevationPath = path.expanduser("~/data/GEBCO_bathymetry_30sec/GEBCO_2014_2D.nc"); #30sec resolution bathymetry
inputRes = int(60/2/1.0); #30 sec resolutioon
targetRes = (180, 360); #output resolution
maxElevation = -500.0; #meters
outputPath = path.join("../../region_masks", "depth_mask_"+str(int(-maxElevation))+"m.nc");

#Read input bathymetry
inputNC = Dataset(gebcoElevationPath, 'r');
inputElevation = inputNC.variables["elevation"][:];

if False: #Debug plotting
    plt.figure();
    plt.imshow(inputElevation);
    plt.colorbar();

#Create grid with minimum depth
minSubgrids = min_of_submatrices(inputElevation, targetRes);
if False: #Debug plotting
    plt.figure();
    plt.imshow(minSubgrids);
    plt.colorbar();

#Apply threshold
depthMask = (minSubgrids < maxElevation).astype(int);
if False: #debug plotting
    plt.figure();
    plt.imshow(depthMask);


#Write netCDF file
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


#Zone 1, equation 1
var = nc.createVariable("depth_mask", int, ("lat", "lon"));
var.long_name = "Depth mask indicating regions where the minimum depth exceeds "+str(int(-maxElevation))+" meters";
var.units = "integer mask";
var[:] = depthMask;

nc.close();
