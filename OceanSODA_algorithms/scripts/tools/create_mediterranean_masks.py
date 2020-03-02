#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:56:01 2020

iho mediterranean basins

@author: tom holding
"""

import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;

#tmp
#ncin = Dataset("Peters_LeeEtAl2006_masks.nc", 'r');
#pregions = ncin.variables["regionNo"][:]; #Peter's region data
#regionNames = ncin.variables["regionName"][:];


outputPath = "../algorithms/algo_data/mediterranean_masks_tmh.nc";
eastBasinMaskPath = "../../misc/med_east_basin.nc";
westBasinMaskPath = "../../misc/med_west_basin.nc";


#create NC file and dimensions
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


#East mediterranean basin
eastnc = Dataset(eastBasinMaskPath, 'r');
eastInput = eastnc.variables["Band1"][:];
ilats, ilons = np.where(eastInput==1);
ilats += int(np.floor(eastnc.variables["lat"][:].min())) + 90;
ilons += int(np.floor(eastnc.variables["lon"][:].min())) + 180;

eastBasin = np.zeros((180, 360), dtype=float);
eastBasin[(ilats, ilons)] = 1;
var = nc.createVariable("east_basin_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = eastBasin;


#West mediterranean basin
westnc = Dataset(westBasinMaskPath, 'r');
westInput = westnc.variables["Band1"][:];
ilats, ilons = np.where(westInput==1);
ilats += int(np.floor(westnc.variables["lat"][:].min())) + 90;
ilons += int(np.floor(westnc.variables["lon"][:].min())) + 180;

westBasin = np.zeros((180, 360), dtype=float);
westBasin[(ilats, ilons)] = 1;
#westBasin = np.flipud(westBasin);
var = nc.createVariable("west_basin_mask", int, ("lat", "lon"));
var.units = "Integer";
var[:] = westBasin;
nc.close();


plt.figure();
plt.imshow(eastBasin+westBasin);
plt.colorbar();



