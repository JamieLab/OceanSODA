#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 12:56:01 2020


@author: tom holding
"""

#import numpy as np;
#from netCDF4 import Dataset;
#import matplotlib.pyplot as plt;
#
##tmp
##ncin = Dataset("Peters_LeeEtAl2006_masks.nc", 'r');
##pregions = ncin.variables["regionNo"][:]; #Peter's region data
##regionNames = ncin.variables["regionName"][:];
#
#
#outputPath = "../algorithms/algo_data/cai2010_masks_tmh.nc";
#oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
#oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);
##plt.figure(); plt.imshow(oceans); plt.colorbar();
#
##ocean mask values
#LAND=0;
#ATLANTIC=30;
#PACIFIC=70;
#INDIAN=50;
#SOUTHERN=90;
#
#
##create NC file and dimensions
#nc = Dataset(outputPath, 'w');
#nc.createDimension("lat", 180);
#nc.createDimension("lon", 360);
#
##dimension variables
#var = nc.createVariable("lat", float, ("lat",));
#var.units = "lat (degrees North)";
#var[:] = np.arange(-90, 90)+0.5;
#var = nc.createVariable("lon", float, ("lon",));
#var.units = "lon (degrees East)";
#var[:] = np.arange(-180, 180)+0.5;
#
#
##Zone 1 atlantic and indian: Equitorial and mid latitudes
#labrador = np.zeros((180, 360), dtype=float);
#labrador[90-30:90+30, :] = 1; #30S - 30N
#labrador[oceans!=ATLANTIC] = 0; #Exclude non-atlantic regions
#var = nc.createVariable("zone1_mask", int, ("lat", "lon"));
#var.units = "Integer";
#var[:] = zone1;
#
#
#
##sanity check for overlaps
##zoneSum = zone1+zone2+zone3+zone4+zone5+zone6a+zone6p;
##var = nc.createVariable("mask_sum_all_zones", int, ("lat", "lon"));
##var.units = "Integer";
##var[:] = zoneSum;
##
##nc.close();
##
##
##plt.figure();
##plt.imshow(zoneSum);
##plt.colorbar();
##plt.ylim(0,180)
#


