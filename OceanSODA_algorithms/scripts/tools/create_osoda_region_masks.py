#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 13:38:52 2020

@author: tom holding
"""

from netCDF4 import Dataset;
from os import path;
import numpy as np;

outputPath = "../../region_masks";

#create the mask matrices as follows (from Land et al 2019):
#global ocean
#the Caribbean (14°N to 30°N, 90°W to 60°W 
#the Amazon plume (2°S to 22°N, 70°W to 32°W)
#the Bay of Bengal (5°N to 24°N, 78°E to 96°E)

globalMask = np.ones((180, 360), dtype=int); #From Pathfinders

pathfinders_caribbeanMask = np.zeros((180, 360), dtype=int);
pathfinders_caribbeanMask[90+14:90+30+1, 180-90:180-60+1] = 1; #From Pathfinders

pathfinders_amazonPlumeMask = np.zeros((180, 360), dtype=int);
pathfinders_amazonPlumeMask[90-2:90+22+1, 180-70:180-32+1] = 1; #From Pathfinders

pathfinders_bayOfBengalMask = np.zeros((180, 360), dtype=int);
pathfinders_bayOfBengalMask[90+5:90+24+1, 180+78:180+96+1] = 1; #From Pathfinders


#######OceanSODA regions
osoda_amazon = pathfinders_amazonPlumeMask;

#Loosely based on Hopkings et al 2013: Hopkins, J., Lucas, M., Dufau, C., Sutton, M., Stum, J., Lauret, O. and Channelliere, C., 2013. Detection and variability of the Congo River plume from satellite derived sea surface temperature, salinity, ocean colour and sea level. Remote sensing of environment, 139, pp.365-385.
#(2W, 15E), (10S, 3N) from various figures (e.g. 14)
osoda_congo = np.zeros((180, 360), dtype=int);
osoda_congo[90-10:90+3+1, 180-2:180+15+1] = 1;

#Manually selection of rough region
#(70W, 56W), (46N, 52N)
osoda_st_lawrence = np.zeros((180, 360), dtype=int);
osoda_st_lawrence[90+46:90+52+1, 180-70:180-55+1] = 1;


#Loosely informed by Fournier et al 2016: Fournier, S., Lee, T. and Gierach, M.M., 2016. Seasonal and interannual variations of sea surface salinity associated with the Mississippi River plume observed by SMOS and Aquarius. Remote sensing of environment, 180, pp.431-439.
#(100W, 82W), (25N, 31N)
osoda_mississippi = np.zeros((180, 360), dtype=int);
osoda_mississippi[90+25:90+31+1, 180-100:180-83+1] = 1;

#Rectangular mask covering enclosed basin, excluding the black sea
#(5W, 36E), (30N, 46N), excluding (27E, 43E), (40N, 50N) (black sea), excluding (12W, 0E), (42N, 50N) (bay of biscay)
osoda_mediterranean = np.zeros((180, 360), dtype=int);
osoda_mediterranean[90+30:90+46+1, 180-5:180+36+1] = 1;
osoda_mediterranean[90+40:90+50+1, 180+27:180+43+1] = 0; #remove black sea
osoda_mediterranean[90+42:90+50+1, 180-12:180+0+1] = 0; #remove Bay of Biscay

#California current system
#28°N – 48°N, #138°W-114°W
osoda_californian_system = np.zeros((180, 360), dtype=int);
osoda_californian_system[90+28:90+48+1, 180-138:180-114+1] = 1;

#Barents sea
osoda_barents = np.zeros((180, 360), dtype=int);
osoda_barents[90+70:90+78+1, 180+21:180+56+1] = 1;

#Canary system
osoda_canary_system = np.zeros((180, 360), dtype=int);
osoda_canary_system[90+20:90+37+1, 180-30:180-7+1] = 1;

#Benguela upwelling system
osoda_benguela_system = np.zeros((180, 360), dtype=int);
osoda_benguela_system[90-25:90-11+1, 180-0:180+20+1] = 1;


osoda_all = np.zeros((180, 360), dtype=int);
osoda_all[np.where(osoda_amazon==1)] = 1;
osoda_all[np.where(osoda_congo==1)] = 2;
osoda_all[np.where(osoda_st_lawrence==1)] = 3;
osoda_all[np.where(osoda_mississippi==1)] = 4;
osoda_all[np.where(osoda_mediterranean==1)] = 5;
osoda_all[np.where(osoda_californian_system==1)] = 6;
osoda_all[np.where(osoda_barents==1)] = 7;
osoda_all[np.where(osoda_canary_system==1)] = 8;
osoda_all[np.where(osoda_benguela_system==1)] = 9;


#create NC file and dimensions
nc = Dataset(path.join(outputPath, "osoda_region_masks_v2.nc"), 'w');
nc.createDimension("lat", 180);
nc.createDimension("lon", 360);

#dimension variables
var = nc.createVariable("lat", float, ("lat",));
var.units = "lat (degrees North)";
var[:] = np.arange(-90, 90)+0.5;
var = nc.createVariable("lon", float, ("lon",));
var.units = "lon (degrees East)";
var[:] = np.arange(-180, 180)+0.5;

###### write masks to NC file
#Pathfinders masks
var = nc.createVariable("global", int, ("lat", "lon"));
var.units = "Integer";
var[:] = globalMask;

var = nc.createVariable("pathfinders_caribbean", int, ("lat", "lon"));
var.units = "Integer";
var[:] = pathfinders_caribbeanMask;

var = nc.createVariable("pathfinders_amazon_plume", int, ("lat", "lon"));
var.units = "Integer";
var[:] = pathfinders_amazonPlumeMask;

var = nc.createVariable("pathfinders_bay_of_bengal", int, ("lat", "lon"));
var.units = "Integer";
var[:] = pathfinders_bayOfBengalMask;


#OceanSODA masks
var = nc.createVariable("oceansoda_amazon_plume", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_amazon;

var = nc.createVariable("oceansoda_congo", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_congo;

var = nc.createVariable("oceansoda_st_lawrence", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_st_lawrence;

var = nc.createVariable("oceansoda_mississippi", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_mississippi;

var = nc.createVariable("oceansoda_mediterranean", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_mediterranean;


var = nc.createVariable("oceansoda_californian_system", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_californian_system;

var = nc.createVariable("oceansoda_barents", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_barents;

var = nc.createVariable("oceansoda_canary_system", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_canary_system;

var = nc.createVariable("oceansoda_benguela_system", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_benguela_system;

var = nc.createVariable("oceansoda_all", int, ("lat", "lon"));
var.units = "Integer";
var[:] = osoda_all;



nc.close();
