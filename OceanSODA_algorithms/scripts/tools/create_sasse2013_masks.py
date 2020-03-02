#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 09:56:28 2020

Zone definitions for AT and DIC algorithms from Sasse et al 2013.
See: Sasse, T.P., McNeil, B.I. and Abramowitz, G., 2013. A novel method for diagnosing seasonal to inter-annual surface ocean carbon dynamics from bottle data using neural networks. Biogeosciences, 10(6), pp.4319-4340.

@author: tom holding
"""

import numpy as np;
from netCDF4 import Dataset;
import matplotlib.pyplot as plt;


outputPath = "../algorithms/algo_data/sasse2013_masks_tmh.nc";

oceanMaskPath = "../../misc/World_Seas-IHO-mask.nc";
oceans = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);

#ocean mask values
LAND=0;
ATLANTIC=30;
PACIFIC=70;
INDIAN=50;
SOUTHERN=90;


#Load Peter's masks and use these to create the mask files
petersMaskPath = "Peters_SasseEtAl2013_masks.nc";
pmasksAT = np.flipud(Dataset(petersMaskPath, 'r').variables["regionNoAt"][:]);
pmasksAT[oceans==LAND]= 0;
pmasksDIC = np.flipud(Dataset(petersMaskPath, 'r').variables["regionNoCt"][:]);
pmasksDIC[oceans==LAND]= 0;

#plt.figure();
#plt.imshow(np.flipud(pmasksAT));
#plt.figure();
#plt.imshow(np.flipud(pmasksDIC));

P_AT_SUBTROPICAL = 1;
P_AT_EQUATORIAL_PACIFIC = 2;
P_AT_NORTH_ATLANTIC = 3;
P_AT_NORTH_PACIFIC = 4;
P_AT_SOUTHERN_OCEAN = 5;

P_DIC_NORTH_PACIFIC = 1;
P_DIC_SOUTHERN_OCEAN = 2;
P_DIC_NORTHWEST_ATLANTIC = 3;
P_DIC_NORTHEAST_ATLANTIC = 4;
P_DIC_EQUATORIAL_PACIFIC = 5;
P_DIC_NORTH_SUBTROPICAL_PACIFIC = 6;
P_DIC_SOUTH_SUBTROPICAL_PACIFIC = 7;
P_DIC_INDIAN_OCEAN = 8;
P_DIC_NORTH_SUBTROPICAL_ATLANTIC = 9;
P_DIC_EQUATORIAL_ATLANTIC = 10;
P_DIC_SOUTH_SUBTROPICAL_ATLANTIC = 11;


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


######## AT zones ######### See Sasse et al 2013 fig 3 and supplemental table 2
#AT_zone 1: Subtropics
at_zone1 = (pmasksAT==P_AT_SUBTROPICAL).astype(int);
var = nc.createVariable("AT_zone1_subtropical", int, ("lat", "lon"));
var.units = "Integer";
var[:] = at_zone1;

#AT_zone 2: Equatorial Pacific
at_zone2 = (pmasksAT==P_AT_EQUATORIAL_PACIFIC).astype(int);
var = nc.createVariable("AT_zone2_equatorial_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = at_zone2;

#AT_zone 3: North Atlantic
at_zone3 = (pmasksAT==P_AT_NORTH_ATLANTIC).astype(int);
var = nc.createVariable("AT_zone3_north_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = at_zone3;

#AT_zone 4: North Pacific
at_zone4 = (pmasksAT==P_AT_NORTH_PACIFIC).astype(int);
var = nc.createVariable("AT_zone4_north_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = at_zone4;

#AT_zone 5: Southern Ocean
at_zone5 = (pmasksAT==P_AT_SOUTHERN_OCEAN).astype(int);
var = nc.createVariable("AT_zone5_southern_ocean", int, ("lat", "lon"));
var.units = "Integer";
var[:] = at_zone5;


######## DIC zones ######### See Sasse et al 2013 fig 3 and supplemental table 1
#DIC_zone 1: North Pacific (summer and winter)
dic_zone1 = (pmasksDIC==P_DIC_NORTH_PACIFIC).astype(int);
var = nc.createVariable("DIC_zone1_north_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone1;

#DIC_zone 2: Southern Ocean (summer and winter)
dic_zone2 = (pmasksDIC==P_DIC_SOUTHERN_OCEAN).astype(int);
var = nc.createVariable("DIC_zone2_southern_ocean", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone2;

#DIC_zone 3: Northwest Atlantic (summer and winter)
dic_zone3 = (pmasksDIC==P_DIC_NORTHWEST_ATLANTIC).astype(int);
var = nc.createVariable("DIC_zone3_northwest_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone3;

#DIC_zone 4: Northeast Atlantic (summer and winter)
dic_zone4 = (pmasksDIC==P_DIC_NORTHEAST_ATLANTIC).astype(int);
var = nc.createVariable("DIC_zone4_northeast_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone4;

#DIC_zone 5: Equatorial Pacific
dic_zone5 = (pmasksDIC==P_DIC_EQUATORIAL_PACIFIC).astype(int);
var = nc.createVariable("DIC_zone5_equatorial_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone5;

#DIC_zone 6: North Subtropical Pacific  (summer and winter)
dic_zone6 = (pmasksDIC==P_DIC_NORTH_SUBTROPICAL_PACIFIC).astype(int);
var = nc.createVariable("DIC_zone6_north_subtropical_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone6;

#DIC_zone 7: South Subtropical Pacific (summer and winter)
dic_zone7 = (pmasksDIC==P_DIC_SOUTH_SUBTROPICAL_PACIFIC).astype(int);
var = nc.createVariable("DIC_zone7_south_subtropical_pacific", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone7;

#DIC_zone 8: Indian Ocean (summer and winter)
dic_zone8 = (pmasksDIC==P_DIC_INDIAN_OCEAN).astype(int);
var = nc.createVariable("DIC_zone8_indian_ocean", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone8;

#DIC_zone 9: North Subtropical Atlantic (summer and winter)
dic_zone9 = (pmasksDIC==P_DIC_NORTH_SUBTROPICAL_ATLANTIC).astype(int);
var = nc.createVariable("DIC_zone9_north_subtropical_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone9;

#DIC_zone 10: Equatorial Atlantc
dic_zone10 = (pmasksDIC==P_DIC_EQUATORIAL_ATLANTIC).astype(int);
var = nc.createVariable("DIC_zone10_equatorial_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone10;

#DIC_zone 11: South Subtropical Atlantic
dic_zone11 = (pmasksDIC==P_DIC_SOUTH_SUBTROPICAL_ATLANTIC).astype(int);
var = nc.createVariable("DIC_zone11_south_subtropical_atlantic", int, ("lat", "lon"));
var.units = "Integer";
var[:] = dic_zone11;

nc.close();



#sanity check for overlaps
zoneSumAT = at_zone1+at_zone2+at_zone3+at_zone4+at_zone5;
plt.figure(); plt.title("AT zone sum");
plt.imshow(zoneSumAT);
plt.colorbar();
plt.ylim(0,180);

zoneSumDIC = dic_zone1+dic_zone2+dic_zone3+dic_zone4+dic_zone5+dic_zone6+dic_zone7+dic_zone8+dic_zone9+dic_zone10+dic_zone11;
plt.figure(); plt.title("DIC zone sum");
plt.imshow(zoneSumDIC);
plt.colorbar();
plt.ylim(0,180);

