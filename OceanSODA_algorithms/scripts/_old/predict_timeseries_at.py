#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:39:41 2020

@author: tom holding
"""

import osoda_global_settings;

import algorithms.at_algorithms;

from string import Template;
from os import path;
import os;
from netCDF4 import Dataset;
from datetime import datetime, timedelta;
import matplotlib.pyplot as plt;


import pandas as pd;
import numpy as np;

settings = osoda_global_settings.get_default_settings();
lonRes = latRes = 1.0;

sstInputPathTemplate = Template(path.join("../../prediction_datasets/OISST_reynolds_SST/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_1.0x1.0.nc"));
sstInputVariableName = "sst_mean";
sssInputPathTemplate = Template(path.join("../../prediction_datasets/WOA_salinity/NOAA_WOA_decav_surface_salinity_month${MM}_1.0x1.0.nc"));
sssInputVariableName = "sal_mean";

ncOutputPathTemplate = Template(path.join("../output/gridded_predictions/${OUTPUTVAR}/gridded_${ALGO}_${LONRES}x${LATRES}.nc"));

algorithmList = osoda_global_settings.get_dic_algorithm_list() + osoda_global_settings.get_at_algorithm_list();
#algorithmList = [algorithms.at_algorithms.Hassoun2015_basins_at,
#              algorithms.at_algorithms.Lefevre2010_at,
#              ];
years = settings["years"];


for AlgorithmFunctor in algorithmList:
    #Create instance of the algorithm
    algorithm = AlgorithmFunctor(settings);
    
    ###########################################
    # Create an output netCDF file to write to
    outputPath = ncOutputPathTemplate.safe_substitute(OUTPUTVAR=algorithm.output_name(), ALGO=type(algorithm).__name__, LATRES=latRes, LONRES=lonRes);
    if path.exists(path.dirname(outputPath)) == False:
        os.makedirs(path.dirname(outputPath));
    
    #create NC file and dimensions
    nc = Dataset(outputPath, 'w');
    nc.createDimension("lat", 180/latRes);
    nc.createDimension("lon", 360/lonRes);
    nc.createDimension("time", len(years)*12);
    
    #dimension variables
    var = nc.createVariable("lat", float, ("lat",));
    var.units = "lat (degrees North)";
    var[:] = np.arange(-90, 90, latRes)+(0.5*latRes);
    var = nc.createVariable("lon", float, ("lon",));
    var.units = "lon (degrees East)";
    var[:] = np.arange(-180, 180, lonRes)+(0.5*lonRes);
    var = nc.createVariable("time", int, ("time",));
    var.units = "seconds since 1980-01-01";
    var[:] = [int((datetime(year, imonth+1, 1)-datetime(1980, 1, 1)).total_seconds()) for year in years for imonth in range(0, 12)];
    
    #create gridded variables
    var = nc.createVariable(algorithm.output_name()+"_pred", float, ("time", "lat", "lon"));
    var.units = "umol kg-1";
    var.long_name = algorithm.output_name()+" predicted by "+type(algorithm).__name__+" ("+var.units+")";
    
    for iyear, year in enumerate(years):
        for imonth in range(0, 12):
            print(AlgorithmFunctor.__name__+":\t", year, format(imonth+1, "02"));
            
            #########################
            # read input netCDF file
            sstnc = Dataset(sstInputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02")));
            sst = sstnc.variables[sstInputVariableName][0,:,:];
            sst += 273.15; #convert from C to K
            sstlon = sstnc.variables["lon"][:];
            sstlat = sstnc.variables["lat"][:];
            
            basedate = datetime(1981, 1, 1);
            curDate = basedate + timedelta(seconds=int(sstnc.variables["time"][0].data));
            curDate = pd.to_datetime(curDate)
            
            sssnc = Dataset(sssInputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02")));
            sss = sssnc.variables[sssInputVariableName][0,:,:];
            if (np.any(sssnc.variables["lon"][:] != sstlon)) | np.any(sssnc.variables["lat"][:] != sstlat):
                raise ValueError("Input matrix dimensions (from gridded netCDF file) do not match for SSS.");
            
            df = pd.DataFrame(sst).stack(dropna=False);
            df = df.rename_axis(["ilat", "ilon"]).reset_index(name='SST');
            df["lon"] = sstlon[df["ilon"]];
            df["lat"] = sstlat[df["ilat"]];
            df["date"] = [curDate]*len(df);
            
            sssdf = pd.DataFrame(sss).stack(dropna=False);
            sssdf = sssdf.rename_axis(["ilat", "ilon"]).reset_index(name='SSS');
            df["SSS"] = sssdf["SSS"];
            
            ################
            # run algorithm
            modelOutput, dataUsed = algorithm(df, predict=True);
            df[algorithm.output_name()+"_pred"] = modelOutput;
            
            
            ###################################
            # write predicted output to netCDF
            griddedOutput = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred"); #unstack into a grid again
            griddedOutput = np.array(griddedOutput);
            var[(iyear*12)+imonth, :, :] = griddedOutput;
    
    
    #After all months and years computers for the current algorith, close the netCDF file.
    nc.close();




#tmp = Dataset("/home/rr/Files/Tasks/20190816_OceanSODA/OceanSODA_algorithms/output/gridded_predictions/_tmp/gridded_Lee2000_dic_1.0x1.0.nc", 'r');
#a = tmp.variables["DIC_pred"][:];
#
#plt.figure();
#plt.imshow(a[0, :, :]);
#plt.figure();
#plt.imshow(a[39, :, :]);






