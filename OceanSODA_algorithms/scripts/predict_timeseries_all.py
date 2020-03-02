#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:39:41 2020

@author: tom holding
"""

import osoda_global_settings;

import algorithms.at_algorithms;
import algorithms.dic_algorithms;

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
sssInputPathTemplate = Template(path.join("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc"));
sssInputVariableName = "salinity_mean";
doInputPathTemplate = Template(path.join("../../prediction_datasets/WOA_dissolved_oxygen/woa18_all_o${MM}_01.nc"));
doInputVariableName = "o_an";
no3InputPathTemplate = Template(path.join("../../prediction_datasets/WOA_nitrate/woa18_all_n${MM}_01.nc"));
no3InputVariableName = "n_an";
sio4InputPathTemplate = Template(path.join("../../prediction_datasets/WOA_silicate/woa18_all_i${MM}_01.nc"));
sio4InputVariableName = "i_an";
po4InputPathTemplate = Template(path.join("../../prediction_datasets/WOA_phosphate/woa18_all_p${MM}_01.nc"));
po4InputVariableName = "p_an";

ncOutputPathTemplate = Template(path.join("../output/gridded_predictions/${OUTPUTVAR}/gridded_${ALGO}_${LONRES}x${LATRES}.nc"));

algorithmList = osoda_global_settings.get_dic_algorithm_list() + osoda_global_settings.get_at_algorithm_list();
#algorithmList = [algorithms.dic_algorithms.Cooley2006a_dic];
#algorithmList = [algorithms.at_algorithms.Hassoun2015_basins_at,
#              algorithms.at_algorithms.Lefevre2010_at,
#


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
    
    var = nc.createVariable("SSS", float, ("time", "lat", "lon"));
    var.units = "PSU";
    var.long_name = "Sea surface salinity used for prediction";
    
    var = nc.createVariable("SST", float, ("time", "lat", "lon"));
    var.units = "Kelvin (k)";
    var.long_name = "Sea surface temperature used for prediction";
    
    for iyear, year in enumerate(years):
        for imonth in range(0, 12):
            print(AlgorithmFunctor.__name__+":\t", year, format(imonth+1, "02"));
            
            #########################
            # read input netCDF file
            try:
                sstnc = Dataset(sstInputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02")));
                sst = sstnc.variables[sstInputVariableName][0,:,:];
                sst += 273.15; #convert from C to K
                sstlon = sstnc.variables["lon"][:];
                sstlat = sstnc.variables["lat"][:];
                
                basedate = datetime(1981, 1, 1);
                curDate = basedate + timedelta(seconds=int(sstnc.variables["time"][0].data));
                curDate = pd.to_datetime(curDate);
                
                sssnc = Dataset(sssInputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02")));
                sss = sssnc.variables[sssInputVariableName][:,:];
                if (np.any(sssnc.variables["lon"][:] != sstlon)) | np.any(sssnc.variables["lat"][:] != sstlat):
                    raise ValueError("Input matrix dimensions (from gridded netCDF file) do not match for SSS.");
                
                donc = Dataset(doInputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02")));
                do = donc.variables[doInputVariableName][0,0,:,:]; #time, depth, lat, lon
                
                no3 = Dataset(no3InputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02"))).variables[no3InputVariableName][0,0,:,:];
                po4 = Dataset(po4InputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02"))).variables[po4InputVariableName][0,0,:,:];
                sio4 = Dataset(sio4InputPathTemplate.safe_substitute(YYYY=year, MM=format(imonth+1, "02"))).variables[sio4InputVariableName][0,0,:,:];
                
                
                df = pd.DataFrame(sst).stack(dropna=False);
                df = df.rename_axis(["ilat", "ilon"]).reset_index(name='SST');
                df["lon"] = sstlon[df["ilon"]];
                df["lat"] = sstlat[df["ilat"]];
                df["date"] = [curDate]*len(df);
                
                sssdf = pd.DataFrame(sss).stack(dropna=False);
                sssdf = sssdf.rename_axis(["ilat", "ilon"]).reset_index(name='SSS');
                df["SSS"] = sssdf["SSS"];
                
                dodf = pd.DataFrame(do).stack(dropna=False);
                dodf = dodf.rename_axis(["ilat", "ilon"]).reset_index(name='DO');
                df["DO"] = dodf["DO"];
                
                no3df = pd.DataFrame(no3).stack(dropna=False);
                no3df = no3df.rename_axis(["ilat", "ilon"]).reset_index(name='NO3');
                df["NO3"] = no3df["NO3"];
                
                po4df = pd.DataFrame(po4).stack(dropna=False);
                po4df = po4df.rename_axis(["ilat", "ilon"]).reset_index(name='PO4');
                df["PO4"] = po4df["PO4"];
                
                sio4df = pd.DataFrame(sio4).stack(dropna=False);
                sio4df = sio4df.rename_axis(["ilat", "ilon"]).reset_index(name='SiO4');
                df["SiO4"] = sio4df["SiO4"];
            except FileNotFoundError as e:
                print("Skipping year", year, format(imonth+1, "02"), "because there was no prediction data for one or more variable(s):", e.args);
                continue;
                
            
            ################
            # run algorithm
            try:
                modelOutput, dataUsed = algorithm(df, predict=True);
                df[algorithm.output_name()+"_pred"] = modelOutput;
            except ValueError:
                print("No data within valid ranges for "+algorithm.__class__.__name__+". No predictions could be made.");
                df[algorithm.output_name()+"_pred"] = [np.nan]*len(df);
            
            ###################################
            # write predicted output to netCDF
            griddedOutput = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred"); #unstack into a grid again
            griddedOutput = np.array(griddedOutput);
            nc.variables[algorithm.output_name()+"_pred"][(iyear*12)+imonth, :, :] = griddedOutput;
            sss[np.where(np.isfinite(griddedOutput)==False)] = np.nan;
            nc.variables["SSS"][(iyear*12)+imonth, :, :] = sss;
            sst[np.where(np.isfinite(griddedOutput)==False)] = np.nan;
            nc.variables["SST"][(iyear*12)+imonth, :, :] = sst;
            
            
    
    
    #After all months and years computers for the current algorith, close the netCDF file.
    nc.close();




#tmp = Dataset("/home/rr/Files/Tasks/20190816_OceanSODA/OceanSODA_algorithms/output/gridded_predictions/_tmp/gridded_Lee2000_dic_1.0x1.0.nc", 'r');
#a = tmp.variables["DIC_pred"][:];
#
#plt.figure();
#plt.imshow(a[0, :, :]);
#plt.figure();
#plt.imshow(a[39, :, :]);






