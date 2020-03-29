#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 12 14:42:11 2020

@author: verwirrt
"""

from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;
from string import Template;
import pyproj;
import matplotlib.pyplot as plt;
from os import path, makedirs;
from datetime import datetime, timedelta;

import sys;
if path.abspath("../../") not in sys.path:
    sys.path.append(path.abspath("../../"));
    sys.path.append(path.abspath("../../OceanSODA_algorithms/scripts"));
#import OceanSODA_algorithms.scripts.utilities as utilities;
import osoda_global_settings;

#Returns area of each grid cell (in meters^2)
def calculate_grid_areas(latRes, lonRes):
    geod = pyproj.Geod(ellps='WGS84');
    gridAreas = np.empty(shape=(int(180/latRes), int(360/lonRes)));
    for y in range(0, int(180/latRes)):
        for x in range(0, int(360/lonRes)):
            lon0 = x-180 / lonRes;
            lon1 = lon0 + lonRes;
            lat0 = y-90 / latRes;
            lat1 = lat0 + latRes;
            area, perimeter = geod.polygon_area_perimeter([lon0, lon1, lon1, lon0], [lat0, lat0, lat1, lat1]);
            gridAreas[y, x] = area;
    return gridAreas;

#Calculates plume volume and estimate the uncertainty associated with the plume volume
#   Uncertainty is estimated by creating SSS samples using a normal distibutions for each grid cell centred on SSS with a stddev of SSS_err.
#   The plume region is extracted by appling the salinity threshold ('threshold') to each sampled SSS.
#   The plume surface area is calculated by summing the area of each grid cell in plume mask
#   Plume depth is calculated by using the relationship in Coles 2013 (useModelledColesRelationship determines whether to use modelled or ANACONDAS derrived relationship).
#   The above steps are repeated numSamples times
#   Plume volume is then calculated for each SSS sample, and the standard deviation of all plume volumes is calculated and returned as the uncertainty
#   Returns time series of plume volume (with uncertainty) and time series of gridded plume volume (with uncertainty)
def calc_plume_volume(sss, sss_err, gridAreas, threshold, numSamples, useModelledColesRelationship=False, verbose=False):
    #Calculate best estimate
    plumeMask = np.zeros(sss.shape, dtype=int);
    plumeMask[np.where(sss<threshold)] = 1;
    
    griddedPlumeSurfaceArea = np.full(sss.shape, np.nan, dtype=float); #m^2
    for t in range(griddedPlumeSurfaceArea.shape[0]):
        griddedPlumeSurfaceArea[t, plumeMask[t,:,:]==1] = gridAreas[plumeMask[t,:,:]==1];
    
    griddedPlumeDepth, _ = plume_thickness_coles2013(sss, plumeMask, useModelled=useModelledColesRelationship, sampleUncertainty=False); #m
    griddedPlumeVolume = griddedPlumeSurfaceArea * griddedPlumeDepth; #m^3
    plumeVolume = np.nansum(griddedPlumeVolume, axis=(1,2));
    
    #Calculate uncertainty
    sampleVolumes = np.full((numSamples,sss.shape[0]), np.nan, dtype=float);
    sampleGriddedVolumes = np.full((numSamples,sss.shape[0],sss.shape[1],sss.shape[2]), np.nan, dtype=float);
    for i in range(0, numSamples):
        if verbose:
            print(i+1, "of", numSamples);
        
        #sample SSS and extract plume mask
        sssSample = np.random.normal(sss, sss_err);
        samplePlumeMask = np.zeros(sssSample.shape, dtype=int);
        samplePlumeMask[np.where(sssSample<threshold)] = 1;
        
        #Calculate plume surface area
        #plumeSurfaceArea = np.array([np.sum(gridAreas[samplePlumeMask[j,:,:]==1]) for j in range(0, samplePlumeMask.shape[0])]); #m^2
        sampleGriddedPlumeSurfaceArea = np.full(sssSample.shape, np.nan, dtype=float);
        for t in range(sampleGriddedPlumeSurfaceArea.shape[0]):
            sampleGriddedPlumeSurfaceArea[t, samplePlumeMask[t,:,:]==1] = gridAreas[samplePlumeMask[t,:,:]==1];
        
        #Calculate plume depth (this is a sample, using the uncertainty)
        sampleGriddedPlumeDepth, griddedPlumeDepthUncertainty = plume_thickness_coles2013(sssSample, samplePlumeMask, useModelled=useModelledColesRelationship, sampleUncertainty=True); #m
        
        #Calculate plume volume
        sampleGriddedPlumeVolume = sampleGriddedPlumeSurfaceArea*sampleGriddedPlumeDepth; #m^3
        samplePlumeVolume = np.nansum(sampleGriddedPlumeVolume, axis=(1,2)); #m^3
        
        sampleVolumes[i,:] = samplePlumeVolume;
        sampleGriddedVolumes[i,:,:,:] = sampleGriddedPlumeVolume;
    
    #calculate stddev in sampled plume volume
    plumeVolumeUncertainty = np.nanstd(sampleVolumes, axis=0);
    griddedPlumeVolumeUncertainty = np.nanstd(sampleGriddedVolumes, axis=(0));
    
    return plumeVolume, plumeVolumeUncertainty, griddedPlumeVolume, griddedPlumeVolumeUncertainty;


#Calculate plume depth estimates (in meters) based on Coles et al 2013 (see table 3): https://doi.org/10.1002/2013JC008981
#useModelled:   When true the salinity-plume depth relationship from the model is used
#salinity:      Sea surface salinity
#plumemaks:     mask describing the horizontal extent of the plume - only calculates plume volume within the mask
def plume_thickness_coles2013(salinity, plumeMask, useModelled=False, sampleUncertainty=False):
    plumeDepth = np.full(salinity.shape, np.nan, float);
    plumeDepthUncertainty = np.full(salinity.shape, np.nan, float);
    
    if useModelled:
        plumeDepth[(salinity>=15) & (salinity<22)] = 5.9;
        plumeDepth[(salinity>=22) & (salinity<31)] = 6.6;
        plumeDepth[(salinity>=31) & (salinity<34)] = 15.1;
        plumeDepth[(salinity>=34) & (salinity<35)] = 25.5;
    else:
        plumeDepth[(salinity>=15) & (salinity<22)] = 5.4;
        plumeDepthUncertainty[(salinity>=15) & (salinity<22)] = 1.3; #estimated from graph in Coles 2013 figure 13
        plumeDepth[(salinity>=22) & (salinity<31)] = 7.3;
        plumeDepthUncertainty[(salinity>=22) & (salinity<31)] = 1.7; #estimated from graph in Coles 2013 figure 13
        plumeDepth[(salinity>=31) & (salinity<34)] = 15.3;
        plumeDepthUncertainty[(salinity>=31) & (salinity<34)] = 2.0; #estimated from graph in Coles 2013 figure 13
        plumeDepth[(salinity>=34) & (salinity<35)] = 14.0;
        plumeDepthUncertainty[(salinity>=34) & (salinity<35)] = 5.5; #estimated from graph in Coles 2013 figure 13
    
    #remove non-plume areas
    plumeDepth[plumeMask==False] = np.nan;
    plumeDepthUncertainty[plumeMask==False] = np.nan;
    
    #If we want to sample plume thickness rather than propagate the uncertainty
    if sampleUncertainty == True:
        plumeDepth += np.random.normal(0, plumeDepthUncertainty);
    
    return plumeDepth;

#Uses relationship between mean plume salinity and sss in Hu 2004 (https://doi.org/10.1016/j.dsr2.2004.04.001, see section 3.4)
#Assumes DIC is perfectly conservative with salinity (no biological activity etc.)
#Uses the uses proportion mean plume salinity to salinity and assumes same relationship, e.g.:
#   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
def calc_mean_plume_dic(surfaceSalinity, surfaceDIC, interceptUncertaintyRatio=0.1, slopeUncertaintyRatio=0.1):
    intercept = np.random.normal(4.352, interceptUncertaintyRatio);
    slope = np.random.normal(0.881, slopeUncertaintyRatio);
    meanPlumeSalinity = intercept + slope*surfaceSalinity;
    
    meanPlumeProportion = meanPlumeSalinity/surfaceSalinity;
    meanPlumeDIC = surfaceDIC * meanPlumeProportion;
    return meanPlumeDIC;


#Calculates an estimate of the monthly DIC outflow from a river, given DIC plume content, river discharge, and plume volume.
def calculate_dic_outflow(outflowMonthlyTotalDischarge, plumeVolume, totalPlumeDIC, outflowMonthlyTotalDischargeUncertainty=None):
    if outflowMonthlyTotalDischargeUncertainty is None:
        outflowMonthlyTotalDischargeUncertainty = [0]*len(outflowMonthlyTotalDischarge);
    
    #discharge volume as a fraction of plume volume
    sampledOutflowMonthlyTotalDischarge = np.random.normal(outflowMonthlyTotalDischarge, outflowMonthlyTotalDischargeUncertainty);
    
    dischargeVolumeProportion = sampledOutflowMonthlyTotalDischarge / plumeVolume;
    dischargeDIC = totalPlumeDIC.data * dischargeVolumeProportion;
    dischargeDIC = dischargeDIC*12.0107 / 1000000000000.0; #convert from mols to TgC
    
    return dischargeDIC;


#Write output variables to a netCDF file.
def write_netCDF(outputPath, carbonateParameterNC, gridAreas, plumeVolume, plumeVolumeUncertainty, griddedPlumeVolume, griddedPlumeVolumeUncertainty,
                 totalPlumeDIC, totalPlumeDICUncertainty, griddedPlumeDIC, griddedPlumeDICUncertainty,
                 dischargeDIC, dischargeDICUncertainty,
                 dischargeTotal, dischargeTotalUncertainty,
                 dischargeDICBetweenYearMean, dischargeDICBetweenYearUncertainty,
                 regionMask, plumeMask, numSamples):
    ncout = Dataset(outputPath, 'w');
    ncout.createDimension("time", len(carbonateParameterNC.variables["time"]));
    ncout.createDimension("lat", len(carbonateParameterNC.variables["lat"]));
    ncout.createDimension("lon", len(carbonateParameterNC.variables["lon"]));
    ncout.createDimension("month", 12);
    ncout.DICAlgorithmName = carbonateParameterNC.getncattr("algorithmName");
    ncout.DICAlgorithmRMSDe = carbonateParameterNC.getncattr("algorithmRMSDe");
    ncout.DICAlgorithmSampleSize = carbonateParameterNC.getncattr("sampleSize");
    ncout.DICAlgorithmNumAlgosCompared = carbonateParameterNC.getncattr("num_algos_compared");
    ncout.DICAlgorithmInputCombinationName = carbonateParameterNC.getncattr("input_combination_name");
    ncout.region = carbonateParameterNC.getncattr("region");
    
    ncout.numSamplesForUncertaintyAnalysis = numSamples;
    
    #dimension variables
    var = ncout.createVariable("lat", float, ("lat",));
    var.units = "lat (degrees North)";
    var[:] = carbonateParameterNC.variables["lat"][:];
    var = ncout.createVariable("lon", float, ("lon",));
    var.units = "lon (degrees East)";
    var[:] = carbonateParameterNC.variables["lon"][:];
    var = ncout.createVariable("time", int, ("time",));
    var.units = "seconds since 1980-01-01";
    var[:] = carbonateParameterNC.variables["time"][:];
    var = ncout.createVariable("month", int, ("month",));
    var.units = "month (1=JAN, 12=DEC)";
    var[:] = range(1,13);
    
    #Write data to netCDF file
    var = ncout.createVariable("DIC_pred", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["DIC_pred"].units;
    var.long_name = carbonateParameterNC.variables["DIC_pred"].long_name;
    var[:] = carbonateParameterNC.variables["DIC_pred"][:]
    
    #Write data to netCDF file
    var = ncout.createVariable("SSS", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["SSS"].units;
    var.long_name = carbonateParameterNC.variables["SSS"].long_name;
    var[:] = carbonateParameterNC.variables["SSS"][:]
    
    #Write data to netCDF file
    var = ncout.createVariable("SSS_err", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["SSS_err"].units;
    var.long_name = carbonateParameterNC.variables["SSS_err"].long_name;
    var[:] = carbonateParameterNC.variables["SSS_err"][:]
    
    var = ncout.createVariable("grid_areas", float, ("lat", "lon"));
    var.units = "m^3";
    var.long_name = "Grid cell surface area";
    var[:] = gridAreas;
    
    var = ncout.createVariable("plume_volume_total", float, ("time"));
    var.units = "m^3";
    var.long_name = "Total plume volume at each point in time";
    var[:] = plumeVolume;
    
    var = ncout.createVariable("plume_volume_total_stddev", float, ("time"));
    var.units = "m^3";
    var.long_name = "Standard deviation in total plume volume at each point in time";
    var[:] = plumeVolumeUncertainty;
    
    var = ncout.createVariable("plume_volume_gridded", float, ("time", "lat", "lon"));
    var.units = "m^3";
    var.long_name = "Plume volume in each grid cell";
    var[:] = griddedPlumeVolume;
    
    var = ncout.createVariable("plume_volume_gridded_stddev", float, ("time", "lat", "lon"));
    var.units = "m^3";
    var.long_name = "tandard deviation in plume volume in each grid cell";
    var[:] = griddedPlumeVolumeUncertainty;
    
    var = ncout.createVariable("plume_dic_total", float, ("time"));
    var.units = "mol";
    var.long_name = "Total DIC contained in the plume over time";
    var[:] = totalPlumeDIC;
    
    var = ncout.createVariable("plume_dic_total_stddev", float, ("time"));
    var.units = "mol";
    var.long_name = "Standard deviation in the total DIC contained in the plume over time";
    var[:] = totalPlumeDICUncertainty;
    
    var = ncout.createVariable("plume_dic_total_gridded", float, ("time", "lat", "lon"));
    var.units = "mol";
    var.long_name = "Total DIC contained in the plume for each grid cell";
    var[:] = griddedPlumeDIC;
    
    var = ncout.createVariable("plume_dic_total_gridded_stddev", float, ("time", "lat", "lon"));
    var.units = "mol";
    var.long_name = "Standard deviation in the total DIC contained in the plume for each grid cell";
    var[:] = griddedPlumeDICUncertainty;
    
    var = ncout.var = ncout.createVariable("discharge_dic", float, ("time"));
    var.units = "TgC";
    var.long_name = "Discharge of DIC into the ocean in teragrams of carbon";
    var[:] = dischargeDIC;
    
    var = ncout.var = ncout.createVariable("discharge_dic_stddev", float, ("time"));
    var.units = "TgC";
    var.long_name = "Standard deviation in the discharge of DIC into the ocean in teragrams of carbon";
    var[:] = dischargeDICUncertainty;
    
    var = ncout.createVariable("discharge", float, ("time"));
    var.units = "m^3";
    var.long_name = "Total monthly river discharge estimated from daily discharge rate from gauging station";
    var[:] = dischargeTotal;
    
    var = ncout.createVariable("discharge_dic_between_year_mean", float, ("month"));
    var.units = "TgC";
    var.long_name = "Mean total monthly river discharge of DIC into the ocean in teragrams of carbon";
    var[:] = dischargeDICBetweenYearMean;
    
    var = ncout.createVariable("discharge_dic_between_year_uncertainty", float, ("month"));
    var.units = "TgC";
    var.long_name = "Standard deviation in the between year monthly mean discharge of DIC into the ocean in teragrams of carbon";
    var[:] = dischargeDICBetweenYearUncertainty;
    
    var = ncout.createVariable("discharge_stddev", float, ("time"));
    var.units = "m^3";
    var.long_name = "Standard deviation of the total monthly river discharge estimated from daily discharge rate from gauging station";
    var[:] = dischargeTotalUncertainty;

    var = ncout.createVariable("region_mask", float, ("lat", "lon"));
    var.units = "integer";
    var.long_name = "Mask defining the region (1=keep, 0=discard)";
    var[:] = regionMask;
    
    var = ncout.createVariable("plume_mask", float, ("time", "lat", "lon"));
    var.units = "integer";
    var.long_name = "Mask defining region defined as the plume, over time (1=keep, 0=discard)";
    var[:] = plumeMask;
    
    ncout.close();



#Settings and constants
carbonateParametersTemplate = Template("../../OceanSODA_algorithms/output/gridded_predictions/gridded_${REGION}_1.0x1.0_${OUTPUTVAR}.nc");
regionMaskPath = "../../OceanSODA_algorithms/region_masks/osoda_region_masks_v2.nc";
regions = ["oceansoda_amazon_plume", "oceansoda_congo", "oceansoda_mississippi", "oceansoda_st_lawrence", "oceansoda_mediterranean"];
amazonFlowPath = "/home/verwirrt/Projects/Work/20190816_OceanSODA/OceanSODA_riverine_outflow/data/obidos/17050001_debits.csv";


settings = osoda_global_settings.get_default_settings();
lonRes = latRes = 1.0;
plumeSalinityThreshold = 35.0; #Below this salinity is considered plume
numSamples = 30; #How many SSS and DIC samples to use
verbose = True;

regions = ["oceansoda_amazon_plume"];
region = "oceansoda_amazon_plume";

carbonateParameters = Dataset(carbonateParametersTemplate.safe_substitute(REGION=region, OUTPUTVAR="DIC"));
regionMaskNC = Dataset(regionMaskPath, 'r');
amazonRegionMask = regionMaskNC.variables["oceansoda_amazon_plume"][:];
amazonRegionMask3D = np.broadcast_to(amazonRegionMask, (len(carbonateParameters.variables["time"]), len(carbonateParameters.variables["lat"]), len(carbonateParameters.variables["lon"])));


#Extract carbonate parameter data
dic = carbonateParameters.variables["DIC_pred"][:];
dic[amazonRegionMask3D != 1] = np.nan;
dic[dic.mask] = np.nan;
dic_err = carbonateParameters.getncattr("algorithmRMSDe");
sss = carbonateParameters.variables["SSS"][:];
sss[amazonRegionMask3D != 1] = np.nan;
sss[sss.mask] = np.nan;
sss_err = carbonateParameters.variables["SSS_err"][:];
sss_err[amazonRegionMask3D != 1] = np.nan;
sss_err[sss_err.mask] = np.nan;
carbonateParamsTimeRange = (datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][0])), datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][-1])));


#Read and pre-process amazon discharge data
amazonOutflow = pd.read_csv(amazonFlowPath, parse_dates=["date"], dayfirst=True);
#carbonate data is monthly averages, so calculate monthly means
amazonOutflow.index = amazonOutflow["date"];
amazonOutflowMonthlyTotal = amazonOutflow.resample("M").sum() * 60*60*24; #convert from per second to per day and sum days to give total monthly discharge (m^3)
amazonOutflowMonthlyTotalSD = amazonOutflow.resample("M").std()["discharge"] * 60*60*24;

amazonOutflowMonthlyTotal.index = [pd.datetime(d.year, d.month, 1) for d in amazonOutflowMonthlyTotal.index]; #Set index date to the start of the month
amazonOutflowMonthlyTotalSD.index = amazonOutflowMonthlyTotal.index;

#subset to just the time range we're interested in (note, this is done after monthly resampling so that whole months are preserved)
toKeep = (amazonOutflowMonthlyTotal.index >= carbonateParamsTimeRange[0]) & (amazonOutflowMonthlyTotal.index <= carbonateParamsTimeRange[1]);
amazonOutflowMonthlyTotal = amazonOutflowMonthlyTotal[toKeep];
amazonOutflowMonthlyTotalSD = amazonOutflowMonthlyTotalSD[toKeep];


#Calculate grid cell surface areas
gridAreas = calculate_grid_areas(latRes, lonRes);



for region in regions:
    #####################
    ### Estimate DIC outflow from the river. Uncertainty is estimated by resampling SSS, DIC etc. according to RMSD
    #####################
    ### Best estimates calculated first
    plumeMask = np.zeros(sss.shape, dtype=int);
    plumeMask[np.where(sss<plumeSalinityThreshold)] = 1;
    
    griddedPlumeSurfaceArea = np.full(sss.shape, np.nan, dtype=float); #m^2
    for t in range(griddedPlumeSurfaceArea.shape[0]):
        griddedPlumeSurfaceArea[t, plumeMask[t,:,:]==1] = gridAreas[plumeMask[t,:,:]==1];
    
    griddedPlumeDepth = plume_thickness_coles2013(sss, plumeMask, useModelled=False, sampleUncertainty=False); #m
    griddedPlumeVolume = griddedPlumeSurfaceArea * griddedPlumeDepth; #m^3
    plumeVolume = np.nansum(griddedPlumeVolume, axis=(1,2));
    
    #Calculate gridded mean plume DIC 
    #This eses relationship between mean plume salinity and sss in Hu 2004 (PS = 0.881*SSS+4.352, see section 3.4 of https://doi.org/10.1016/j.dsr2.2004.04.001)
    #and assumes DIC is perfectly conservative with salinity
    #I.e. the proportion mean plume salinity to surface salinity is assumed to equal mean plume DIC / surface DIC, thus rearranging gives:
    #   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
    meanGriddedPlumeDIC = calc_mean_plume_dic(sss, dic, interceptUncertaintyRatio=0.0, slopeUncertaintyRatio=0.0);
    
    #Estimate total DIC in the plume
    gridPlumeDIC = meanGriddedPlumeDIC*griddedPlumeVolume;
    gridPlumeDIC = gridPlumeDIC * 1029 / 1000000.0; #Convert from concentration (umol kg-1) to quantity (mol). 1m^3 seawater weights approx. 1029 kg.
    totalPlumeDIC = np.nansum(gridPlumeDIC, axis=(1,2));
    
    #calculate DIC outflow
    dischargeDIC = calculate_dic_outflow(amazonOutflowMonthlyTotal["discharge"], plumeVolume, totalPlumeDIC, outflowMonthlyTotalDischargeUncertainty=amazonOutflowMonthlyTotalSD);
    
    #################
    ### Now calculate uncertainty
    sampleVolumes = np.full((numSamples,sss.shape[0]), np.nan, dtype=float);
    sampleGriddedVolumes = np.full((numSamples,sss.shape[0],sss.shape[1],sss.shape[2]), np.nan, dtype=float);
    sampleTotalPlumeDICs = np.full((numSamples,sss.shape[0]), np.nan, dtype=float);
    sampleMeanGriddedPlumeDICs = np.full((numSamples,sss.shape[0],sss.shape[1],sss.shape[2]), np.nan, dtype=float);
    sampleDischargeDICs = np.full((numSamples,sss.shape[0]), np.nan, dtype=float);
    for i in range(0, numSamples):
        if verbose:
            print(i+1, "of", numSamples);
        
        #sample SSS and dic
        #extract plume mask
        sssSample = np.random.normal(sss, sss_err);
        dicSample = np.random.normal(dic, dic_err);
        
        samplePlumeMask = np.zeros(sssSample.shape, dtype=int);
        samplePlumeMask[np.where(sssSample<plumeSalinityThreshold)] = 1;
        
        #Calculate plume surface area
        #plumeSurfaceArea = np.array([np.sum(gridAreas[samplePlumeMask[j,:,:]==1]) for j in range(0, samplePlumeMask.shape[0])]); #m^2
        sampleGriddedPlumeSurfaceArea = np.full(sssSample.shape, np.nan, dtype=float);
        for t in range(sampleGriddedPlumeSurfaceArea.shape[0]):
            sampleGriddedPlumeSurfaceArea[t, samplePlumeMask[t,:,:]==1] = gridAreas[samplePlumeMask[t,:,:]==1];
        
        #Calculate plume depth (this is a sample, using the uncertainty)
        sampleGriddedPlumeDepth = plume_thickness_coles2013(sssSample, samplePlumeMask, useModelled=False, sampleUncertainty=True); #m
        
        #Calculate plume volume
        sampleGriddedPlumeVolume = sampleGriddedPlumeSurfaceArea*sampleGriddedPlumeDepth; #m^3
        samplePlumeVolume = np.nansum(sampleGriddedPlumeVolume, axis=(1,2)); #m^3
        #store volume data for each sample (to calculate and output uncertainties later)
        sampleVolumes[i,:] = samplePlumeVolume;
        sampleGriddedVolumes[i,:,:,:] = sampleGriddedPlumeVolume;
        
        #Calculate gridded mean plume DIC of the sample
        #This eses relationship between mean plume salinity and sss in Hu 2004 (PS = 0.881*SSS+4.352, see section 3.4 of https://doi.org/10.1016/j.dsr2.2004.04.001)
        #and assumes DIC is perfectly conservative with salinity
        #I.e. the proportion mean plume salinity to surface salinity is assumed to equal mean plume DIC / surface DIC, thus rearranging gives:
        #   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
        sampleMeanGriddedPlumeDIC = calc_mean_plume_dic(sssSample, dicSample);
        
        #Estimate total DIC in the plume for this sample
        sampleGridPlumeDIC = sampleMeanGriddedPlumeDIC*sampleGriddedPlumeVolume;
        sampleGridPlumeDIC = sampleGridPlumeDIC * 1029 / 1000000.0; #Convert from concentration (umol kg-1) to quantity (mol). 1m^3 seawater weights approx. 1029 kg.
        sampleTotalPlumeDIC = np.nansum(sampleGridPlumeDIC, axis=(1,2));
        
        #Store plume DIC data
        sampleTotalPlumeDICs[i, :] = sampleTotalPlumeDIC;
        sampleMeanGriddedPlumeDICs[i,:,:,:] = sampleMeanGriddedPlumeDIC
        
        #Calculate DIC outflow
        sampledDischargeDIC = calculate_dic_outflow(amazonOutflowMonthlyTotal["discharge"], samplePlumeVolume, sampleTotalPlumeDIC, outflowMonthlyTotalDischargeUncertainty=amazonOutflowMonthlyTotalSD);
        sampleDischargeDICs[i,:] = sampledDischargeDIC;
    
    #calculate stddev in sampled plume volume
    plumeVolumeUncertainty = np.nanstd(sampleVolumes, axis=0);
    griddedPlumeVolumeUncertainty = np.nanstd(sampleGriddedVolumes, axis=(0));
    totalPlumeDICUncertainty = np.nanstd(sampleTotalPlumeDICs, axis=0);
    meanGriddedPlumeDICUncertainty = np.nanstd(sampleMeanGriddedPlumeDICs, axis=(0));
    dischargeDICUncertainty = np.nanstd(sampleDischargeDICs, axis=0);
    
    if False: #Visualise plume volume area over time (with estimated uncertainty)
        plt.figure();
        plt.fill_between(range(0,len(plumeVolume)), plumeVolume+plumeVolumeUncertainty, plumeVolume-plumeVolumeUncertainty, alpha=0.75);
        plt.plot(plumeVolume, 'k', alpha=0.5);
    if False: #Visualise total plume DIC over time (with estimated uncertainty)
        plt.figure();
        plt.fill_between(range(0,len(totalPlumeDIC)), totalPlumeDIC+totalPlumeDICUncertainty, totalPlumeDIC-totalPlumeDICUncertainty, alpha=0.75);
        plt.plot(totalPlumeDIC, 'k', alpha=0.5);
    if False: #Visualise DIC discharge over time (with estimated uncertainty)
        plt.figure();
        plt.fill_between(range(0,len(dischargeDIC)), dischargeDIC+dischargeDICUncertainty, dischargeDIC-dischargeDICUncertainty, alpha=0.75);
        plt.plot(dischargeDIC, 'k', alpha=0.5);
        plt.title("DIC discharge from the Amazon");
        plt.xlabel("time (months)");
        plt.ylabel("Monthly discharge of DIC (TgC)");
    
    
    #Summary stats (whole period)
    dischargeDICMean = np.mean(dischargeDIC);
    dischargeDICSD = np.std(dischargeDIC);
    amazonOutflowMonthlyTotal["dischargeDIC"] = dischargeDIC;
    amazonOutflowMonthlyTotal["dischargeDICUncertainty"] = dischargeDICUncertainty;
    amazonOutflowMonthlyTotal["year"] = [d.year for d in amazonOutflowMonthlyTotal.index];
    amazonOutflowMonthlyTotal["month"] = [d.month for d in amazonOutflowMonthlyTotal.index];
    
    dischargeDICBetweenYearGroup = amazonOutflowMonthlyTotal.groupby(["month"]);
    dischargeDICBetweenYearMean = dischargeDICBetweenYearGroup["dischargeDIC"].mean();
    #dischargeDICBetweenYearSD = dischargeDICBetweenYearGroup["dischargeDIC"].std();
    
    dischargeDICBetweenYearUncertainty = np.zeros((12,), dtype=float);
    for month in range(1,13):
        monthUncertainties = amazonOutflowMonthlyTotal[amazonOutflowMonthlyTotal["month"]==month]["dischargeDICUncertainty"];
        #Propagate error, sandard deviation of the mean, see Taylor section 4.4
        dischargeDICBetweenYearUncertainty[month-1] = np.nansum(monthUncertainties) / np.sqrt(np.sum(np.isfinite(monthUncertainties))); #sum reuslt of isfinite instead of length as to exclude nans
    
    if False:
        plt.figure();
        plt.fill_between(range(1, 13), dischargeDICBetweenYearMean-dischargeDICBetweenYearUncertainty, dischargeDICBetweenYearMean+dischargeDICBetweenYearUncertainty, alpha=0.5);
        plt.plot(range(1, 13), dischargeDICBetweenYearMean);
        plt.xlabel("month");
        plt.ylabel("Amazon DIC source to the ocean (TgC)");
        plt.title("Between year monthly mean");

    
    ### Write to netCDF
    outputPath = "../output/dic_outflow_"+region+".nc";
    if path.exists(path.dirname(outputPath)) == False:
        makedirs(path.dirname(outputPath));
    write_netCDF(outputPath, carbonateParameters, gridAreas, plumeVolume, plumeVolumeUncertainty, griddedPlumeVolume, griddedPlumeVolumeUncertainty,
                 totalPlumeDIC, totalPlumeDICUncertainty, meanGriddedPlumeDIC, meanGriddedPlumeDICUncertainty,
                 dischargeDIC, dischargeDICUncertainty, amazonOutflowMonthlyTotal["discharge"].as_matrix(), amazonOutflowMonthlyTotalSD.as_matrix(),
                 dischargeDICBetweenYearMean.as_matrix(), dischargeDICBetweenYearUncertainty,
                 regionMaskNC.variables[region][:], plumeMask, numSamples);
    
    #Summary stats (whole period)
    dischargeAnnual = np.sum(dischargeDICBetweenYearMean);
    dischargeAnnualUncertainty = np.sqrt(np.nansum([uncertainty**2 for uncertainty in dischargeDICBetweenYearUncertainty]));
    
    #Table output
    print("DIC annual discharge (PgC)", dischargeAnnual/1000, dischargeAnnualUncertainty/1000, sep=",");





