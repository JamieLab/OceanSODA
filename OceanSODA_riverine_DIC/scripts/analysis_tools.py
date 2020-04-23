#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  7 13:37:27 2020

@author: verwirrt
"""

import pyproj;
import numpy as np;
from datetime import datetime;
import pandas as pd;
from os import path, makedirs;
from netCDF4 import Dataset;


####################
# Simple utilities #
####################

#Convert from concentration (umol kg-1) to quantity (mol).
def concentration_to_moles(data):
    return data * 1029 * 1000000.0; #1m^3 seawater weights approx. 1029 kg.

#moles C to TgC conversion
def mol_to_TgC(valueInMols):
    return valueInMols * 12.0107 / 1000000000000.0; #convert from mols to TgC

#return an array of datatime objects for the start of each month for an inclusive date range
def create_date_month_array(startDate, endDateInclusive):
    dates = [];
    for year in range(startDate.year, endDateInclusive.year+1):
        for month in range(1, 13):
            newDatetime = datetime.strptime("01-%s-%s"%(format(month, "02d"), str(year)), "%d-%m-%Y");
            if (newDatetime >= startDate) and (newDatetime <= endDateInclusive):
                dates.append(newDatetime);
    return dates;



########################################
# Read, process or generate input data #
########################################

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


#Extracts DIC and SSS data for a specific region and time point
def extract_data(carbonateParameters, regionMask, t):
    dic = carbonateParameters.variables["DIC_pred"][t,:,:];
    dic[regionMask != 1] = np.nan;
    dic_err = carbonateParameters.getncattr("algorithmRMSDe");
    
    sss = carbonateParameters.variables["SSS"][t,:,:];
    sss[regionMask != 1] = np.nan;
    sss[sss.mask] = np.nan;
    sss_err = carbonateParameters.variables["SSS_err"][t,:,:];
    sss_err[regionMask != 1] = np.nan;
    sss_err[sss_err.mask] = np.nan;
    
    return dic, dic_err, sss, sss_err;


#resample river discharge input to monthly totals. Removes non-overlapping time points
#scaleFactor allows conversion e.g. from discharge per second to discharge per day (must be set correctly to ensure each original discharge amount of correct when scaled to the full original time period)
def resample_discharge_monthly(riverOutflow, carbonateDates, scaleFactor=60*60*24):
    monthlyDF = pd.DataFrame();
    monthlyDF["date"] = carbonateDates;
    monthlyDF.index = carbonateDates;
    
    #Calculate monthly discharge
    monthlyOutflow = riverOutflow.resample("M", label="left", loffset="1D", on="date").sum();
    monthlyOutflow["discharge"] = monthlyOutflow["discharge"] * scaleFactor; #convert from per second to per day (daily values already summed)
    inDateRange = (monthlyOutflow.index >= carbonateDates[0]) & (monthlyOutflow.index <= carbonateDates[-1]); #remove dates that have no carbonate data for
    monthlyDF["monthly_discharge"] = monthlyOutflow["discharge"][inDateRange];
    
    
    #Calculate monthly discharge uncertainty
    monthlyOutflowSD = riverOutflow.resample("M", label="left", loffset="1D", on="date").std();
    monthlyOutflowSD["discharge"] = monthlyOutflowSD["discharge"] * scaleFactor; #convert from per second to per day (daily values already summed)
    monthlyDF["monthly_discharge_sd"] = monthlyOutflowSD["discharge"][inDateRange];
    
    return monthlyDF;



#############################################
# Components of the DIC outflow calculation #
#############################################

#Given SSS data, or a set of sampled SSS data, calculate the 2D plume mask/s
def calculate_plume_mask(sssData, plumeSalinityThreshold=35.0):
    plumeMask = np.zeros(sssData.shape, dtype=int);
    plumeMask[np.where(sssData<plumeSalinityThreshold)] = 1;
    return plumeMask;


#Returns the total plume surface area and gridded areas for each grid cell within the plume
#If a set of plume masks (3D array) is provided, set sampleDimension=True and the total plume surface areas will be calculated for each sample.
#Sample dimension is assumed to be the first.
def calculate_plume_surface_area(plumeMask, gridAreas, calculateForSamples=False):
    griddedSurfaceArea = np.full(plumeMask.shape, np.nan, dtype=float);
    
    if calculateForSamples == False:
        griddedSurfaceArea[plumeMask==1] = gridAreas[plumeMask==1];
        totalSurfaceArea = np.nansum(griddedSurfaceArea);
        return totalSurfaceArea, griddedSurfaceArea;
    else:
        gridAreas3D = np.broadcast_to(gridAreas, plumeMask.shape);
        griddedSurfaceArea[plumeMask==1] = gridAreas3D[plumeMask==1];
        totalSurfaceArea = np.nansum(griddedSurfaceArea, axis=(1,2));
        return totalSurfaceArea, griddedSurfaceArea;


#Calculate plume depth estimates (in meters) based on Coles et al 2013 (see table 3): https://doi.org/10.1002/2013JC008981
#useModelled:   When true the salinity-plume depth relationship from the model is used. Uncertainty can't be used when this is set to True
#sampleUncertainty:    When True, the model coefficients are sampled based in uncertainty estimates
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
    
    #calculate the mean of the gridded depths
    if len(salinity.shape) == 2: #if not a collection of samples
        meanPlumeDepth = np.nanmean(plumeDepth);
    else: #first dimension is treated as sample index
        meanPlumeDepth = np.nanmean(plumeDepth, axis=(1,2));
    
    return meanPlumeDepth, plumeDepth;


#Calculates plume volume from gridded surface area and gridded plume thickness.
#calculateForSamples: When true, interprets inputs as 3D arrays for multiple samples with sample index in the first dimension
def calculate_plume_volume(griddedSurfaceArea, griddedPlumeThickness, calculateForSamples=False):
    griddedPlumeVolume = griddedSurfaceArea * griddedPlumeThickness;
    if calculateForSamples == False:
        plumeVolume = np.nansum(griddedPlumeVolume);
    else:
        plumeVolume = np.nansum(griddedPlumeVolume, axis=(1,2));
    return plumeVolume, griddedPlumeVolume;


#Uses relationship between mean plume salinity and sss in Hu 2004 (https://doi.org/10.1016/j.dsr2.2004.04.001, see section 3.4)
#Assumes DIC is perfectly conservative with salinity (no biological activity etc.)
#Uses the uses proportion mean plume salinity to salinity and assumes same relationship, e.g.:
#   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
#   calculateForSamples: When true, interprets inputs as 3D arrays for multiple samples with sample index in the first dimension
def calculate_mean_dic_sss(surfaceSalinity, surfaceDIC, plumeMask, interceptUncertaintyRatio=0.1, slopeUncertaintyRatio=0.1, calculateForSamples=False):
    if calculateForSamples==False:
        samplingSize = (1,);
    else:
        samplingSize = (surfaceSalinity.shape[0], 1, 1);
    
    intercept = np.random.normal(4.352, interceptUncertaintyRatio, size=samplingSize);
    slope = np.random.normal(0.881, slopeUncertaintyRatio, size=samplingSize);
    meanPlumeSalinity = intercept + slope*surfaceSalinity;
    meanPlumeSalinity[plumeMask!=1] = np.nan;
        
    meanPlumeProportion = meanPlumeSalinity/surfaceSalinity;
    meanPlumeDIC = surfaceDIC * meanPlumeProportion;
    return meanPlumeDIC, meanPlumeSalinity;


#Calculate the total amount of DIC in the plume, and in each plume grid cell.
#This takes griddedMeanDICConc (mean concentration of DIC, umol kg-1) and grdded plume volumes (m^3)
#calculateForSamples: When true, interprets inputs as 3D arrays for multiple samples with sample index in the first dimension
def calculate_total_plume_dic(griddedMeanDIC, griddedVolume, calculateForSamples=False):
    griddedTotalDIC = griddedMeanDIC*griddedVolume;
    
    griddedTotalDIC = concentration_to_moles(griddedTotalDIC);
    
    if calculateForSamples==True:
        totalDIC = np.nansum(griddedTotalDIC, axis=(1,2));
    else:
        totalDIC = np.nansum(griddedTotalDIC);
    
    return totalDIC, griddedTotalDIC;


#Calculates an estimate of the monthly DIC outflow from a river, given DIC plume content, river discharge, and plume volume.
#If outflowMonthlyTotalDischargeUncertainty is not None, samples are assumed with number of samples equal to the length of the first dimension
def calculate_dic_outflow(outflowMonthlyTotalDischarge, plumeVolume, totalPlumeDIC, outflowMonthlyTotalDischargeUncertainty=None):
    #discharge volume as a fraction of plume volume
    if outflowMonthlyTotalDischargeUncertainty is not None:
        sampledOutflowMonthlyTotalDischarge = np.random.normal(outflowMonthlyTotalDischarge, outflowMonthlyTotalDischargeUncertainty, size=plumeVolume.shape[0]);
    else:
        sampledOutflowMonthlyTotalDischarge = outflowMonthlyTotalDischarge;
    
    dischargeVolumeProportion = sampledOutflowMonthlyTotalDischarge / plumeVolume;
    dischargeDIC = totalPlumeDIC * dischargeVolumeProportion;
    
    return dischargeDIC;



##############################
# Summary table calculations #
##############################

#returns a data frame containing the inter-year monthly values
def calculate_inter_year_monthly_means(monthlyDF):
    def inter_year_sd(monthlyDF, colName):
        arr = [];
        for month in range(1, 13):
            monthData = monthlyDF[colName][monthlyDF["date"].dt.month == month];
            n = np.sum(np.isfinite(monthData));
            val = np.sqrt(np.nansum(monthData**2)) / n; #mean
            arr.append(val);
        return arr;
    
    def inter_year_mean(monthlyDF, colName):
        arr = [];
        for month in range(1, 13):
            monthData = monthlyDF[colName][monthlyDF["date"].dt.month == month];
            val = np.nanmean(monthData);
            arr.append(val);
        return arr;
 
    interyearDF = pd.DataFrame();
    interyearDF["month"] = range(1,13);
    interyearDF.index = interyearDF["month"];
    
    for colName in monthlyDF.keys():
        if colName == "date":
            continue; #ignore date
        if colName[-3:] == "_sd": #standard deviations
            interyearDF[colName] = inter_year_sd(monthlyDF, colName);
        else: #values
            interyearDF[colName] = inter_year_mean(monthlyDF, colName);
    
    return interyearDF;


#returns a data frame containing mean annual values
def calculate_annual_values(interyearDF):
    def calc_sum_sd(interyearDF, colName):
        return np.sqrt(np.sum(interyearDF[colName]**2));
    
    def calc_mean_sd(interyearDF, colName):
        return calc_sum_sd(interyearDF, colName) / 12;
    
    columnsToSum = ["plume_total_dic", "plume_total_dic_sd", "dic_outflow", "dic_outflow_sd"];
    columnsToAverage = ["plume_surface_area", "plume_surface_area_sd", "plume_mean_thickness", "plume_mean_thickness_sd", "plume_volume", "plume_volume_sd"];
    
    output = {};
    
    #Some columns sum, some average, to give annual means
    #handle them separately
    for colName in columnsToSum:
        if colName[-3:] == "_sd": #if a standard deviation column
            output[colName+"_annual"] = calc_sum_sd(interyearDF, colName);
        else:
            output[colName+"_annual"] = np.sum(interyearDF[colName]);
    
    for colName in columnsToAverage:
        if colName[-3:] == "_sd": #if a standard deviation column
            output[colName+"_annual"] = calc_mean_sd(interyearDF, colName);
        else:
            output[colName+"_annual"] = np.mean(interyearDF[colName]);
    
    #convert to df
    df = pd.DataFrame();
    df["parameter"] = columnsToAverage + columnsToSum;
    df["value"] = [output[key] for key in output.keys()];
    return df;


##############
# outputting #
##############

#Create an empty netCDF file ready to write results of the uncertainty analysis into
def create_netCDF_file(outputPath, carbonateParameterNC, numSamples, regionMask, gridAreas):
    if path.exists(path.dirname(outputPath)) == False:
        makedirs(path.dirname(outputPath));
    
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
    
    #data variables
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
    
    var = ncout.createVariable("plume_volume_total", float, ("time"));
    var.units = "m^3";
    var.long_name = "Total plume volume at each point in time";
    
    var = ncout.createVariable("plume_volume_total_stddev", float, ("time"));
    var.units = "m^3";
    var.long_name = "Standard deviation in total plume volume at each point in time";
    
    var = ncout.createVariable("plume_volume_gridded", float, ("time", "lat", "lon"));
    var.units = "m^3";
    var.long_name = "Plume volume in each grid cell";
    
    var = ncout.createVariable("plume_volume_gridded_stddev", float, ("time", "lat", "lon"));
    var.units = "m^3";
    var.long_name = "tandard deviation in plume volume in each grid cell";
    
    var = ncout.createVariable("plume_dic_total", float, ("time"));
    var.units = "mol";
    var.long_name = "Total DIC contained in the plume over time";
    
    var = ncout.createVariable("plume_dic_total_stddev", float, ("time"));
    var.units = "mol";
    var.long_name = "Standard deviation in the total DIC contained in the plume over time";
    
    var = ncout.createVariable("plume_dic_total_gridded", float, ("time", "lat", "lon"));
    var.units = "mol";
    var.long_name = "Total DIC contained in the plume for each grid cell";
    
    var = ncout.createVariable("plume_dic_total_gridded_stddev", float, ("time", "lat", "lon"));
    var.units = "mol";
    var.long_name = "Standard deviation in the total DIC contained in the plume for each grid cell";
    
    var = ncout.var = ncout.createVariable("discharge_dic", float, ("time"));
    var.units = "TgC";
    var.long_name = "Discharge of DIC into the ocean in teragrams of carbon";
    
    var = ncout.var = ncout.createVariable("discharge_dic_stddev", float, ("time"));
    var.units = "TgC";
    var.long_name = "Standard deviation in the discharge of DIC into the ocean in teragrams of carbon";
    
    var = ncout.createVariable("discharge", float, ("time"));
    var.units = "m^3";
    var.long_name = "Total monthly river discharge estimated from daily discharge rate from gauging station";
    
    var = ncout.createVariable("discharge_dic_between_year_mean", float, ("month"));
    var.units = "TgC";
    var.long_name = "Mean total monthly river discharge of DIC into the ocean in teragrams of carbon";
    
    var = ncout.createVariable("discharge_dic_between_year_uncertainty", float, ("month"));
    var.units = "TgC";
    var.long_name = "Standard deviation in the between year monthly mean discharge of DIC into the ocean in teragrams of carbon";
    
    var = ncout.createVariable("discharge_stddev", float, ("time"));
    var.units = "m^3";
    var.long_name = "Standard deviation of the total monthly river discharge estimated from daily discharge rate from gauging station";

    var = ncout.createVariable("region_mask", float, ("lat", "lon"));
    var.units = "integer";
    var.long_name = "Mask defining the region (1=keep, 0=discard)";
    
    var = ncout.createVariable("plume_mask", float, ("time", "lat", "lon"));
    var.units = "integer";
    var.long_name = "Mask defining region defined as the plume, over time (1=keep, 0=discard)";
    
    #Write some one off 2D data
    ncout.variables["region_mask"][:] = regionMask;
    ncout.variables["grid_areas"][:] = gridAreas;
    
    return ncout;
    

#Update an open netCDF file with gridded outputs for a single time point
def update_gridded_time_point_netCDF(ncout, timeIndex, griddedPlumeVolume, griddedPlumeVolumeUncertainty, griddedPlumeDIC, griddedPlumeDICUncertainty, plumeMask):
    ncout.variables["plume_volume_gridded"][timeIndex,:,:] = griddedPlumeVolume;
    ncout.variables["plume_volume_gridded_stddev"][timeIndex,:,:] = griddedPlumeVolumeUncertainty;
    ncout.variables["plume_dic_total_gridded"][timeIndex,:,:] = griddedPlumeDIC;
    ncout.variables["plume_dic_total_gridded_stddev"][timeIndex,:,:] = griddedPlumeDICUncertainty;
    ncout.variables["plume_mask"][timeIndex,:,:] = plumeMask;


#Write non-gridded time series outputs (e.g. best estimate time series and estimated uncertainty time series)
def write_timeseries_to_netCDF(ncout, plumeVolume, plumeVolumeUncertainty, totalPlumeDIC, totalPlumeDICUncertainty, dischargeDIC, dischargeDICUncertainty,
                               dischargeTotal, dischargeTotalUncertainty, dischargeDICBetweenYearMean, dischargeDICBetweenYearMeanUncertainty):
    ncout.variables["plume_volume_total"][:] = plumeVolume;
    ncout.variables["plume_volume_total_stddev"][:] = plumeVolumeUncertainty;
    ncout.variables["plume_dic_total"][:] = totalPlumeDIC;
    ncout.variables["plume_dic_total_stddev"][:] = totalPlumeDICUncertainty;
    
    ncout.variables["discharge_dic"][:] = dischargeDIC;
    ncout.variables["discharge_dic_stddev"][:] = dischargeDICUncertainty;
    
    ncout.variables["discharge"][:] = dischargeTotal;
    ncout.variables["discharge_stddev"][:] = dischargeTotalUncertainty;
    
    ncout.variables["discharge_dic_between_year_mean"][:] = dischargeDICBetweenYearMean;
    ncout.variables["discharge_dic_between_year_uncertainty"][:] = dischargeDICBetweenYearMeanUncertainty;




######################
# plotting functions #
######################




########## Misc / old


def average_month_plumes(plumeMask, dates):
    monthlyPlumes = np.empty((12, plumeMask.shape[1], plumeMask.shape[2]), dtype=float);
    for month in range(1, 13):
        inMonth = dates.month==month;
        n = np.sum(inMonth);
        monthlyPlumes[month-1, :, :] = np.sum(plumeMask[inMonth,:,:], axis=0) / n;
    return monthlyPlumes


def plot_average_month_plumes(monthlyPlumes):
    pass;