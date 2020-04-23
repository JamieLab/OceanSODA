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
import matplotlib.pyplot as plt;
from os import path, makedirs;
from datetime import datetime, timedelta;

import analysis_tools;

import sys;
if path.abspath("../../") not in sys.path:
    sys.path.append(path.abspath("../../"));
    sys.path.append(path.abspath("../../OceanSODA_algorithms/scripts"));
#import OceanSODA_algorithms.scripts.utilities as utilities;
import osoda_global_settings;


##### Global settings
writeNetCDFOutput = True;
writeCSVOutput = True;
baseOutputPath = "../output/";
lonRes = latRes = 1.0;
plumeSalinityThreshold = 35.0; #Below this salinity is considered plume
numSamples = 3;#25; #How many SSS and DIC samples to use
verbose = True;
#regions = ["oceansoda_amazon_plume", "oceansoda_congo", "oceansoda_mississippi", "oceansoda_st_lawrence", "oceansoda_mediterranean"];
regions = ["oceansoda_amazon_plume"];
region = "oceansoda_amazon_plume";

carbonateParametersTemplate = Template("../../OceanSODA_algorithms/output/gridded_predictions/gridded_${REGION}_1.0x1.0_${OUTPUTVAR}.nc");
regionMaskPath = "../../OceanSODA_algorithms/region_masks/osoda_region_masks_v2.nc";
riverDischargePaths = {"oceansoda_amazon_plume": "/home/verwirrt/Projects/Work/20190816_OceanSODA/OceanSODA_riverine_outflow/data/obidos/17050001_debits.csv",
                       };

settings = osoda_global_settings.get_default_settings();
gridAreas = analysis_tools.calculate_grid_areas(latRes, lonRes);






#Read data sets
carbonateParameters = Dataset(carbonateParametersTemplate.safe_substitute(REGION=region, OUTPUTVAR="DIC"));
regionMaskNC = Dataset(regionMaskPath, 'r');


#extract mask data
regionMask = regionMaskNC.variables["oceansoda_amazon_plume"][:];
#regionMask3D = np.broadcast_to(amazonRegionMask, (len(carbonateParameters.variables["time"]), len(carbonateParameters.variables["lat"]), len(carbonateParameters.variables["lon"])));

#Time data
numTimePoints = len(carbonateParameters.variables["time"]);
carbonateParamsTimeRange = (datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][0])), datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][-1])));
carbonateDates = analysis_tools.create_date_month_array(carbonateParamsTimeRange[0], carbonateParamsTimeRange[1]);


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
    

def update_gridded_time_point_netCDF(ncout, timeIndex, griddedPlumeVolume, griddedPlumeVolumeUncertainty, griddedPlumeDIC, griddedPlumeDICUncertainty, plumeMask):
    ncout.variables["plume_volume_gridded"][timeIndex,:,:] = griddedPlumeVolume;
    ncout.variables["plume_volume_gridded_stddev"][timeIndex,:,:] = griddedPlumeVolumeUncertainty;
    ncout.variables["plume_dic_total_gridded"][timeIndex,:,:] = griddedPlumeDIC;
    ncout.variables["plume_dic_total_gridded_stddev"][timeIndex,:,:] = griddedPlumeDICUncertainty;
    ncout.variables["plume_mask"][timeIndex,:,:] = plumeMask;


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
    


for region in regions:
    ### Create netCDF file to store results - this will be updated with each iteration
    if writeNetCDFOutput == True:
        outputPathNetCDF = path.join(baseOutputPath, "dic_outflow_"+region+".nc");
        ncout = create_netCDF_file(outputPathNetCDF, carbonateParameters, numSamples, regionMask, gridAreas);
        
        
    #Read and pre-process amazon discharge data
    riverDischarge = pd.read_csv(riverDischargePaths[region], parse_dates=["date"], dayfirst=True);
    monthlyDischarge = analysis_tools.resample_discharge_monthly(riverDischarge, carbonateDates, scaleFactor=60*60*24); #Scale factor to scale per second discharge to per day
    
    
    #Create a data frame to contain all the values used in the calculation
    #River discharge. carbonate data is monthly averages, so calculate monthly means of river discharge.
    monthlyDF = pd.DataFrame();
    monthlyDF["date"] = carbonateDates;
    monthlyDF.index = carbonateDates;
    monthlyDF["plume_surface_area"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_surface_area_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_mean_thickness"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_mean_thickness_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_volume"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_volume_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_total_dic"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["plume_total_dic_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["dic_outflow"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["dic_outflow_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
    monthlyDF["river_discharge"] = monthlyDischarge["monthly_discharge"];
    monthlyDF["river_discharge_sd"] = monthlyDischarge["monthly_discharge_sd"];
    
       
    
    
    #Split calculations by time point because memory usage is too large to store many samples for all time points at once
    for t in range(0, numTimePoints):
        currentDatetime = carbonateDates[t];
        #extract data for this time point
        dic, dic_err, sss, sss_err = analysis_tools.extract_data(carbonateParameters, regionMask, t=t);
        
        #Don't perform computations if everything is nan
        if np.all(np.isfinite(dic)==False):
            continue;
        
        #Generate DIC and SSS samplestime point
        sssSamples = np.random.normal(sss, sss_err, size=(numSamples,)+sss.shape);
        dicSamples = np.random.normal(dic, dic_err, size=(numSamples,)+dic.shape);
        
        #Calculate plume masks (1==inside plume)
        plumeMask = analysis_tools.calculate_plume_mask(sss, plumeSalinityThreshold=plumeSalinityThreshold);
        plumeMaskSamples = analysis_tools.calculate_plume_mask(sssSamples, plumeSalinityThreshold=plumeSalinityThreshold);
        
        
        #Calculate plume surface area (m^2)
        surfaceArea, griddedSurfaceArea = analysis_tools.calculate_plume_surface_area(plumeMask, gridAreas);
        surfaceAreaSamples, griddedSurfaceAreaSamples = analysis_tools.calculate_plume_surface_area(plumeMaskSamples, gridAreas, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_surface_area"] = surfaceArea;
        monthlyDF.loc[currentDatetime, "plume_surface_area_sd"] = np.std(surfaceAreaSamples);
        
        
        #Calculate plume thickness(m)
        meanPlumeThickness, griddedPlumeThickness = analysis_tools.plume_thickness_coles2013(sss, plumeMask, useModelled=False, sampleUncertainty=False);
        meanPlumeThicknessSamples, griddedPlumeThicknessSamples = analysis_tools.plume_thickness_coles2013(sssSamples, plumeMaskSamples, useModelled=False, sampleUncertainty=True);
        monthlyDF.loc[currentDatetime, "plume_mean_thickness"] = meanPlumeThickness;
        monthlyDF.loc[currentDatetime, "plume_mean_thickness_sd"] = np.std(meanPlumeThicknessSamples);
        
        
        #Calculate plume volume (m^3)
        plumeVolume, griddedPlumeVolume = analysis_tools.calculate_plume_volume(griddedSurfaceArea, griddedPlumeThickness, calculateForSamples=False);
        plumeVolumeSamples, griddedPlumeVolumeSamples = analysis_tools.calculate_plume_volume(griddedSurfaceAreaSamples, griddedPlumeThicknessSamples, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_volume"] = plumeVolume;
        monthlyDF.loc[currentDatetime, "plume_volume_sd"] = np.std(plumeVolumeSamples);
        
        
        #Calculate gridded mean plume DIC, in the process also calculate meanplume salinity (umol kg-1)
        #This uses relationship between mean plume salinity and sss in Hu 2004 (PS = 0.881*SSS+4.352, see section 3.4 of https://doi.org/10.1016/j.dsr2.2004.04.001)
        #and assumes DIC is perfectly conservative with salinity
        #I.e. the proportion mean plume salinity to surface salinity is assumed to equal mean plume DIC / surface DIC, thus rearranging gives:
        #   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
        griddedMeanDICConc, griddedMeanSSS = analysis_tools.calculate_mean_dic_sss(sss, dic, plumeMask, interceptUncertaintyRatio=0.0, slopeUncertaintyRatio=0.0, calculateForSamples=False);
        griddedMeanDICConcSamples, griddedMeanSSSSamples = analysis_tools.calculate_mean_dic_sss(sssSamples, dicSamples, plumeMaskSamples, interceptUncertaintyRatio=0.1, slopeUncertaintyRatio=0.1, calculateForSamples=True);
        
        
        #calculate total plume DIC (mols C)
        totalDIC, griddedTotalDIC = analysis_tools.calculate_total_plume_dic(griddedMeanDICConc, griddedPlumeVolume, calculateForSamples=False);
        totalDICSamples, griddedTotalDICSamples = analysis_tools.calculate_total_plume_dic(griddedMeanDICConcSamples, griddedPlumeVolumeSamples, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_total_dic"] = totalDIC;
        monthlyDF.loc[currentDatetime, "plume_total_dic_sd"] = np.std(totalDICSamples);
        
        
        #Calculate DIC outflow (moles, then converted to TgC)
        dicOutflow = analysis_tools.calculate_dic_outflow(monthlyDischarge.loc[currentDatetime]["monthly_discharge"], plumeVolume, totalDIC, outflowMonthlyTotalDischargeUncertainty=None);
        dicOutflow = analysis_tools.mol_to_TgC(dicOutflow); #convert from mols to TgC
        dicOutflowSamples = analysis_tools.calculate_dic_outflow(monthlyDischarge.loc[currentDatetime]["monthly_discharge"], plumeVolumeSamples, totalDICSamples, outflowMonthlyTotalDischargeUncertainty=monthlyDischarge.loc[currentDatetime]["monthly_discharge_sd"]);
        dicOutflowSamples = analysis_tools.mol_to_TgC(dicOutflowSamples); #convert from mols to TgC
        monthlyDF.loc[currentDatetime, "dic_outflow"] = dicOutflow;
        monthlyDF.loc[currentDatetime, "dic_outflow_sd"] = np.std(dicOutflowSamples);
        
        
        #Update output netCDF file for this time point
        if writeNetCDFOutput == True:
            griddedPlumeVolumeSD = np.nanstd(griddedPlumeVolumeSamples, axis=0);
            griddedTotalDICSD = np.nanstd(griddedTotalDICSamples, axis=0);
            update_gridded_time_point_netCDF(ncout, t, griddedPlumeVolume, griddedPlumeVolumeSD, griddedTotalDIC, griddedTotalDICSD, plumeMask);
        
    

    #Calculate inter-annual values for each month        
    interyearDF = analysis_tools.calculate_inter_year_monthly_means(monthlyDF);
    
    #calculate mean annual values
    annualDF = analysis_tools.calculate_annual_values(interyearDF);
    
    
    #Write time series to netCDF file
    if writeNetCDFOutput == True:
        write_timeseries_to_netCDF(ncout, monthlyDF["plume_volume"].values, monthlyDF["plume_volume_sd"].values,
                                   monthlyDF["plume_total_dic"].values, monthlyDF["plume_total_dic_sd"].values,
                                   monthlyDF["dic_outflow"].values, monthlyDF["dic_outflow_sd"].values,
                                   monthlyDF["river_discharge"].values, monthlyDF["river_discharge_sd"].values,
                                   interyearDF["dic_outflow"].values, interyearDF["dic_outflow_sd"].values);
        print("Written netCDF file output to: " + baseOutputPath);
        
        ncout.close();
    
    #Write csv files
    if writeCSVOutput == True:
        if path.exists(baseOutputPath) == False:
            makedirs(baseOutputPath);
            
        monthlyPath = path.join(baseOutputPath, "monthly_timeseries_"+region+".csv");
        monthlyDF.to_csv(monthlyPath, sep=",");
        
        interyearPath = path.join(baseOutputPath, "interyear_timeseries_"+region+".csv");
        interyearDF.to_csv(interyearPath, sep=",");
        
        annualPath = path.join(baseOutputPath, "annual_means_"+region+".csv");
        annualDF.to_csv(annualPath, sep=",");
        
        print("Written csv file output to:" + baseOutputPath);

    
    

    
    















