#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;
from string import Template;
from os import path, makedirs;
from datetime import datetime, timedelta;
from preprocess_data_sources import load_discharge_data;

import analysis_tools;

import sys;
if path.abspath("../../") not in sys.path:
    sys.path.append(path.abspath("../../"));
    sys.path.append(path.abspath("../../OceanSODA_algorithms/scripts"));
#import OceanSODA_algorithms.scripts.utilities as utilities;
import osoda_global_settings;


##### Global settings
useMinRangeCarbonateParameters = True; #Use gridded prediction data derived from algorithms and input data sets which meet a minimum temporal range requirement (e.g. 8 year range)
if useMinRangeCarbonateParameters == True:
    CARBONATE_VERSION_STR = "_min_range"
else:
    CARBONATE_VERSION_STR = "";
writeNetCDFOutput = True;
writeCSVOutput = True;
baseOutputTemplate = Template("../output/${REGION}"+CARBONATE_VERSION_STR+"/");
baseDataPath = "../data/";
makePlots = True;
closePlots = True; #Close plots automatically?
basePlotTemplate = Template("../plots/${REGION}"+CARBONATE_VERSION_STR+"/");
lonRes = latRes = 1.0;
plumeSalinityThreshold = 35.0; #Below this salinity is considered plume
numSamples = 100; #How many SSS and DIC samples to use
verbose = True;
regions = ["oceansoda_amazon_plume", "oceansoda_congo", "oceansoda_mississippi", "oceansoda_st_lawrence"];

carbonateParametersTemplate = Template("../../OceanSODA_algorithms/output/gridded_predictions/gridded_${REGION}_1.0x1.0_${OUTPUTVAR}.nc");
carbonateParametersMinRangeTemplate = Template("../../OceanSODA_algorithms/output/gridded_predictions_min_year_range/gridded_${REGION}_1.0x1.0_${OUTPUTVAR}.nc"); #Gridded predictions using the algorithm/input data sets which have a minimum range
regionMaskPath = "../../OceanSODA_algorithms/region_masks/osoda_region_masks_v2.nc";

settings = osoda_global_settings.get_default_settings();
gridAreas = np.genfromtxt(path.join(baseDataPath, "grid_areas_"+str(lonRes)+"x"+str(latRes)+".csv"), delimiter=","); #manually read as newer pyproj isn't compatible with qt5 dependencies


zeroPlumes=[];

skipped = [];
for region in regions:
    #Setup region specific output directories
    baseOutputPath = baseOutputTemplate.safe_substitute(REGION=region);
    if path.exists(baseOutputPath) == False:
            makedirs(baseOutputPath);
    basePlotPath = basePlotTemplate.safe_substitute(REGION=region);
    if path.exists(basePlotPath) == False:
            makedirs(basePlotPath);
    
    #Read input data sets
    #Carbonate parameters
    if useMinRangeCarbonateParameters == True: #Use the predictions that cover a minimum temporal range, even if it isn't the absolute best
        carbonatePath = carbonateParametersMinRangeTemplate.safe_substitute(REGION=region, OUTPUTVAR="DIC")
    else: #useMinRangeCarbonateParameters is False, so use overall best performing algorithm / input combination
        carbonatePath = carbonateParametersTemplate.safe_substitute(REGION=region, OUTPUTVAR="DIC");
    try:
        carbonateParameters = Dataset(carbonatePath);
        regionMaskNC = Dataset(regionMaskPath, 'r');
    except FileNotFoundError:
        print("Skipping", region, "because no gridded carbonate parameter file found for it.");
        skipped.append(region, carbonatePath);
        continue;
    
    #mask data
    regionMask = regionMaskNC.variables[region][:];
    #regionMask3D = np.broadcast_to(amazonRegionMask, (len(carbonateParameters.variables["time"]), len(carbonateParameters.variables["lat"]), len(carbonateParameters.variables["lon"])));
    
    #Time data
    numTimePoints = len(carbonateParameters.variables["time"]);
    carbonateParamsTimeRange = (datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][0])), datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][-1])));
    carbonateDates = analysis_tools.create_date_month_array(carbonateParamsTimeRange[0], carbonateParamsTimeRange[1]);
    
    #River discharge data
    monthlyDischarge = load_discharge_data(carbonateDates, region);

    
    ### Create netCDF file to store results - this will be updated with each iteration
    if writeNetCDFOutput == True:
        outputPathNetCDF = path.join(baseOutputPath, "dic_outflow_"+region+".nc");
        ncout = analysis_tools.create_netCDF_file(outputPathNetCDF, carbonateParameters, numSamples, regionMask, gridAreas);
    
    
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
        if verbose:
            print(region, t+1, "of", numTimePoints, "("+str(currentDatetime)+")");
        
        #extract data for this time point #SSS in PSU, DIC in umol kg-1
        dic, dic_err, sss, sss_err = analysis_tools.extract_data(carbonateParameters, regionMask, t=t);
        
        #Don't perform computations if everything is nan
        if np.all(np.isfinite(dic)==False):
            continue;
        
        #Generate DIC and SSS samplestime point
        sssSamples = np.random.normal(sss, sss_err, size=(numSamples,)+sss.shape);
        dicSamples = np.random.normal(dic, dic_err, size=(numSamples,)+dic.shape);
        
        #dicV = dic[np.where(np.isnan(dic)==False)]
        #dic_errV = dic_err[np.where(np.isnan(dic_err)==False)]
        
        #Calculate plume masks (1==inside plume)
        plumeMask = analysis_tools.calculate_plume_mask(sss, plumeSalinityThreshold=plumeSalinityThreshold);
        plumeMaskSamples = analysis_tools.calculate_plume_mask(sssSamples, plumeSalinityThreshold=plumeSalinityThreshold);
        
        if np.sum(plumeMask) == 0:
            zeroPlumes.append(carbonateDates[t]);
        
        #Calculate plume surface area (m^2)
        surfaceArea, griddedSurfaceArea = analysis_tools.calculate_plume_surface_area(plumeMask, gridAreas); #surfaceArea in m^2
        surfaceAreaSamples, griddedSurfaceAreaSamples = analysis_tools.calculate_plume_surface_area(plumeMaskSamples, gridAreas, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_surface_area"] = surfaceArea;
        monthlyDF.loc[currentDatetime, "plume_surface_area_sd"] = np.std(surfaceAreaSamples);
        
        
        #Calculate plume thickness(m)
        meanPlumeThickness, griddedPlumeThickness = analysis_tools.plume_thickness_coles2013(sss, plumeMask, useModelled=False, sampleUncertainty=False); #thickness in metres (m)
        meanPlumeThicknessSamples, griddedPlumeThicknessSamples = analysis_tools.plume_thickness_coles2013(sssSamples, plumeMaskSamples, useModelled=False, sampleUncertainty=True);
        monthlyDF.loc[currentDatetime, "plume_mean_thickness"] = meanPlumeThickness;
        monthlyDF.loc[currentDatetime, "plume_mean_thickness_sd"] = np.std(meanPlumeThicknessSamples);
        
        
        #Calculate plume volume (m^3)
        plumeVolume, griddedPlumeVolume = analysis_tools.calculate_plume_volume(griddedSurfaceArea, griddedPlumeThickness, calculateForSamples=False); #volumes in m^3
        plumeVolumeSamples, griddedPlumeVolumeSamples = analysis_tools.calculate_plume_volume(griddedSurfaceAreaSamples, griddedPlumeThicknessSamples, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_volume"] = plumeVolume;
        monthlyDF.loc[currentDatetime, "plume_volume_sd"] = np.std(plumeVolumeSamples);
        
        
        #Calculate gridded mean plume DIC, in the process also calculate meanplume salinity (umol kg-1)
        #This uses relationship between mean plume salinity and sss in Hu 2004 (PS = 0.881*SSS+4.352, see section 3.4 of https://doi.org/10.1016/j.dsr2.2004.04.001)
        #and assumes DIC is perfectly conservative with salinity
        #I.e. the proportion mean plume salinity to surface salinity is assumed to equal mean plume DIC / surface DIC, thus rearranging gives:
        #   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
        #mean plume DIC given in umol kg-1, mean plume SSS given in PSU
        griddedMeanDICConc, griddedMeanSSS = analysis_tools.calculate_mean_dic_sss(sss, dic, plumeMask, interceptUncertaintyRatio=0.0, slopeUncertaintyRatio=0.0, calculateForSamples=False); #DIC in umol kg-1, salinity in PSU
        griddedMeanDICConcSamples, griddedMeanSSSSamples = analysis_tools.calculate_mean_dic_sss(sssSamples, dicSamples, plumeMaskSamples, interceptUncertaintyRatio=0.1, slopeUncertaintyRatio=0.1, calculateForSamples=True);
        
        
        #calculate total plume DIC (mols C)
        totalDIC, griddedTotalDIC = analysis_tools.calculate_total_plume_dic(griddedMeanDICConc, griddedPlumeVolume, calculateForSamples=False); #DIC in mol
        totalDICSamples, griddedTotalDICSamples = analysis_tools.calculate_total_plume_dic(griddedMeanDICConcSamples, griddedPlumeVolumeSamples, calculateForSamples=True);
        monthlyDF.loc[currentDatetime, "plume_total_dic"] = totalDIC;
        monthlyDF.loc[currentDatetime, "plume_total_dic_sd"] = np.std(totalDICSamples);
        
        
        #Calculate DIC outflow (moles, then converted to TgC)
        dicOutflow = analysis_tools.calculate_dic_outflow(monthlyDischarge.loc[currentDatetime]["monthly_discharge"], plumeVolume, totalDIC, outflowMonthlyTotalDischargeUncertainty=None); #monthly DIC in mols
        dicOutflow = analysis_tools.mol_to_TgC(dicOutflow); #convert from mols to TgC
        dicOutflowSamples = analysis_tools.calculate_dic_outflow(monthlyDischarge.loc[currentDatetime]["monthly_discharge"], plumeVolumeSamples, totalDICSamples, outflowMonthlyTotalDischargeUncertainty=monthlyDischarge.loc[currentDatetime]["monthly_discharge_sd"]);
        dicOutflowSamples = analysis_tools.mol_to_TgC(dicOutflowSamples); #convert from mols to TgC
        monthlyDF.loc[currentDatetime, "dic_outflow"] = dicOutflow;
        monthlyDF.loc[currentDatetime, "dic_outflow_sd"] = np.std(dicOutflowSamples);
        
        
        #Update output netCDF file for this time point
        if writeNetCDFOutput == True:
            griddedPlumeVolumeSD = np.nanstd(griddedPlumeVolumeSamples, axis=0);
            griddedTotalDICSD = np.nanstd(griddedTotalDICSamples, axis=0);
            analysis_tools.update_gridded_time_point_netCDF(ncout, t, griddedPlumeVolume, griddedPlumeVolumeSD, griddedTotalDIC, griddedTotalDICSD, plumeMask);

    

    #Calculate inter-annual values for each month        
    interyearDF = analysis_tools.calculate_inter_year_monthly_means(monthlyDF);
    
    #calculate mean annual values
    annualDF = analysis_tools.calculate_annual_values(interyearDF);
    
    
    #Write time series to netCDF file
    if writeNetCDFOutput == True:
        analysis_tools.write_timeseries_to_netCDF(ncout, monthlyDF["plume_volume"].values, monthlyDF["plume_volume_sd"].values,
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
        annualDF.insert(0, "region", region);
        annualDF.to_csv(annualPath, sep=",");
        
        print("Written csv file output to:" + baseOutputPath);

    
    
    #plotting
    if makePlots == True:
        #Visualise plume volume area over time
        outputPath = path.join(basePlotPath, "monthly_plume_volume.pdf");
        analysis_tools.plot_with_uncertainty(monthlyDF, "plume_volume", "date", r"plume volume ($m^3$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
        
        #Visualise total plume DIC content over time
        outputPath = path.join(basePlotPath, "monthly_plume_dic.pdf");
        analysis_tools.plot_with_uncertainty(monthlyDF, "plume_total_dic", "date", r"DIC content of plume ($mols C$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
        
        #Visualise DIC discharge over time
        outputPath = path.join(basePlotPath, "monthly_dic_discharge.pdf");
        analysis_tools.plot_with_uncertainty(monthlyDF, "dic_outflow", "date", r"riverine DIC discharge ($TgC$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
        
        #Visualise mean (inter-year) monthly plume volume
        outputPath = path.join(basePlotPath, "interyear_plume_volume.pdf");
        analysis_tools.plot_with_uncertainty(interyearDF, "plume_volume", "", r"mean plume volume ($m^3$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
        
        #Visualise mean (inter-year) monthly total plume DIC
        outputPath = path.join(basePlotPath, "interyear_plume_dic.pdf");
        analysis_tools.plot_with_uncertainty(interyearDF, "plume_total_dic", "", r"mean DIC content of plume ($mols C$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
        
        #Visualise mean (inter-year) monthly DIC discharge
        outputPath = path.join(basePlotPath, "interyear_dic_discharge.pdf");
        analysis_tools.plot_with_uncertainty(interyearDF, "dic_outflow", "", r"mean riverine DIC discharge ($TgC$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
        

if len(skipped) != 0:
    print("Some regions were skipped due to missing gridded carbonate parameter files:");
    for tup in skipped:
        print("\t"+tup[0]+":", tup[1]);


#from netCDF4 import Dataset;
#import matplotlib.pyplot as plt;
#import numpy as np;
#nc = Dataset("/home/verwirrt/Projects/Work/20190816_OceanSODA/OceanSODA_algorithms/output/gridded_predictions/gridded_oceansoda_amazon_plume_1.0x1.0_DIC.nc", 'r');
#
#
#dic = nc.variables["DIC_pred"][:];
#dic.shape
#means = np.nanmean(dic, axis=(1,2));
#means.shape
#
#plt.figure(); plt.plot(means);
#










