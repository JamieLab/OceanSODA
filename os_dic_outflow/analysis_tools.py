#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import pyproj;
import numpy as np;
from datetime import datetime;
import pandas as pd;
from os import path, makedirs;
from netCDF4 import Dataset;
import calendar;
import matplotlib.pyplot as plt;


######################
### Simple utilities #
######################

#Convert from concentration (umol kg-1) to quantity (mol).
#def concentration_to_moles(data):
#    return data * 1029 * 1000000.0; #1m^3 seawater weights approx. 1029 kg.

#moles C to TgC conversion
def mol_to_TgC(valueInMols):
    return valueInMols * 12.0107 / (10.0**12); #convert from mols to TgC. 12.0107 is the mean molecular mass of carbon

#return an array of datatime objects for the start of each month for an inclusive date range
def create_date_month_array(startDate, endDateInclusive):
    dates = [];
    for year in range(startDate.year, endDateInclusive.year+1):
        for month in range(1, 13):
            newDatetime = datetime.strptime("01-%s-%s"%(format(month, "02d"), str(year)), "%d-%m-%Y");
            if (newDatetime >= startDate) and (newDatetime <= endDateInclusive):
                dates.append(newDatetime);
    return dates;



##########################################
### Read, process or generate input data #
##########################################

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

# mannually create this as newer version of pyproj isn't compatible with gt5 dependencies
# areas = calculate_grid_areas(1.0, 1.0);
# np.savetxt("grid_areas_1.0x1.0.csv", areas, delimiter=",");


#Extracts DIC and SSS data for a specific region and time point
def extract_data(carbonateParameters, regionMask, t):
    dic = carbonateParameters.variables["DIC"][t,:,:];
    dic[regionMask != 1] = np.nan;
    #dic_uncertainty = carbonateParameters.getncattr("algorithmRMSDe");
    dic_uncertainty = carbonateParameters.variables["DIC_Combined_uncertainty_dueto_RMSD_and_Bias"][t,:,:];
    dic_uncertainty[regionMask != 1] = np.nan;
    
    sss = carbonateParameters.variables["SSS"][t,:,:];
    sss[regionMask != 1] = np.nan;
    sss[sss.mask] = np.nan;
    sss_uncertainty = carbonateParameters.variables["SSS_uncertainty"][t,:,:];
    sss_uncertainty[regionMask != 1] = np.nan;
    sss_uncertainty[sss_uncertainty.mask] = np.nan;
    
    return dic, dic_uncertainty, sss, sss_uncertainty;


###############################################
### Components of the DIC outflow calculation #
###############################################

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
#griddedVolume in m^3
#griddedMeanDIC in umol kg-1
#outputs are in mols
def calculate_total_plume_dic(griddedMeanDIC, griddedVolume, calculateForSamples=False):
    griddedMass = griddedVolume*1020.0; #1020kg in 1 m^3 sea water. Results is in kg
    griddedTotalDIC_umols = griddedMass*griddedMeanDIC; #umol kg-1 / kg: Result is umol
    griddedTotalDIC_mols = griddedTotalDIC_umols / (10**6); #scale umol to mol
    
#    griddedTotalDIC = griddedMeanDIC*griddedVolume;
#    #import pdb; pdb.set_trace()
#    griddedTotalDIC = concentration_to_moles(griddedTotalDIC);
    
    if calculateForSamples==True:
        totalDIC = np.nansum(griddedTotalDIC_mols, axis=(1,2));
    else:
        totalDIC = np.nansum(griddedTotalDIC_mols);
    
    return totalDIC, griddedTotalDIC_mols;


#Calculates the total and gridded DIC estimated to originate from the river. This assumes conservative mixing with saliniyt.
#Uses a linear interpolation between 0 salinity (100% riverine origin) and 35 salinity (0% riverine origin).
#returns units in the same griddedTotalDIC (e.g. can also be used for concentrations)
def interpolate_riverine_dic(griddedTotalDIC, griddedMeanSSS, calculateForSamples=False):
    #TODO: Check this logic against conservative mixing model in Cooley et al: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2005JC002954
    proportionOriginatingRiver = (griddedMeanSSS / -35.0)+1.0;
    proportionOriginatingRiver[np.where(proportionOriginatingRiver<0.0)] = 0.0;
    proportionOriginatingRiver[np.where(proportionOriginatingRiver>1.0)] = 1.0; #sanity check, but should be impossible (barring errors originating from floating point representation)
    
    griddedTotalRiverineDIC = griddedTotalDIC * proportionOriginatingRiver;
    if calculateForSamples==True: #If the calculation is being carried out on lots of samples instead of a single matrix
        totalRiverineDIC = np.nansum(griddedTotalRiverineDIC, axis=(1,2));
    else:
        totalRiverineDIC = np.nansum(griddedTotalRiverineDIC);
    
    return totalRiverineDIC, griddedTotalRiverineDIC; #Returns DIC in the same units (e.g. mols)


#Calculates an estimate of the monthly DIC outflow from a river, given DIC plume content, river discharge, and plume volume.
#If outflowMonthlyTotalDischargeUncertainty is not None, samples are assumed with number of samples equal to the length of the first dimension
#dicharge in m^3 month-1
#volume in m^3
#DIC in mols
#output in mols
def calculate_dic_outflow(outflowMonthlyTotalDischarge, plumeVolume, totalPlumeDIC, outflowMonthlyTotalDischargeUncertainty=None):
    #discharge volume as a fraction of plume volume
    if outflowMonthlyTotalDischargeUncertainty is not None:
        sampledOutflowMonthlyTotalDischarge = np.random.normal(outflowMonthlyTotalDischarge, outflowMonthlyTotalDischargeUncertainty, size=plumeVolume.shape[0]);
    else:
        sampledOutflowMonthlyTotalDischarge = outflowMonthlyTotalDischarge;
    
    dischargeVolumeProportion = sampledOutflowMonthlyTotalDischarge / plumeVolume;
    dischargeDIC = totalPlumeDIC * dischargeVolumeProportion;
    
    return dischargeDIC;



################################
### Summary table calculations #
################################

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
    interyearDF["month_name"] = [calendar.month_abbr[i] for i in range(1, 13)];
    interyearDF.index = interyearDF["month"];
    
    for colName in monthlyDF.keys():
        if colName == "date":
            continue; #ignore date
        if colName[-3:] == "_sd": #standard deviations
            interyearDF[colName] = inter_year_sd(monthlyDF, colName);
        else: #values
            interyearDF[colName] = inter_year_mean(monthlyDF, colName);
    
    return interyearDF;


def format_annual_dataframe(dataDict):
    unitsDict = {"plume_volume": "m^3",
                 "plume_surface_area": "m^2",
                 "plume_mean_thickness": "m",
                 "plume_total_dic": "g C",
                 "dic_outflow": "Tg C",
                 "river_discharge": "m^3 year-1",
                 };
    
    names = dataDict.keys();
    parameters = [name for name in names if "_sd" not in name];
    values = [dataDict[param] for param in parameters];
    uncertainties = [dataDict[param+"_sd"] for param in parameters];
    units = [unitsDict[param] for param in parameters];
    
    df = pd.DataFrame();
    df["parameters"] = parameters;
    df["annual_mean"] = values;
    df["uncertainty"] = uncertainties;
    df["units"] = units;
    
    return df;

#returns a data frame containing mean annual values
def calculate_annual_values(interyearDF):
    def calc_sum_sd(interyearDF, colName):
        return np.sqrt(np.sum(interyearDF[colName]**2));
    
    def calc_mean_sd(interyearDF, colName):
        return calc_sum_sd(interyearDF, colName) / 12;
    
    columnsToSum = ["dic_outflow", "dic_outflow_sd", "river_discharge", "river_discharge_sd"];
    columnsToAverage = ["plume_total_dic", "plume_total_dic_sd", "plume_surface_area", "plume_surface_area_sd", "plume_mean_thickness", "plume_mean_thickness_sd", "plume_volume", "plume_volume_sd"];
    
    output = {};
    
    #Some columns sum, some average, to give annual means
    #handle them separately
    for colName in columnsToSum:
        if colName[-3:] == "_sd": #if a standard deviation column
            output[colName] = calc_sum_sd(interyearDF, colName);
        else:
            output[colName] = np.sum(interyearDF[colName]);
    
    for colName in columnsToAverage:
        if colName[-3:] == "_sd": #if a standard deviation column
            output[colName] = calc_mean_sd(interyearDF, colName);
        else:
            output[colName] = np.mean(interyearDF[colName]);
    
    #convert to df
    df = format_annual_dataframe(output);
    return df;


################
### outputting #
################

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
    var = ncout.createVariable("DIC", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["DIC"].units;
    var.long_name = carbonateParameterNC.variables["DIC"].long_name;
    var[:] = carbonateParameterNC.variables["DIC"][:]
    
    #Write data to netCDF file
    var = ncout.createVariable("SSS", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["SSS"].units;
    var.long_name = carbonateParameterNC.variables["SSS"].long_name;
    var[:] = carbonateParameterNC.variables["SSS"][:]
    
    #Write data to netCDF file
    var = ncout.createVariable("SSS_uncertainty", float, ("time", "lat", "lon"));
    var.units = carbonateParameterNC.variables["SSS_uncertainty"].units;
    var.long_name = carbonateParameterNC.variables["SSS_uncertainty"].long_name;
    var[:] = carbonateParameterNC.variables["SSS_uncertainty"][:]
    
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
    
    var = ncout.createVariable("plume_riverine_dic_total", float, ("time"));
    var.units = "mol";
    var.long_name = "Total DIC originating from the river, contained in the plume over time";
    
    var = ncout.createVariable("plume_riverine_dic_total_stddev", float, ("time"));
    var.units = "mol";
    var.long_name = "Standard deviation in the total DIC originating from the river, contained in the plume over time";
    
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
def write_timeseries_to_netCDF(ncout, plumeVolume, plumeVolumeUncertainty, totalPlumeDIC, totalPlumeDICUncertainty,
                               totalPlumeRiverineDIC, totalPlumeRiverineDICUnvertainty, dischargeDIC, dischargeDICUncertainty,
                               dischargeTotal, dischargeTotalUncertainty, dischargeDICBetweenYearMean, dischargeDICBetweenYearMeanUncertainty):
    ncout.variables["plume_volume_total"][:] = plumeVolume;
    ncout.variables["plume_volume_total_stddev"][:] = plumeVolumeUncertainty;
    ncout.variables["plume_dic_total"][:] = totalPlumeDIC;
    ncout.variables["plume_dic_total_stddev"][:] = totalPlumeDICUncertainty;
    ncout.variables["plume_riverine_dic_total"][:] = totalPlumeRiverineDIC;
    ncout.variables["plume_riverine_dic_total_stddev"][:] = totalPlumeRiverineDICUnvertainty;
    
    ncout.variables["discharge_dic"][:] = dischargeDIC;
    ncout.variables["discharge_dic_stddev"][:] = dischargeDICUncertainty;
    
    ncout.variables["discharge"][:] = dischargeTotal;
    ncout.variables["discharge_stddev"][:] = dischargeTotalUncertainty;
    
    ncout.variables["discharge_dic_between_year_mean"][:] = dischargeDICBetweenYearMean;
    ncout.variables["discharge_dic_between_year_uncertainty"][:] = dischargeDICBetweenYearMeanUncertainty;




########################
### plotting functions #
########################
#plots colName and colName+"_sd"
def plot_with_uncertainty(df, colName, xlabel, ylabel, x=None, title=None, outputPath=None, closeAfter=False):
    if x is None:
        x = range(0, len(df));

    value = np.array(df[colName]);
    uncertainty = np.array(df[colName+"_sd"]);
    
    plt.figure();
    plt.fill_between(x, value+uncertainty, value-uncertainty, alpha=0.75);
    plt.plot(x, value, 'k', alpha=0.5);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    if title is not None:
        plt.title(title);
    plt.tight_layout();
    
    #save fig
    if outputPath is not None:
        if path.exists(path.dirname(outputPath)) == False:
            makedirs(path.dirname(outputPath));
        plt.savefig(outputPath);
    
    if closeAfter == True:
        plt.close();



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