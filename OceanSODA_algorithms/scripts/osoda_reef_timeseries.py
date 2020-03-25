#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:41:50 2020

Uses reef base 'reef location' data to produce summaries of number of reefs in each OSODA region.

@author: tom holding
"""

import pandas as pd;
import numpy as np;
from os import path;
import os;
from netCDF4 import Dataset;
from string import Template;
import datetime;
#import matplotlib.pyplot as plt;
#import matplotlib.ticker as plticker

import osoda_global_settings;
import utilities;

def convert_time(time, baseline = datetime.datetime(1980,1,1)):
    return np.array([baseline+datetime.timedelta(seconds=int(t)) for t in time]);

extractReefData = True;

settings = osoda_global_settings.get_default_settings();
reefLocationsPath = "../reefbase/ReefLocations.csv";
finalScoresTemplatePath = Template(path.join(settings["outputPathMetrics"], "${INPUTCOMBINATION}/${OUTPUTVAR}/${REGION}/final_scores.csv"));
griddedPredictionResolution = 1.0; #resolution of the gridded predictions
reefIndividualOutputPathTemplate = Template(path.join(settings["outputPathRoot"], "reef_outputs", "individual_${REGION}/reef_${REEFID}.csv"));
reefSummaryOutputPathTemplate = Template(path.join(settings["outputPathRoot"], "reef_outputs", "all_reef_${SUMMARYVAR}_summary_metrics.csv"));


#Sort reef locations by region. Returns a dictionary of region:reefsInLocationDF
def sort_reefs_by_region(settings, reefLocationsPath):
    reefsByRegion = {};
    #Load reef data
    reefsDF = pd.read_csv(reefLocationsPath, sep=",", encoding="ISO-8859-1");
    reefsDF.columns = [key.lower() for key in reefsDF.keys()];
    reefsDF = reefsDF[(np.isfinite(reefsDF["lon"])) & (np.isfinite(reefsDF["lat"]))];
    
    #Select reefs in each region
    regionMaskNC = Dataset(settings["regionMasksPath"], 'r');
    for region in settings["regions"]:
        reefsByRegion[region] = utilities.subset_from_mask(reefsDF, regionMaskNC, region);
    
    return reefsByRegion;

def get_grid_indices_from_latlon(lon, lat, resolution):
    ilat = int(lat / resolution) + 90;
    ilon = int(lon / resolution) + 180;
    return ilat, ilon;
    
#Returns a dictionary containing summary metrics for a particular variable in a reef time series
#colToSummarise: string column name in the reefTimeseriesDF dataframe
def calculate_reef_summary_metrics(reefTimeseriesDF, colToSummarise):
    reefData = reefTimeseriesDF[colToSummarise];
    
    #Check it isn't an all nan reef
    if np.all(np.isfinite(reefData)==False): #If all the values are nan fill the rest with nans and return
        return None;
    
    row = {};
    
    #mean, standard deviation, #number of time points (nt)
    row["mean"] = np.nanmean(reefData);
    row["median"] = np.nanmedian(reefData);
    row["sd"] = np.nanstd(reefData);
    row["nt"] = np.sum(np.isfinite(reefData)); #number of time points (nt)
    
    #max, min, range
    row["max"] = np.nanmax(reefData);
    row["min"] = np.nanmin(reefData);
    row["range"] = row["max"] - row["min"];
    
    #annual statistics (assumes whole years in monthly resolution)
    if len(reefData) % 12 != 0:
        raise ValueError("Must be monthly time points and whole years!");
    byYear = reefData.values.reshape((len(reefData)//12, 12));
    
    #mean annual maxmimum, mean annual minimum and their standard deviations
    row["annual_max_mean"] = np.nanmean(np.nanmax(byYear, axis=1));
    row["annual_max_sd"] = np.nanstd(np.nanmax(byYear, axis=1));
    row["annual_min_mean"] = np.nanmean(np.nanmin(byYear, axis=1));
    row["annual_min_sd"] = np.nanstd(np.nanmin(byYear, axis=1));
    
    #mean annual range
    row["annual_range_mean"] = np.nanmean(np.nanmax(byYear, axis=1) - np.nanmin(byYear, axis=1));
    row["annual_range_sd"] = np.nanstd(np.nanmax(byYear, axis=1) - np.nanmin(byYear, axis=1));
    
    return row;


### Sort reefs by region. Returns a dictionary of reef dataframes with region as key
reefsByRegion = sort_reefs_by_region(settings, reefLocationsPath);

#Extract ocean carbonate predictions for each reef
if extractReefData:
    #Create a dataframe to store the summary metrics for each reef
    summaryColNames = ["region", "algorithm", "reef_id", "mean", "median", "sd", "nt", "max", "min", "annual_max_mean", "annual_min_mean",
                       "annual_max_sd", "annual_min_sd", "annual_range_mean", "annual_range_sd"]; #Defines order of columns. The specific names must match those used in calculate_reef_summary_metrics
    reefSummaryMetricsDF_AT = pd.DataFrame(columns=summaryColNames);
    reefSummaryMetricsDF_DIC = pd.DataFrame(columns=summaryColNames);
    
    for region in settings["regions"]:
        if len(reefsByRegion[region]) == 0: #No reefs in this region, so move onto the next region
            continue;
        
        #load the gridded carbonate parameter predictions for this input/region combination
        griddedPredictionsPath = settings["griddedPredictionOutputTemplate"].safe_substitute(REGION=region, LATRES=griddedPredictionResolution, LONRES=griddedPredictionResolution, OUTPUTVAR="DIC");
        griddedPredictionsPathAT = settings["griddedPredictionOutputTemplate"].safe_substitute(REGION=region, LATRES=griddedPredictionResolution, LONRES=griddedPredictionResolution, OUTPUTVAR="AT");
        try:
            griddedPredictionNC = Dataset(griddedPredictionsPath, 'r');
            griddedPredictionNC_at = Dataset(griddedPredictionsPathAT, 'r');
        except FileNotFoundError: #This can occur if there was no best algorithm (e.g. because no matchup data for this region and inpu combination)
            continue; #Ignore the region
        
        #For each reef, extract a time series and output to file
        for r, reefRow in reefsByRegion[region].iterrows():
            print(region, "reef_"+str(r));
            
            #grid indices corresponding to the reef location
            ilat, ilon = get_grid_indices_from_latlon(reefRow["lon"], reefRow["lat"], griddedPredictionResolution);
            
            #Copy time series for the single grid point into a data frame
            reefDF = pd.DataFrame();
            reefDF["time_s_since_1980"] = griddedPredictionNC.variables["time"][:];
            for varName in list(griddedPredictionNC.variables.keys()) + ["AT_pred"]:
                if varName in ["time", "lat", "lon"]: #Depending on the algorithm, some inputs may not be present
                    continue;
                
                if varName == "AT_pred": #Special case, AT comes from a differentfile. This is awkward but avoids duplicate input layers
                    var = griddedPredictionNC_at.variables[varName][:, ilat, ilon]
                    var[var.mask] = np.nan;
                    reefDF[varName] = var;
                else: #All other variables come from DIC netCDF file
                    var = griddedPredictionNC.variables[varName][:, ilat, ilon];
                    var[var.mask] = np.nan; #Replace default missing value with nan
                    reefDF[varName] = var;
            
            #Write reef time series to csv file
            reefOutputPath = reefIndividualOutputPathTemplate.safe_substitute(REGION=region, REEFID=reefRow["id"]);
            if path.exists(path.dirname(reefOutputPath)) == False:
                os.makedirs(path.dirname(reefOutputPath))
            reefDF.to_csv(reefOutputPath, index=False, sep=",");
            
            ### Summary metrics for each reef
            #Construct summary data for AT
            summaryMetricsAT = calculate_reef_summary_metrics(reefDF, "AT_pred");
            if summaryMetricsAT is not None:
                summaryMetricsAT["reef_id"] = reefRow["id"];
                summaryMetricsAT["region"] = region;
                summaryMetricsAT["algorithm"] = griddedPredictionNC_at.getncattr("algorithmName");
            
                #append this reef's summary metrics to the end of the summary dataframe (with the columns in the prescribed order)
                reefSummaryMetricsDF_AT.loc[len(reefSummaryMetricsDF_AT)] = [summaryMetricsAT[key] for key in summaryColNames];
            else: #None, because no data was available to summarise for the selected variable, fill row with nan
                reefSummaryMetricsDF_AT.loc[len(reefSummaryMetricsDF_AT)] = [np.nan for key in summaryColNames];
            
            #Construct summary data for DIC
            summaryMetricsDIC = calculate_reef_summary_metrics(reefDF, "DIC_pred");
            if summaryMetricsDIC is not None:
                summaryMetricsDIC["reef_id"] = reefRow["id"];
                summaryMetricsDIC["region"] = region;
                summaryMetricsDIC["algorithm"] = griddedPredictionNC.getncattr("algorithmName");
            
                #append this reef's summary metrics to the end of the summary dataframe (with the columns in the prescribed order)
                reefSummaryMetricsDF_DIC.loc[len(reefSummaryMetricsDF_DIC)] = [summaryMetricsDIC[key] for key in summaryColNames];
            else: #None, because no data was available to summarise for the selected variable, fill row with nan
                reefSummaryMetricsDF_DIC.loc[len(reefSummaryMetricsDF_DIC)] = [np.nan for key in summaryColNames];
            
            
        #write reef summary metrics dataframe to file
        reefSummaryOutputPathAT = reefSummaryOutputPathTemplate.safe_substitute(SUMMARYVAR="AT");
        if path.exists(path.dirname(reefSummaryOutputPathAT)) == False:
            os.makedirs(path.dirname(reefSummaryOutputPathAT));
        reefSummaryMetricsDF_AT.to_csv(reefSummaryOutputPathAT);
        reefSummaryOutputPathDIC = reefSummaryOutputPathTemplate.safe_substitute(SUMMARYVAR="DIC");
        reefSummaryMetricsDF_DIC.to_csv(reefSummaryOutputPathDIC);

            


#if calculateMetaMetrics:
#    #meta metrics: summary of number of reefs fitting into different categories
#    #??? Number of reefs where the mean annual range > some threshold?
#    for outputVar in ["AT", "DIC"]:
#        for region in settings["regions"]:
#            try:
#                reefMetrics = pd.read_csv(reefTimeSeriesMetricsPath.safe_substitute(OUTPUTVAR=outputVar, REGION=region, sep=","));
#            except FileNotFoundError:
#                continue;
#            meanAnnualRange = reefMetrics["annual_range_mean"];
#            meanAnnualRange = meanAnnualRange[meanAnnualRange<10000];
#            plt.figure();
#            plt.hist(meanAnnualRange);
#            plt.xlabel(r"mean annual range in "+outputVar+" ($\mu$mol/kg)");
#            plt.ylabel("frequency");
#            plt.title(region+" "+outputVar);
#            plt.savefig(path.join(reefTimeSeriesPlotsPath, "hist_mean_annual_range_"+region+"_"+outputVar+".png"));
#    
#
#

















