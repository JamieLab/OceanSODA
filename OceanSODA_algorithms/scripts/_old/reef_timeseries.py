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
import matplotlib.pyplot as plt;
import matplotlib.ticker as plticker

import osoda_global_settings;
from utilities import subset_from_mask;

def convert_time(time, baseline = datetime.datetime(1980,1,1)):
    return np.array([baseline+datetime.timedelta(seconds=int(t)) for t in time]);

extractReefData = True;
makePlots = False;
calculateReefMetrics = False;
calculateMetaMetrics = False;

settings = osoda_global_settings.get_default_settings();
regionMaskPath = settings["regionMasksPath"];#path.join("../region_masks/osoda_region_masks_v2.nc");
reefLocationsPath = "../reefbase/ReefLocations.csv";

finalScoresTemplatePath = Template(path.join(settings["outputPathMetrics"], "${OUTPUTVAR}/${REGION}/final_scores.csv"));
predictedTimeSeriesTemplatePath = Template("../output/gridded_predictions/${OUTPUTVAR}/gridded_${ALGONAME}_1.0x1.0.nc");
reefTimeSeriesByRegionPath = Template("../output/reef_time_series/reef_time_series_${REGION}_${OUTPUTVAR}.csv"); #output
reefTimeSeriesByReefPath = Template("../output/reef_time_series/individual/reef_${ID}_${REGION}.csv"); #output
reefTimeSeriesPlotsPath = "../output/reef_time_series/plots/"; #output plots
reefTimeSeriesMetricsPath = Template("../output/reef_time_series/metrics/reef_time_series_metrics_${REGION}_${OUTPUTVAR}.csv"); #output metrics

predictedResolution = 1.0; #resolution of the gridded predictions


#Load reef position data
if extractReefData:
    bestAlgoInfo = [];
    reefsDF = pd.read_csv(reefLocationsPath, sep=",", encoding="ISO-8859-1");
    reefsDF.columns = [key.lower() for key in reefsDF.keys()];
    reefsDF = reefsDF[(np.isfinite(reefsDF["lon"])) & (np.isfinite(reefsDF["lat"]))]
    
    regionMaskNC = Dataset(regionMaskPath, 'r');
    for region in settings["regions"]:
        reefsInRegion = subset_from_mask(reefsDF, regionMaskNC, region);
        
        if len(reefsInRegion) == 0:
            print("No reefs in region bounds for "+region+".");
            continue;
        
        #Find best AT and DIC algorithms for this region, extracted predicted data for each reef
        ### First for AT
        finalScoresAT = pd.read_csv(finalScoresTemplatePath.safe_substitute(OUTPUTVAR="AT", REGION=region));
        try:
            ibestAlgoAT = np.nanargmin(finalScoresAT["final_wrmsd"]);
            bestAlgoNameAT = finalScoresAT["algorithm"][ibestAlgoAT];
            print("Best algorithm: AT", region, bestAlgoNameAT);
            bestAlgoInfo.append(("AT", region, bestAlgoNameAT, finalScoresAT["final_wrmsd"][ibestAlgoAT]));
            
            #Read predicted time series data for AT and DIC
            predictedAT = Dataset(predictedTimeSeriesTemplatePath.safe_substitute(OUTPUTVAR="AT", ALGONAME=bestAlgoNameAT), 'r');
            
            #Extract reef positions
            reefTimeSeriesAT = pd.DataFrame();
            reefTimeSeriesAT["time"] = convert_time(predictedAT.variables["time"][:]);
            
            #Debug plot, why is there no data for some reefs?
            if False:
                plt.figure();
                plt.imshow(predictedAT.variables["AT_pred"][0,:,:]);
                convLat = (reefsInRegion["lat"] / predictedResolution) + 90;
                convLon = (reefsInRegion["lon"] / predictedResolution) + 180;
                plt.scatter(convLon, convLat, color='r');
                plt.ylim(0, 180);
            
            
            for r, row in reefsInRegion.iterrows():
                print("AT", region, r, "of", len(reefsInRegion));
                ilat = int(row["lat"] / predictedResolution) + 90;
                ilon = int(row["lon"] / predictedResolution) + 180;
                
                pointTimeSeriesAT = predictedAT.variables["AT_pred"][:, ilat, ilon];
                pointTimeSeriesAT[pointTimeSeriesAT.mask] = np.nan;
                reefTimeSeriesAT[str(row["id"])+"_"+bestAlgoNameAT] = pointTimeSeriesAT;
            
            #Write to file / plot
            outputPathAT = reefTimeSeriesByRegionPath.safe_substitute(OUTPUTVAR="AT", REGION=region);
            if path.exists(path.dirname(outputPathAT)) == False:
                os.makedirs(path.dirname(outputPathAT))
            reefTimeSeriesAT.to_csv(outputPathAT, index=False);
            
        except ValueError:
            print("*** All NaN encountered in finalScores.csv wRMSD row for", region, "AT");
            print("    \tThis means no pairwise weighted metrics for this region could be calculated (e.g. because there were no spatially overlapping algorithms or no algorithms reported their RMSD.", region);
        
        ### Repeat for DIC
        finalScoresDIC = pd.read_csv(finalScoresTemplatePath.safe_substitute(OUTPUTVAR="DIC", REGION=region));
        
        try:
            ibestAlgoDIC = np.nanargmin(finalScoresDIC["final_wrmsd"]);
            bestAlgoNameDIC = finalScoresDIC["algorithm"][ibestAlgoDIC];
            print("Best algorithm: DIC", region, bestAlgoNameDIC);
            bestAlgoInfo.append(("DIC", region, bestAlgoNameDIC, finalScoresDIC["final_wrmsd"][ibestAlgoDIC]));

            #Read predicted time series data for AT and DIC
            predictedDIC = Dataset(predictedTimeSeriesTemplatePath.safe_substitute(OUTPUTVAR="DIC", ALGONAME=bestAlgoNameDIC), 'r');
            
            #Extract reef positions
            reefTimeSeriesDIC = pd.DataFrame();
            reefTimeSeriesDIC["time"] = convert_time(predictedDIC.variables["time"][:]);
            
            #Debug plot, why is there no data for some reefs?
            if False:
                plt.figure();
                plt.imshow(predictedDIC.variables["DIC_pred"][0,:,:]);
                convLon = (reefsInRegion["lon"] / predictedResolution) + 180;
                convLat = (reefsInRegion["lat"] / predictedResolution) + 90;
                plt.scatter(convLon, convLat, color='r');
                plt.ylim(0, 180);
            
            for r, row in reefsInRegion.iterrows():
                print("DIC", region, r, "of", len(reefsInRegion));
                ilat = int(row["lat"] / predictedResolution) + 90;
                ilon = int(row["lon"] / predictedResolution) + 180;
                
                pointTimeSeriesDIC = predictedDIC.variables["DIC_pred"][:, ilat, ilon];
                pointTimeSeriesDIC[pointTimeSeriesDIC.mask] = np.nan;
                reefTimeSeriesDIC[str(row["id"])+"_"+bestAlgoNameDIC] = pointTimeSeriesDIC;
                
                if np.any(np.isfinite(predictedDIC.variables["DIC_pred"][:, ilat, ilon])):
                    print(row["id"]);
            
            #Write to file / plot
            outputPathDIC = reefTimeSeriesByRegionPath.safe_substitute(OUTPUTVAR="DIC", REGION=region);
            if path.exists(path.dirname(outputPathDIC)) == False:
                os.makedirs(path.dirname(outputPathDIC))
            reefTimeSeriesDIC.to_csv(outputPathDIC, index=False);
            
        except ValueError:
            print("*** All NaN encountered in finalScores.csv wRMSD row for", region, "DIC");
            print("    \tThis means no pairwise weighted metrics for this region could be calculated (e.g. because there were no spatially overlapping algorithms or no algorithms reported their RMSD.", region);
        
        ### Individual reef data
        reefTimeSeries = pd.DataFrame();
        reefTimeSeries["time"] = convert_time(predictedAT.variables["time"][:]);
        
        for r, row in reefsInRegion.iterrows():
            print("AT", region, r, "of", len(reefsInRegion));
            ilat = int(row["lat"] / predictedResolution) + 90;
            ilon = int(row["lon"] / predictedResolution) + 180;
            
            pointTimeSeriesSSS = predictedAT.variables["SSS"][:, ilat, ilon];
            pointTimeSeriesSSS[pointTimeSeriesSSS.mask] = np.nan;
            reefTimeSeries["SSS"] = pointTimeSeriesSSS;
            
            pointTimeSeriesSST = predictedAT.variables["SST"][:, ilat, ilon];
            pointTimeSeriesSST[pointTimeSeriesSST.mask] = np.nan;
            reefTimeSeries["SST"] = pointTimeSeriesSST;
            
            pointTimeSeriesAT = predictedAT.variables["AT_pred"][:, ilat, ilon];
            pointTimeSeriesAT[pointTimeSeriesAT.mask] = np.nan;
            reefTimeSeries["AT_pred"] = pointTimeSeriesAT;
            
            pointTimeSeriesDIC = predictedDIC.variables["DIC_pred"][:, ilat, ilon];
            pointTimeSeriesDIC[pointTimeSeriesDIC.mask] = np.nan;
            reefTimeSeries["DIC_pred"] = pointTimeSeriesDIC;
            
            outputPathIndividual = reefTimeSeriesByReefPath.safe_substitute(ID=str(row["id"]), REGION=region);
            if path.exists(path.dirname(outputPathIndividual)) == False:
                os.makedirs(path.dirname(outputPathIndividual))
            reefTimeSeries.to_csv(outputPathIndividual, index=False);
        
    
    for line in bestAlgoInfo:
        print(line);



#### Plotting
if makePlots:
    if path.exists(reefTimeSeriesPlotsPath) == False:
        os.makedirs(reefTimeSeriesPlotsPath);
    
    for outputVar in ["AT", "DIC"]:
        for region in settings["regions"]:
            #Read data
            try:
                df = pd.read_csv(reefTimeSeriesByRegionPath.safe_substitute(OUTPUTVAR=outputVar, REGION=region, sep=","));
            except FileNotFoundError:
                print("Reef time series .csv file not found for", outputVar, region, "so this will be skipped.");
                continue;
            
            #Extract the algorithm name (i.e. the best found algorithm that was used to extract the time series predictions)
            selectedAlgorithmName = "_".join(df.keys()[1].split("_")[1:]);
            
            #do the plotting
            data = df.values; #Convert to numpy array for ease of plotting
            
            ############### TMP remove very high values. Need to investigate why this is. ######################
            #data[:,1:][data[:,1:]>10000] = np.nan;
            #print("Warning, very large values have been removed. Check why this is...");
            ############### !!!!!!!!!!!!!!!!!!!!!
            
            
            fig, ax = plt.subplots(1,1, figsize=(12, 7));
            plt.plot(data[:,0], data[:,1:]);
            loc = plticker.MultipleLocator(base=24.0); #ticks every two years
            ax.xaxis.set_major_locator(loc);
            plt.xlabel("time");
            plt.ylabel(outputVar+r" ($\mu$mol/kg)");
            plt.title(region+" using "+selectedAlgorithmName);
            plt.savefig(path.join(reefTimeSeriesPlotsPath, "reef_time_series_"+region+"_"+outputVar+".png"));
            


if calculateReefMetrics:
    #for each output variable and regions
    for outputVar in ["AT", "DIC"]:
        for region in settings["regions"]:
            #create empty dataframe to store summary stats in
            columnNames = ["region", "algo_name", "reef_id", "mean", "sd", "nt", "max", "min", "annual_max_mean", "annual_min_mean",
                       "annual_max_sd", "annual_min_sd", "annual_range_mean", "annual_range_sd"];
            reefMetrics = pd.DataFrame(columns=columnNames);
            
            #read reef data
            try:
                reefTimeSeries = pd.read_csv(reefTimeSeriesByRegionPath.safe_substitute(OUTPUTVAR=outputVar, REGION=region, sep=","));
            except FileNotFoundError:
                print("Reef time series .csv file not found for", outputVar, region, "so this will be skipped.");
                continue;
            
            time = reefTimeSeries["time"];
            for reefCol in reefTimeSeries.columns[1:]: #1 to skip time
                reefData = reefTimeSeries[reefCol];
                row = {};
                
                #basic identifying info
                algoIdTokens = reefCol.split("_");
                algoName = "_".join(algoIdTokens[1:]);
                reefId = int(algoIdTokens[0]);
                row["region"] = region;
                row["algo_name"] = algoName;
                row["reef_id"] = reefId;
                
                #Check it isn't an all nan reef
                if np.all(np.isfinite(reefData)==False): #If all the values are nan fill the rest with nans
                    reefMetrics.loc[len(reefMetrics)] = [row["region"], row["algo_name"], row["reef_id"]]+[np.nan]*(len(columnNames)-3);
                    continue;
                
                #mean, standard deviation, #number of time points (nt)
                row["mean"] = np.nanmean(reefData);
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
                
                #append to data frame
                reefMetrics.loc[len(reefMetrics)] = [row[key] for key in columnNames];
            
            #write reef metrics to file
            reefMetricsOutputPath = reefTimeSeriesMetricsPath.safe_substitute(REGION=region, OUTPUTVAR=outputVar);
            if path.exists(path.dirname(reefMetricsOutputPath)) == False:
                os.makedirs(path.dirname(reefMetricsOutputPath));
            reefMetrics.to_csv(reefMetricsOutputPath);


if calculateMetaMetrics:
    #meta metrics: summary of number of reefs fitting into different categories
    #??? Number of reefs where the mean annual range > some threshold?
    for outputVar in ["AT", "DIC"]:
        for region in settings["regions"]:
            try:
                reefMetrics = pd.read_csv(reefTimeSeriesMetricsPath.safe_substitute(OUTPUTVAR=outputVar, REGION=region, sep=","));
            except FileNotFoundError:
                continue;
            meanAnnualRange = reefMetrics["annual_range_mean"];
            meanAnnualRange = meanAnnualRange[meanAnnualRange<10000];
            plt.figure();
            plt.hist(meanAnnualRange);
            plt.xlabel(r"mean annual range in "+outputVar+" ($\mu$mol/kg)");
            plt.ylabel("frequency");
            plt.title(region+" "+outputVar);
            plt.savefig(path.join(reefTimeSeriesPlotsPath, "hist_mean_annual_range_"+region+"_"+outputVar+".png"));
    



















