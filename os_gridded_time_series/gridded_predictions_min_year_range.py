#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:39:41 2020
A version of gridded predictions driver script which allows a minimum temporal range to be specified. All input data combinations which have shorter temporal ranges will be excluded from consideration for 'best algorithm/input data combination'.
@author: tom holding
"""


#from string import Template;
from os import path;
import os;
import pandas as pd;

import osoda_global_settings;
import osoda_gridded_predictions as gp;


if __name__ == "__main__":
    settings = osoda_global_settings.get_default_settings();
    lonRes = latRes = 1.0;
    years = settings["years"];
    
    overallBestAlgosOutputPath = path.join(settings["outputPathMetrics"], "overall_best_algos_min_years=8.csv");
    bestAlgoTable = pd.read_csv(overallBestAlgosOutputPath);
    
    #make gridded time series predictions for each region, using the best input combination and algorithms
    for region in settings["regions"]:
        #Find the best input combination and algorithm for DIC and AT
        atAlgoInfo = bestAlgoTable[(bestAlgoTable["region"]==region) & (bestAlgoTable["output_var"]=="AT")];
        if len(atAlgoInfo) != 1:
            raise ValueError("Error: 0 or more than 1 entries returned for the best AT algorithm in "+region);
        else:
            atAlgoInfo = atAlgoInfo.iloc[0]; #get row from dataframe with only one row
        
        dicAlgoInfo = bestAlgoTable[(bestAlgoTable["region"]==region) & (bestAlgoTable["output_var"]=="DIC")];
        if len(dicAlgoInfo) != 1:
            raise ValueError("Error: 0 or more than 1 entries returned for the best DIC algorithm in "+region);
        else:
            dicAlgoInfo = dicAlgoInfo.iloc[0]; #get row from dataframe with only one row
        
        #create instances of the algorithm functors by searching the available algorithms using the algorithm name
        if type(atAlgoInfo.algo_name) == str:
            atAlgo = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == atAlgoInfo["algo_name"]][0];
        else:
            atAlgo = None;
        if type(dicAlgoInfo.algo_name) == str:
            dicAlgo = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == dicAlgoInfo["algo_name"]][0];
        else:
            dicAlgo = None;

        
        #Create output file path
        griddedPredictionOutputPathAT = settings["griddedPredictionMinYearsOutputTemplate"].safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="AT");
        griddedPredictionOutputPathDIC = settings["griddedPredictionMinYearsOutputTemplate"].safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="DIC");
        if path.exists(path.dirname(griddedPredictionOutputPathAT)) == False:
            os.makedirs(path.dirname(griddedPredictionOutputPathAT));
        #calculate the gridded time series and write to file for this input / region combination
        gp.calculate_gridded_timeseries_driver(griddedPredictionOutputPathAT, griddedPredictionOutputPathDIC, atAlgo, atAlgoInfo, dicAlgo, dicAlgoInfo, settings, years, latRes, lonRes, verbose=True);



