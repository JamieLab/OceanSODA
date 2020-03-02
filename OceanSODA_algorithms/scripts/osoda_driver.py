#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:36:57 2019

@author: tom holding
"""

from os import makedirs;
from os import path;
import numpy as np;
#import pandas as pd;
from netCDF4 import Dataset;
import pickle;
import logging;
import datetime;

import metrics;
import utilities;
import osoda_global_settings;

#tmp flags
reloadData = True;

#Setup logger
logger = logging.getLogger('osoda_driver');
logger.setLevel(logging.DEBUG);
loggerFileHandle = logging.FileHandler('logs/osoda_driver.log', mode='w');
loggerFileHandle.setLevel(logging.DEBUG);
logger.addHandler(loggerFileHandle);

logger.info("Started execution at: "+str(datetime.datetime.now()));


####################
# Define constants #
####################
settings = osoda_global_settings.get_default_settings(); #TODO: Hannah - Change this to call your settings function
outputVariableSelector = "AT"; #TODO: set 'AT' or 'DIC' depending on whether you want to use AT or DIC algorithms

#################
# Read datasets #
#################
if reloadData:
    matchupData = utilities.load_matchup_to_dataframe(settings); #each year concatinated is concatinated to create a single dataframe
    
#    if False: #Debug plots
#        import matplotlib.pyplot as plt;
#        
#        at = matchupData[np.isfinite(matchupData["AT"])];
#        plt.figure();
#        plt.scatter(at["lon"], at["lat"], c=at["AT"]);
#        plt.title("AT");
#        
#        dic = matchupData[np.isfinite(matchupData["DIC"])];
#        plt.figure();
#        plt.scatter(at["lon"], at["lat"], c=at["DIC"])
#        plt.title("DIC");


##################
# Run algorithms #
##################
for region in settings["regions"]:
    print("Running for region:", region);
    logger.info("Beginning new region: "+region);
    algorithmFunctorList = [];
    modelOutputList = [];
    matchupRowsUsedList = [];
    #for AlgorithmFunctor in settings["dic_algorithms"]: #+settings["at_algorithms"]:
    for AlgorithmFunctor in settings[outputVariableSelector.lower()+"_algorithms"]:
        outputVariable = AlgorithmFunctor.output_name(); #The common name of the output variable (i.e. 'DIC' or 'AT')
        algorithm = AlgorithmFunctor(settings); #Create an instance of the algorithm functor, this is the handle used to access the algorithm
        
        ####Subset matchup dataset based on region
        regionMaskNC = Dataset(settings["regionMasksPath"], 'r');
        subsetData = utilities.subset_from_mask(matchupData, regionMaskNC, region);
        
        ####Subset based on depth and distance to coast.
        if settings["depthMaskPath"] != None:
            depthMaskNC = Dataset(settings["depthMaskPath"], 'r');
            subsetData = utilities.subset_from_mask(subsetData, depthMaskNC, settings["depthMaskVar"]);
        if settings["distToCoastMaskPath"] != None:
            distToCoastMaskNC = Dataset(settings["distToCoastMaskPath"], 'r');
            subsetData = utilities.subset_from_mask(subsetData, distToCoastMaskNC, settings["distToCoastMaskVar"]);
        
        ####Make predictions
        try:
            modelOutput, dataUsed = algorithm(subsetData);
            logger.info("Output calculated (region:"+region+", algo: "+algorithm.__class__.__name__+")");
        except ValueError as e: #After subsetting there's no data left to make predictions with.
            print(algorithm, e.args[0]);
            logger.info("No matchup data left after subsettings (region:"+region+", algo: "+algorithm.__class__.__name__+")");
            continue;
        
        #Store the data/objects which we'll used later
        algorithmFunctorList.append(algorithm);
        modelOutputList.append(modelOutput);
        matchupRowsUsedList.append(dataUsed.index);
    
    #####################
    # Calculate metrics #
    #####################
    if len(algorithmFunctorList) == 0:
        print("*** WARNING: No algorithms were executed for region", region, " possibly because there was no matchup data in this region.");
        logger.error("No algorithms were executed for region '"+region+"' possibly because there was no matchup data in this region.");
        continue; #Skip calculating metrics for this region
    
    basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores = \
        metrics.calc_all_metrics(algorithmFunctorList, modelOutputList, matchupRowsUsedList, matchupData, settings);
    
    #########################
    # Write outputs to file #
    #########################
    outputDirectory = path.join(settings["outputPathMetrics"], outputVariableSelector.upper(), region);
    if path.exists(outputDirectory) == False:
        makedirs(outputDirectory);
    print("Writing metrics to: ", outputDirectory);
    pickle.dump(basicMetrics, open(path.join(outputDirectory, "basic_metrics.json"), 'wb'));
    np.savetxt(path.join(outputDirectory, "n_intersect_matrix.csv"), nIntersectMatrix, delimiter=",", fmt='%i');
    np.savetxt(path.join(outputDirectory, "paired_score_matrix.csv"), pairedScoreMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_wscore_matrix.csv"), pairedWScoreMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_rmsd_matrix.csv"), pairedRmsdMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_wrmsd_matrix.csv"), pairedWRmsdMatrix, delimiter=",");
    finalScores.to_csv(path.join(outputDirectory, "final_scores.csv"), sep=",", na_rep="nan", index=False);
    
    #Construct new matchup datasets with additional predicted column and store output
    for i in range(len(algorithmFunctorList)):
        algorithmName = type(algorithmFunctorList[i]).__name__;
        outputVariable = algorithmFunctorList[i].output_name();
        reconstructedData = matchupData.loc[matchupRowsUsedList[i]];
        algorithmOutput = basicMetrics[i]["model_output"];
        reconstructedData[outputVariableSelector.upper()+"_pred"] = algorithmOutput;
        reconstructedData.to_csv(path.join(outputDirectory, "matchup_appended_"+algorithmName+".csv"), sep=",");
        
        #tmp sanity checking
        if region == "global":
            from osoda_analysis_tools import prediction_accuracy_plot;
            prediction_accuracy_plot(reconstructedData[outputVariable], reconstructedData[outputVariable+"_pred"], algorithmName, outputVariable)

    #write algorithm list/order to file
    algoNameOrder = np.array([type(algoFunctor).__name__ for algoFunctor in algorithmFunctorList]);
    with open(path.join(outputDirectory, "algorithms_used.csv"), 'w') as file:
        for algoName in algoNameOrder:
            file.write(algoName+",");

#Shutdown logger
loggerFileHandle.close();












