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
from string import Template;
import pandas as pd;

import metrics;
import utilities;
import osoda_global_settings;

#control flags (turn parts of the script on/off)
runAlgorithmsAndCalculateMetrics = True;
diagnosticPlots = False; #Diagnostic plots to help see if algorithms are making sensible predictions


#Run an algorithm on the matchup dataset and returns the predicted (model) output and the rows of the matchup database which those predictions were made from
#algorithmInstance is an object of a class derived from BaseAlgorithm
#matchupData is a pandas DataFrame containing the complete matchup data set
#regionMaskPath, depthMaskPath, distToCoastMaskpath:    file paths to spatial mask netCDF files
#region, depthMaskVar, distToCoastMaskVar:    string variable names for accessing the variable containing gridded mask (0=filter, 1=keep) for each of the mask netCDFs
def run_algorithm(algorithmInstance, matchupData, regionMaskPath=None, region=None, depthMaskPath=None, depthMaskVar=None, distToCoastMaskPath=None, distToCoastMaskVar=None):
    
    subsetData = matchupData; #new view of the matchup data
    
    #Subset matchup dataset based on region
    if regionMaskPath != None:
        regionMaskNC = Dataset(regionMaskPath, 'r');
        subsetData = utilities.subset_from_mask(matchupData, regionMaskNC, region);
        
    #Subset based on depth
    if depthMaskPath != None:
        depthMaskNC = Dataset(depthMaskPath, 'r');
        subsetData = utilities.subset_from_mask(subsetData, depthMaskNC, depthMaskVar);
    
    #Subset based on distance to coast
    if distToCoastMaskPath != None:
        distToCoastMaskNC = Dataset(distToCoastMaskPath, 'r');
        subsetData = utilities.subset_from_mask(subsetData, distToCoastMaskNC, distToCoastMaskVar);
    
    #Run the algorithm to make predictions using the subsetted matchup dataset as input
    try:
        modelOutput, dataUsed = algorithm(subsetData); #model output is the predictions, dataUsed indicates a further subset of the matchup dataset (e.g. some rows may not contain values for all required inputs, or may exceed the valid ranges for the algorithm)
    except ValueError: #After subsetting there's no data left to make predictions with.
        raise; #Pass exception up the stack
    
    return modelOutput, dataUsed;
  
#Writes basic metrics objects, paired metrics matrices, final score dataframe, appended matchup dataset subsets (with predicted values), and algorithm name ordered list to file
def write_metrics_to_file(outputDirectory, basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores, algorithmFunctorList, matchupRowsUsedList):
    if path.exists(outputDirectory) == False:
        makedirs(outputDirectory);
    
    #Write basicMetrics object and paired metric matrices
    pickle.dump(basicMetrics, open(path.join(outputDirectory, "basic_metrics.json"), 'wb'));
    np.savetxt(path.join(outputDirectory, "n_intersect_matrix.csv"), nIntersectMatrix, delimiter=",", fmt='%i');
    np.savetxt(path.join(outputDirectory, "paired_score_matrix.csv"), pairedScoreMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_wscore_matrix.csv"), pairedWScoreMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_rmsd_matrix.csv"), pairedRmsdMatrix, delimiter=",");
    np.savetxt(path.join(outputDirectory, "paired_wrmsd_matrix.csv"), pairedWRmsdMatrix, delimiter=",");
    finalScores.to_csv(path.join(outputDirectory, "final_scores.csv"), sep=",", na_rep="nan", index=False);
    
    #Construct new matchup datasets with additional predicted column and output this to file
    for i in range(len(algorithmFunctorList)):
        algorithmName = type(algorithmFunctorList[i]).__name__;
        outputVariable = algorithmFunctorList[i].output_name(); #The common name of the output variable (i.e. 'AT' or 'DIC')
        reconstructedData = matchupData.loc[matchupRowsUsedList[i]];
        algorithmOutput = basicMetrics[i]["model_output"];
        reconstructedData[outputVariable.upper()+"_pred"] = algorithmOutput;
        reconstructedData.to_csv(path.join(outputDirectory, "matchup_appended_"+algorithmName+".csv"), sep=",");

    #write algorithm list/order to file
    algoNameOrder = np.array([type(algoFunctor).__name__ for algoFunctor in algorithmFunctorList]);
    with open(path.join(outputDirectory, "algorithms_used.csv"), 'w') as file:
        for algoName in algoNameOrder:
            file.write(algoName+",");


if __name__ == "__main__":
    #Setup logger
    logger = logging.getLogger('osoda_driver');
    logger.setLevel(logging.DEBUG);
    loggerFileHandle = logging.FileHandler('logs/osoda_driver.log', mode='w');
    loggerFileHandle.setLevel(logging.DEBUG);
    logger.addHandler(loggerFileHandle);
    
    logger.info("Started execution at: "+str(datetime.datetime.now()));
    
    #Load the settings for this run
    settings = osoda_global_settings.get_default_settings(); #TODO: Hannah - Change this to call your settings function
    
    #We'll run each algorithm for every possible combination of input variables.
    #This gets a list of dictionaries, each dictionary represents the mapping of the ocean parameter names to database variable names for one specific combination input variables
    #Also returns a string representation of the names for these (which is used later for creating output directories and distinguishing combinations from one another)
    specificVariableToDatabaseMaps, specificVariableToDatabaseMapNames = utilities.get_dataset_variable_map_combinations(settings);
    
    
    if runAlgorithmsAndCalculateMetrics == True:
        ### Iterate through each specific combination of input variables (aka specific variable to database mappings)
        for inputCombinationName, inputCombination in zip(specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps):
            logger.info("Begining new combination of input variables: "+inputCombinationName);
            print("Begining new combination of input variables: "+inputCombinationName);
            
            #Create a directory to contain the metrics for this combination
            currentCombinationOutputDirectory = path.join(settings["outputPathMetrics"], inputCombinationName);
            if path.exists(currentCombinationOutputDirectory) == False:
                makedirs(currentCombinationOutputDirectory);
            
            #Write information about the combination of input datasets used:
            utilities.write_specific_variable_to_database_mapping(inputCombination, path.join(currentCombinationOutputDirectory, "inputs_used.txt"), inputCombinationName);
            
            ### Read matchup database using the current set of input data
            years = utilities.calculate_years_for_input_combination(settings, inputCombination);
            matchupData = utilities.load_matchup_to_dataframe(settings, inputCombination, years=years); #each year concatinated is concatinated to create a single dataframe
    
    
            ### For each region, run all the algorithms that are relevant for that region
            for region in settings["regions"]:
                print("Running for region:", region);
                logger.info("Beginning new region: "+region);
                algorithmFunctorList = [];
                modelOutputList = [];
                matchupRowsUsedList = [];
                
                
                ######################################
                ### For each algorithm in the region, run the algorithm using relevent matchup database rows and store model output
                for AlgorithmFunctor in settings["algorithmRegionMapping"][region]:
                    algorithm = AlgorithmFunctor(settings); #Create an instance of the algorithm functor, this is the handle used to access the algorithm
                    
                    try:
                        modelOutput, dataUsed = run_algorithm(algorithm, matchupData,
                                                              regionMaskPath=settings["regionMasksPath"], region=region,
                                                              depthMaskPath=settings["depthMaskPath"], depthMaskVar=settings["depthMaskVar"],
                                                              distToCoastMaskPath=settings["distToCoastMaskPath"], distToCoastMaskVar=settings["distToCoastMaskVar"]);
                        logger.info("Output calculated (region:"+region+", algo: "+algorithm.__class__.__name__+")");
                    except ValueError as e: #Raised if there are no matchup data rows left after spatial mask and algorithm internal subsetting has taken place
                        print(algorithm, e.args[0]);
                        logger.info("No matchup data left after subsettings (region:"+region+", algo: "+algorithm.__class__.__name__+")");
                        continue;
                    
                    #Store the data/objects which we'll used later
                    algorithmFunctorList.append(algorithm);
                    modelOutputList.append(modelOutput);
                    matchupRowsUsedList.append(dataUsed.index);
                    
                    ######################
                    ### Diagnostic plots #
                    #Diagnostic plot to compare model and reference output
                    if diagnosticPlots == True:
                        from osoda_analysis_tools import prediction_accuracy_plot;
                        outputVariable = algorithm.output_name();
                        prediction_accuracy_plot(dataUsed[outputVariable], modelOutput, algorithm.__class__.__name__, outputVariable);
                    
                
                #######################################
                ### Calculate metrics for each region #
                if len(algorithmFunctorList) == 0:
                    print("*** WARNING: No algorithms were executed for region", region, " possibly because there was no matchup data in this region.");
                    logger.error("No algorithms were executed for region '"+region+"' possibly because there was no matchup data in this region.");
                    continue; #Skip calculating metrics for this region
                
                #Only compare algorithms which predict the same output variable, so we need a list of all the output variables
                uniqueOutputVars = np.unique([algorithm.output_name() for algorithm in algorithmFunctorList]); #Get a list of all the output variables
                for currentOutputVar in uniqueOutputVars:
                    #filter algorithm list, output list etc. by the current output variable (in order to compare like-for-like)
                    toKeep = [i for i in range(len(algorithmFunctorList)) if algorithmFunctorList[i].output_name()==currentOutputVar];
                    currentAlgorithmFunctorList = [algorithmFunctorList[i] for i in toKeep];
                    currentModelOutputList = [modelOutputList[i] for i in toKeep];
                    currentMatchupRowsUsedList = [matchupRowsUsedList[i] for i in toKeep];
                    
                    basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores = \
                        metrics.calc_all_metrics(currentAlgorithmFunctorList, currentModelOutputList, currentMatchupRowsUsedList, matchupData, settings);
                    
                    ###########################
                    ### Write outputs to file #
                    outputDirectory = path.join(currentCombinationOutputDirectory, currentOutputVar, region);
                    print("Writing metrics to: ", outputDirectory);
                    write_metrics_to_file(outputDirectory, basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores, currentAlgorithmFunctorList, currentMatchupRowsUsedList);
        
        
    ##############################
    ### Calculating best metrics #
    #FIXME: appending columns to dataframe in two different ways, this is confusing and annoying...
    bestAlgosAT = []; bestAlgosDIC = []; bestRMSDesAT = []; bestRMSDesDIC = []; regionNamesList = []; combinationNamesList = []; numAlgosComparedAT = []; numAlgosComparedDIC = [];
    for region in settings["regions"]:
        for inputCombinationName in specificVariableToDatabaseMapNames:
            metricsRootDirectory = path.join(settings["outputPathMetrics"], inputCombinationName);
            bestAlgorithmInfo = utilities.find_best_algorithm(metricsRootDirectory, region, useWeightedRMSDe=settings["assessUsingWeightedRMSDe"], verbose=False);
            
            #Append to columns
            combinationNamesList.append(inputCombinationName);
            regionNamesList.append(region);
            if bestAlgorithmInfo["AT"] != None: #If no paired metrics could be calculated (e.g. no overlapping algorithms, no RMSDs reported for algorithms, no matchup data for a region) None is returned instead of a tuple
                bestAlgosAT.append(bestAlgorithmInfo["AT"][0]);
                bestRMSDesAT.append(bestAlgorithmInfo["AT"][1]);
                numAlgosComparedAT.append(bestAlgorithmInfo["AT"][2]);
            else: #No best algorithm, fill with default data
                bestAlgosAT.append(np.nan);
                bestRMSDesAT.append(np.nan);
                numAlgosComparedAT.append(np.nan);
            if bestAlgorithmInfo["DIC"] != None: #If no paired metrics could be calculated (e.g. no overlapping algorithms, no RMSDs reported for algorithms, no matchup data for a region) None is returned instead of a tuple
                bestAlgosDIC.append(bestAlgorithmInfo["DIC"][0]);
                bestRMSDesDIC.append(bestAlgorithmInfo["DIC"][1]);
                numAlgosComparedDIC.append(bestAlgorithmInfo["DIC"][2]);
            else: #No best algorithm, fill with default data
                bestAlgosDIC.append(np.nan);
                bestRMSDesDIC.append(np.nan);
                numAlgosComparedDIC.append(np.nan);
    
    summaryTable = pd.DataFrame();
    summaryTable["region"] = regionNamesList;
    summaryTable["input_combination"] = combinationNamesList;
    summaryTable["AT_best_in_region"] = [False]*len(summaryTable);
    summaryTable["DIC_best_in_region"] = [False]*len(summaryTable);
    summaryTable["overall_best_in_region"] = [False]*len(summaryTable);
    summaryTable["AT_best_algorithm"] = bestAlgosAT;
    summaryTable["AT_RMSDe"] = bestRMSDesAT;
    summaryTable["AT_n"] = [np.nan]*len(summaryTable);
    summaryTable["AT_algos_compared"] = numAlgosComparedAT;
    summaryTable["DIC_best_algorithm"] = bestAlgosDIC;
    summaryTable["DIC_RMSDe"] = bestRMSDesDIC;
    summaryTable["DIC_n"] = [np.nan]*len(summaryTable);
    summaryTable["DIC_algos_compared"] = numAlgosComparedDIC;
    summaryTable["n_years"] = [0]*len(summaryTable);
    summaryTable["min_year"] = [0]*len(summaryTable);
    summaryTable["max_year"] = [0]*len(summaryTable);
    
    #Some meta metrics
    #Which input combination is the best for AT, DIC and both, for each region?
    for region in settings["regions"]:
        #Mark the input combination with the RMSDe in each region (first for AT)
        rmsdesInRegionAT = summaryTable[summaryTable["region"]==region]["AT_RMSDe"];
        if rmsdesInRegionAT.isnull().all() == False: #If any are real numbers
            bestRmsdAT = rmsdesInRegionAT[rmsdesInRegionAT==np.nanmin(rmsdesInRegionAT)];
            summaryTable.loc[bestRmsdAT.index, "AT_best_in_region"] = True;
            #summaryTable["AT_best_in_region"][(summaryTable["region"]==region) & (summaryTable["AT_RMSDe"]==np.nanmin(rmsdesInRegionAT))] = True;
        
        #Now again with DIC
        rmsdesInRegionDIC = summaryTable[summaryTable["region"]==region]["DIC_RMSDe"];
        if rmsdesInRegionDIC.isnull().all() == False: #If any are real numbers
            bestRmsdDIC = rmsdesInRegionDIC[rmsdesInRegionDIC==np.nanmin(rmsdesInRegionDIC)];
            summaryTable.loc[bestRmsdDIC.index, "DIC_best_in_region"] = True;
            #summaryTable["DIC_best_in_region"][(summaryTable["region"]==region) & (summaryTable["DIC_RMSDe"]==np.nanmin(rmsdesInRegionDIC))] = True;
        
        #Now the sum of both
        rmsdesInRegionSum = rmsdesInRegionAT+rmsdesInRegionDIC;
        if rmsdesInRegionSum.isnull().all() == False: #If any are real numbers
            bestRmsdSum = rmsdesInRegionSum[rmsdesInRegionSum==np.nanmin(rmsdesInRegionSum)];
            summaryTable.loc[bestRmsdSum.index, "overall_best_in_region"] = True;
    
    #How many years, and what was the min and max time range, for each combination?
    for (inputCombinationName, inputCombination) in zip(specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps):
        yearsUsed = utilities.calculate_years_for_input_combination(settings, inputCombination);
        currentCombinationRows = summaryTable[summaryTable["input_combination"]==inputCombinationName].index;
        summaryTable.loc[currentCombinationRows, "n_years"] = len(yearsUsed);
        summaryTable.loc[currentCombinationRows, "min_year"] = min(yearsUsed);
        summaryTable.loc[currentCombinationRows, "max_year"] = max(yearsUsed);
        
    
    #number of data points ran for
    for r, row in summaryTable.iterrows():
        if isinstance(row["AT_best_algorithm"], str):
            predictedDataPathAT = path.join(settings["outputPathMetrics"], row["input_combination"], "AT", row["region"], "matchup_appended_"+row["AT_best_algorithm"]+".csv");
            summaryTable.loc[r, "AT_n"] = len(pd.read_csv(predictedDataPathAT));
        
        if isinstance(row["DIC_best_algorithm"], str):
            predictedDataPathDIC = path.join(settings["outputPathMetrics"], row["input_combination"], "DIC", row["region"], "matchup_appended_"+row["DIC_best_algorithm"]+".csv");
            summaryTable.loc[r, "DIC_n"] = len(pd.read_csv(predictedDataPathDIC));
    
    #Write to file
    summaryTableOutputPath = path.join(settings["outputPathMetrics"], "summary_best_algos.csv");
    summaryTable.to_csv(summaryTableOutputPath, sep=",", index=False);     
    print("Full summary table written to:", path.abspath(summaryTableOutputPath));
    
#    print("Combination name variable key as follows:");
#    utilities.print_combination_name_keys(settings);
    
    ### Calculate overall best algorithms / input combinations for each egion
    overallBestAlgos = pd.DataFrame(columns=["region", "output_var", "input_combination", "algo_name", "RMSDe", "n", "algos_compared"]);
    for region in settings["regions"]:
        #find best AT algorithm info by sorting a subset of the summary table for just this region
        regionTable = summaryTable[summaryTable["region"] == region];
        regionTable = regionTable.sort_values(by=["AT_RMSDe", "AT_n", "AT_algos_compared", "DIC_RMSDe"], ascending=[True, False, False, True]);
        overallBestAlgos.loc[len(overallBestAlgos)] = [region, "AT", regionTable.iloc[0]["input_combination"], regionTable.iloc[0]["AT_best_algorithm"], regionTable.iloc[0]["AT_RMSDe"], regionTable.iloc[0]["AT_n"], regionTable.iloc[0]["AT_algos_compared"]];
        
        #find best DIC algorithm info by sorting a subset of the summary table for just this region
        regionTable = summaryTable[summaryTable["region"] == region];
        regionTable = regionTable.sort_values(by=["DIC_RMSDe", "DIC_n", "DIC_algos_compared", "AT_RMSDe"], ascending=[True, False, False, True]);
        overallBestAlgos.loc[len(overallBestAlgos)] = [region, "DIC", regionTable.iloc[0]["input_combination"], regionTable.iloc[0]["DIC_best_algorithm"], regionTable.iloc[0]["DIC_RMSDe"], regionTable.iloc[0]["DIC_n"], regionTable.iloc[0]["DIC_algos_compared"]];
    
    overallBestAlgosOutputPath = path.join(settings["outputPathMetrics"], "overall_best_algos.csv");
    overallBestAlgos.to_csv(overallBestAlgosOutputPath, sep=",", index=False);     
    print("Overall best algorithm table written to:", path.abspath(overallBestAlgosOutputPath));
    
    #Shutdown logger
    loggerFileHandle.close();












