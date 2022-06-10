#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 14:36:57 2019

This utility performs the algorithm comparison step, and computes weighted and unweighted metrics for each region/input data combination.
It uses the methodology of Land et al 2019 ( https://www.sciencedirect.com/science/article/pii/S0034425719304882 ).
It's output directory and settings are controlled by the global settings file (osoda_global_settings.py:get_default_settings)

By default, the following output files are created:
    overall_best_algos.csv, overall_best_algos_min_years=8.csv   These give a summary of the best performing algorithm and input data combination for each region
                                                                     with and without the minimum 8 year temporal range constraint
    summary_best_algos.csv, summary_best_algos_unweighted.csv    These give the best algorithms for each input data combination (to assess difference between input data combinations)
                                                                     for weighted and unweighted metrics.
    combination<N>_<input_datasets> directories                  These are separate directories for each input data combination, where N is the combination enumeration and <input_datasets> is
                                                                     the named SST and SSS data sets used. Inside are separate files for each region and output variable (AT/DIC), which contain
                                                                     all the detailed info for the runs for that region/input combination. This includes:
                                                                            * a list of the algorithms assessed (algorithms_used.csv)
                                                                            * the basic (non-paired) metrics for each algorithm (basic_metrics.json) - note you can read this using pickle in python
                                                                            * subset matchup data set used for each algorithm, including the predicted DIC/AT from the algorithm (matchup_appended_<algorithm>.csv)
                                                                            * data files for the paired metrics. These are .csv files containing a matrix, where the rows/columns correspond to each algorithm
                                                                                 and the order of algorithms is as defined by algorithms_used.csv. Included are:
                                                                                     * number of intersecting matchup database rows (n_intersect_matrix.csv)
                                                                                     * paired_wrmsd_matrix.csv, paired_rmsd_matrix.csv (the RMSDe for weighted and unweighted calculations)
                                                                                     * paired_wscore_matrix.csv, paired_score_matrix.csv (the final scores for each algorithm pair, weighted and unweighted respectively)
    (development, use with care)

@author: tom holding
"""

#TODO: This whole script needs refactoring. Lots of duplicate code, different patterns to do the same thing, confusing layout...

from os import path, makedirs;
import numpy as np;
from netCDF4 import Dataset;
import pickle;
import json;
import logging;
import datetime;
import pandas as pd;
import csv

import os_algorithms.metrics as metrics;
import os_algorithms.utilities as utilities;
import osoda_global_settings;
from os_algorithms.diagnostic_plotting import prediction_accuracy_plot;

#control flags (turn parts of the script on/off)
runAlgorithmsAndCalculateMetrics = True;
diagnosticPlots = True; #Diagnostic plots to help see if algorithms are making sensible predictions


#Run an algorithm on the matchup dataset and returns the predicted (model) output and the rows of the matchup database which those predictions were made from
#algorithmInstance is an object of a class derived from BaseAlgorithm
#matchupData is a pandas DataFrame containing the complete matchup data set
#regionMaskPath, depthMaskPath, distToCoastMaskpath:    file paths to spatial mask netCDF files
#region, depthMaskVar, distToCoastMaskVar:    string variable names for accessing the variable containing gridded mask (0=filter, 1=keep) for each of the mask netCDFs
def run_algorithm(algorithmInstance, matchupData, regionMaskPath=None, region=None, useDepthMask=False, depthMaskPath=None, depthMaskVar=None, useDistToCoastMask=False, distToCoastMaskPath=None, distToCoastMaskVar=None):
    
    subsetData = matchupData; #new view of the matchup data
    
    #Subset matchup dataset based on region
    if regionMaskPath != None:
        regionMaskNC = Dataset(regionMaskPath, 'r');
        subsetData = utilities.subset_from_mask(matchupData, regionMaskNC, region);
        regionMaskNC.close();     
    
    #Subset based on depth
    if useDepthMask == True:
        depthMaskNC = Dataset(depthMaskPath, 'r');
        subsetData = utilities.subset_from_mask(subsetData, depthMaskNC, depthMaskVar);
        depthMaskNC.close();
    
    #Subset based on distance to coast
    if useDistToCoastMask == True:
        distToCoastMaskNC = Dataset(distToCoastMaskPath, 'r');
        subsetData = utilities.subset_from_mask(subsetData, distToCoastMaskNC, distToCoastMaskVar);
        distToCoastMaskNC.close();
    
    #Run the algorithm to make predictions using the subsetted matchup dataset as input
    
    try:
        algorithmOutputTuple = algorithmInstance(subsetData); #model output is the predictions, dataUsed indicates a further subset of the matchup dataset (e.g. some rows may not contain values for all required inputs, or may exceed the valid ranges for the algorithm)
        modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed = algorithmOutputTuple;
    except ValueError: #After subsetting there's no data left to make predictions with.
        raise; #Pass exception up the stack
    
    return modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed,subsetData;


#Writes basic metrics objects, paired metrics matrices, final score dataframe, appended matchup dataset subsets (with predicted values), and algorithm name ordered list to file
def write_metrics_to_file(outputDirectory, matchupData, basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores, currentAlgorithmOutputs):
    if path.exists(outputDirectory) == False:
        makedirs(outputDirectory);
    
    #write algorithm list/order to file
    algoNameOrder = np.array([algorithmOutput["name"] for algorithmOutput in currentAlgorithmOutputs]);
    with open(path.join(outputDirectory, "algorithms_used.csv"), 'w') as file:
        for algoName in algoNameOrder:
            file.write(algoName+",");
    
    
    #Write basicMetrics object
    pickle.dump(basicMetrics, open(path.join(outputDirectory, "basic_metrics.pickle"), 'wb'));
    #Do not output long vectors of data in the json summary file, so copy the full basic metrics then delete the vectors
    basicMetricsJson = {algoNameOrder[i]: basicMetrics[i].copy() for i in range(len(basicMetrics))}
    for key in basicMetricsJson.keys():
        del basicMetricsJson[key]["model_output"];
        del basicMetricsJson[key]["reference_output_uncertainty"];
        del basicMetricsJson[key]["weights"];
        del basicMetricsJson[key]["model_uncertainty"];
        del basicMetricsJson[key]["model_propagated_input_uncertainty"];
        del basicMetricsJson[key]["model_combined_output_uncertainty"];
    with open(path.join(outputDirectory, "basic_metrics.json"), 'w') as file:
        json.dump(basicMetricsJson, file, indent=4);
    

    #Write paired metrics
    nIntersectMatrixDF = pd.DataFrame(nIntersectMatrix, dtype=int, columns=algoNameOrder, index=algoNameOrder);
    nIntersectMatrixDF.to_csv(path.join(outputDirectory, "n_intersect_matrix.csv"));
    pairedScoreMatrixDF = pd.DataFrame(pairedScoreMatrix, dtype=float, columns=algoNameOrder, index=algoNameOrder);
    pairedScoreMatrixDF.to_csv(path.join(outputDirectory, "paired_score_matrix.csv"), na_rep="nan");
    
    pairedWScoreMatrixDF = pd.DataFrame(pairedWScoreMatrix, dtype=float, columns=algoNameOrder, index=algoNameOrder);
    pairedWScoreMatrixDF.to_csv(path.join(outputDirectory, "paired_wscore_matrix.csv"), na_rep="nan");
    
    pairedRmsdMatrixDF = pd.DataFrame(pairedRmsdMatrix, dtype=float, columns=algoNameOrder, index=algoNameOrder);
    pairedRmsdMatrixDF.to_csv(path.join(outputDirectory, "paired_rmsd_matrix.csv"), na_rep="nan");
    
    pairedWRmsdMatrixDF = pd.DataFrame(pairedWRmsdMatrix, dtype=float, columns=algoNameOrder, index=algoNameOrder);
    pairedWRmsdMatrixDF.to_csv(path.join(outputDirectory, "paired_wrmsd_matrix.csv"), na_rep="nan");
    
    finalScores.to_csv(path.join(outputDirectory, "final_scores.csv"), sep=",", na_rep="nan", index=False);
    
    #Construct new matchup datasets with additional predicted column and output this to file
    for i in range(len(currentAlgorithmOutputs)):
        #commented to debug
        # if(currentAlgorithmOutputs[i]["instance"])==None:
        #     algorithmName=currentAlgorithmOutputs[i]["name"]
        # else:
        algorithmName = type(currentAlgorithmOutputs[i]["instance"]).__name__;
        outputVariable = currentAlgorithmOutputs[i]["outputVar"]; #The common name of the output variable (i.e. 'AT' or 'DIC')
        reconstructedData = matchupData.loc[currentAlgorithmOutputs[i]["dataUsedIndices"]];
        reconstructedData[outputVariable.upper()+"_pred"] = basicMetrics[i]["model_output"];
        reconstructedData[outputVariable.upper()+"_pred_model_uncertainty"] = basicMetrics[i]["model_uncertainty"]; #RMSD of the model
        reconstructedData[outputVariable.upper()+"_pred_propagated_input_uncertainty"] = basicMetrics[i]["model_propagated_input_uncertainty"]; #Uncertainty due purely to inputs used by the model
        reconstructedData[outputVariable.upper()+"_pred_combined_output_uncertainty"] = basicMetrics[i]["model_combined_output_uncertainty"]; #combined RMSD of model and input uncertainty
        reconstructedData[outputVariable.upper()+"_reference_output_uncertainty"] = basicMetrics[i]["reference_output_uncertainty"]; #uncertainty of the reference (in situ) outputs
        
        
        #reconstructedData[outputVariable.upper()+"_pred_err"] = basicMetrics[i]["model_output_err"];
        reconstructedData.to_csv(path.join(outputDirectory, "matchup_appended_"+algorithmName+".csv"), sep=",");


#Extracts the best algorithm for each region and input data combination.
#Compiles information into a dataframe and returns
def create_summary_table(n_threshold,settings, inputCombinationNames, inputCombinationVariableMaps, useWeighted, regionOverload=None):
    #FIXME: appending columns to dataframe in two different ways, this is confusing and annoying...
    bestAlgosAT = [];
    bestAlgosDIC = [];
    bestRMSDesAT = [];
    bestRMSDesDIC = [];
    regionNamesList = [];
    combinationNamesList = [];
    numAlgosComparedAT = [];
    numAlgosComparedDIC = [];
    bestbiasAT=[];
    bestbiasDIC=[];
    bestuncendtoendAT=[];
    bestRMSD_DIC=[]; 
    bestRMSD_AT=[];
    bestuncendtoendDIC=[]; 
    bestwRMSD_DIC=[]; 
    bestwRMSD_AT=[];
    
    if regionOverload is not None:
        regions = regionOverload;
    else:
        regions = settings["regions"];

    for region in regions:
        for inputCombinationName in inputCombinationNames:
            metricsRootDirectory = path.join(settings["outputPathMetrics"], inputCombinationName);
            bestAlgorithmInfo = utilities.find_best_algorithm(n_threshold,metricsRootDirectory, region, useWeightedRMSDe=useWeighted, verbose=False);
            
            #Append to columns
            regionNamesList.append(region);
            combinationNamesList.append(inputCombinationName);
            if bestAlgorithmInfo["AT"] != None: #If no paired metrics could be calculated (e.g. no overlapping algorithms, no RMSDs reported for algorithms, no matchup data for a region) None is returned instead of a tuple
                bestAlgosAT.append(bestAlgorithmInfo["AT"][0]);
                bestRMSDesAT.append(bestAlgorithmInfo["AT"][1]);
                numAlgosComparedAT.append(bestAlgorithmInfo["AT"][2]);
                bestbiasAT.append(bestAlgorithmInfo["AT"][3]);
                bestuncendtoendAT.append(bestAlgorithmInfo["AT"][4]);
                bestRMSD_AT.append(bestAlgorithmInfo["AT"][5]);
                bestwRMSD_AT.append(bestAlgorithmInfo["AT"][6]);


            else: #No best algorithm, fill with default data
                bestAlgosAT.append(np.nan);
                bestRMSDesAT.append(np.nan);
                numAlgosComparedAT.append(np.nan);
                bestbiasAT.append(np.nan); 
                bestuncendtoendAT.append(np.nan);
                bestRMSD_AT.append(np.nan);
                bestwRMSD_AT.append(np.nan);


            if bestAlgorithmInfo["DIC"] != None: #If no paired metrics could be calculated (e.g. no overlapping algorithms, no RMSDs reported for algorithms, no matchup data for a region) None is returned instead of a tuple
                bestAlgosDIC.append(bestAlgorithmInfo["DIC"][0]);
                bestRMSDesDIC.append(bestAlgorithmInfo["DIC"][1]);
                numAlgosComparedDIC.append(bestAlgorithmInfo["DIC"][2]);
                bestbiasDIC.append(bestAlgorithmInfo["DIC"][3]);
                bestuncendtoendDIC.append(bestAlgorithmInfo["DIC"][4]);
                bestRMSD_DIC.append(bestAlgorithmInfo["DIC"][5]);
                bestwRMSD_DIC.append(bestAlgorithmInfo["DIC"][6]);

            else: #No best algorithm, fill with default data
                bestAlgosDIC.append(np.nan);
                bestRMSDesDIC.append(np.nan);
                numAlgosComparedDIC.append(np.nan);
                bestbiasDIC.append(np.nan);
                bestuncendtoendDIC.append(np.nan);
                bestRMSD_DIC.append(np.nan);
                bestwRMSD_DIC.append(np.nan);

    #Concatinate to a dataframe
    summaryTable = pd.DataFrame();
    summaryTable["region"] = regionNamesList;
    summaryTable["input_combination"] = combinationNamesList;
    summaryTable["AT_best_in_region"] = [False]*len(summaryTable);
    summaryTable["DIC_best_in_region"] = [False]*len(summaryTable);
    
    summaryTable["AT_best_algorithm"] = bestAlgosAT;
    summaryTable["AT_RMSDe"] = bestRMSDesAT;
    summaryTable["AT_n"] = [np.nan]*len(summaryTable);
    summaryTable["AT_algos_compared"] = numAlgosComparedAT;
    summaryTable["AT_bias"] = bestbiasAT;
    summaryTable["AT_uncendtoend"] = bestuncendtoendAT;
    summaryTable["AT_RMSD"] = bestRMSD_AT;
    summaryTable["AT_wRMSD"] = bestwRMSD_AT;

    summaryTable["DIC_best_algorithm"] = bestAlgosDIC;
    summaryTable["DIC_RMSDe"] = bestRMSDesDIC;
    summaryTable["DIC_n"] = [np.nan]*len(summaryTable);
    summaryTable["DIC_algos_compared"] = numAlgosComparedDIC;
    summaryTable["DIC_bias"] = bestbiasDIC;
    summaryTable["DIC_uncendtoend"] = bestuncendtoendDIC;
    summaryTable["DIC_RMSD"] = bestRMSD_DIC;
    summaryTable["DIC_wRMSD"] = bestwRMSD_DIC;

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
    
    #How many years, and what was the min and max time range, for each combination?
    for (inputCombinationName, inputCombination) in zip(inputCombinationNames, inputCombinationVariableMaps):
        yearsUsed = utilities.calculate_years_for_input_combination(settings, inputCombination);
        currentCombinationRows = summaryTable[summaryTable["input_combination"]==inputCombinationName].index;
        summaryTable.loc[currentCombinationRows, "n_years"] = len(yearsUsed);
        summaryTable.loc[currentCombinationRows, "min_year"] = min(yearsUsed);
        summaryTable.loc[currentCombinationRows, "max_year"] = max(yearsUsed);
        
    
    #number of data points ran for
    for r, row in summaryTable.iterrows():
        if isinstance(row["AT_best_algorithm"], str):
            predictedDataPathAT = path.join(settings["outputPathMetrics"], row["input_combination"], "AT", row["region"], "matchup_appended_"+row["AT_best_algorithm"].split(":")[0]+".csv");
            summaryTable.loc[r, "AT_n"] = len(pd.read_csv(predictedDataPathAT));
        
        if isinstance(row["DIC_best_algorithm"], str):
            predictedDataPathDIC = path.join(settings["outputPathMetrics"], row["input_combination"], "DIC", row["region"], "matchup_appended_"+row["DIC_best_algorithm"].split(":")[0]+".csv");
            summaryTable.loc[r, "DIC_n"] = len(pd.read_csv(predictedDataPathDIC));
    
    #Return summary table
    return summaryTable;


#Returns a data frame comtaining key metrics for the best algorithm and input data combination for each region.
#   summaryTable: the metric summary data for all algorithms and input combinations. If a string this is assumed to be a path to the table and read into a dataframe
#   regions: ocean soda region strings
#   minTimeSpanYears: minimum number of years for which the input data combination covers. Combinations below this number will be excluded from consideration
def get_best_algorithms(_summaryTable, regions, minTimeSpanYears = 0, minMatchupsTA = 0, minMatchupsDIC =0):
    if isinstance(_summaryTable, str):
        summaryTable = pd.read_csv(_summaryTable, sep=",");
    else:
        summaryTable = _summaryTable;
    
    overallBestAlgos = pd.DataFrame(columns=["region", "output_var", "input_combination", "algo_name", "RMSDe", "n", "algos_compared","bias","uncendtoend","RMSD","wRMSD", "n_years", "min_year", "max_year"]);
    for region in regions:
        #find best AT algorithm info by sorting a subset of the summary table for just this region, and only include entries where the time range meets the criteria
        regionTable = summaryTable[(summaryTable["region"] == region) & (summaryTable["n_years"] >= minTimeSpanYears) & (summaryTable["AT_n"] >= minMatchupsTA)];
        if len(regionTable) > 0: #if there are any entries left, sort and select the best
            regionTable = regionTable.sort_values(by=["AT_RMSDe", "AT_n", "AT_algos_compared","AT_bias","AT_uncendtoend","AT_RMSD","AT_wRMSD", "DIC_RMSDe"], ascending=[True, False, False,True,True, True,True, True]);
            overallBestAlgos.loc[len(overallBestAlgos)] = [region, "AT", regionTable.iloc[0]["input_combination"], regionTable.iloc[0]["AT_best_algorithm"], regionTable.iloc[0]["AT_RMSDe"], regionTable.iloc[0]["AT_n"],  regionTable.iloc[0]["AT_algos_compared"],regionTable.iloc[0]["AT_bias"], regionTable.iloc[0]["AT_uncendtoend"], regionTable.iloc[0]["AT_RMSD"],regionTable.iloc[0]["AT_wRMSD"],regionTable.iloc[0]["n_years"], regionTable.iloc[0]["min_year"], regionTable.iloc[0]["max_year"]];
        else: #no entries, so inserts nans
            overallBestAlgos.loc[len(overallBestAlgos)] = [region, "AT", np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan,np.nan,np.nan, 0, 0, 0];
        
        #find best DIC algorithm info by sorting a subset of the summary table for just this region
        regionTable = summaryTable[(summaryTable["region"] == region) & (summaryTable["n_years"] >= minTimeSpanYears) & (summaryTable["DIC_n"] >= minMatchupsDIC)];
        if len(regionTable) > 0: #if there are any entries left, sort and select the best
            regionTable = regionTable.sort_values(by=["DIC_RMSDe", "DIC_n", "DIC_algos_compared","DIC_bias","DIC_uncendtoend","DIC_RMSD","DIC_wRMSD","AT_RMSDe"], ascending=[True, False, False,True,True,True,True, True]);
            overallBestAlgos.loc[len(overallBestAlgos)] = [region, "DIC", regionTable.iloc[0]["input_combination"], regionTable.iloc[0]["DIC_best_algorithm"], regionTable.iloc[0]["DIC_RMSDe"], regionTable.iloc[0]["DIC_n"], regionTable.iloc[0]["DIC_algos_compared"],regionTable.iloc[0]["DIC_bias"],regionTable.iloc[0]["DIC_uncendtoend"],regionTable.iloc[0]["DIC_RMSD"],regionTable.iloc[0]["DIC_wRMSD"],   regionTable.iloc[0]["n_years"], regionTable.iloc[0]["min_year"], regionTable.iloc[0]["max_year"]];
        else: #no entries, so insert nans
            overallBestAlgos.loc[len(overallBestAlgos)] = [region, "DIC", np.nan, np.nan, np.nan, np.nan,np.nan, np.nan, np.nan,np.nan,np.nan, 0, 0, 0];
        
    return overallBestAlgos;

def main(settings, extraAlgosToTest=[]):
    #Setup logger
    logPath = path.join(settings["logDirectoryRoot"], "osoda_algorithm_comparison.log");
    if path.exists(path.dirname(logPath)) == False:
        makedirs(path.dirname(logPath));
    logger = logging.getLogger("osoda_algorithm_comparison");
    logger.setLevel(logging.DEBUG);
    loggerFileHandle = logging.FileHandler(logPath, mode='w');
    loggerFileHandle.setLevel(logging.DEBUG);
    logger.addHandler(loggerFileHandle);
    logger.info("Started execution at: "+str(datetime.datetime.now()));
    
    # these are output dictionarys so we can get the final scores output 
    # in a consolidated format which is very useful  for plotting
    final_scores_allcombos_oceansoda_amazon_plume_AT={};
    final_scores_allcombos_oceansoda_amazon_plume_DIC={};
    
    final_scores_allcombos_oceansoda_congo_AT={};
    final_scores_allcombos_oceansoda_congo_DIC={};
    
    final_scores_allcombos_oceansoda_mississippi_AT={};
    final_scores_allcombos_oceansoda_mississippi_DIC={};
    
    final_scores_allcombos_oceansoda_st_lawrence_AT={};
    final_scores_allcombos_oceansoda_st_lawrence_DIC={};
    
    final_scores_allcombos_oceansoda_mediterranean_AT={};
    final_scores_allcombos_oceansoda_mediterranean_DIC={};
    
    #We'll run each algorithm for every possible combination of input variables.
    #This gets a list of dictionaries, each dictionary represents the mapping of the ocean parameter names to database variable names for one specific combination input variables
    #Also returns a string representation of the names for these (which is used later for creating output directories and distinguishing combinations from one another)
    specificVariableToDatabaseMaps, specificVariableToDatabaseMapNames = utilities.get_dataset_variable_map_combinations(settings);
    
    if runAlgorithmsAndCalculateMetrics == True:
        #Grab some settings
        if settings["subsetWithDepthMask"] == True:
            depthMaskPath=settings["depthMaskPath"];
        else:
            depthMaskPath=None;
        if settings["subsetWithDistToCoast"] == True:
            distToCoastMaskPath=settings["distToCoastMaskVar"]
        else:
            distToCoastMaskPath=None;
            
        earthobs_dataset_iteration_number=1;
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
            matchupData = utilities.load_matchup_to_dataframe(settings, inputCombination, years=years); #each year is concatinated to create a single dataframe
            
            #checks on the Matchup databse files to check that they are
            #realistic and fall within the expected range
            
            SST_max=settings["MDB_flags"]["SST_max"];
            SST_min=settings["MDB_flags"]["SST_min"];
            SSS_max=settings["MDB_flags"]["SSS_max"];
            SSS_min=settings["MDB_flags"]["SSS_min"];      
            DIC_max=settings["MDB_flags"]["DIC_max"];
            DIC_min=settings["MDB_flags"]["DIC_min"]; 
            pH_max=settings["MDB_flags"]["pH_max"];
            pH_min=settings["MDB_flags"]["pH_min"];
            # pCO2_max=settings["MDB_flags"]["pCO2_max"];
            # pCO2_min=settings["MDB_flags"]["pCO2_min"];
            TA_max=settings["MDB_flags"]["TA_max"];
            TA_min=settings["MDB_flags"]["TA_min"];
            
            #SST
            mdb_SST = matchupData["SST"];
            
            if np.nanmean(mdb_SST.values) > 200.0: #Convert SST from C to K, if required
                mdb_SST_numeric=mdb_SST.values-273.15;#unit conversion
            else:
                mdb_SST_numeric=mdb_SST.values
            #Upper realistic limit 
            states=mdb_SST_numeric>SST_max;
            index_temp_exceed=np.where(states)[0]
            #Lower realistic limit -10 degrees
            states2=mdb_SST_numeric<SST_min;
            index_temp_below=np.where(states2)[0]
            
            #Salinity
            mdb_SSS = matchupData["SSS"];
            mdb_SSS_numeric=mdb_SSS.values;
            #Upper realistic limit 50 PSU
            states3=mdb_SSS_numeric>SSS_max;
            index_sal_exceed=np.where(states3)[0]
            #Lower realistic limit <0 PSU
            states4=mdb_SSS_numeric<SSS_min;
            index_sal_below=np.where(states4)[0]        
            
            #DIC
            mdb_DIC = matchupData["DIC"];
            mdb_DIC_numeric=mdb_DIC.values;
            #Upper realistic limit 2500 UMOLKG?
            states5=mdb_DIC_numeric>DIC_max;
            index_DIC_exceed=np.where(states5)[0]
            #Lower realistic limit <500 UMOL KG
            states6=mdb_DIC_numeric<DIC_min;
            index_DIC_below=np.where(states6)[0] 
            
            #pH
            mdb_pH = matchupData["region_pH_mean"];
            mdb_pH_numeric=mdb_pH.values;
            #Upper realistic limit 8.5
            states7=mdb_pH_numeric>pH_max;
            index_pH_exceed=np.where(states7)[0]
            #Lower realistic limit <7
            states8=mdb_pH_numeric<pH_min;
            index_pH_below=np.where(states8)[0]
            
            # #pCO2
            # mdb_pco2 = matchupData["region_pco2w_mean"];
            # mdb_pco2_numeric=mdb_pco2.values;
            # #Upper realistic limit 700 ppm
            # states9=mdb_pco2_numeric>pCO2_max;
            # index_pco2_exceed=np.where(states9)[0]
            # #Lower realistic limit <200 ppm
            # states10=mdb_pco2_numeric<pCO2_min;
            # index_pco2_below=np.where(states10)[0]
            
            #TA
            mdb_TA = matchupData["AT"];
            mdb_TA_numeric=mdb_TA.values;
            #Upper realistic limit 3000 umol kg
            states11=mdb_TA_numeric>TA_max;
            index_TA_exceed=np.where(states11)[0]
            #Lower realistic limit <500 umol kg
            states12=mdb_TA_numeric<TA_min;
            index_TA_below=np.where(states12)[0]        
            
            #note that other variables could also have bounds placed on them but not applied 
            #for this code iteration e.g. lat long date OC chla DO NO3 PO4 SiO4
            
            # now produce a file with all of the out of bounds data points from the mdb
            mdb_flag_warnings_list=['SST greater than maximum value of', 'SST less than minimum value of', 'SSS greater than maximum value of', 'SST less than minimum value of'\
                                    , 'DIC greater than maximum value of', 'DIC less than minimum value of','pH greater than maximum value of','pH less than minimum value of'\
                                    ,'TA greater than maximum value of','TA less than minimum value of']
        
                                    #,'pCO2 greater than maximum value of','pCO2 less than minimum value of'
            mdb_flag_index_list=[index_temp_exceed, index_temp_below, index_sal_exceed, index_sal_below, index_DIC_exceed, index_DIC_below,index_pH_exceed,index_pH_below,index_TA_exceed,index_TA_below]
                                             #index_pco2_exceed,index_pco2_below,
                               
                
            mdb_flag_limits_list=[SST_max,SST_min,SSS_max,SSS_min,DIC_max, DIC_min,  pH_max, pH_min,pCO2_max, pCO2_min, TA_max,TA_min]  
            
            for idx, g in enumerate(mdb_flag_index_list):
                print(idx, g)
                if len(g) == 0:
                    print("list is empty")
                else:
                    #this prints a header for what the entries have been flagged for
                    with open('output/mdb_flag.csv','a') as flag_result_file:
                        wr = csv.writer(flag_result_file, dialect='excel')
                        wr.writerow([mdb_flag_warnings_list[idx]]+ [mdb_flag_limits_list[idx]])
                    #this prints the values that have been flagged to the csv    
                    print(matchupData.loc[g].to_csv("output/mdb_flag.csv",mode='a'));
        
            #now filter those mdb from the analysis
            mdb_ind_rmv =np.concatenate(mdb_flag_index_list) #combine all the numpy arrays into a single array of 'bad data point indexes'
            #new matchupdatabase filtered for bad data
            matchupData = matchupData.drop(matchupData.index[mdb_ind_rmv])
                    
            ### For each region, run all the algorithms that are relevant for that region
            for region in settings["regions"]:
                print("Running for region:", region);
                logger.info("Beginning new region: "+region);
                algorithmOutputList = [];
                 
                ######################################
                ### For each algorithm in the region, run the algorithm using relevent matchup database rows and store model output
    
                for AlgorithmFunctor in settings["algorithmRegionMapping"][region]:
                    
                    algorithm = AlgorithmFunctor(settings); #Create an instance of the algorithm functor, this is the handle used to access the algorithm
                    
                    try:
                        modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed,subsetData = \
                                  run_algorithm(algorithm, matchupData,
                                  regionMaskPath=settings["regionMasksPath"], region=region,
                                  useDepthMask=settings["subsetWithDepthMask"],
                                  depthMaskPath=settings["depthMaskPath"], depthMaskVar=settings["depthMaskVar"],
                                  useDistToCoastMask=settings["subsetWithDistToCoast"],
                                  distToCoastMaskPath=settings["distToCoastMaskPath"], distToCoastMaskVar=settings["distToCoastMaskVar"]);
                        algorithmOutput = {};
                        algorithmOutput["instance"] = algorithm;
                        algorithmOutput["name"] = algorithm.__str__();
                        algorithmOutput["outputVar"] = algorithm.output_name();
                        algorithmOutput["modelOutput"] = modelOutput;
                        algorithmOutput["propagatedInputUncertainty"] = propagatedInputUncertainty;
                        algorithmOutput["rmsd"] = rmsd;
                        algorithmOutput["combinedUncertainty"] = combinedUncertainty;
                        algorithmOutput["dataUsedIndices"] = dataUsedIndices;

                        logger.info("Output calculated (region:"+region+", algo: "+algorithm.__class__.__name__+")");
                    except ValueError as e: #Raised if there are no matchup data rows left after spatial mask and algorithm internal subsetting has taken place
                        print(algorithm, e.args[0]);
                        logger.info("No matchup data left after subsettings (region:"+region+", algo: "+algorithm.__class__.__name__+")");
                        continue;
                    
                    
                    #filter here the algorithms below our n threshold
                    if len(algorithmOutput["dataUsedIndices"])<30:
                        print("number of matchups less than n. Dont include this algorithm in paired calculations")
                    else:
                        #Store the data/objects which we'll used later
                        algorithmOutputList.append(algorithmOutput);
  
                    ######################
                    ### Diagnostic plots #
                    #Diagnostic plot to compare model and reference output
                    if diagnosticPlots == True:
                        outputVariable = algorithm.output_name();
                        if path.exists(path.join(currentCombinationOutputDirectory, algorithmOutput["outputVar"], region, "diagnostic_plots")) == False:
                            makedirs(path.join(currentCombinationOutputDirectory, algorithmOutput["outputVar"], region, "diagnostic_plots"));
                        savePath = path.join(currentCombinationOutputDirectory, algorithmOutput["outputVar"], region, "diagnostic_plots",  algorithmOutput["name"].split(":")[0]+".png");
                        prediction_accuracy_plot(dataUsed[outputVariable], modelOutput, algorithm.__class__.__name__, outputVariable, savePath=savePath);
                    #increase iteration
                
                #######################################
                ### Calculate metrics for each region #
                if len(algorithmOutputList) == 0:
                    print("*** WARNING: No algorithms were executed for region", region, " possibly because there was no matchup data in this region.Or the number of matchups n was less then specified in settings file.");
                    logger.error("No algorithms were executed for region '"+region+"' possibly because there was no matchup data in this region. Or the number of matchups n was less then specified in settings file.");
                    continue; #Skip calculating metrics for this region
                
                #Only compare algorithms which predict the same output variable, so we need a list of all the output variables
                uniqueOutputVars = np.unique([algorithmOutput["instance"].output_name() for algorithmOutput in algorithmOutputList]); #Get a list of all the output variables
                for currentOutputVar in uniqueOutputVars:
                    #filter algorithm list, output list etc. by the current output variable (in order to compare like-for-like)
                    algosMatchingCurrentOutputVar = [i for i in range(len(algorithmOutputList)) if algorithmOutputList[i]["instance"].output_name()==currentOutputVar];
                    currentAlgorithmOutputs = [algorithmOutputList[i] for i in algosMatchingCurrentOutputVar];
                    
                    
                    basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores = \
                        metrics.calc_all_metrics(currentAlgorithmOutputs, matchupData, settings,currentOutputVar);
                    

                    #x is the string of the dictionary to add to
                    #add final scores output to the dictionary
                    x="final_scores_allcombos_{0}_{1}".format(region,currentOutputVar)
                    eval(x)[inputCombinationName]=finalScores


                    ###########################
                    ### Write outputs to file #
                    outputDirectory = path.join(currentCombinationOutputDirectory, currentOutputVar, region);
                    print("Writing metrics to: ", outputDirectory);
                    write_metrics_to_file(outputDirectory, matchupData, basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores, currentAlgorithmOutputs);

            earthobs_dataset_iteration_number=earthobs_dataset_iteration_number+1;
    
    #Calculate summary table for weighted metrics and output to file
    n_threshold=0
    summaryTable_weighted_no_n_req = create_summary_table(n_threshold,settings, specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps, useWeighted=True);
    summaryTableOutputPath_no_n_req = path.join(settings["outputPathMetrics"], "summary_best_algos_no_n_req.csv");
    summaryTable_weighted_no_n_req.to_csv(summaryTableOutputPath_no_n_req, sep=",", index=False);     
    print("Full weighted summary table written to:", path.abspath(summaryTableOutputPath_no_n_req));
   
    #Calculate summary table for weighted metrics and output to file
    n_threshold=30
    summaryTable_weighted = create_summary_table(n_threshold,settings, specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps, useWeighted=True);
    summaryTableOutputPath = path.join(settings["outputPathMetrics"], "summary_best_algos.csv");
    summaryTable_weighted.to_csv(summaryTableOutputPath, sep=",", index=False);     
    print("Full weighted summary table written to:", path.abspath(summaryTableOutputPath));
   
    #Calculate summary table for unweighted metrics and output to file
    n_threshold=0
    summaryTable_unweighted_no_n_req = create_summary_table(n_threshold, settings, specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps,useWeighted=False);
    summaryTableOutputPathUnweighted_no_n_req = path.join(settings["outputPathMetrics"], "summary_best_algos_unweighted_no_n_req.csv");
    summaryTable_unweighted_no_n_req.to_csv(summaryTableOutputPathUnweighted_no_n_req, sep=",", index=False);     
    print("Full unweighted summary table written to:", path.abspath(summaryTableOutputPathUnweighted_no_n_req));
    
    #Calculate summary table for unweighted metrics and output to file
    n_threshold=30
    summaryTable_unweighted = create_summary_table(n_threshold,settings, specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps, useWeighted=False);
    summaryTableOutputPathUnweighted = path.join(settings["outputPathMetrics"], "summary_best_algos_unweighted.csv");
    summaryTable_unweighted.to_csv(summaryTableOutputPathUnweighted, sep=",", index=False);     
    print("Full unweighted summary table written to:", path.abspath(summaryTableOutputPathUnweighted));
    
    ##### For weighted metrics
    minMatchupsTArange = 0;
    minMatchupsDICrange = 0;
    minYearRange = 0;
    ### Calculate overall best algorithms / input combinations for each region
    overallBestAlgos = get_best_algorithms(summaryTableOutputPath_no_n_req, settings["regions"], minTimeSpanYears=0, minMatchupsTA = 0, minMatchupsDIC=0);
    overallBestAlgosOutputPath = path.join(settings["outputPathMetrics"], "overall_best_algos.csv");
    overallBestAlgos.to_csv(overallBestAlgosOutputPath, sep=",", index=False);     
    print("Overall best algorithm table written to:", path.abspath(overallBestAlgosOutputPath));
    
    ### Calculate overall best algorithms / input combinations for each region again, but this time ensuring a minimum time series range and number of matchup points.
    minMatchupsTArange = 30;
    minMatchupsDICrange = 30;
    minYearRange = 8;
    overallBestAlgosMinRange = get_best_algorithms(summaryTableOutputPath, settings["regions"], minTimeSpanYears=minYearRange, minMatchupsTA = minMatchupsTArange, minMatchupsDIC=minMatchupsDICrange);
    overallBestAlgosOutputPath = path.join(settings["outputPathMetrics"], "overall_best_algos_min_years="+str(minYearRange)+".csv");
    overallBestAlgosMinRange.to_csv(overallBestAlgosOutputPath, sep=",", index=False);     
    print("Overall best algorithm (with min year range="+str(minYearRange)+") table written to:", path.abspath(overallBestAlgosOutputPath));
    
    ##### For unweighted metrics
    ### Calculate overall best algorithms / input combinations for each egion
    minMatchupsTArange = 0;
    minMatchupsDICrange = 0;
    minYearRange = 0;
    overallBestAlgosUnweighted = get_best_algorithms(summaryTableOutputPathUnweighted_no_n_req, settings["regions"], minTimeSpanYears=0, minMatchupsTA = 0, minMatchupsDIC=0);
    overallBestAlgosOutputPathUnweighted = path.join(settings["outputPathMetrics"], "overall_best_algos_unweighted.csv");
    overallBestAlgosUnweighted.to_csv(overallBestAlgosOutputPathUnweighted, sep=",", index=False);     
    print("Overall best unweighted algorithm table written to:", path.abspath(overallBestAlgosOutputPathUnweighted));
    
    ### Calculate overall best algorithms / input combinations for each region again, but this time ensuring a minimum time series range and number of matchup points.
    minMatchupsTArange = 30;
    minMatchupsDICrange = 30;
    minYearRange = 8;
    overallBestAlgosUnweightedMinRange = get_best_algorithms(summaryTableOutputPathUnweighted, settings["regions"], minTimeSpanYears=minYearRange, minMatchupsTA = minMatchupsTArange, minMatchupsDIC=minMatchupsDICrange);
    overallBestAlgosOutputPathUnweighted = path.join(settings["outputPathMetrics"], "overall_best_algos_unweighted_min_years="+str(minYearRange)+".csv");
    overallBestAlgosUnweightedMinRange.to_csv(overallBestAlgosOutputPathUnweighted, sep=",", index=False);     
    print("Overall best unweighted algorithm (with min year range="+str(minYearRange)+") table written to:", path.abspath(overallBestAlgosOutputPathUnweighted));
    
    
    #Store these dictionaries at pickle txt files for plotting
    #They contain the final scores data as tables
    import pickle
    with open("output/algo_metrics/Amazon_AT_finalscores_pickled.txt", "wb") as myFile:
        pickle.dump(final_scores_allcombos_oceansoda_amazon_plume_AT, myFile)
    
    with open("output/algo_metrics/Amazon_DIC_finalscores_pickled.txt", "wb") as myFile2:
        pickle.dump(final_scores_allcombos_oceansoda_amazon_plume_DIC, myFile2) 
        
    with open("output/algo_metrics/Congo_AT_finalscores_pickled.txt", "wb") as myFile3:
        pickle.dump(final_scores_allcombos_oceansoda_congo_AT, myFile3) 
        
    with open("output/algo_metrics/Congo_DIC_finalscores_pickled.txt", "wb") as myFile4:
        pickle.dump(final_scores_allcombos_oceansoda_congo_DIC, myFile4)
        
    #Shutdown logger
    loggerFileHandle.close();

if __name__ == "__main__":
    #Load the settings for this run
    settings = osoda_global_settings.get_default_settings(); #TODO: Hannah - Change this to call your settings function
    
    main(settings);












