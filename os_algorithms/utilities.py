#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:49:55 2020

Misc. utilities used for loading and processing data
Includes:
    Functions for subsetting the matchup database in different ways
    Wrapper for seacarb
    Functions

@author: tom holding
"""

from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;
from os import path, makedirs;
from string import Template;


#Creates a list of dictionaries containing variable to database mappings for every unique combination of input variables
#Also generates unique names for each combination and returns a list of these
def get_dataset_variable_map_combinations2(settings):
    datasetInfoMap = settings["datasetInfoMap"];
    commonVariableNames = datasetInfoMap.keys();

    #Make sure all the mappings are lists, even if there is only one possibility.
    for variableName in commonVariableNames:
        if isinstance(datasetInfoMap[variableName], list) == False:
            datasetInfoMap[variableName] = [datasetInfoMap[variableName]];

    #Build a new variable name to database name maps which is are specific to each possible combination
    #Create a list combination names and a dictionary mapping common variable name to datasetInfo objects for each
    combinations = [{}];
    labels = ["_"];
    for variableName in commonVariableNames:
        newCombinations = [];
        newLabels = [];
        for j, datasetInfo in enumerate(datasetInfoMap[variableName]):
            for i, combination in enumerate(combinations): #For each combination of the previous iteration (not yet including the current variable name), append a combination with also includes this dataset
                newCombination = combination.copy();
                newCombination[variableName] = datasetInfo;
                newCombinations.append(newCombination);

                #If this variable is one we're trying different combinations of (e.g. more than one dataset listed for it), then include it in the name of the combination
                if len(datasetInfoMap[variableName]) > 1:
                    newLabel = labels[i] + "_" + datasetInfo.datasetName; #New label/name is the previous label plus the current dataset name appended to it
                else: #Not trying different combinations of this variable, so it doesn't need to be added to the combination name
                    newLabel = labels[i];
                newLabels.append(newLabel);
        combinations = newCombinations;
        labels = newLabels;

    #append "combination" to the start of each label to give the final combination name
    labels = ["combination"+str(i)+labels[i] for i in range(len(labels))];

    #Return a list of specific variable to database mappings
    return combinations, labels;


#Creates a list of dictionaries containing variable to database mappings for every unique combination of input variables
#Also generates unique names for each combination and returns a list of these
def get_dataset_variable_map_combinations(settings):
    alwaysIncludeInLabel = ["SSS", "SST"]; #Dataset names for these variables will always be included in the combination label, even if there is just one of them.

    datasetInfoMap = settings["datasetInfoMap"];
    commonVariableNames = datasetInfoMap.keys();

    #Make sure all the mappings are lists, even if there is only one possibility.
    for variableName in commonVariableNames:
        if isinstance(datasetInfoMap[variableName], list) == False:
            datasetInfoMap[variableName] = [datasetInfoMap[variableName]];

    #Build a new variable name to database name maps which is are specific to each possible combination
    #Create a list combination names and a dictionary mapping common variable name to datasetInfo objects for each
    combinations = [{}];
    labels = ["_"];
    for variableName in commonVariableNames:
        newCombinations = [];
        newLabels = [];
        for j, datasetInfo in enumerate(datasetInfoMap[variableName]):
            for i, combination in enumerate(combinations): #For each combination of the previous iteration (not yet including the current variable name), append a combination with also includes this dataset
                newCombination = combination.copy();
                newCombination[variableName] = datasetInfo;
                newCombinations.append(newCombination);

                #If this variable is one we're trying different combinations of (e.g. more than one dataset listed for it), then include it in the name of the combination
                if (len(datasetInfoMap[variableName]) > 1) or (variableName in alwaysIncludeInLabel):
                    newLabel = labels[i] + "_" + datasetInfo.datasetName; #New label/name is the previous label plus the current dataset name appended to it
                else: #Not trying different combinations of this variable, so it doesn't need to be added to the combination name
                    newLabel = labels[i];
                newLabels.append(newLabel);
        combinations = newCombinations;
        labels = newLabels;

    #append "combination" to the start of each label to give the final combination name
    labels = ["combination"+str(i)+labels[i] for i in range(len(labels))];

    #Return a list of specific variable to database mappings
    return combinations, labels;


#Write the specific combination of database inputs to a file
def write_specific_variable_to_database_mapping(specificVariableToDatabaseMap, outputPath, combinationName=None):
    if path.exists(path.dirname(outputPath)) == False:
        makedirs(path.dirname(outputPath));
    with open(outputPath, 'w') as file:
        if combinationName==None:
            combinationName = "Unnamed combination\n\n";
        file.write(combinationName+"\n");
        for variableName in specificVariableToDatabaseMap.keys():
            file.write(variableName+":"+specificVariableToDatabaseMap[variableName].datasetName+"\n");

##prints the matchup database inputs that correspond to each combination name code
#def print_combination_name_keys(settings):
#    settingsMapAll = settings["variableToDatabaseMap"];
#    variableNames = settingsMapAll.keys();
#
#    #Make sure all the mappings are lists, even if there is only one possibility.
#    for variableName in variableNames:
#        if isinstance(settingsMapAll[variableName], list) == False:
#            settingsMapAll[variableName] = [settingsMapAll[variableName]];
#
#        if len(settingsMapAll[variableName]) > 1:
#            for i, datasetName in enumerate(settingsMapAll[variableName]):
#                print(variableName+str(i)+":\t"+datasetName);


#Given a set of inputs, this will return a list of years/time points for which the matchup dataset contains all of these inputs
def calculate_years_for_input_combination(settings, inputCombination, minYear=1900, maxYear=2050):
    years = [];

    for year in range(minYear, maxYear+1):
        filePath = settings["matchupDatasetTemplate"].safe_substitute(YYYY=year);
        if path.exists(filePath):
            try:
                nc = Dataset(filePath, 'r');
            # except FileNotFoundError:
            #     continue; #No file for this, continue as before.
            except Exception as e:
                print("SOMETHING WENT WRONG!");
                print(e);

            #Make sure all input variables exist in this file
            allExist = True;
            for key, datasetInfo in inputCombination.items():
                if datasetInfo.matchupVariableName not in nc.variables.keys():
                    allExist = False;
                    break;
            nc.close();

            #include the current year if all input variables exist
            if allExist == True:
                years.append(year);

    return years;


#Reads yearly matchup database netCDF files and concatinates them into a single dataframe.
#   years: and iterable of integer years
#   settings: the global settings dictionary
#   datasetInfoMap: dictionary mapping common variable names to DatasetInfo objects which contain all the information needed to load a dataset
def load_matchup_to_dataframe(settings, datasetInfoMap, years=None, commonNames=None):
    #Concatinate each year
    dfList = [];
    print(years)
    if commonNames == None:
        commonNames = datasetInfoMap.keys()
    if years == None:
        years = settings["years"];

    for year in years:
        # We need to add a check that the matchup_file exists
        matchup_file = settings["matchupDatasetTemplate"].safe_substitute(YYYY=year)
        if path.exists(matchup_file):
            print(matchup_file + ' - File exists...')
            matchupNC = Dataset(matchup_file, 'r');

            #Create a pandas dataframe from the netCDF file
            df = pd.DataFrame();
            for commonName in commonNames:
                try:
                    df[commonName] = matchupNC[datasetInfoMap[commonName].matchupVariableName][:];
                    try: #If there is a missing data value, filter out missing values and replace with nans
                        missingValue = matchupNC[datasetInfoMap[commonName].matchupVariableName]._FillValue;
                        if np.isnan(missingValue) == False:
                            df.loc[df[commonName]==missingValue, commonName] = np.nan;
                    except:
                        pass;

                except IndexError:
                    print("Missing data: ", year, commonName, datasetInfoMap[commonName].datasetName, datasetInfoMap[commonName].matchupVariableName);
                if datasetInfoMap[commonName].matchupDatabaseError is not None:
                    try:
                        df[commonName+"_err"] = matchupNC[datasetInfoMap[commonName].matchupDatabaseError][:];
                        try: #If there is a missing data value, filter out missing values and replace with nans
                            missingValue = matchupNC[datasetInfoMap[commonName].matchupDatabaseError]._FillValue;
                            if np.isnan(missingValue) == False:
                                df.loc[df[commonName+"_err"]==missingValue, commonName+"_err"] = np.nan;
                        except:
                            pass;
                    except IndexError:
                        print("Missing uncertainty data: ", year, commonName, datasetInfoMap[commonName].datasetName, datasetInfoMap[commonName].matchupDatabaseError);
            #Convert any C temperature units to K
            #if 'SST' in df.values:




            dfList.append(df);
            matchupNC.close();
    #print(dfList)

    matchupData = pd.concat(dfList, ignore_index=True);
    # We move this out of the loop above as for some years we have no data for the SST dataset (i.e CCI-SST in 1957 doesn't exist).
    # This means we extract all the data and then check :-)

    if np.nanmean(matchupData["SST"]) < 200.0: #Convert SST from C to K, if required
        matchupData["SST"]=matchupData["SST"]+ 273.15 ;

    #print(list(matchupData.columns.values))
    #Convert date from time in seconds since 1980-01-01 to a pd.datetime object
    matchupData["date"] = convert_time_to_date(matchupData["date"]);
    return matchupData;


#returns data frame containing data from the matchup database for each variable in 'cols'
#Differs from load_matchup_to_dataframe in that it doesn't require DatasetInfo objects, and therefore doesn't try to load uncertainty variables automatically
def read_matchup_cols(matchupTemplate, cols, years):

    dfList = [];
    for year in years:
        print(year)
        if path.isfile(matchupTemplate.safe_substitute(YYYY=year)):
            matchupNC = Dataset(matchupTemplate.safe_substitute(YYYY=year), 'r');
            df = pd.DataFrame();
            for col in cols:
                try:
                    if col == "date": #some variables have different names in the matchup. This is hacky, but manually matches the commonName and matchupdatabase names. This is probably ok, as these are unlikely to change. Should really look into the variable:matchup mapping in the global settings file
                        df[col] = matchupNC.variables["time"][:];
                        try: #If there is a missing data value, filter out missing values and replace with nans
                            missingValue = matchupNC.variables["time"]._FillValue;
                            if np.isnan(missingValue) == False:
                                df.loc[df["time"]==missingValue, "time"] = np.nan;
                        except:
                            pass;
                    elif col == "AT":
                        df[col] = matchupNC.variables["insitu_ta_mean"][:];
                        try: #If there is a missing data value, filter out missing values and replace with nans
                            missingValue = matchupNC.variables["insitu_ta_mean"]._FillValue;
                            if np.isnan(missingValue) == False:
                                df.loc[df["AT"]==missingValue, "AT"] = np.nan;
                        except:
                            pass;

                    elif col == "DIC":
                        df[col] = matchupNC.variables["insitu_dic_mean"][:];
                        try: #If there is a missing data value, filter out missing values and replace with nans
                            missingValue = matchupNC.variables["insitu_dic_mean"]._FillValue;
                            if np.isnan(missingValue) == False:
                                df.loc[df["DIC"]==missingValue, "DIC"] = np.nan;
                        except:
                            pass;
                    # elif col in ["AT", "DIC"]: #some variables have different names in the matchup. This is hacky, but manually matches the commonName and matchupdatabase names. This is probably ok, as these are unlikely to change. Should really look into the variable:matchup mapping in the global settings file
                    #     df[col] = matchupNC.variables[col+"_mean"][:];
                    else:
                        #list variables and check it exists
                        matchup_var_list=matchupNC.variables.keys()
                        if col in matchup_var_list:
                            print("this will execute")
                            df[col] = matchupNC.variables[col][:];
                            try: #If there is a missing data value, filter out missing values and replace with nans
                                missingValue = matchupNC.variables[col]._FillValue;
                                if np.isnan(missingValue) == False:
                                    df.loc[df[col]==missingValue, col] = np.nan;
                            except:
                                pass;
                        try: #If there is a missing data value, filter out missing values and replace with nans
                            missingValue = matchupNC.variables[col]._FillValue;
                            if np.isnan(missingValue) == False:
                                df.loc[df[col]==missingValue, col] = np.nan;
                        except:
                            pass;

                except IndexError:
                    print("Couldn't find {0} in the matchup database for year {1}. Continuing without this variable.".format(col, year));
            try:
                print(np.nanmax(df[cols[-1]]))
            except:
                pass;
            dfList.append(df);
            matchupNC.close();
        else:
            print('File '+ matchupTemplate.safe_substitute(YYYY=year)+ ' does not exist...')
    #merge dataframes from each year
    data = pd.concat(dfList, ignore_index=True);

    #Convert date from time in seconds since 1980-01-01 to a pd.datetime object
    data["date"] = convert_time_to_date(data["date"]);
    return data;


#Converts time in seconds to data
#   dt: time delta in seconds from the base date (pandas series)
#   baseDate: base date from which the change in time (dt) is from
def convert_time_to_date(dt, baseDate=pd.to_datetime("1950/1/1")):
    dates = baseDate + pd.to_timedelta(dt, 'S');
    return dates;


#Subset a dataframe by including any points which fall inside any boxes formed by a list of lon and lat ranges
def subset_from_inclusive_coord_list(includedLons, includedLats, data):
    if len(includedLons) != len(includedLats):
        raise ValueError("There must be an equal number of lon and lat ranges");

    #If empty, default to global (no regional subsetting needed)
    if len(includedLons) == 0:
        return data;

    toKeep = np.full((len(data,)), False, dtype=bool);
    for i in range(len(includedLons)): #for each pair of lon/lat ranges, keep any datapoints that fall into this range.
        lonMin = includedLons[i][0];
        lonMax = includedLons[i][1]
        latMin = includedLats[i][0];
        latMax = includedLats[i][1];
        toKeep = toKeep | ((data["lon"] >= lonMin) & (data["lon"] <= lonMax) & (data["lat"] >= latMin) & (data["lat"] <= latMax));

    subset = data.loc[toKeep];
    return subset;


#Given a dataframe, returns a subset containing only complete rows.
def subset_complete_rows(df):
    selectedRows = np.all(df.notna(), axis=1)
    subset = df.loc[selectedRows];
    return subset;


#subset from a gridded netCDF mask. Makes no assumptions as to shape of mask
def subset_from_mask(data, maskNC, maskName, maskValue=1):
    #Read the mask and lon/lat dimension values
    mask = maskNC.variables[maskName][:];
    lats = maskNC.variables["lat"][:];
    lons = maskNC.variables["lon"][:];

    #find index of closest mask lats value to each row
    a, b = np.meshgrid(lats, data["lat"]);
    latIndices = np.abs(a-b).argmin(axis=1);
    a, b = np.meshgrid(lons, data["lon"]); #and again for longitudes
    lonIndices = np.abs(a-b).argmin(axis=1);

    subset = data.iloc[mask[(latIndices, lonIndices)] == maskValue];
    return subset;



#Returns the name of the best AT and DIC algorithm for a particular input combination and region
#Intended to be used after the metrics have been computed for each algorithm and region
#metricsRootDirectory should include the input combination directory, if using input combinations
def find_best_algorithm(n_threshold,metricsRootDirectory, region, outputVars=["AT", "DIC"], useWeightedRMSDe=True, verbose=False):
    finalScoresTemplatePath = Template(path.join(metricsRootDirectory, "${OUTPUTVAR}/${REGION}/final_scores.csv"));
    rmsdeCol="final_wrmsd" if useWeightedRMSDe else "final_rmsd";

    bestAlgorithms = {}; #Store the names of the best algorithms for each output variable
    for outputVar in outputVars:
        finalScoresPath = finalScoresTemplatePath.safe_substitute(OUTPUTVAR=outputVar, REGION=region);
        try:
            finalScores = pd.read_csv(finalScoresPath);
            # below threshold
            bool_remove=(finalScores["n"] >= n_threshold);
            finalScores=finalScores.loc[bool_remove,];
            finalScores.reset_index(drop=True, inplace=True);
        except:
            if verbose:
                print("No output file found at:", finalScoresPath);
            bestAlgorithms[outputVar] = None;
            continue;



        if np.all(np.isfinite(finalScores[rmsdeCol])==False):
            if verbose:
                print("*** All NaN encountered in finalScores.csv", rmsdeCol, "row at", finalScoresPath);
                print("    \tThis probably means no pairwise weighted metrics for this region could be calculated (e.g. because there were no spatially overlapping algorithms or no algorithms reported their RMSD.");
            bestAlgorithms[outputVar] = None;
            continue;

        #Now we know there is at least one non-NaN value, find the best algorithm and store its name
        ibestAlgo = np.nanargmin(finalScores[rmsdeCol]);

        bestAlgobias = finalScores["wbias"][ibestAlgo];
        bestAlgo_unc_end_end = finalScores["unc_end_end"][ibestAlgo];
        bestAlgo_RMSD = finalScores["final_rmsd"][ibestAlgo];
        bestAlgo_wRMSD = finalScores["final_wrmsd"][ibestAlgo];

        bestAlgoName = finalScores["algorithm"][ibestAlgo];
        numAlgosCompared = sum(finalScores[rmsdeCol].isna()==False);
        bestAlgorithms[outputVar] = (bestAlgoName, finalScores[rmsdeCol][ibestAlgo], numAlgosCompared,bestAlgobias,bestAlgo_unc_end_end,bestAlgo_RMSD,bestAlgo_wRMSD); #store tuple of algorithm name and selected RMSDe

        if verbose:
            print("Best algorithm:", outputVar, region, bestAlgoName);

    return bestAlgorithms;
