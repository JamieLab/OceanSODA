#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:39:41 2020

@author: tom holding
"""


#from string import Template;
from os import path;
import os;
from netCDF4 import Dataset;
from datetime import datetime;
import pandas as pd;
import numpy as np;
#import matplotlib.pyplot as plt;

import osoda_global_settings;
import utilities;



#Returns the file handle to an open netCDF file suitable for storing the gridded output predictions and associated inputs
#   outputPath: specific file path to the .nc file that will be created
#   algorithmFunctorsToSupport: a list of algorithm classes that you want to create a file to store the output in
#   latRes, lonRes, years: describe the dimension sizes of the file
def create_gridded_timeseries_output_netCDF_file(outputPath, algorithmAT, algorithmDIC, latRes, lonRes, years, includeSeaCarbVariables=True):
    ### Create netCDF file to store output
    ncout = Dataset(outputPath, 'w');
    ncout.createDimension("lat", 180/latRes);
    ncout.createDimension("lon", 360/lonRes);
    ncout.createDimension("time", len(years)*12);
    ncout.algorithmNameAT = algorithmAT.__class__.__name__;
    ncout.algorithmNameDIC = algorithmDIC.__class__.__name__;
    
    #dimension variables
    var = ncout.createVariable("lat", float, ("lat",));
    var.units = "lat (degrees North)";
    var[:] = np.arange(-90, 90, latRes)+(0.5*latRes);
    var = ncout.createVariable("lon", float, ("lon",));
    var.units = "lon (degrees East)";
    var[:] = np.arange(-180, 180, lonRes)+(0.5*lonRes);
    var = ncout.createVariable("time", int, ("time",));
    var.units = "seconds since 1980-01-01";
    var[:] = [int((datetime(year, imonth+1, 1)-datetime(1980, 1, 1)).total_seconds()) for year in years for imonth in range(0, 12)];

    
    #create output/predicted variables, and find which input variables will be needed
    inputVariablesToCreate = [];
    for algorithm in [algorithmAT, algorithmDIC]:
        if algorithm is not None:
            var = ncout.createVariable(algorithm.output_name()+"_pred", float, ("time", "lat", "lon"));
            var.units = "umol kg-1";
            var.long_name = algorithm.output_name()+" predicted by "+type(algorithm).__name__+" ("+var.units+")";
            inputVariablesToCreate += algorithm.input_names();
    
    inputVariablesToCreate = np.unique(inputVariablesToCreate+["SSS", "SST"]); #Always add SSS and SST as these are used by seacarb
    
    #always add salinity and sst because this will be used by seacarb
    var = ncout.createVariable("SSS", float, ("time", "lat", "lon"));
    var.units = "PSU";
    var.long_name = "Sea surface salinity used for prediction";
    
    var = ncout.createVariable("SST", float, ("time", "lat", "lon"));
    var.units = "Kelvin (k)";
    var.long_name = "Sea surface temperature used for prediction";
    
    #Only add these fields if the algorithm uses them
    if "DO" in algorithm.input_names():
        var = ncout.createVariable("DO", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Dissolved oxygen";
    
    if "NO3" in algorithm.input_names():
        var = ncout.createVariable("NO3", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Nitrate concentration";
    
    if "PO4" in algorithm.input_names():
        var = ncout.createVariable("PO4", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Phosphate concentration";
    
    if "SiO4" in algorithm.input_names():
        var = ncout.createVariable("SiO4", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Silicate concentration";
    
    if includeSeaCarbVariables:
        var = ncout.createVariable("pH", float, ("time", "lat", "lon"));
        var.units = "pH";
        var.long_name = "pH calculated from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("CO2", float, ("time", "lat", "lon"));
        var.units = "mol kg-1";
        var.long_name = "CO2 concentration calculated from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("fCO2", float, ("time", "lat", "lon"));
        var.units = "uatm";
        var.long_name = "'standard' CO2 fugacity computed from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("HCO3", float, ("time", "lat", "lon"));
        var.units = "mol kg-1";
        var.long_name = "HCO3 concentration calculated from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("CO3", float, ("time", "lat", "lon"));
        var.units = "mol kg-1";
        var.long_name = "CO3 concentration calculated from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("OmegaAragonite", float, ("time", "lat", "lon"));
        var.units = "none";
        var.long_name = "Aragonite saturation state calculated from AT and DIC using SeaCarb";
        
        var = ncout.createVariable("OmegaCalcite", float, ("time", "lat", "lon"));
        var.units = "none";
        var.long_name = "Calcite saturation state calculated from AT and DIC using SeaCarb";
    return ncout;


#loads gridded input data for a given year and month
#returns a dictionary of variableName:griddedData, and arrays for latitude and longitude
#returns (None, None, None) if any input data is missing
def load_input_data(inputDataPathInfo, year, monthStr, verbose=False):
    loadedInputData = {};
    #Use SST input to get the resolution and lon/lat values
    try:
        sstnc = Dataset(inputDataPathInfo["SST"][1].safe_substitute(YYYY=year, MM=monthStr));
        sst = sstnc.variables[inputDataPathInfo["SST"][0]][0,:,:];
        if np.nanmean(sst) < 200.0: #Convert SST from C to K, if required
            sst[sst.mask==False] += 273.15;
        loadedInputData["SST"] = sst;
        sstLons = sstnc.variables["lon"][:];
        sstLats = sstnc.variables["lat"][:];
    except FileNotFoundError as e:
        if verbose:
            print("Skipping month", year, monthStr, "because prediction data for SST was not found:", e.args);
        return None, None, None;
    
    #Now that the dataframe is setup add all the other inputs, as required
    try:
        for inputVariable in inputDataPathInfo.keys():
            if inputVariable == "SST": #Already added SST when creating the data frame
                continue;
                
            griddedInputNC = Dataset(inputDataPathInfo[inputVariable][1].safe_substitute(YYYY=year, MM=monthStr));
            if len(griddedInputNC.variables[inputDataPathInfo[inputVariable][0]].shape)==2: #i.e.lat, lon
                griddedInput = griddedInputNC.variables[inputDataPathInfo[inputVariable][0]][:];
            elif len(griddedInputNC.variables[inputDataPathInfo[inputVariable][0]].shape)==3: #i.e. time, lat, lon
                griddedInput = griddedInputNC.variables[inputDataPathInfo[inputVariable][0]][0,:,:];
            elif len(griddedInputNC.variables[inputDataPathInfo[inputVariable][0]].shape)==4: #i.e. depth, time, lat, lont
                griddedInput = griddedInputNC.variables[inputDataPathInfo[inputVariable][0]][0,0,:,:];
            else:
                raise ValueError("Unsupported number of dimensions ("+griddedInputNC.variables[inputDataPathInfo[inputVariable][0]].shape+") in input variable: "+inputVariable);
                
            loadedInputData[inputVariable] = griddedInput; #Store for later use
    except FileNotFoundError as e:
        if verbose:
            print("Skipping month", year, monthStr, "because prediction data for", inputVariable, "was not found:", e.args);
        return None, None, None;
    
    return loadedInputData, sstLats, sstLons;


#Given a dictionary of gridded data matrices, on a lon lat grid, apply an algorithm and return gridded predicted output
def calculate_gridded_output_from_inputs(algorithm, inputVariables, lon, lat, curDate):
    #Convert gridded SST data into a dataframe, and then just add more columns for the other inputs
    df = pd.DataFrame(inputVariables["SST"]).stack(dropna=False);
    df = df.rename_axis(["ilat", "ilon"]).reset_index(name='SST');
    df["lon"] = lon[df["ilon"]];
    df["lat"] = lat[df["ilat"]];
    df["date"] = [curDate]*len(df);
    
    #convert each of the other input matrices into a pandas dataframe which can be used with the algorithm.
    for variableName in inputVariables.keys():
        if variableName in algorithm.input_names():
            if variableName == "SST": #SST is always already added
                continue;
            #Append as a column in the dataframe
            varDF = pd.DataFrame(inputVariables[variableName]).stack(dropna=False); #temporary dataframe to convert from gridded to table format
            varDF = varDF.rename_axis(["ilat", "ilon"]).reset_index(name=variableName); #Keep the same indices as the main dataframe
            df[variableName] = varDF[variableName]; #Add to the main dataframe (which will be used by the algorithm)
    
    try:
        modelOutput, dataUsed = algorithm(df, predict=True);
        df[algorithm.output_name()+"_pred"] = modelOutput;
    except ValueError:
        print("No data within valid ranges for "+algorithm.__class__.__name__+". No predictions could be made.");
        df[algorithm.output_name()+"_pred"] = [np.nan]*len(df);
    
    griddedOutput = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred"); #unstack into a grid again
    griddedOutput = np.array(griddedOutput);
    
    return griddedOutput, df; #return the gridded output and the dataframe that was used by the algorithm to make the predictions

#Writes predicted algorithm output and input data used to make predictions for a single year and month to an already open netCDF file
def write_yearmonth_to_netCDF(ncFileHandle, iyear, imonth, griddedOutputAT, griddedOutputDIC, loadedInputData, carbonateParameters):
    if griddedOutputAT is not None: ncFileHandle.variables["AT_pred"][(iyear*12)+imonth, :, :] = griddedOutputAT;
    if griddedOutputDIC is not None: ncFileHandle.variables["DIC_pred"][(iyear*12)+imonth, :, :] = griddedOutputDIC;
    
    for variableName in loadedInputData.keys():
        #loadedInputData[variableName][np.where(np.isfinite(griddedOutput)==False)] = np.nan; #Where there is no predicted output variable, set the input to nan for consistency.
        ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = loadedInputData[variableName];
    
    #Write carbonate parameters, if supplied
    if carbonateParameters is not None:
        for variableName in carbonateParameters:
            if variableName in ncFileHandle.variables.keys():
                ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
    
    return;

def calculate_gridded_timeseries_driver(outputPath, algorithmAT, algorithmDIC, years, inputDataPathInfo, latRes=1.0, lonRes=1.0, verbose=False):
    #requiredInputs = np.unique(["SSS", "SST"]+algorithmAT.input_names()+algorithmDIC.input_names()); #Always require SSS and SST
    
    #Create a netCDF file to store the output
    if (algorithmAT is None) and (algorithmDIC is None):
        return;
    
    ncout = create_gridded_timeseries_output_netCDF_file(outputPath, algorithmAT, algorithmDIC, latRes, lonRes, years);
    
    #For each month, run the prediction algorithms and append outputs (and inputs) to the output netCDF file
    for iyear, year in enumerate(years):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02");
            if verbose:
                print("Beginning year month:", year, monthStr);
            
            #read input netCDF files and create a DataFrame which can be used by the algorithms
            loadedInputData, inputLats, inputLons = load_input_data(inputDataPathInfo, year, monthStr);
            if loadedInputData == None: #If there was missing input data, move to the next month
                continue;
            
            #run algorithm to get the predicted gridded outputs
            griddedOutputAT = griddedOutputDIC = None;
            curDate = pd.to_datetime(datetime(year, imonth+1, 1));
            if algorithmAT is not None:
                griddedOutputAT, dfUsedByAlgorithmAT = calculate_gridded_output_from_inputs(algorithmAT, loadedInputData, inputLons, inputLats, curDate);
            if algorithmDIC is not None:
                griddedOutputDIC, dfUsedByAlgorithmDIC = calculate_gridded_output_from_inputs(algorithmDIC, loadedInputData, inputLons, inputLats, curDate);
            
            #########################################
            ### calculate other carbonate parameters
            if verbose:
                print("Calculate carbonate parameters using SeaCarb...");
            flag = 15; #This is the SeaCarb flag for calculating carbonate system from AT and DIC
            at = dfUsedByAlgorithmAT["AT_pred"].values/1000.0; #AT values converted to mol kg-1
            dic = dfUsedByAlgorithmDIC["DIC_pred"].values/1000.0; #DIC values converted to mol kg-1
            sss = dfUsedByAlgorithmAT["SSS"].values; #SSS values
            sst = dfUsedByAlgorithmAT["SST"].values - 273.15; #SST values converted to C
            pAtm = 1.0; #Pressure at sea surface, in atm
            pHydrostatic = 0.0; #Hydrostatic pressure, 0=surface
            carbonateParameters = utilities.calculate_carbonate_parameters(flag, at, dic, sssData=sss, sstData=sst, pAtm=pAtm, pHydrostatic=pHydrostatic, k1k2="x");
            carbonateParameters["lat"] = dfUsedByAlgorithmAT["lat"]; #Add lon lat information
            carbonateParameters["lon"] = dfUsedByAlgorithmAT["lon"]; #Add lon lat information
            carbonateParameters = utilities.convert_dataframe_to_gridded_list(carbonateParameters, "lat", "lon"); #convert to a list of 2D matrices ready for writing to netCDF
            
            #write predicted output for this year/month to the netCDF file
            write_yearmonth_to_netCDF(ncout, iyear, imonth, griddedOutputAT, griddedOutputDIC, loadedInputData, carbonateParameters);
        
    #All months and years have been computed so close the netCDF file.
    ncout.close();



if __name__ == "__main__":
    settings = osoda_global_settings.get_default_settings();
    inputDataPathInfo = settings["predictionDataPaths"];
    lonRes = latRes = 1.0;
    years = settings["years"];
    
    #For each input data combination, extract the best DIC and AT algorithm for each region
    for inputCombination in utilities.get_dataset_variable_map_combinations(settings)[1]: #just the combination names
        currentMetricsRootDirectory = path.join(settings["outputPathMetrics"], inputCombination);
        for region in settings["regions"]:
            print("Starting", inputCombination, region);
            #Get an instance of the 'best' AT and DIC algorithms for this region/input combo.
            #This is done by comparing the best algorithm name to each name of the algorithms used in the current region
            bestAlgorithms = utilities.find_best_algorithm(currentMetricsRootDirectory, region, useWeightedRMSDe=settings["assessUsingWeightedRMSDe"], verbose=False)
            if bestAlgorithms["AT"] is not None:
                AlgorithmFunctorAT = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == bestAlgorithms["AT"][0]][0];
                bestAlgorithmAT = AlgorithmFunctorAT(settings);
            else:
                bestAlgorithmAT = None;
            if bestAlgorithms["DIC"] is not None:
                AlgorithmFunctorDIC = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == bestAlgorithms["DIC"][0]][0];
                bestAlgorithmDIC = AlgorithmFunctorDIC(settings);
            else:
                bestAlgorithmDIC = None;
            
            #Create output file path
            griddedPredictionOutputPath = settings["griddedPredictionOutputTemplate"].safe_substitute(INPUTCOMBINATION=inputCombination, REGION=region, LATRES=latRes, LONRES=lonRes);
            if path.exists(path.dirname(griddedPredictionOutputPath)) == False:
                os.makedirs(path.dirname(griddedPredictionOutputPath));

            #find the union of the inputs required for these algorithms
            inputsRequired = ["SST", "SSS"]; #Always include SSS and SST
            if bestAlgorithmAT is not None:
                inputsRequired += bestAlgorithmAT.input_names();
            if bestAlgorithmDIC is not None:
                inputsRequired += bestAlgorithmDIC.input_names();
            inputsRequired = np.unique(inputsRequired);
            inputDataPathInfo = {inputName:settings["predictionDataPaths"][inputName] for inputName in inputsRequired};
            
            #calculate the gridded time series and write to file for this input / region combination
            calculate_gridded_timeseries_driver(griddedPredictionOutputPath, bestAlgorithmAT, bestAlgorithmDIC, years, inputDataPathInfo, latRes, lonRes, verbose=True);






