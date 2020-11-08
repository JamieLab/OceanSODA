#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:39:41 2020

@author: tom holding
"""


from string import Template;
from os import path;
import os;
from netCDF4 import Dataset;
from datetime import datetime;
import pandas as pd;
import numpy as np;

import osoda_global_settings;


#Returns the file handle to an open netCDF file suitable for storing the gridded output predictions and associated inputs
#   outputPath: specific file path to the .nc file that will be created
#   algorithmFunctorsToSupport: a list of algorithm classes that you want to create a file to store the output in
#   latRes, lonRes, years: describe the dimension sizes of the file
def create_gridded_timeseries_output_netCDF_file(outputPath, algoInfo, datasetInfoMap, latRes, lonRes, years, includeSeaCarbVariables=True):
    ### Create netCDF file to store output
    ncout = Dataset(outputPath, 'w');
    ncout.createDimension("lat", 180/latRes);
    ncout.createDimension("lon", 360/lonRes);
    ncout.createDimension("time", len(years)*12);
    
    #write attributes / best algorithm information
    ncout.algorithmName = algoInfo["algo_name"];
    ncout.algorithmRMSDe = algoInfo["RMSDe"];
    ncout.sampleSize = algoInfo["n"];
    ncout.input_combination_name = algoInfo["input_combination"];
    ncout.num_algos_compared = algoInfo["algos_compared"];
    ncout.region = algoInfo["region"];
    
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
    variableList = list(datasetInfoMap.keys()) + [algoInfo["output_var"]+"_pred"];
    outputVar = algoInfo["output_var"];
    algoName = algoInfo["algo_name"];
    
    var = ncout.createVariable(outputVar+"_pred", float, ("time", "lat", "lon"));
    var.units = "umol kg-1";
    var.long_name = outputVar+" predicted by "+algoName+" ("+var.units+")";
    var.RMSDe = algoInfo["RMSDe"];
    
    var = ncout.createVariable(outputVar+"_pred_uncertainty_due_to_algorithm", float, ("time", "lat", "lon"));
    var.units = "umol kg-1";
    var.long_name = "uncertainty in "+outputVar+" originating from the algorithm ('"+algoName+"') uncertainty";
    
    var = ncout.createVariable(outputVar+"_pred_uncertainty_due_to_input_uncertainty", float, ("time", "lat", "lon"));
    var.units = "umol kg-1";
    var.long_name = "uncertainty in "+outputVar+" originating from the combined input data uncertainty";
    
    var = ncout.createVariable(outputVar+"_pred_combined_uncertainty", float, ("time", "lat", "lon"));
    var.units = "umol kg-1";
    var.long_name = "combined uncertainty in "+outputVar+" (i.e. combining input data uncertainty with algorithm uncertainty for algorithm '"+algoName+"')";
    
    #always add salinity and sst because this will be used by seacarb
    var = ncout.createVariable("SSS", float, ("time", "lat", "lon"));
    var.units = "PSU";
    var.long_name = "Sea surface salinity used for prediction of "+outputVar;
    if datasetInfoMap["SSS"].predictionDatasetError is not None:
        var = ncout.createVariable("SSS_err", float, ("time", "lat", "lon"));
        var.units = "PSU";
        var.long_name = "Uncertainty associated with the sea surface salinity used for prediction of "+outputVar;
    
    var = ncout.createVariable("SST", float, ("time", "lat", "lon"));
    var.units = "Kelvin (k)";
    var.long_name = "Sea surface temperature used for prediction of "+outputVar;
    if datasetInfoMap["SST"].predictionDatasetError is not None:
        var = ncout.createVariable("SST_err", float, ("time", "lat", "lon"));
        var.units = "Kelvin (k)";
        var.long_name = "Uncertainty associated with the sea surface temperature used for prediction of "+outputVar;
    
    #Only add these fields if the algorithm uses them
    if "DO" in variableList:
        var = ncout.createVariable("DO", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Dissolved oxygen used to predict "+outputVar;
        if datasetInfoMap["DO"].predictionDatasetError is not None:
            var = ncout.createVariable("DO_err", float, ("time", "lat", "lon"));
            var.units = "umol kg-1";
            var.long_name = "Uncertainty associated with the dissolved oxygen used for prediction of "+outputVar;
    
    if "NO3" in variableList:
        var = ncout.createVariable("NO3", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Nitrate concentration used to predict "+outputVar;
        if datasetInfoMap["NO3"].predictionDatasetError is not None:
            var = ncout.createVariable("NO3_err", float, ("time", "lat", "lon"));
            var.units = "umol kg-1";
            var.long_name = "Uncertainty associated with the nitrate concentration used for prediction of "+outputVar;
    
    if "PO4" in variableList:
        var = ncout.createVariable("PO4", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Phosphate concentration used to predict "+outputVar;
        if datasetInfoMap["PO4"].predictionDatasetError is not None:
            var = ncout.createVariable("PO4_err", float, ("time", "lat", "lon"));
            var.units = "umol kg-1";
            var.long_name = "Uncertainty associated with the phosphate concentration used for prediction of "+outputVar;
    
    if "SiO4" in variableList:
        var = ncout.createVariable("SiO4", float, ("time", "lat", "lon"));
        var.units = "umol kg-1";
        var.long_name = "Silicate concentration used to predict "+outputVar;
        if datasetInfoMap["SiO4"].predictionDatasetError is not None:
            var = ncout.createVariable("SiO4_err", float, ("time", "lat", "lon"));
            var.units = "umol kg-1";
            var.long_name = "Uncertainty associated with the silicate concentration used for prediction of "+outputVar;
    
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
def load_input_data(datasetInfoMap, year, monthStr, verbose=False):
    #Internal function which extracts 2D gridded data from 2, 3 and 4 dimensional netCDF files (always assuming the last two dimensions are lat and lon)
    def do_variable_extraction(datasetNC, datasetVariable):
        if len(datasetNC.variables[datasetVariable].shape)==2: #i.e.lat, lon
            extractedVariable = datasetNC.variables[datasetVariable][:];
        elif len(datasetNC.variables[datasetVariable].shape)==3: #i.e. time, lat, lon
            extractedVariable = datasetNC.variables[datasetVariable][0,:,:];
        elif len(datasetNC.variables[datasetVariable].shape)==4: #i.e. depth, time, lat, lont
            extractedVariable = datasetNC.variables[datasetVariable][0,0,:,:];
        else:
            raise ValueError("Unsupported number of dimensions ("+datasetNC.variables[datasetVariable].shape+") in input variable: "+inputVariable+" "+datasetInfo.datasetName);
        return extractedVariable;
    
    loadedInputData = {};
    #Use SST input to get the resolution and lon/lat values
    try:
        sstnc = Dataset(datasetInfoMap["SST"].predictionDatasetTemplate.safe_substitute(YYYY=year, MM=monthStr));
        #sst = sstnc.variables[datasetInfoMap["SST"].predictionDatasetVariable][0,:,:];
        sst = do_variable_extraction(sstnc, datasetInfoMap["SST"].predictionDatasetVariable);
        if np.nanmean(sst) < 200.0: #Convert SST from C to K, if required
            sst[sst.mask==False] += 273.15;
        loadedInputData["SST"] = sst;
        sstLons = sstnc.variables["lon"][:];
        sstLats = sstnc.variables["lat"][:];
        #If an error variable is specified, also load this
        if datasetInfoMap["SST"].predictionDatasetError is not None:
            loadedInputData["SST_err"] = do_variable_extraction(sstnc, datasetInfoMap["SST"].predictionDatasetError);
    except FileNotFoundError as e:
        if verbose:
            print("Skipping month", year, monthStr, "because prediction data for SST was not found:", e.args);
        return None, None, None;
    
    #Now that the dataframe is setup add all the other inputs, as required
    try:
        for commonName in datasetInfoMap:
            datasetInfo = datasetInfoMap[commonName];
            inputVariable = datasetInfo.commonName;
            if inputVariable == "SST": #Already added SST when creating the data frame
                continue;
            
            griddedInputNC = Dataset(datasetInfo.predictionDatasetTemplate.safe_substitute(YYYY=year, MM=monthStr));
            griddedInput = do_variable_extraction(griddedInputNC, datasetInfo.predictionDatasetVariable);
            loadedInputData[inputVariable] = griddedInput; #Store for later use
            
            #If an error variable is specified, also load this
            if datasetInfo.predictionDatasetError is not None:
                griddedError = do_variable_extraction(griddedInputNC, datasetInfo.predictionDatasetError);
                loadedInputData[inputVariable+"_err"] = griddedError;
    
    except FileNotFoundError as e:
        if verbose:
            print("Skipping month", year, monthStr, "because prediction data for", datasetInfo.datasetName, "was not found:", e.args);
        return None, None, None;
    
    return loadedInputData, sstLats, sstLons;


#Given a dictionary of gridded data matrices, on a lon lat grid, apply an algorithm and return gridded predicted output
def calculate_gridded_output_from_inputs(AlgorithmClass, inputVariables, lon, lat, curDate, settings):
    #Convert gridded SST data into a dataframe, and then just add more columns for the other inputs
    df = pd.DataFrame(inputVariables["SST"]).stack(dropna=False);
    df = df.rename_axis(["ilat", "ilon"]).reset_index(name='SST');
    df["lon"] = lon[df["ilon"]];
    df["lat"] = lat[df["ilat"]];
    df["date"] = [curDate]*len(df);
    
    algorithm = AlgorithmClass(settings);
    
    #convert each of the other input matrices into a pandas dataframe which can be used with the algorithm.
    for variableName in inputVariables.keys():
        if variableName in algorithm.input_names()+[name+"_err" for name in algorithm.input_names()]:
            if variableName == "SST": #SST is always already added
                continue;
            
            #Append as a column in the dataframe
            varDF = pd.DataFrame(inputVariables[variableName]).stack(dropna=False); #temporary dataframe to convert from gridded to table format
            varDF = varDF.rename_axis(["ilat", "ilon"]).reset_index(name=variableName); #Keep the same indices as the main dataframe
            df[variableName] = varDF[variableName]; #Add to the main dataframe (which will be used by the algorithm)
    
    try:
        algorithmOutputTuple = algorithm(df, predict=True);
        modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed = algorithmOutputTuple;
        df[algorithm.output_name()+"_pred"] = modelOutput;
        df[algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"] = pd.Series(propagatedInputUncertainty, index=modelOutput.index);
        df[algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"] = pd.Series(rmsd, index=modelOutput.index);
        df[algorithm.output_name()+"_pred_combined_uncertainty"] = pd.Series(combinedUncertainty, index=modelOutput.index);
    except ValueError as e:
        print("No data within valid ranges for "+algorithm.__class__.__name__+". No predictions could be made.\n");
        print(e);
        #Fill with nans
        df[algorithm.output_name()+"_pred"] = [np.nan]*len(df);
        df[algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"] = [np.nan]*len(df);
        df[algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"] = [np.nan]*len(df);
        df[algorithm.output_name()+"_pred_combined_uncertainty"] = [np.nan]*len(df);
    
    griddedOutput = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred"); #unstack into a grid again
    griddedOutput = np.array(griddedOutput);
    
    griddedRMSD = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"); #unstack into a grid again
    griddedRMSD = np.array(griddedRMSD);
    griddedInputUncertainty = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"); #unstack into a grid again
    griddedInputUncertainty = np.array(griddedInputUncertainty);
    griddedCombinedUncertainty = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_combined_uncertainty"); #unstack into a grid again
    griddedCombinedUncertainty = np.array(griddedCombinedUncertainty);
    
    #Only return data with both model output and combined uncertainty 
    inconsistentCells = np.where(np.isnan(griddedOutput) | np.isnan(griddedCombinedUncertainty));
    griddedOutput[inconsistentCells] = np.nan;
    griddedRMSD[inconsistentCells] = np.nan;
    griddedInputUncertainty[inconsistentCells] = np.nan;
    griddedCombinedUncertainty[inconsistentCells] = np.nan;
    
    return griddedOutput, griddedRMSD, griddedInputUncertainty, griddedCombinedUncertainty, df; #return the gridded output and the dataframe that was used by the algorithm to make the predictions


#Calculates other carbonate parameters from AT and DIC using the SeaCarb package for R
#Requires numpy objects as inputs
#See for details of SeaCarb's 'carb' function:
#https://cran.r-project.org/web/packages/seacarb/seacarb.pdf
def calculate_carbonate_parameters(tflag, atData, dicData, sssData, sstData, pAtm, pHydrostatic, k1k2="x"):
    import rpy2.robjects as ro;
    import rpy2.robjects.numpy2ri; #Convert between R objects and numpy objects
    import rpy2.robjects.pandas2ri; #Convert between R objects and pandas objects
    import rpy2;
    from rpy2.robjects.packages import importr;
    
    #Get a handle to the SeaCarb R library
    try:
        seacarb = importr("seacarb");
    except rpy2.rinterface.RRuntimeError: #If the library isn't installed, install it for the rpy2 version of r
        print("Installing R package: seacarb");
        utils = importr('utils');
        utils.install_packages('seacarb', repos='https://cloud.r-project.org');
        seacarb = importr("seacarb");
    
    rpy2.robjects.numpy2ri.activate();
    output = seacarb.carb(tflag, atData, dicData, S=35.0, T=20.0, Patm=pAtm, P=pHydrostatic, k1k2=k1k2);
    #output = seacarb.carb(15, at, dic, S=35.0, T=20.0, Patm=1.0, Pt=0.0, k1k2="x");
    
    #No longer needed in rpy2 v3.3.6
    #output = rpy2.robjects.pandas2ri.ri2py(output);
    output = pd.DataFrame(output); #convert form recarray to pandas dataframe
    
#    from rpy2.robjects.conversion import localconverter;
#    with localconverter(ro.default_converter + pandas2ri.converter):
#        tmp = ro.conversion.ri2py(output)
#    
#    output = pd.DataFrame(output);
    
    return output;


#Converts a pandas DataFrame into a list of 2D numpy arrays (gridded spatially according to lonlatIndices)
def convert_dataframe_to_gridded_list(df, latColName, lonColName):
    output = {};
    for variable in df.keys():
        if variable not in [latColName, lonColName]:
            griddedOutput = df.pivot(index=latColName, columns=lonColName, values=variable); #unstack into a grid again
            output[variable] = np.array(griddedOutput);
    return output;


#Writes predicted algorithm output and input data used to make predictions for a single year and month to an already open netCDF file
def write_yearmonth_to_netCDF(ncFileHandle, iyear, imonth, outputVar, griddedModelOutput, griddedRMSD, griddedInputUncertainty, griddedCombinedUncertainty, loadedInputData, carbonateParameters):
    if griddedModelOutput is not None: ncFileHandle.variables[outputVar+"_pred"][(iyear*12)+imonth, :, :] = griddedModelOutput;
    if griddedRMSD is not None: ncFileHandle.variables[outputVar+"_pred_uncertainty_due_to_algorithm"][(iyear*12)+imonth, :, :] = griddedRMSD;
    if griddedInputUncertainty is not None: ncFileHandle.variables[outputVar+"_pred_uncertainty_due_to_input_uncertainty"][(iyear*12)+imonth, :, :] = griddedInputUncertainty;
    if griddedCombinedUncertainty is not None: ncFileHandle.variables[outputVar+"_pred_combined_uncertainty"][(iyear*12)+imonth, :, :] = griddedCombinedUncertainty;
    
    for variableName in loadedInputData.keys():
        #loadedInputData[variableName][np.where(np.isfinite(griddedOutput)==False)] = np.nan; #Where there is no predicted output variable, set the input to nan for consistency.
        ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = loadedInputData[variableName];
    
    #Write carbonate parameters, if supplied
    if carbonateParameters is not None:
        for variableName in carbonateParameters:
            if variableName in ncFileHandle.variables.keys():
                ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
    
    return;


#Return a dictionary mapping common variable names to DatasetInfo objects, for a given algorithm and input combination
def get_combination_dataset_info(settings, algoObj, inputCombinationName, alwaysRequired=["SSS", "SST"]):
    datasetInfoMapAll = settings["datasetInfoMap"];
    datasetInfoMapAlgo = {}
    
    requiredInputs = np.unique(alwaysRequired + algoObj.input_names());
    for commonName in requiredInputs:
        if isinstance(datasetInfoMapAll[commonName], list):
            for algoInfo in datasetInfoMapAll[commonName]:
                if algoInfo.datasetName in inputCombinationName.split("_"):
                    datasetInfoMapAlgo[commonName] = algoInfo;
                    break;
        else:
            datasetInfoMapAlgo[commonName] = datasetInfoMapAll[commonName];
    
    return datasetInfoMapAlgo;

def calculate_gridded_timeseries_single_region(outputPathAT, outputPathDIC, atAlgo, atAlgoInfo, dicAlgo, dicAlgoInfo, settings, years, latRes=1.0, lonRes=1.0, verbose=False):
    if (atAlgo is None) and (dicAlgo is None):
        return;
    
    if atAlgo is not None:
        datasetInfoMapAT = get_combination_dataset_info(settings, atAlgo, atAlgoInfo["input_combination"]);
        #Create a netCDF file to store the output
        ncoutAT = create_gridded_timeseries_output_netCDF_file(outputPathAT, atAlgoInfo, datasetInfoMapAT, latRes, lonRes, years);
    if dicAlgo is not None:
        datasetInfoMapDIC = get_combination_dataset_info(settings, dicAlgo, dicAlgoInfo["input_combination"]);
        #Create a netCDF file to store the output
        ncoutDIC = create_gridded_timeseries_output_netCDF_file(outputPathDIC, dicAlgoInfo, datasetInfoMapDIC, latRes, lonRes, years);
    
    #For each month, run the prediction algorithms and append outputs (and inputs) to the output netCDF file
    for iyear, year in enumerate(years):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02");
            if verbose:
                print("Beginning year month:", year, monthStr);
            
            #read input netCDF files and create a DataFrame which can be used by the algorithms
            if atAlgo is not None:
                loadedInputDataAT, inputLatsAT, inputLonsAT = load_input_data(datasetInfoMapAT, year, monthStr);
            if dicAlgo is not None:
                loadedInputDataDIC, inputLatsDIC, inputLonsDIC = load_input_data(datasetInfoMapDIC, year, monthStr);
            if (atAlgo is not None) and (dicAlgo is not None):
                if (loadedInputDataAT is None) & (loadedInputDataDIC is None): #No input data for this year/month so move to the next
                    continue;
            curDate = pd.to_datetime(datetime(year, imonth+1, 1));
            griddedOutputAT = griddedOutputDIC = None;
            if atAlgo is not None:
                if loadedInputDataAT != None: #If there was missing input data, move to the next month
                    #run algorithm to get the predicted gridded outputs
                    if atAlgo is not None:
                        griddedOutputAT, griddedRMSDAT, griddedInputUncertaintyAT, griddedCombinedUncertaintyAT, dfUsedByAlgorithmAT = calculate_gridded_output_from_inputs(atAlgo, loadedInputDataAT, inputLonsAT, inputLatsAT, curDate, settings);
            if dicAlgo is not None:
                if loadedInputDataDIC != None: #If there was missing input data, move to the next month
                    #run algorithm to get the predicted gridded outputs
                    if dicAlgo is not None:
                        griddedOutputDIC, griddedRMSDDIC, griddedInputUncertaintyDIC, griddedCombinedUncertaintyDIC, dfUsedByAlgorithmDIC = calculate_gridded_output_from_inputs(dicAlgo, loadedInputDataDIC, inputLonsDIC, inputLatsDIC, curDate, settings);
            #TODO: Worst coding ever. Refactor to reduce repetition in handling of None
            
            #########################################
            ### calculate other carbonate parameters
            carbonateParameters = None;
            if (griddedOutputAT is not None) & (griddedOutputDIC is not None):
                if verbose:
                    print("Calculate carbonate parameters using SeaCarb...");
                flag = 15; #This is the SeaCarb flag for calculating carbonate system from AT and DIC
                at = dfUsedByAlgorithmAT["AT_pred"].values/1000.0; #AT values converted to mol kg-1
                dic = dfUsedByAlgorithmDIC["DIC_pred"].values/1000.0; #DIC values converted to mol kg-1
                sss = dfUsedByAlgorithmAT["SSS"].values; #SSS values
                sst = dfUsedByAlgorithmAT["SST"].values - 273.15; #SST values converted to C
                pAtm = 1.0; #Pressure at sea surface, in atm
                pHydrostatic = 0.0; #Hydrostatic pressure, 0=surface
                carbonateParameters = calculate_carbonate_parameters(flag, at, dic, sssData=sss, sstData=sst, pAtm=pAtm, pHydrostatic=pHydrostatic, k1k2="x");
                carbonateParameters["lat"] = dfUsedByAlgorithmAT["lat"]; #Add lon lat information
                carbonateParameters["lon"] = dfUsedByAlgorithmAT["lon"]; #Add lon lat information
                carbonateParameters = convert_dataframe_to_gridded_list(carbonateParameters, "lat", "lon"); #convert to a list of 2D matrices ready for writing to netCDF
            
            #write predicted output for this year/month to the netCDF file
            if griddedOutputAT is not None:
                write_yearmonth_to_netCDF(ncoutAT, iyear, imonth, "AT", griddedOutputAT, griddedRMSDAT, griddedInputUncertaintyAT, griddedCombinedUncertaintyAT, loadedInputDataAT, carbonateParameters);
            if griddedOutputDIC is not None:
                write_yearmonth_to_netCDF(ncoutDIC, iyear, imonth, "DIC", griddedOutputDIC, griddedRMSDDIC, griddedInputUncertaintyDIC, griddedCombinedUncertaintyDIC, loadedInputDataDIC, carbonateParameters);
        
    #All months and years have been computed so close the netCDF file.
    if atAlgo is not None:
        ncoutAT.close();
    if dicAlgo is not None:
        ncoutDIC.close();
    
    
def calculate_gridded_timeseries_all_regions(tablePath, outputPathTemplate, years, regions=None):
    settings = osoda_global_settings.get_default_settings();    
    if regions is None:
        regions = settings["regions"];
    
    bestAlgoTable = pd.read_csv(tablePath);
    
    latRes = lonRes = 1.0;
    
    if isinstance(outputPathTemplate, str):
        outputPathTemplate = Template(outputPathTemplate);
    
    #make gridded time series predictions for each region, using the best input combination and algorithms
    for region in regions:
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
        griddedPredictionOutputPathAT = outputPathTemplate.safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="AT");
        griddedPredictionOutputPathDIC = outputPathTemplate.safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="DIC");
        #griddedPredictionOutputPathAT = settings["griddedPredictionOutputTemplate"].safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="AT");
        #griddedPredictionOutputPathDIC = settings["griddedPredictionOutputTemplate"].safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="DIC");
        if path.exists(path.dirname(griddedPredictionOutputPathAT)) == False:
            os.makedirs(path.dirname(griddedPredictionOutputPathAT));
        #calculate the gridded time series and write to file for this input / region combination
        calculate_gridded_timeseries_single_region(griddedPredictionOutputPathAT, griddedPredictionOutputPathDIC, atAlgo, atAlgoInfo, dicAlgo, dicAlgoInfo, settings, years, latRes, lonRes, verbose=True);



