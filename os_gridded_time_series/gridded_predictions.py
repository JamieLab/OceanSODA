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
import PyCO2SYS as pyco2
import osoda_global_settings;


#Returns the file handle to an open netCDF file suitable for storing the gridded output predictions and associated inputs
#   outputPath: specific file path to the .nc file that will be created
#   algorithmFunctorsToSupport: a list of algorithm classes that you want to create a file to store the output in
#   latRes, lonRes, years: describe the dimension sizes of the file
def create_gridded_timeseries_output_netCDF_file(outputPath, algoInfo, datasetInfoMap, latRes, lonRes, years, includepyco2sysVariables=True):
    ### Create netCDF file to store output
    
    # create dataset and provide dimensions
    ncout = Dataset(outputPath, 'w');
    ncout.createDimension("lat", 180/latRes);
    ncout.createDimension("lon", 360/lonRes);
    ncout.createDimension("time", len(years)*12);

    #create output/predicted variables, and find which input variables will be needed
    variableList = list(datasetInfoMap.keys()) + [algoInfo["output_var"]];
    algoName = algoInfo["algo_name"];
    outputVar = algoInfo["output_var"];

    # variable is AT but should be TA in all the files
    if outputVar == "AT":
        outputVar="TA"
    # Need long names defining the variables in the netCDFs
    if outputVar == "TA":
        outputVar_longname="total alkalinity (TA)"
        outputVar_other="DIC"
        outputVar_other_longname="dissolved inorganic carbon(DIC)"
    elif outputVar == "DIC":
        outputVar_longname="dissolved inorganic carbon(DIC)"
        outputVar_other="TA"
        outputVar_other_longname="total alkalinity (TA)"
    
    #meta data that will require updating for future revisions
    CodeAccessTime="February 2022"
    DatasetAccessTime="November 2021" # note that all the datasets nuts and SST/SSS should be downloaded together with get_datasets script
    PO4_ref="WOA volume 4 (Garcia et al., 2013)"
    SiO4_ref="WOA volume 4 (Garcia et al., 2013)"
    NO3_ref="WOA volume 4 (Garcia et al., 2013)"
    DO_ref="WOA volume 3 (Garcia et al., 2013)"

    #need to get the references for the SST and SSS datasets, this is messy but works
    if algoInfo["input_combination"]== "combination0__SST-ESACCI_SSS-ESACCI":
        SST_ref="European Space Agency Climate Change Initiative (ESACCI) SST v2.1 (Merchant et al., 2019;Good et al., 2019)"
        SSS_ref="European Space Agency Climate Change Initiative (ESACCI) SSS v2.31(Boutin et al., 2021;Boutin et al., 2020)"
    
    elif algoInfo["input_combination"]== "combination1__SST-CORA_SSS-ESACCI":
        SST_ref="Coriolis Ocean database for ReAnalysis (CORA) v5.2 (Szekely et al., 2019)"
        SSS_ref="European Space Agency Climate Change Initiative (ESACCI) SSS v2.31(Boutin et al., 2021;Boutin et al., 2020)"
    
    elif algoInfo["input_combination"]== "combination2__SST-OISST_SSS-ESACCI":
        SST_ref="Optimum interpolation seas surface temperature (OISST) v2.1 (Huang et al., 2021;Banzon et al., 2016)"
        SSS_ref="European Space Agency Climate Change Initiative (ESACCI) SSS v2.31(Boutin et al., 2021;Boutin et al., 2020)"
                  
    elif algoInfo["input_combination"]== "combination3__SST-ESACCI_SSS-CORA":
        SST_ref="European Space Agency Climate Change Initiative (ESACCI) SST v2.1 (Merchant et al., 2019;Good et al., 2019)"
        SSS_ref="Coriolis Ocean database for ReAnalysis (CORA) SSS v5.2 (Szekely et al., 2019)"
        
    elif algoInfo["input_combination"]== "combination4__SST-CORA_SSS-CORA":
        SST_ref="Coriolis Ocean database for ReAnalysis (CORA) SST v5.2  (Szekely et al., 2019)"
        SSS_ref="Coriolis Ocean database for ReAnalysis (CORA) SSS v5.2 (Szekely et al., 2019)"
        
    elif algoInfo["input_combination"]== "combination5__SST-OISST_SSS-CORA":
        SST_ref="Optimum interpolation seas surface temperature (OISST) v2.1 (Huang et al., 2021;Banzon et al., 2016)"
        SSS_ref="Coriolis Ocean database for ReAnalysis (CORA) SSS v5.2 (Szekely et al., 2019)"   
        
    elif algoInfo["input_combination"]== "combination6__SST-ESACCI_SSS-RSS-SMAP":
        SST_ref="European Space Agency Climate Change Initiative (ESACCI) SST v2.1 (Merchant et al., 2019;Good et al., 2019)"
        SSS_ref="Remote sensing systems - Soil Moisture Active Passive (RSS-SMAP) SSS v4.0 (Meissner et al., 2018;Meissner et al., 2019)"
        
    elif algoInfo["input_combination"]== "combination7__SST-CORA_SSS-RSS-SMAP":
        SST_ref="Coriolis Ocean database for ReAnalysis (CORA) SST v5.2  (Szekely et al., 2019)"
        SSS_ref="Remote sensing systems - Soil Moisture Active Passive (RSS-SMAP) SSS v4.0 (Meissner et al., 2018;Meissner et al., 2019)"
        
    elif algoInfo["input_combination"]== "combination8__SST-OISST_SSS-RSS-SMAP":
        SST_ref="Optimum interpolation seas surface temperature (OISST) v2.1 (Huang et al., 2021;Banzon et al., 2016)"
        SSS_ref="Remote sensing systems - Soil Moisture Active Passive (RSS-SMAP) SSS v4.0 (Meissner et al., 2018;Meissner et al., 2019)"   
        
    elif algoInfo["input_combination"]== "combination9__SST-ESACCI_SSS-ISAS":
        SST_ref="European Space Agency Climate Change Initiative (ESACCI) SST v2.1 (Merchant et al., 2019;Good et al., 2019)"
        SSS_ref="In Situ Analysis System (ISAS) ISAS-15 SSS (Kolodziejczyk et al., 2017;Gaillard et al., 2016)"
        
    elif algoInfo["input_combination"]== "combination10__SST-CORA_SSS-ISAS":
        SST_ref="Coriolis Ocean database for ReAnalysis (CORA) SST v5.2  (Szekely et al., 2019)"
        SSS_ref="In Situ Analysis System (ISAS) ISAS-15 SSS (Kolodziejczyk et al., 2017;Gaillard et al., 2016)"
        
    elif algoInfo["input_combination"]== "combination11__SST-OISST_SSS-ISAS":
        SST_ref="Optimum interpolation seas surface temperature (OISST) v2.1 (Huang et al., 2021;Banzon et al., 2016)"
        SSS_ref="In Situ Analysis System (ISAS) ISAS-15 SSS (Kolodziejczyk et al., 2017;Gaillard et al., 2016)"  
             

    today=datetime.today()
    #write global attributes 
    ncout.Title = "OceanSODA-UNEXE surface carbonate system dataset file for the " +algoInfo["region"] + " using SST and SSS from "+outputVar_longname+"algorithm evaluation"
    ncout.Contact = "For dataset enquiries contact Richard Sims r.sims2@exeter.ac.uk or Jamie Shutler j.d.shutler@exeter.ac.uk" ;
    ncout.Metadatafile = "The metadata file for this dataset is hosted on www.pangaea.de"
    ncout.Reference = "This dataset is described in the journal Earth System Science Data (Sims et.al. 2022)"
    ncout.ProcessingCode = "The processing code can be found at https://github.com/Richard-Sims/OceanSODA this dataset was run on the "+CodeAccessTime+" code, the original code is at https://github.com/JamieLab/OceanSODA"
    ncout.ProjectWebsite = "The project website is https://esa-oceansoda.org/"
    ncout.CarbonSystemEquations = "The carbonate system equations were performed with PyCO2SYS v1.7 (Humphreys et.al 2022) , the website is https://pyco2sys.readthedocs.io/en/latest/"
    ncout.Conventions = "COARDS" ;
    ncout.Filename = "gridded_"+algoInfo["region"] +"_1degx1deg_"+outputVar+".nc"
    ncout.History = "File generated on "+today.strftime("%d/%m/%y %H:%M:%S") ;
    ncout.ProductionDateTime = "File generated on "+today.strftime("%d/%m/%y %H:%M:%S");
    #ncout.ModificationDateTime = "File generated on: Mon Mar 17 16:18:09 2014 GMT" ;
    ncout.VersionID = "1.0" ;
    ncout.Format = "NetCDF-4" ;
    ncout.Grid = "World Geodetic System (WGS84)" ;
    ncout.Delta_Lon = "1.f" ;
    ncout.Delta_Lat = "1.f" ;
    ncout.SpatialCoverage = "Global" ;
    ncout.NLayers = 1;            
    ncout.Start_Date = 19570101 ;
    ncout.Start_Time = "00:00:00.0" ;
    ncout.End_Date = 20211201 ;
    ncout.End_Time = "00:00:00.0" ;
    
    #write best algorithm information
    ncout.DIC_algorithm_evaluation_number_algos_compared = algoInfo["algos_compared"];
    ncout.DIC_algorithm_evaluation_sampleSize_in_matchupdatabase = algoInfo["n"];
    ncout.DIC_algorithm_evaluation_region = algoInfo["region"];
    ncout.DIC_algorithm_Name = algoInfo["algo_name"];
    ncout.DIC_algorithm_SST_dataset_ref = "The dataset reference is "+SST_ref+". The dataset was accessed on "+DatasetAccessTime;
    ncout.DIC_algorithm_SSS_dataset_ref = "The dataset reference is "+SSS_ref+". The dataset was accessed on "+DatasetAccessTime;
    ncout.DIC_algorithm_input_combination_usedincode = algoInfo["input_combination"];
    ncout.DIC_algorithm_RMSDe = algoInfo["RMSDe"];
    ncout.DIC_algorithm_RMSD = algoInfo["RMSD"];
    ncout.DIC_algorithm_Bias=algoInfo["bias"];
    ncout.DIC_algorithm_uncertainty =algoInfo["uncendtoend"];
    
    #dimension variables - (latitude, longitude and time)
    var = ncout.createVariable("lat", float, ("lat",));
    var.long_name = "Latitude";
    var.units = "degrees_north";
    var[:] = np.arange(-90, 90, latRes)+(0.5*latRes);
    
    var = ncout.createVariable("lon", float, ("lon",));
    var.long_name = "Longitude";
    var.units = "degrees_east";
    var[:] = np.arange(-180, 180, lonRes)+(0.5*lonRes);
    
    var = ncout.createVariable("time", int, ("time",));
    var.long_name = "Time";
    var.calendar = "standard" ;
    var.units = "months since 1957-01-01 00:00:00";
    var[:] = [int((datetime(year, imonth+1, 1)-datetime(1957, 1, 1)).total_seconds()) for year in years for imonth in range(0, 12)];

    

    var = ncout.createVariable(outputVar, float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = outputVar_longname+" predicted by "+algoName+"with input variables (SST,SSS,nutrients etc) provided in this netCDF file";
    
    var = ncout.createVariable(outputVar_other, float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = outputVar_other_longname+" extracted from "+ outputVar_other+ " netCDF dataset";
    
    #these uncertainties use the RMSDe and bias generated during the Algo evaluation 
    #to create an end to end uncertainty analysis estimate
    #uncomment this to define extra TA/DIC Uncertainties

    # #RMSDe
    # var = ncout.createVariable(outputVar+"_RMSD", float, ("time", "lat", "lon"), zlib=True);
    # var.units = "umol kg-1";
    # var.long_name = "end to end - uncertainty in "+outputVar+" originating from RMSDe component from algorithm evaluation ('"+algoName+"') uncertainty";

    # #Bias
    # var = ncout.createVariable(outputVar+"_Bias", float, ("time", "lat", "lon"), zlib=True);
    # var.units = "umol kg-1";
    # var.long_name = "end to end - uncertainty in "+outputVar+" originating from bias component from algorithm evaluation ('"+algoName+"') uncertainty";

    #Combined
    var = ncout.createVariable(outputVar+"_uncertainty", float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = "Combined standard uncertainty in "+outputVar+", combining measurement bias and RMSD from algorithm evaluation ('"+algoName+"')";


    var = ncout.createVariable(outputVar_other+"_uncertainty", float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = "Combined standard uncertainty in "+outputVar_other+ " extracted from "+ outputVar_other+ " netCDF dataset";


    #uncomment this to define extra TA/DIC Uncertainties
    # #these uncertainties are for the bottom up estimate
    # var = ncout.createVariable(outputVar+"_pred_uncertainty_due_to_algorithm", float, ("time", "lat", "lon"), zlib=True);
    # var.units = "umol kg-1";
    # var.long_name = "Bottom up - uncertainty in "+outputVar+" originating from the literature algorithm ('"+algoName+"') uncertainty";
    
    # var = ncout.createVariable(outputVar+"_pred_uncertainty_due_to_input_uncertainty", float, ("time", "lat", "lon"), zlib=True);
    # var.units = "umol kg-1";
    # var.long_name = "Bottom up - uncertainty in "+outputVar+" originating from the combined input data uncertainty";
    
    # var = ncout.createVariable(outputVar+"_pred_combined_uncertainty", float, ("time", "lat", "lon"), zlib=True);
    # var.units = "umol kg-1";
    # var.long_name = "Bottom up - combined uncertainty in "+outputVar+" (i.e. combining input data uncertainty with literature algorithm uncertainty for algorithm '"+algoName+"')";
    
    
    #always add salinity and sst because this will be used by pyco2sys
    var = ncout.createVariable("SSS", float, ("time", "lat", "lon"), zlib=True);
    var.units = "salinity";
    var.long_name = "Sea surface salinity (SSS), dataset reference is "+SSS_ref+". The dataset was accessed on "+DatasetAccessTime;
    if datasetInfoMap["SSS"].predictionDatasetError is not None:
        var = ncout.createVariable("SSS_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "salinity";
        var.long_name = "Uncertainty in sea surface salinity (SSS), dataset reference is "+SSS_ref+".The dataset was accessed on "+DatasetAccessTime;
    
    var = ncout.createVariable("SST", float, ("time", "lat", "lon"), zlib=True);
    var.units = "kelvin";
    var.long_name = "Sea surface temperature (SST), dataset reference is  "+SST_ref+". The dataset was accessed on "+DatasetAccessTime;
    if datasetInfoMap["SST"].predictionDatasetError is not None:
        var = ncout.createVariable("SST_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "kelvin";
        var.long_name = "Uncertainty in sea surface temperature (SST), dataset reference is "+SST_ref+".The dataset was accessed on "+DatasetAccessTime ;
    
    #Only add these fields if the algorithm uses them
    if "DO" in variableList:
        var = ncout.createVariable("DO", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "Dissolved oxygen (DO). The dataset reference is"+DO_ref+". The dataset was accessed on "+DatasetAccessTime;
        if datasetInfoMap["DO"].predictionDatasetError is not None:
            var = ncout.createVariable("DO_err", float, ("time", "lat", "lon"), zlib=True);
            var.units = "umol kg-1";
            var.long_name = "Uncertainty associated with the dissolved oxygen (DO).The dataset reference is"+DO_ref+". The dataset was accessed on "+DatasetAccessTime;
    
    if "NO3" in variableList:
        var = ncout.createVariable("NO3", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "Nitrate  (NO3) content.The dataset reference is"+NO3_ref+". The dataset was accessed on "+DatasetAccessTime;
        if datasetInfoMap["NO3"].predictionDatasetError is not None:
            var = ncout.createVariable("NO3_err", float, ("time", "lat", "lon"), zlib=True);
            var.units = "umol kg-1";
            var.long_name = "Uncertainty in nitrate (NO3) content. The dataset reference is"+NO3_ref+". The dataset was accessed on "+DatasetAccessTime;
    
    #Always save as used by PyCO2SYS
    var = ncout.createVariable("PO4", float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = "Phosphate (PO4) content. The dataset reference is"+PO4_ref+". The dataset was accessed on "+DatasetAccessTime;
    if datasetInfoMap["PO4"].predictionDatasetError is not None:
        var = ncout.createVariable("PO4_err", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "Uncertainty in phosphate (PO4) content. The dataset reference is"+PO4_ref+". The dataset was accessed on "+DatasetAccessTime;
    
    #Always save as used by PyCO2SYS
    var = ncout.createVariable("SiO4", float, ("time", "lat", "lon"), zlib=True);
    var.units = "umol kg-1";
    var.long_name = "Silicate(SiO4) content. The dataset reference is"+SiO4_ref+". The dataset was accessed on "+DatasetAccessTime;
    if datasetInfoMap["SiO4"].predictionDatasetError is not None:
        var = ncout.createVariable("SiO4_err", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "Uncertainty in silicate (SiO4) content.The dataset reference is"+SiO4_ref+". The dataset was accessed on "+DatasetAccessTime;

    #these are variable names which match up to the pyco2sys outputs
    if includepyco2sysVariables:
        
        var = ncout.createVariable("pH_seawater_scale", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "pH on the seawater scale calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("pH_seawater_scale_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "propogated uncertainty in pH on the seawater scale calculated using pyco2syS";
        
        
        
        
        var = ncout.createVariable("pH_total_scale", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "pH on the total scale calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("pH_total_scale_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "propogated uncertainty in on the seawater scale calculated using pyco2sys";
        
        
        
        
        var = ncout.createVariable("pH_free_scale", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "pH on the free scale calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("pH_free_scale_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "pH";
        var.long_name = "propogated uncertainty in pH on the free scale calculated using pyco2sys";
        
        
        
        
        var = ncout.createVariable("H+", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "hydrogen free ions(H+) calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("H+_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "propogated uncertainty in hydrogen free ions (H+) calculated using pyco2sys";
         
        
        
        
        var = ncout.createVariable("pCO2", float, ("time", "lat", "lon"), zlib=True);
        var.units = "ppm";
        var.long_name = "partial pressure of carbon dioxide (pCO2) calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("pCO2_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "ppm";
        var.long_name = "propogated uncertainty in partial pressure of carbon dioxide (pCO2) calculated using pyco2sys";
        
        
        
        
        var = ncout.createVariable("fCO2", float, ("time", "lat", "lon"), zlib=True);
        var.units = "uatm";
        var.long_name = "fugacity of carbon dioxide (fCO2) calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("fCO2_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "uatm";
        var.long_name = "fugacity of carbon dioxide (fCO2) calculated using pyco2sys";
        
        
        
        
        var = ncout.createVariable("HCO3-", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "bicarbonate ion(HCO3-) calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("HCO3-_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "propogated uncertainty in bicarbonate ion (HCO3-) calculated using pyco2sys.";




        var = ncout.createVariable("CO3-2", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "carbonate ion (CO3-2) calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("CO3-2_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "umol kg-1";
        var.long_name = "propogated uncertainty in the carbonate ion (CO3-2) calculated using pyco2sys.";
        
        
        
        var = ncout.createVariable("omega_aragonite", float, ("time", "lat", "lon"), zlib=True);
        var.units = "omega_aragonite";
        var.long_name = "aragonite saturation state calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("omega_aragonite_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "omega_aragonite";
        var.long_name = "propogated uncertainty in aragonite saturation state calculated using pyco2sys.";
        
        
        
        var = ncout.createVariable("omega_calcite", float, ("time", "lat", "lon"), zlib=True);
        var.units = "omega_calcite";
        var.long_name = "calcite saturation state calculated from TA and DIC using pyco2sys";
        
        var = ncout.createVariable("omega_calcite_uncertainty", float, ("time", "lat", "lon"), zlib=True);
        var.units = "omega_calcite";
        var.long_name = "propogated uncertainty in calcite saturation state calculated using pyco2sys.";
    
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
def calculate_gridded_output_from_inputs(AlgorithmClass,AlgoInfo, inputVariables, lon, lat, curDate, region, regionMaskNC, settings):
    #
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
    
    #Convert gridded SST data into a dataframe, and then just add more columns for the other inputs
    df = pd.DataFrame(inputVariables["SST"]).stack(dropna=False);
    df = df.rename_axis(["ilat", "ilon"]).reset_index(name='SST');
    df["lon"] = lon[df["ilon"]];
    df["lat"] = lat[df["ilat"]];
    df["date"] = [curDate]*len(df);

    algoRMSDe=AlgoInfo["RMSDe"];
    algoRMSD=AlgoInfo["RMSD"];
    algoBias=AlgoInfo["bias"];
    algouncendtoend=AlgoInfo["uncendtoend"];
    algorithm = AlgorithmClass(settings);
    
    #convert each of the other input matrices into a pandas dataframe which can be used with the algorithm.
    for variableName in inputVariables.keys():
        #create a set of variables that need adding. Not that SST is already added when creating the initial data frame
        inputsToAdd = set(algorithm.input_names()+[name+"_err" for name in algorithm.input_names()] + ["SST_err", "SSS", "SSS_err","SiO4","SiO4_err","PO4","PO4_err"]);
        #add variables here
        if variableName in inputsToAdd:
            if variableName == "SST": #SST is always already added
                continue;
            
            #Append as a column in the dataframe
            varDF = pd.DataFrame(inputVariables[variableName]).stack(dropna=False); #temporary dataframe to convert from gridded to table format
            varDF = varDF.rename_axis(["ilat", "ilon"]).reset_index(name=variableName); #Keep the same indices as the main dataframe
            df[variableName] = varDF[variableName]; #Add to the main dataframe (which will be used by the algorithm)
            
    
    
    #apply region mask to input dataframe
    df = subset_from_mask(df, regionMaskNC, region);
    
    #Run the algorithm using the input dataframe (df)
    try:
        algorithmOutputTuple = algorithm(df, predict=True);
        modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed = algorithmOutputTuple;
        df[algorithm.output_name()] = modelOutput;
        
        #uncomment this to define extra TA/DIC Uncertainties

        # df[algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"] = pd.Series(propagatedInputUncertainty, index=modelOutput.index);
        # df[algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"] = pd.Series(rmsd, index=modelOutput.index);
        # df[algorithm.output_name()+"_pred_combined_uncertainty"] = pd.Series(combinedUncertainty, index=modelOutput.index);
    
        # df[algorithm.output_name()+"_RMSD"] = pd.Series(algoRMSD, index=modelOutput.index);
        # df[algorithm.output_name()+"_Bias"] = pd.Series(algoBias, index=modelOutput.index);
        df[algorithm.output_name()+"_uncertainty"] = pd.Series(algouncendtoend, index=modelOutput.index);
    
    
    except ValueError as e:
        print("No data within valid ranges for "+algorithm.__class__.__name__+". No predictions could be made.\n");
        print(e);
        #Fill with nans
        df[algorithm.output_name()] = [np.nan]*len(df);
        
        #uncomment this to define extra TA/DIC Uncertainties

        # df[algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"] = [np.nan]*len(df);
        # df[algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"] = [np.nan]*len(df);
        # df[algorithm.output_name()+"_pred_combined_uncertainty"] = [np.nan]*len(df);
    
        # df[algorithm.output_name()+"_RMSD"] = [np.nan]*len(df);
        # df[algorithm.output_name()+"_Bias"] = [np.nan]*len(df);
        df[algorithm.output_name()+"_uncertainty"] = [np.nan]*len(df);
    
    def centre_df_data_to_grid(df, varCol, latCol="lat", lonCol="lon"):
        #Get min lon and lat used for centring gridded output on 180 by 360 grid
        ilatOffset = int((df[latCol].min() + 89.5));
        ilonOffset = int((df[lonCol].min() + 179.5));
        
        pivoted = df.pivot(index=latCol, columns=lonCol, values=varCol); #unstack into a grid again
        gridded = np.full((180, 360), np.nan);
        gridded[ilatOffset:ilatOffset+pivoted.shape[0], ilonOffset:ilonOffset+pivoted.shape[1]] = np.array(pivoted);
        gridded = np.array(gridded);
        
        return gridded;
        
    
#    pivoted = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred"); #unstack into a grid again
#    griddedOutput = np.full((180, 360), np.nan);
#    griddedOutput[ilatOffset:ilatOffset+pivoted.shape[0], ilonOffset:ilonOffset+pivoted.shape[1]] = np.array(pivoted);
#    griddedOutput = np.array(griddedOutput);
    
    griddedOutput = centre_df_data_to_grid(df, algorithm.output_name());
    # griddedRMSD = centre_df_data_to_grid(df, algorithm.output_name()+"_pred_uncertainty_due_to_algorithm");
    # griddedInputUncertainty = centre_df_data_to_grid(df, algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty");
    # griddedCombinedUncertainty = centre_df_data_to_grid(df, algorithm.output_name()+"_pred_combined_uncertainty");
    
    # griddedRMSD = centre_df_data_to_grid(df, algorithm.output_name()+"_RMSD");
    # griddedBias = centre_df_data_to_grid(df, algorithm.output_name()+"_Bias");
    gridded_uncertainty = centre_df_data_to_grid(df, algorithm.output_name()+"_uncertainty");
    
    
        
    
#    griddedRMSD = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_uncertainty_due_to_algorithm"); #unstack into a grid again
#    griddedRMSD = np.array(griddedRMSD);
#    griddedInputUncertainty = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_uncertainty_due_to_input_uncertainty"); #unstack into a grid again
#    griddedInputUncertainty = np.array(griddedInputUncertainty);
#    griddedCombinedUncertainty = df.pivot(index="lat", columns="lon", values=algorithm.output_name()+"_pred_combined_uncertainty"); #unstack into a grid again
#    griddedCombinedUncertainty = np.array(griddedCombinedUncertainty);
    
    #Only return data with both model output and combined uncertainty 
    inconsistentCells = np.where(np.isnan(griddedOutput) | np.isnan(gridded_uncertainty));
    griddedOutput[inconsistentCells] = np.nan;
    # griddedRMSD[inconsistentCells] = np.nan;
    # griddedInputUncertainty[inconsistentCells] = np.nan;
    # griddedCombinedUncertainty[inconsistentCells] = np.nan;
    # griddedRMSD[inconsistentCells] = np.nan;
    # griddedBias[inconsistentCells] = np.nan;
    gridded_uncertainty[inconsistentCells] = np.nan;
    return griddedOutput, gridded_uncertainty, df; #return the gridded output and the dataframe that was used by the algorithm to make the predictions



#Writes predicted algorithm output and input data used to make predictions for a single year and month to an already open netCDF file
def write_outputvars_to_netCDF(ncFileHandle, iyear, imonth, outputVar,griddedModelOutput, gridded_uncertainty, loadedInputData, carbonateParameters,griddedModelOutput_other, gridded_uncertainty_other):
    
    #this saves either TA or DIC and uncertainties to the dataset
    if outputVar == "AT":
        outputVar="TA"
    
    if outputVar == "TA":
        outputVar_other="DIC"
    elif outputVar == "DIC":
        outputVar_other="TA"
    
    #Save the main variable TA or DIC
    if griddedModelOutput is not None: ncFileHandle.variables[outputVar][(iyear*12)+imonth, :, :] = griddedModelOutput;
    
    #Save the missing carbonate variable DIC or TA
    if griddedModelOutput_other is not None: ncFileHandle.variables[outputVar_other][(iyear*12)+imonth, :, :] = griddedModelOutput_other;

    #uncomment this to define extra TA/DIC Uncertainties

    # if griddedRMSD is not None: ncFileHandle.variables[outputVar+"_pred_uncertainty_due_to_algorithm"][(iyear*12)+imonth, :, :] = griddedRMSD;
    # if griddedInputUncertainty is not None: ncFileHandle.variables[outputVar+"_pred_uncertainty_due_to_input_uncertainty"][(iyear*12)+imonth, :, :] = griddedInputUncertainty;
    # if griddedCombinedUncertainty is not None: ncFileHandle.variables[outputVar+"_pred_combined_uncertainty"][(iyear*12)+imonth, :, :] = griddedCombinedUncertainty;
    
    # if griddedRMSD is not None: ncFileHandle.variables[outputVar+"_RMSD"][(iyear*12)+imonth, :, :] = griddedRMSD;
    # if griddedBias is not None: ncFileHandle.variables[outputVar+"_Bias"][(iyear*12)+imonth, :, :] = griddedBias;
    
    #Save the uncertainty of the main variable TA or DIC
    if gridded_uncertainty is not None: ncFileHandle.variables[outputVar+"_uncertainty"][(iyear*12)+imonth, :, :] = gridded_uncertainty;
    
    #Save the uncertainty of the missing carbonate variable DIC or TA
    if gridded_uncertainty_other is not None: ncFileHandle.variables[outputVar_other+"_uncertainty"][(iyear*12)+imonth, :, :] = gridded_uncertainty_other;


    for variableName in loadedInputData.keys():
        #loadedInputData[variableName][np.where(np.isfinite(griddedOutput)==False)] = np.nan; #Where there is no predicted output variable, set the input to nan for consistency.
        if variableName=="SST_err":
            ncFileHandle.variables["SST_uncertainty"][(iyear*12)+imonth, :, :] = loadedInputData[variableName];
        elif variableName=="SSS_err":
            ncFileHandle.variables["SSS_uncertainty"][(iyear*12)+imonth, :, :] = loadedInputData[variableName];
        else:
            ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = loadedInputData[variableName];

    #Write carbonate parameters, if supplied
    if carbonateParameters is not None:
        for variableName in carbonateParameters:
            # pyco2sys naming convention for uncertainties is a bit short, so these are changed here. 
            if variableName=="u_CO3":
                ncFileHandle.variables["CO3-2_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_fCO2":
                ncFileHandle.variables["fCO2_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_HCO3":
                ncFileHandle.variables["HCO3-_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_hydrogen_free":
                ncFileHandle.variables["H+_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_pCO2":
                ncFileHandle.variables["pCO2_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_pH_sws":
                ncFileHandle.variables["pH_seawater_scale_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_pH_free":
                ncFileHandle.variables["pH_free_scale_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_pH_total":
                ncFileHandle.variables["pH_total_scale_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_saturation_aragonite":
                ncFileHandle.variables["omega_aragonite_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="u_saturation_calcite":
                ncFileHandle.variables["omega_calcite_uncertainty"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="pH_total":
                ncFileHandle.variables["pH_total_scale"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="pH_free":
                ncFileHandle.variables["pH_free_scale"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="pH_sws":
                ncFileHandle.variables["pH_seawater_scale"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="hydrogen_free":
                ncFileHandle.variables["H+"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="HCO3":
                ncFileHandle.variables["HCO3-"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="CO3":
                ncFileHandle.variables["CO3-2"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="saturation_aragonite":
                ncFileHandle.variables["omega_aragonite"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
            elif variableName=="saturation_calcite":
                ncFileHandle.variables["omega_calcite"][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];

            #for all other variable CO2SYS names are output variable names
            elif variableName in ncFileHandle.variables.keys():
                ncFileHandle.variables[variableName][(iyear*12)+imonth, :, :] = carbonateParameters[variableName];
    
    return;


#Return a dictionary mapping common variable names to DatasetInfo objects, for a given algorithm and input combination
#Note silicate and phosphate are always required for PYCO2SYS
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
    
    #,"SiO4","PO4" are always used by PyCO2SYS so add them here everytime
    datasetInfoMapAlgo["SiO4"] = datasetInfoMapAll["SiO4"];
    datasetInfoMapAlgo["PO4"] = datasetInfoMapAll["PO4"];

    return datasetInfoMapAlgo;


#Converts a pandas DataFrame into a list of 2D numpy arrays (gridded spatially according to lonlatIndices)
#Assumes centred grid cells and 1x1 degree spatial resolution
def convert_dataframe_to_gridded_list(df, latColName, lonColName):
    output = {};
    for variable in df.keys():
        if variable not in [latColName, lonColName]:
            griddedOutput = df.pivot(index=latColName, columns=lonColName, values=variable); #unstack into a grid again
            griddedOutput = np.array(griddedOutput);
            
            #place in a 180x360 grid, choosing position based on existing lat values
            #find ilat and ilon index offsets
            ilatOffset = int((df[latColName].min() + 89.5));
            ilonOffset = int((df[lonColName].min() + 179.5));
            centredGriddedOutput = np.full((180,360), np.nan);
            centredGriddedOutput[ilatOffset:ilatOffset+griddedOutput.shape[0], ilonOffset:ilonOffset+griddedOutput.shape[1]] = griddedOutput;
            
            output[variable] = np.array(centredGriddedOutput);
    return output;

def calculate_gridded_timeseries_single_region(outputPathAT, outputPathDIC, atAlgo, atAlgoInfo, dicAlgo, dicAlgoInfo, settings, years, region, regionMaskNC, latRes=1.0, lonRes=1.0, verbose=False):
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
                print("test1:");
                loadedInputDataAT, inputLatsAT, inputLonsAT = load_input_data(datasetInfoMapAT, year, monthStr);
            
            if dicAlgo is not None:
                loadedInputDataDIC, inputLatsDIC, inputLonsDIC = load_input_data(datasetInfoMapDIC, year, monthStr);
        
            if (atAlgo is not None) and (dicAlgo is not None):
                if (loadedInputDataAT is None) & (loadedInputDataDIC is None): #No input data for this year/month so move to the next
                    print("test3:");
                    continue;
            curDate = pd.to_datetime(datetime(year, imonth+1, 1));
            griddedOutputAT = griddedOutputDIC = gridded_uncertaintyAT =gridded_uncertaintyDIC= None;
            if atAlgo is not None:
                if loadedInputDataAT != None: #If there was missing input data, move to the next month
                    #run algorithm to get the predicted gridded outputs
                    print("test4:");
                    if atAlgo is not None:
                        print("test5:");
                        algoinfotofunction=atAlgoInfo
                        griddedOutputAT,gridded_uncertaintyAT, dfUsedByAlgorithmAT = calculate_gridded_output_from_inputs(atAlgo,algoinfotofunction, loadedInputDataAT, inputLonsAT, inputLatsAT, curDate, region, regionMaskNC, settings);
            if dicAlgo is not None:
                if loadedInputDataDIC != None: #If there was missing input data, move to the next month
                    #run algorithm to get the predicted gridded outputs
                    if dicAlgo is not None:
                        algoinfotofunction=dicAlgoInfo
                        griddedOutputDIC, gridded_uncertaintyDIC, dfUsedByAlgorithmDIC = calculate_gridded_output_from_inputs(dicAlgo,algoinfotofunction, loadedInputDataDIC, inputLonsDIC, inputLatsDIC, curDate, region, regionMaskNC, settings);
            #TODO: Worst coding ever. Refactor to reduce repetition in handling of None


            
            #########################################
            ### calculate other carbonate parameters
            carbonateParametersAT = None;
            carbonateParametersDIC = None;
            if (griddedOutputAT is not None) & (griddedOutputDIC is not None):
                if verbose:
                    print("Calculate carbonate parameters using pyco2...");
                
                #CO2SYS on the data - USE TA best temp and sal setup
                kwargs = dict(
                par1 = dfUsedByAlgorithmAT["AT"].values,  # Value of the first parameter
                par2 = dfUsedByAlgorithmDIC["DIC"].values,  # Value of the second parameter
                par1_type = 1,  # The first parameter supplied is of type "1", which is "alkalinity"
                par2_type = 2,  # The second parameter supplied is of type "2", which is "DIC"
                salinity = dfUsedByAlgorithmAT["SSS"].values,  # Salinity of the sample
                temperature = dfUsedByAlgorithmAT["SST"].values - 273.15,  # Temperature at input conditions
                pressure = 1.0,  # Pressure    at input conditions#Pressure at sea surface, in atm
                total_silicate = dfUsedByAlgorithmAT["SiO4"].values,  # Concentration of silicate  in the sample (in umol/kg)
                total_phosphate = dfUsedByAlgorithmAT["PO4"].values,  # Concentration of phosphate in the sample (in umol/kg)
                opt_k_carbonic = 4,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
                opt_k_bisulfate = 1,);  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
                    
                
                carbonateParametersAT = pyco2.sys(**kwargs,
                                                uncertainty_into=["pCO2", "pH","pH_total","pH_sws","pH_free","hydrogen_free","fCO2","CO3","HCO3","saturation_aragonite","saturation_calcite"],
                                                uncertainty_from={"par1": dfUsedByAlgorithmAT["AT_uncertainty"].values,
                                                                  "par2": dfUsedByAlgorithmDIC["DIC_uncertainty"].values,
                                                                  "total_silicate": dfUsedByAlgorithmAT["SiO4_err"].values,
                                                                  "total_phosphate": dfUsedByAlgorithmAT["PO4_err"].values,
                                                                  "temperature": dfUsedByAlgorithmAT["SST_err"].values,
                                                                  "salinity": dfUsedByAlgorithmAT["SSS_err"].values});
                
                #CO2SYS on the data - USE DIC best temp and sal setup
                kwargs2 = dict(
                par1 = dfUsedByAlgorithmAT["AT"].values,  # Value of the first parameter
                par2 = dfUsedByAlgorithmDIC["DIC"].values,  # Value of the second parameter
                par1_type = 1,  # The first parameter supplied is of type "1", which is "alkalinity"
                par2_type = 2,  # The second parameter supplied is of type "2", which is "DIC"
                salinity = dfUsedByAlgorithmDIC["SSS"].values,  # Salinity of the sample
                temperature = dfUsedByAlgorithmDIC["SST"].values - 273.15,  # Temperature at input conditions
                pressure = 1.0,  # Pressure    at input conditions#Pressure at sea surface, in atm
                total_silicate = dfUsedByAlgorithmDIC["SiO4"].values,  # Concentration of silicate  in the sample (in umol/kg)
                total_phosphate = dfUsedByAlgorithmDIC["SiO4"].values,  # Concentration of phosphate in the sample (in umol/kg)
                opt_k_carbonic = 4,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
                opt_k_bisulfate = 1,);  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
                    
                
                carbonateParametersDIC = pyco2.sys(**kwargs2,
                                                uncertainty_into=["pCO2", "pH","pH_total","pH_sws","pH_free","hydrogen_free","fCO2","CO3","HCO3","saturation_aragonite","saturation_calcite"],
                                                uncertainty_from={"par1": dfUsedByAlgorithmAT["AT_uncertainty"].values,
                                                                  "par2": dfUsedByAlgorithmDIC["DIC_uncertainty"].values,
                                                                  "total_silicate": dfUsedByAlgorithmDIC["SiO4_err"].values,
                                                                  "total_phosphate": dfUsedByAlgorithmDIC["PO4_err"].values,
                                                                  "temperature": dfUsedByAlgorithmDIC["SST_err"].values,
                                                                  "salinity": dfUsedByAlgorithmDIC["SSS_err"].values});
                
                


                #this variable in the output is 'auto' which is a string, donta ctually need this
                #so just delete the element
                del carbonateParametersAT["buffers_mode"];
                del carbonateParametersDIC["buffers_mode"];

                #carbonateParametersAT = calculate_carbonate_parameters(flag, at[wNoNanAT], dic[wNoNanAT], sssData=sss_at[wNoNanAT], sstData=sst_at[wNoNanAT], atErr=atErr[wNoNanAT], dicErr=dicErr[wNoNanAT], SErr=sssErr_at[wNoNanAT], TErr=sstErr_at[wNoNanAT], pAtm=pAtm, pHydrostatic=pHydrostatic, k1k2="x");
                carbonateParametersAT["lat"] = dfUsedByAlgorithmAT["lat"].values; #Add lon lat information
                carbonateParametersAT["lon"] = dfUsedByAlgorithmAT["lon"].values; #Add lon lat information
                carbonateParametersAT=pd.DataFrame.from_dict(carbonateParametersAT);
                carbonateParametersAT = convert_dataframe_to_gridded_list(carbonateParametersAT, "lat", "lon"); #convert to a list of 2D matrices ready for writing to netCDF
                
                #carbonateParametersDIC = calculate_carbonate_parameters(flag, at[wNoNanDIC], dic[wNoNanDIC], sssData=sss_dic[wNoNanDIC], sstData=sst_dic[wNoNanDIC], atErr=atErr[wNoNanDIC], dicErr=dicErr[wNoNanDIC], SErr=sssErr_dic[wNoNanDIC], TErr=sstErr_dic[wNoNanDIC], pAtm=pAtm, pHydrostatic=pHydrostatic, k1k2="x");
                carbonateParametersDIC["lat"] = dfUsedByAlgorithmDIC["lat"].values; #Add lon lat information
                carbonateParametersDIC["lon"] = dfUsedByAlgorithmDIC["lon"].values; #Add lon lat information
                carbonateParametersDIC=pd.DataFrame.from_dict(carbonateParametersDIC);
                carbonateParametersDIC = convert_dataframe_to_gridded_list(carbonateParametersDIC, "lat", "lon"); #convert to a list of 2D matrices ready for writing to netCDF
                
            #write predicted output for this year/month to the netCDF file
            if griddedOutputAT is not None:
                write_outputvars_to_netCDF(ncoutAT, iyear, imonth, "AT", griddedOutputAT, gridded_uncertaintyAT, loadedInputDataAT, carbonateParametersAT,griddedOutputAT, gridded_uncertaintyAT);
            if griddedOutputDIC is not None:
                write_outputvars_to_netCDF(ncoutDIC, iyear, imonth, "DIC", griddedOutputDIC,gridded_uncertaintyDIC, loadedInputDataDIC, carbonateParametersDIC, griddedOutputDIC,gridded_uncertaintyDIC);
        
    #All months and years have been computed so close the netCDF file.
    if atAlgo is not None:
        ncoutAT.close();
    if dicAlgo is not None:
        ncoutDIC.close();
    
    
def calculate_gridded_timeseries_all_regions(tablePath, outputPathTemplate, years, regions=None, regionMaskPath=None):
    settings = osoda_global_settings.get_default_settings();    
    if regions is None:
        regions = settings["regions"];
    if regionMaskPath is None:
        regionMaskPath = settings["regionMasksPath"];
    
    bestAlgoTable = pd.read_csv(tablePath);
    
    latRes = lonRes = 1.0;
    
    if isinstance(outputPathTemplate, str):
        outputPathTemplate = Template(outputPathTemplate);
        
    regionMaskNC = Dataset(regionMaskPath, 'r');
    
    #make gridded time series predictions for each region, using the best input combination and algorithms
    for region in regions:
        if region =="oceansoda_st_lawrence":  
            pass
        elif region =="oceansoda_mississippi":
            pass
        else:
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
                atAlgo = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == atAlgoInfo["algo_name"].split(":")[0]][0];
            else:
                atAlgo = None;
            if type(dicAlgoInfo.algo_name) == str:
                dicAlgo = [algorithm for algorithm in settings["algorithmRegionMapping"][region] if algorithm.__name__ == dicAlgoInfo["algo_name"].split(":")[0]][0];
            else:
                dicAlgo = None;
            
            
            #Create output file path
            griddedPredictionOutputPathAT = outputPathTemplate.safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="TA");
            griddedPredictionOutputPathDIC = outputPathTemplate.safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="DIC");
            if path.exists(path.dirname(griddedPredictionOutputPathAT)) == False:
                os.makedirs(path.dirname(griddedPredictionOutputPathAT));
            #calculate the gridded time series and write to file for this input / region combination
            calculate_gridded_timeseries_single_region(griddedPredictionOutputPathAT, griddedPredictionOutputPathDIC, atAlgo, atAlgoInfo, dicAlgo, dicAlgoInfo, settings, years, region, regionMaskNC, latRes, lonRes, verbose=True);




