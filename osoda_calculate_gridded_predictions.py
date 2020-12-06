#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  4 01:10:09 2020
Commandline tool to calculate gridded DIC and AT time series using a set of algorithms
determined to be 'optimal'.

@author: tom holding
"""

#Test arguments:
# "output/algo_metrics/overall_best_algos_min_years=8.csv" "output/gridded_predictions_min_year_range/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc"

from os import path;
import numpy as np;
import argparse; #parsing command line arguments

import osoda_global_settings;
from os_gridded_time_series.gridded_predictions import calculate_gridded_timeseries_all_regions;

#Return a dictionary mapping common variable names to DatasetInfo objects, for a given algorithm and input combination
#TODO: remove? - Is this functionality is in utilities now?
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


def main(overallBestAlgosOutputPath, outputPathTemplate, years, regions, regionMaskPath):
    calculate_gridded_timeseries_all_regions(overallBestAlgosOutputPath, outputPathTemplate, years, regions=regions, regionMaskPath=regionMaskPath);


if __name__ == "__main__":
    settings = osoda_global_settings.get_default_settings();
    lonRes = latRes = 1.0;
    
    #Setup command line parser, and parse arguments.
    description = """Utility for calculating gridded surface DIC and AT time series for the OceanSODA project.
        This utility will download and pre-process input files into 1x1 degree monthly files. Then it uses the algorithms identified as 'optimal' to calculate gridded carbonate system data and output them.
        Carbonate chemistry data is outputted to netCDF files along with input data and propagated uncertainties.
        (development, use with care)
        """;
    clParser = argparse.ArgumentParser(description=description);
    clParser.add_argument("best_algo_table_path", help="Path to the CSV file containing the overall best algorithms data (containing the OceanSODA algorithm metrics)");
    clParser.add_argument("output_path_template", help=r"Template path to output netCDF file. This will generate a unique file path for different regions and target variable (e.g. DIC / AT) combinations. Use '${REGION}' and '${OUTPUTVAR}' to represent the region and output variable names in generated output file paths.");
    clParser.add_argument("--input_data_root", "-i", help="Path to root directory containing gridded input data. If no data exists here, it will be downloaded, processed, and stored in this directory. Defaults to the prediction dataset path provided by the global settings file.", default=settings["predictionDatasetsRoot"]);
    clParser.add_argument("--start_year", "-s", help="Gridded time series will be calculated from this date forward (inclusive). Whole years only (e.g. 2010).", type=int, default=settings["years"][0]);
    clParser.add_argument("--end_year", "-e", help="Gridded time series will be calculated up to and including this date. Whole years only (e.g. 2020).", type=int, default=settings["years"][-1]);
    clParser.add_argument("--regions", "-r", nargs="+", help="List of region names separated by spaces. These must correspond to regions in the global settings file (osoda_global_settings.py). Defaults to the osoda_global_settings.py:get_default_setting()['regions'] value.", default=settings["regions"]);
    clParser.add_argument("--mask_path", '-m', help="Path to the netCDF file containing region masks. Each mask variable should be named after the region name (set using --regions) and contain 1s to indicate inclusion in the region and 0s to indicated exclusion. Defaults to the osoda_global_settings.py:get_default_setting()['regionMasksPath'] value.", default=settings["regionMasksPath"]);
    clArgs = clParser.parse_args();
    
    #Collect values for the main function call
    years = range(clArgs.start_year, clArgs.end_year+1);
    overallBestAlgosOutputPath = path.join(clArgs.best_algo_table_path);
    regions = clArgs.regions;
    regionMaskPath = clArgs.mask_path
    outputPathTemplate = clArgs.output_path_template;
    
    main(overallBestAlgosOutputPath, outputPathTemplate, years, regions, regionMaskPath);
    