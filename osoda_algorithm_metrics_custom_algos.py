#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on Sun Jun 20 07:14:07 2021

# Calculate algorithm metrics including algorithms for which there are no python implementations (i.e. only model output and std dev data are available)

# @author: tom holding


import osoda_global_settings;
import osoda_algorithm_comparison;
import os_algorithms.utilities as utilities;
import os_algorithms.metrics as metrics;
from os_algorithms.diagnostic_plotting import prediction_accuracy_plot;
from os import path, makedirs;
import csv
#import PyCO2SYS as pyco2
import math
import numpy as np;
import pandas as pd;
pd.set_option("mode.chained_assignment", None)
import pickle
from string import Template; # Added DJF 29/11/2024 - So we can update the input directory for the MMDB (for testing...)
import geopandas as gpd
from shapely.geometry import Point
from netCDF4 import Dataset

import faulthandler;
faulthandler.enable();

runPairedMetrics = False;
runBasicMetrics = True;


from math import cos, sin, asin, sqrt, radians
def calc_distance(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees):
    from: https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points/4913653#4913653
    """
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) # convert decimal degrees to radians
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2  #haversine formula
    c = 2 * asin(sqrt(a))
    km = 6371 * c
    return km

def calc_distance_to_coastline(longitude,latitude ):
    target_coordinate=Point(longitude,latitude )
    return coastline.distance(target_coordinate).values[0]

def distance_degrees_to_kilometers(distance,coord=[0,0]):
    coord_plus=[c+distance for c in coord]
    coord_minus=[c-distance for c in coord]
    return (calc_distance(*coord,*coord_plus)+calc_distance(*coord,*coord_minus))*0.5

def calc_distance_to_coastline_km(longitude,latitude ):
    target_coordinate=Point(longitude,latitude )
    return distance_degrees_to_kilometers(coastline.distance(target_coordinate).values[0],[longitude,latitude])

#outputVar is AT or DIC etc. definiing the parameter which is the target output of the algorithm
#matchupVariableName is the string name indicating where to find it in the matchup database
#algoRMSD is the goodness of fit from the original algorithm fitting process
#inputUncertainty is the uncertainty associated with the input data
"""
Old MMDB setup - DJF: 29/11/2024 - Setting up the code for the OceanHealth version of the dataset.
"""
# customAlgorithmInfo = [{"name": "ethz_at", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar": "AT", #DIC or AT
#                        "matchupVariableName": "ethz_ta_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
#                        "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                       # {"name": "ethz_dic", #Human readable name, this can be set to anything and is only used as a label
#                       #  "outputVar": "DIC", #DIC or AT
#                       #  "matchupVariableName": "ethz_dic_mean", #netCDF variable name of the model output (algorithm prediction)
#                       #  "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       #  "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
#                       #  "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       #  },
#                       {"name": "ethz_ph", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_pH_mean", # "insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "ethz_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                        "combinedUncertainty": 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                       {"name": "ethz_pco2", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar": "region_pco2w_mean", #DIC or AT
#                        "matchupVariableName": "ethz_pco2_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                        {"name": "cmems_ph", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar":"insitu_pH_mean", #"insitu_ph_mean", #DIC or AT
#                         "matchupVariableName": "cmems_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD": 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #propagated input data uncertainty
#                         "combinedUncertainty": 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         },
#                        {"name": "cmems_pco2", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar": "region_pco2w_mean", #DIC or AT
#                         "matchupVariableName": "cmems_pco2_mean", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #propagated input data uncertainty
#                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         },
#                        # {"name": "pml_at", #Human readable name, this can be set to anything and is only used as a label
#                        #  "outputVar": "AT", #DIC or AT
#                        #  "matchupVariableName": "pml_ta_mu", #netCDF variable name of the model output (algorithm prediction)
#                        #  "algoRMSD": 17.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        #  "inputUncertaintyName": None, #propagated input data uncertainty
#                        #  "combinedUncertainty": 22.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        #  },
#                     ];



settings = osoda_global_settings.get_default_settings();
settings["matchupDatasetTemplate"] = Template(path.join(settings["dataPathRoot"], "matchup_datasets/oceansoda/", "OCEANSODA-MMDB-${YYYY}-fv01.nc")); #DJF Added 29/11/2024 - Added to allow the matchup dataset location to be changed for testing :-)
outputRoot=path.join("output/new_example_test_custom_algo_metrics/");
if path.exists(outputRoot) == False:
    makedirs(outputRoot);
    makedirs(path.join(outputRoot,'diagnostic_plots'))
    makedirs(path.join(outputRoot,'algorithm_outputs'))
diagnosticPlots = True;
regions = list(settings["regions"]);
sstDatasetName = "SST-ESACCI"; #This is the matchup database variable name corresponding to the SST dataset used as input for the Ethz data
sssDatasetName = "SSS-ESACCI"; #This is the matchup database variable name corresponding to the SSS dataset used as input for the Ethz data

name_amazon={};
out_amazon={};
name_congo={};
out_congo={};
name_med={};
out_med={};

#First find the combination that corresponds to the Ethz input data.
#We'll then run all algorithms using this input data combination, in order to compare like-for-like.
specificVariableToDatabaseMaps, specificVariableToDatabaseMapNames = utilities.get_dataset_variable_map_combinations(settings);
combinationMap = combinationName = None;
for i, combiName in enumerate(specificVariableToDatabaseMapNames):
    if (sstDatasetName in combiName) and (sssDatasetName in combiName):
        combinationMap = specificVariableToDatabaseMaps[i];
        combinationName = combiName;
        break;
if combinationMap is None:
    raise ValueError("Couldn't find input data combination for the specified SST ({0}) and SSS({1}) data sets.".format(sstDatasetName, sssDatasetName));

if runPairedMetrics == True:
    buildInAlgorithmsMap = settings["algorithmRegionMapping"];
    for region in regions:
        print("Starting regions:", region);
        builtInAlgorithms = settings["algorithmRegionMapping"][region]; #[at_algorithms.Sasse2013_global_at, dic_algorithms.Sasse2013_global_dic];
        customAlgorithms = [cai for cai in customAlgorithmInfo if cai["outputVar"] in ["AT", "DIC"]];


        #### Create the directory to contain outputs. This is structured to be analogous to the main() output, even though it only contains a single input combination run
        currentCombinationOutputDirectory = path.join(outputRoot, combinationName);
        if path.exists(currentCombinationOutputDirectory) == False:
            makedirs(currentCombinationOutputDirectory);

        #### Write information about the combination of input datasets used:
        utilities.write_specific_variable_to_database_mapping(combinationMap, path.join(currentCombinationOutputDirectory, "inputs_used.txt"), combinationName);


        ### #### The non-custom algorithms need to be ran on the matchup database's input data for the selected input combination
        #### So this data must be extracted
        years = utilities.calculate_years_for_input_combination(settings, combinationMap); #Find the years where there is overlap in the selected input data combinations
        matchupData = utilities.load_matchup_to_dataframe(settings, combinationMap, years=years); #each year is concatinated to create a single dataframe

        # perform checks on the Matchup databse files to check that they are
        # #realistic and fall within the expected range

        SST_max=settings["MDB_flags"]["SST_max"];
        SST_min=settings["MDB_flags"]["SST_min"];
        SSS_max=settings["MDB_flags"]["SSS_max"];
        SSS_min=settings["MDB_flags"]["SSS_min"];
        DIC_max=settings["MDB_flags"]["DIC_max"];
        DIC_min=settings["MDB_flags"]["DIC_min"];
        pH_max=settings["MDB_flags"]["pH_max"];
        pH_min=settings["MDB_flags"]["pH_min"];
        pCO2_max=settings["MDB_flags"]["pCO2_max"];
        pCO2_min=settings["MDB_flags"]["pCO2_min"];
        TA_max=settings["MDB_flags"]["TA_max"];
        TA_min=settings["MDB_flags"]["TA_min"];

        #SST
        mdb_SST = matchupData["SST"];
        mdb_SST_numeric=mdb_SST.values-273.15;
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

        #pCO2
        mdb_pco2 = matchupData["region_pco2w_mean"];
        mdb_pco2_numeric=mdb_pco2.values;
        #Upper realistic limit 700 ppm
        states9=mdb_pco2_numeric>pCO2_max;
        index_pco2_exceed=np.where(states9)[0]
        #Lower realistic limit <200 ppm
        states10=mdb_pco2_numeric<pCO2_min;
        index_pco2_below=np.where(states10)[0]

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
                                ,'pCO2 greater than maximum value of','pCO2 less than minimum value of','TA greater than maximum value of','TA less than minimum value of']

        mdb_flag_index_list=[index_temp_exceed, index_temp_below, index_sal_exceed, index_sal_below, index_DIC_exceed, index_DIC_below,index_pH_exceed,index_pH_below,index_pco2_exceed,index_pco2_below,index_TA_exceed,index_TA_below]

        mdb_flag_limits_list=[SST_max,SST_min,SSS_max,SSS_min,DIC_max, DIC_min,  pH_max, pH_min,pCO2_max, pCO2_min, TA_max,TA_min]

        for idx, g in enumerate(mdb_flag_index_list):
            print(idx, g)
            if len(g) == 0:
                print("list is empty")
            else:
                #this prints a header for what the entries have been flagged for
                with open('mdb_flag.csv','a') as result_file:
                    wr = csv.writer(result_file, dialect='excel')
                    wr.writerow([mdb_flag_warnings_list[idx]]+ [mdb_flag_limits_list[idx]])
                #this prints the values that have been flagged to the csv
                print(matchupData.loc[g].to_csv("mdb_flag.csv",mode='a'));

        #now filter those mdb from the analysis
        mdb_ind_rmv =np.concatenate(mdb_flag_index_list) #combine all the numpy arrays into a single array of 'bad data point indexes'
        matchupData = matchupData.drop(matchupData.index[mdb_ind_rmv])

        del states, states2, states3 ,states4, states5, states6, states7 ,states8, states9, states10, states11, states12
        del mdb_DIC,mdb_DIC_numeric,mdb_SSS,mdb_SSS_numeric,mdb_SST,mdb_SST_numeric,mdb_TA,mdb_TA_numeric,mdb_pH,mdb_pH_numeric,mdb_pco2,mdb_pco2_numeric
        #Rich_edits end

        ### Run each of the build in algorithms to compute model output, extract RMSD etc
        allAlgoOutputInfo = [];
        for ialgorithm, AlgorithmFunctor in enumerate(builtInAlgorithms):
            print("Calculating model outputs for implemented algorithms ({0}/{1}): {2}".format(ialgorithm+1, len(builtInAlgorithms), str(AlgorithmFunctor)));
            algorithm = AlgorithmFunctor(settings);
            try:
                modelOutput, propagatedInputUncertainty, rmsd, combinedUncertainty, dataUsedIndices, dataUsed,subsetData = \
                          osoda_algorithm_comparison.run_algorithm(algorithm, matchupData,
                                                                   regionMaskPath=settings["regionMasksPath"], region=region,
                                                                   useDepthMask=settings["subsetWithDepthMask"],
                                                                   depthMaskPath=settings["depthMaskPath"],
                                                                   depthMaskVar=settings["depthMaskVar"],
                                                                   useDistToCoastMask=settings["subsetWithDistToCoast"],
                                                                   distToCoastMaskPath=settings["distToCoastMaskPath"],
                                                                   distToCoastMaskVar=settings["distToCoastMaskVar"]
                                                                   );

                ### Store the algorithm output information in a dictionary, and then append it to the list
                algorithmOutput = {};
                algorithmOutput["instance"] = algorithm;
                algorithmOutput["name"] = algorithm.__str__();
                algorithmOutput["outputVar"] = algorithm.output_name();
                algorithmOutput["modelOutput"] = modelOutput;
                algorithmOutput["propagatedInputUncertainty"] = propagatedInputUncertainty;
                algorithmOutput["rmsd"] = rmsd;
                algorithmOutput["combinedUncertainty"] = combinedUncertainty;
                algorithmOutput["dataUsedIndices"] = dataUsedIndices;
                allAlgoOutputInfo.append(algorithmOutput);

                if region=="oceansoda_amazon_plume":
                    name_amazon[AlgorithmFunctor]=algorithmOutput["name"];
                    out_amazon[AlgorithmFunctor]=algorithmOutput["modelOutput"];
                elif region=="oceansoda_congo":
                    name_congo[AlgorithmFunctor]=algorithmOutput["name"];
                    out_congo[AlgorithmFunctor]=algorithmOutput["modelOutput"];
                elif region=="oceansoda_mediterranean":
                    name_med[AlgorithmFunctor]=algorithmOutput["name"];
                    out_med[AlgorithmFunctor]=algorithmOutput["modelOutput"];
                else:
                    pass

                if diagnosticPlots == True:
                    outputVariable = algorithm.output_name();
                    savePath = path.join(outputRoot, "diagnostic_plots", algorithm.__str__().split(":")[0]+".png");
                    prediction_accuracy_plot(dataUsed[outputVariable], modelOutput, algorithm.__class__.__name__, outputVariable, savePath=savePath);


                print("Output calculated (region:"+region+", algo: "+algorithm.__class__.__name__+")");
            except ValueError as e: #Raised if there are no matchup data rows left after spatial mask and algorithm internal subsetting has taken place
                print(algorithm, e.args[0]);
                print("No matchup data left after subsettings (region:"+region+", algo: "+algorithm.__class__.__name__+")");
                continue;


        ###### Now add extract the custom algorithm data (e.g. rmsd, output), construct output dictionaries, and append them to the algo output info list.
        for customAlgo in customAlgorithms:

            #Extract the model output from the matchup data
            colsToExtract = ["date", customAlgo["outputVar"], customAlgo["matchupVariableName"]];
            ##### TODO: Missing: combined uncertainty, matchupRMSD, matchupBias?
            customAlgoData = utilities.read_matchup_cols(settings["matchupDatasetTemplate"], colsToExtract, years); #returns data frame containing data from the matchup database for each variable in 'cols'

            #Extract the model output from the matchup data
            colsToExtract = ["date", customAlgo["outputVar"], customAlgo["matchupVariableName"]];
            ##### TODO:
            ##### Missing: combined uncertainty, matchupRMSD, matchupBias?
            customAlgoData = utilities.read_matchup_cols(settings["matchupDatasetTemplate"], colsToExtract, years); #returns data frame containing data from the matchup database for each variable in 'cols'

            #this takes the subset information e.g. region/depth/distance coast filters etc
            #and applies it to the custom algos as well
            customAlgoData=customAlgoData.loc[subsetData.index,]

            # data are loaded in here again, delete rows that are removed by QC
            A=customAlgoData.index#these are the rows in the subset
            B=mdb_ind_rmv#these are the bad rows to removed
            #this finds the indexes
            B_unique_sorted, B_idx = np.unique(B, return_index=True)
            B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
            inthesubset=B_idx[B_in_A_bool]
            C=B[inthesubset]
            #this line then drops the bad rows
            customAlgoData = customAlgoData.drop(customAlgoData.index[inthesubset]); #remove where there is no reference outputVar data

            ###Subset to remove where model data is NaN
            customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["matchupVariableName"]])]; #remove where there is no model predictions
            customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["outputVar"]])]; #remove where there is no reference outputVar data

            if customAlgo["name"] == "cmems_pco2": #unit conversion
                customAlgoData["cmems_pco2_mean"] = customAlgoData["cmems_pco2_mean"]*0.00000986923 * 1000000;

            algorithmOutput = {};
            algorithmOutput["instance"] = None;
            algorithmOutput["name"] = customAlgo["name"];
            algorithmOutput["outputVar"] = customAlgo["outputVar"]; #e.g. AT or DIC
            algorithmOutput["modelOutput"] = customAlgoData[customAlgo["matchupVariableName"]]; #The modelled output DIC or AT from matchup database
            algorithmOutput["propagatedInputUncertainty"] = None;
            algorithmOutput["rmsd"] = customAlgo["algoRMSD"]; #Goodness of fit from the original algorithm fit
            algorithmOutput["combinedUncertainty"] = customAlgo["combinedUncertainty"];
            algorithmOutput["dataUsedIndices"] = customAlgoData.index;
            allAlgoOutputInfo.append(algorithmOutput);

            if diagnosticPlots == True:
                savePath = path.join(outputRoot, "diagnostic_plots", customAlgo["name"]+".png");
                prediction_accuracy_plot(customAlgoData[customAlgo["outputVar"]], customAlgoData[customAlgo["matchupVariableName"]], customAlgo["name"], customAlgo["outputVar"], savePath=savePath);

        ####Calculate metrics
        for currentOutputVar in ["AT", "DIC"]:
            ##Subset the algorithm output data by the target output variable
            algorithmOutputGroup = [algorithmOutput for algorithmOutput in allAlgoOutputInfo if algorithmOutput["outputVar"]==currentOutputVar];


            basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores = \
                            metrics.calc_all_metrics(algorithmOutputGroup, matchupData, settings,currentOutputVar);
            print("Completed metrics for", currentOutputVar);

            ###########################
            ### Write outputs to file #
            outputDirectory = path.join(currentCombinationOutputDirectory, currentOutputVar, region);
            print("Writing metrics to: ", outputDirectory);
            osoda_algorithm_comparison.write_metrics_to_file(outputDirectory, matchupData, basicMetrics, nIntersectMatrix, pairedScoreMatrix, pairedWScoreMatrix, pairedRmsdMatrix, pairedWRmsdMatrix, finalScores, algorithmOutputGroup);

        n_threshold=30
        #Calculate summary table for weighted metrics and output to file
        summaryTable_weighted = osoda_algorithm_comparison.create_summary_table(n_threshold,settings, [combinationName], [combinationMap], useWeighted=True, regionOverload = regions);
        summaryTableOutputPath = path.join(outputRoot, combinationName, "summary_best_algos.csv");
        summaryTable_weighted.to_csv(summaryTableOutputPath, sep=",", index=False);
        print("Full weighted summary table written to:", path.abspath(summaryTableOutputPath));

        #Calculate summary table for unweighted metrics and output to fill
        summaryTable_unweighted = osoda_algorithm_comparison.create_summary_table(n_threshold,settings, [combinationName], [combinationMap], useWeighted=False, regionOverload = regions);
        summaryTableOutputPathUnweighted = path.join(outputRoot, combinationName, "summary_best_algos_unweighted.csv");
        summaryTable_unweighted.to_csv(summaryTableOutputPathUnweighted, sep=",", index=False);
        print("Full unweighted summary table written to:", path.abspath(summaryTableOutputPathUnweighted));


"""
Old MMBD statisitcs extraction
"""
# customAlgorithmInfo = [{"name": "ethz_at", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar": "insitu_ta_mean", #DIC or AT
#                        "matchupVariableName": "ethz_ta_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 13.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_ta_stddev", #propagated input data uncertainty
#                        "combinedUncertainty": 21.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                       # {"name": "ethz_dic", #Human readable name, this can be set to anything and is only used as a label
#                       #  "outputVar": "DIC", #DIC or AT
#                       #  "matchupVariableName": "ethz_dic_mean", #netCDF variable name of the model output (algorithm prediction)
#                       #  "algoRMSD": 16.3, #netCDF variable name of the RMSD of the (original) algorithm fit
#                       #  "inputUncertaintyName": None, #"ethz_dic_stddev", #propagated input data uncertainty
#                       #  "combinedUncertainty": 16.3, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                       #  },
#                       {"name": "ethz_ph", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar":"insitu_ph_mean", # "insitu_ph_mean", #DIC or AT
#                        "matchupVariableName": "ethz_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD":None, # value from algo - 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                        "combinedUncertainty":None, # value from algo - 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                        {"name": "ethz_ph_tadic", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar":"insitu_ph_from_ta_dic_mean", # "insitu_ph_mean", #DIC or AT
#                         "matchupVariableName": "ethz_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD":None, # value from algo - 0.024, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #"ethz_ph_stddev", #propagated input data uncertainty
#                         "combinedUncertainty":None, # value from algo - 0.024, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         },
#                       {"name": "ethz_pco2", #Human readable name, this can be set to anything and is only used as a label
#                        "outputVar": "insitu_pco2w_mean", #DIC or AT
#                        "matchupVariableName": "ethz_pco2_mean", #netCDF variable name of the model output (algorithm prediction)
#                        "algoRMSD": 12.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        "inputUncertaintyName": None, #"ethz_pco2_stddev", #propagated input data uncertainty
#                        "combinedUncertainty": 14.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        },
#                        {"name": "cmems_ph", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar":"insitu_ph_mean", #"insitu_ph_mean", #DIC or AT
#                         "matchupVariableName": "cmems_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD":None,# 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #propagated input data uncertainty
#                         "combinedUncertainty":None, #value from algo - 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         },
#                         {"name": "cmems_ph_tadic", #Human readable name, this can be set to anything and is only used as a label
#                          "outputVar":"insitu_ph_from_ta_dic_mean", #"insitu_ph_mean", #DIC or AT
#                          "matchupVariableName": "cmems_ph_mean", #netCDF variable name of the model output (algorithm prediction)
#                          "algoRMSD":None,# 0.03, #netCDF variable name of the RMSD of the (original) algorithm fit
#                          "inputUncertaintyName": None, #propagated input data uncertainty
#                          "combinedUncertainty":None, #value from algo - 0.03, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                          },
#                        {"name": "cmems_pco2", #Human readable name, this can be set to anything and is only used as a label
#                         "outputVar": "insitu_pco2w_mean", #DIC or AT
#                         "matchupVariableName": "cmems_pco2_mean", #netCDF variable name of the model output (algorithm prediction)
#                         "algoRMSD": 14.8, #netCDF variable name of the RMSD of the (original) algorithm fit
#                         "inputUncertaintyName": None, #propagated input data uncertainty
#                         "combinedUncertainty": 17.97, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                         }
#                        # {"name": "pml_at", #Human readable name, this can be set to anything and is only used as a label
#                        #  "outputVar": "AT", #DIC or AT
#                        #  "matchupVariableName": "pml_ta_mu", #netCDF variable name of the model output (algorithm prediction)
#                        #  "algoRMSD": 17.0, #netCDF variable name of the RMSD of the (original) algorithm fit
#                        #  "inputUncertaintyName": None, #propagated input data uncertainty
#                        #  "combinedUncertainty": 22.0, #netCDF variable name of the propagated uncertainty combining input data uncertainty with algorithm fit uncertainty
#                        #  },
#                     ];
"""
New MMDB statistics extraction - Some names for the products have changes (I.e ETHZ OHOA is now the high res version)
"""

coastal_distance_from_coast = 300
natural_earth_coastline_file = './data/coastline/ne_10m_coastline/ne_10m_coastline.shp'
gebco_bathymetry_file = 'F:/Data/Bathymetry/GEBCO_2023.nc'
coastal_depth_definition = 1000

# c = Dataset(gebco_bathymetry_file,'r')
# lon = np.array(c['lon'])
# lat = np.array(c['lat'])
# depth = np.array(c['elevation'])
# c.close()

import osoda_algorithm_metrics_custom_algos_config
customAlgorithmInfo= osoda_algorithm_metrics_custom_algos_config.customAlgorithmInfo

#note that in this configuration the basic metrics have not been subset at all meaning
#they are comparing globally
if runBasicMetrics == True:
    ##########
    #calculate individual level metrics for the non-AT/DIC parameters
    import os_algorithms.metrics as metrics;
    import pickle;
    import json;


    # years = utilities.calculate_years_for_input_combination(settings, combinationMap);
    # print(years)
    years = range(1980,2024,1)
    # matchupData = utilities.load_matchup_to_dataframe(settings, combinationMap, years); #each year is concatinated to create a single dataframe
    # print(list(matchupData.columns.values))
    # #note that the mdb ph is at 25 so we need to convert it to the SST in the SST column
    # #This should be fixed in future versions of the MDB 3.4 onwards!
    #
    # #for every ph entry in matchupdatabase with pH
    # #look for ph - ta pairing
    # #look for ph - pco2 pairing
    # #look for ph dic pairing
    # # if no ph pairing remove that entry.
    #
    # #create a new variable in matchupdata for ph that is temp adjusted and populate with nans
    # matchupData["ph_corr_insitu_temp"] = ""
    # NaN = np.nan
    # matchupData["ph_corr_insitu_temp"]  = NaN
    #
    # matchupData["ph_corr_insitu_temp_err"] = ""
    # NaN = np.nan
    # matchupData["ph_corr_insitu_temp_err"]  = NaN
    #
    # matchupData["hydrogen_free"] = ""
    # NaN = np.nan
    # matchupData["hydrogen_free"]  = NaN
    #
    # loopsize=matchupData.region_pH_mean.size
    # for ph_loop in range(0,loopsize):
    #     print(ph_loop)
    #     #if the ph value is nan make new variable nan
    #     if math.isnan(matchupData.region_pH_mean[ph_loop])==True:
    #         #print("skipping at nan step")
    #         pass
    #         #do nothing - no ph data to temp adjust
    #         #look for ph - ta pairing
    #     elif math.isnan(matchupData.region_pH_mean[ph_loop])==False and math.isnan(matchupData.AT[ph_loop])==False:
    #         #CO2SYS tHE DATA
    #         #print("use ta and ph")
    #         kwargs = dict(
    #         par1 = matchupData.AT[ph_loop],  # Value of the first parameter
    #         par2 = matchupData.region_pH_mean[ph_loop],  # Value of the second parameter
    #         par1_type = 1,  # The first parameter supplied is of type "1", which is "alkalinity"
    #         par2_type = 3,  # The second parameter supplied is of type "2", which is "pH"
    #         salinity = matchupData.SSS[ph_loop],  # Salinity of the sample
    #         temperature = 25,  # Temperature at input conditions
    #         temperature_out = matchupData.SST[ph_loop]-273.15,  # Temperature at output conditions
    #         pressure = 0,  # Pressure    at input conditions
    #         pressure_out = 0,  # Pressure    at output conditions
    #         total_silicate = matchupData.SiO4[ph_loop],  # Concentration of silicate  in the sample (in umol/kg)
    #         total_phosphate = matchupData.PO4[ph_loop],  # Concentration of phosphate in the sample (in umol/kg)
    #         opt_k_carbonic = 4,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    #         opt_k_bisulfate = 1,);  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    #         results = pyco2.sys(**kwargs);
    #         matchupData["ph_corr_insitu_temp"][ph_loop]=results["pH_out"];
    #         matchupData["ph_corr_insitu_temp_err"][ph_loop]=matchupData["region_pH_mean_err"][ph_loop]
    #         matchupData["hydrogen_free"][ph_loop]=results["hydrogen_free_out"]*1e-6;
    #     #look for ph - pco2 pairing
    #     elif math.isnan(matchupData.region_pH_mean[ph_loop])==False and math.isnan(matchupData.region_pco2w_mean[ph_loop])==False:
    #         #print("use ta and pco2")
    #         #CO2SYS tHE DATA
    #         kwargs = dict(
    #         par1 = matchupData.region_pco2w_mean[ph_loop],  # Value of the first parameter
    #         par2 = matchupData.region_pH_mean[ph_loop],  # Value of the second parameter
    #         par1_type = 4,  # The first parameter supplied is of type "1", which is "PCO2"
    #         par2_type = 3,  # The second parameter supplied is of type "2", which is "pH"
    #         salinity = matchupData.SSS[ph_loop],  # Salinity of the sample
    #         temperature = 25,  # Temperature at input conditions
    #         temperature_out = matchupData.SST[ph_loop]-273.15,  # Temperature at output conditions
    #         pressure = 0,  # Pressure    at input conditions
    #         pressure_out = 0,  # Pressure    at output conditions
    #         total_silicate = matchupData.SiO4[ph_loop],  # Concentration of silicate  in the sample (in umol/kg)
    #         total_phosphate = matchupData.PO4[ph_loop],  # Concentration of phosphate in the sample (in umol/kg)
    #         opt_k_carbonic = 4,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    #         opt_k_bisulfate = 1,)  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    #         results = pyco2.sys(**kwargs)
    #         matchupData["ph_corr_insitu_temp"][ph_loop]=results["pH_out"]
    #         matchupData["ph_corr_insitu_temp_err"][ph_loop]=matchupData["region_pH_mean_err"][ph_loop]
    #         matchupData["hydrogen_free"][ph_loop]=results["hydrogen_free_out"]*1e-6;
    #     #look for ph - dic pairing
    #     elif math.isnan(matchupData.region_pH_mean[ph_loop])==False and math.isnan(matchupData.DIC[ph_loop])==False:
    #         #print("use ta and dic")
    #         #CO2SYS tHE DATA
    #         kwargs = dict(
    #         par1 = matchupData.DIC[ph_loop],  # Value of the first parameter
    #         par2 = matchupData.region_pH_mean[ph_loop],  # Value of the second parameter
    #         par1_type = 2,  # The first parameter supplied is of type "1", which is "DIC"
    #         par2_type = 3,  # The second parameter supplied is of type "2", which is "pH"
    #         salinity = matchupData.SSS[ph_loop],  # Salinity of the sample
    #         temperature = 25,  # Temperature at input conditions
    #         temperature_out = matchupData.SST[ph_loop]-273.15,  # Temperature at output conditions
    #         pressure = 0,  # Pressure    at input conditions
    #         pressure_out = 0,  # Pressure    at output conditions
    #         total_silicate = matchupData.SiO4[ph_loop],  # Concentration of silicate  in the sample (in umol/kg)
    #         total_phosphate = matchupData.PO4[ph_loop],  # Concentration of phosphate in the sample (in umol/kg)
    #         opt_k_carbonic = 4,  # Choice of H2CO3 and HCO3- dissociation constants K1 and K2 ("4" means "Mehrbach refit")
    #         opt_k_bisulfate = 1,)  # Choice of HSO4- dissociation constants KSO4 ("1" means "Dickson")
    #         results = pyco2.sys(**kwargs)
    #         matchupData["ph_corr_insitu_temp"][ph_loop]=results["pH_out"]
    #         matchupData["ph_corr_insitu_temp_err"][ph_loop]=matchupData["region_pH_mean_err"][ph_loop]
    #         matchupData["hydrogen_free"][ph_loop]=results["hydrogen_free_out"]*1e-6;
    #     #if there is no secondary carbonate variable the correction is not possible
    #     else:
    #         pass
    #         #print('Congratulations! You guessed it.')
    #         #do nothing pass
    #
    # #Rich_edits start_ perform checks on the Matchup databse files to check that they are
    # #realistic and fall within the expected range
    #
    # SST_max=settings["MDB_flags"]["SST_max"];
    # SST_min=settings["MDB_flags"]["SST_min"];
    # SSS_max=settings["MDB_flags"]["SSS_max"];
    # SSS_min=settings["MDB_flags"]["SSS_min"];
    # DIC_max=settings["MDB_flags"]["DIC_max"];
    # DIC_min=settings["MDB_flags"]["DIC_min"];
    # pH_max=settings["MDB_flags"]["pH_max"];
    # pH_min=settings["MDB_flags"]["pH_min"];
    # pCO2_max=settings["MDB_flags"]["pCO2_max"];
    # pCO2_min=settings["MDB_flags"]["pCO2_min"];
    # TA_max=settings["MDB_flags"]["TA_max"];
    # TA_min=settings["MDB_flags"]["TA_min"];
    # hfree_max=1e-7;
    # hfree_min=0;
    #
    # #SST
    # mdb_SST = matchupData["SST"];
    # mdb_SST_numeric=mdb_SST.values-273.15;
    # #Upper realistic limit
    # states=mdb_SST_numeric>SST_max;
    # index_temp_exceed=np.where(states)[0]
    # #Lower realistic limit -10 degrees
    # states2=mdb_SST_numeric<SST_min;
    # index_temp_below=np.where(states2)[0]
    #
    # #Salinity
    # mdb_SSS = matchupData["SSS"];
    # mdb_SSS_numeric=mdb_SSS.values;
    # #Upper realistic limit 50 PSU
    # states3=mdb_SSS_numeric>SSS_max;
    # index_sal_exceed=np.where(states3)[0]
    # #Lower realistic limit <0 PSU
    # states4=mdb_SSS_numeric<SSS_min;
    # index_sal_below=np.where(states4)[0]
    #
    # #DIC
    # mdb_DIC = matchupData["DIC"];
    # mdb_DIC_numeric=mdb_DIC.values;
    # #Upper realistic limit 2500 UMOLKG?
    # states5=mdb_DIC_numeric>DIC_max;
    # index_DIC_exceed=np.where(states5)[0]
    # #Lower realistic limit <500 UMOL KG
    # states6=mdb_DIC_numeric<DIC_min;
    # index_DIC_below=np.where(states6)[0]
    #
    # #pH
    # mdb_pH = matchupData["ph_corr_insitu_temp"];
    # mdb_pH_numeric=mdb_pH.values;
    # #Upper realistic limit 8.5
    # states7=mdb_pH_numeric>pH_max;
    # index_pH_exceed=np.where(states7)[0]
    # #Lower realistic limit <7
    # states8=mdb_pH_numeric<pH_min;
    # index_pH_below=np.where(states8)[0]
    #
    # #pH - as H+
    # mdb_hfree = matchupData["hydrogen_free"];
    # mdb_hfree_numeric=mdb_hfree.values;
    # #Upper realistic limit 8.5
    # states13=mdb_hfree_numeric>hfree_max;
    # index_hfree_exceed=np.where(states13)[0]
    # #Lower realistic limit <7
    # states14=mdb_hfree_numeric<hfree_min;
    # index_hfree_below=np.where(states14)[0]
    #
    # #pCO2
    # mdb_pco2 = matchupData["region_pco2w_mean"];
    # mdb_pco2_numeric=mdb_pco2.values;
    # #Upper realistic limit 700 ppm
    # states9=mdb_pco2_numeric>pCO2_max;
    # index_pco2_exceed=np.where(states9)[0]
    # #Lower realistic limit <200 ppm
    # states10=mdb_pco2_numeric<pCO2_min;
    # index_pco2_below=np.where(states10)[0]
    #
    # #TA
    # mdb_TA = matchupData["AT"];
    # mdb_TA_numeric=mdb_TA.values;
    # #Upper realistic limit 3000 umol kg
    # states11=mdb_TA_numeric>TA_max;
    # index_TA_exceed=np.where(states11)[0]
    # #Lower realistic limit <500 umol kg
    # states12=mdb_TA_numeric<TA_min;
    # index_TA_below=np.where(states12)[0]
    #
    #
    # #these variables could also have bounds placed on them but not applied
    # #for this iteration
    # # lat long date OC chla DO NO3 PO4 SiO4
    #
    # # now produce a file with all of the out of bounds data points from the mdb
    #
    # mdb_flag_warnings_list=['SST greater than maximum value of', 'SST less than minimum value of', 'SSS greater than maximum value of', 'SST less than minimum value of'\
    #                         , 'DIC greater than maximum value of', 'DIC less than minimum value of','pH greater than maximum value of','pH less than minimum value of'\
    #                         ,'pCO2 greater than maximum value of','pCO2 less than minimum value of','TA greater than maximum value of','TA less than minimum value of'\
    #                             ,'H+ greater than maximum value of','H+ less than minimum value of']
    #
    # mdb_flag_index_list=[index_temp_exceed, index_temp_below, index_sal_exceed, index_sal_below, index_DIC_exceed, index_DIC_below\
    #                      ,index_pH_exceed,index_pH_below,index_pco2_exceed,index_pco2_below,index_TA_exceed,index_TA_below,index_hfree_exceed,index_hfree_below]
    #
    #
    # mdb_flag_limits_list=[SST_max,SST_min,SSS_max,SSS_min,DIC_max, DIC_min,  pH_max, pH_min,pCO2_max, pCO2_min, TA_max,TA_min,hfree_max,hfree_min]
    #
    # for idx, g in enumerate(mdb_flag_index_list):
    #     print(idx, g)
    #     if len(g) == 0:
    #         print("list is empty")
    #     else:
    #         #this prints a header for what the entries have been flagged for
    #         with open('mdb_flag.csv','a') as result_file:
    #             wr = csv.writer(result_file, dialect='excel')
    #             wr.writerow([mdb_flag_warnings_list[idx]]+ [mdb_flag_limits_list[idx]])
    #         #this prints the values that have been flagged to the csv
    #         print(matchupData.loc[g].to_csv("mdb_flag.csv",mode='a'));
    #
    # #now filter those mdb from the analysis
    # mdb_ind_rmv =np.concatenate(mdb_flag_index_list) #combine all the numpy arrays into a single array of 'bad data point indexes'
    # del states, states2, states3 ,states4, states5, states6, states7 ,states8, states9, states10, states11, states12
    # del mdb_DIC,mdb_DIC_numeric,mdb_SSS,mdb_SSS_numeric,mdb_SST,mdb_SST_numeric,mdb_TA,mdb_TA_numeric,mdb_pH,mdb_pH_numeric,mdb_pco2,mdb_pco2_numeric
    #
    #
    #
    # matchupData = matchupData.drop(matchupData.index[mdb_ind_rmv])
    # print(list(matchupData.columns.values))

    #Rich_edits end


    #Extract algorithm output for each custom algorithm
    algorithmOutputs = [];
    dataUsedList = [];
    for customAlgo in customAlgorithmInfo:
        print("Extracting data for custom algorithm: {0}".format(customAlgo["name"]));

        #Extract the model output from the matchup data
        if customAlgo['thrid'] != None:
            colsToExtract = ["date", "lat","lon", customAlgo["outputVar"], customAlgo["matchupVariableName"],customAlgo['secondary'],customAlgo['thrid']];
        elif customAlgo['secondary'] != None:
            colsToExtract = ["date", "lat","lon", customAlgo["outputVar"], customAlgo["matchupVariableName"],customAlgo['secondary']];
        else:
            colsToExtract = ["date", "lat","lon", customAlgo["outputVar"], customAlgo["matchupVariableName"]];

        #if (coastal == 'coastal') | (coastal == 'open'):
        colsToExtract.append('elevation')
        colsToExtract.append('distance_to_coast')
        # colsToExtract.append('insitu_waterdepth_mean')
        # colsToExtract.append('insitu_dist2coast_mean')
        ##### TODO:
        ##### Missing: combined uncertainty, matchupRMSD, matchupBias?
        customAlgoData = utilities.read_matchup_cols(settings["matchupDatasetTemplate"], colsToExtract, years); #returns data frame containing data from the matchup database for each variable in 'cols'

        if customAlgo["matchupVariableName"] == 'cmems_biocarbon_mean_spco2':
            customAlgoData[customAlgo["matchupVariableName"]] = customAlgoData[customAlgo["matchupVariableName"]]/101325 * 10**6
        #customAlgoData = matchupData.copy()
        ###Subset to remove where model data is NaN
        customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["matchupVariableName"]])]; #remove where there is no model predictions
        customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["outputVar"]])]; #remove where there is no reference outputVar data
        if customAlgo['secondary'] != None:
            customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["secondary"]])];
        if customAlgo['thrid'] != None:
            customAlgoData = customAlgoData.loc[np.isfinite(customAlgoData[customAlgo["thrid"]])];
        """
        This adds spatial region splitting using the lat and lon definitions
        """
        if customAlgo['lat_limits'] != None:
            customAlgoData = customAlgoData.loc[(customAlgoData["lat"] >= customAlgo['lat_limits'][0]) & (customAlgoData["lat"] <= customAlgo['lat_limits'][1])]
        if customAlgo['lon_limits'] != None:
            if customAlgo['lon_or'] == True:
                customAlgoData = customAlgoData.loc[((customAlgoData["lon"] >= customAlgo['lon_limits'][0][0]) & (customAlgoData["lon"] <= customAlgo['lon_limits'][0][1])) | ((customAlgoData["lon"] >= customAlgo['lon_limits'][1][0]) & (customAlgoData["lon"] <= customAlgo['lon_limits'][1][1]))]
            else:
                customAlgoData = customAlgoData.loc[((customAlgoData["lon"] >= customAlgo['lon_limits'][0]) & (customAlgoData["lon"] <= customAlgo['lon_limits'][1]))]

        if customAlgo['exclude_box'] != None:
            for v in customAlgo['exclude_box']:
                customAlgoData = customAlgoData.loc[((customAlgoData["lon"] <= v[0]) | (customAlgoData["lon"] >= v[1]) | (customAlgoData["lat"] <= v[2]) | (customAlgoData["lat"] >= v[3]))]
        """
        This adds splitting based upon whether its "coastal" as defined by a distance from a coastline and bathymetry data
        """
        if customAlgo['coastal'] == 'coastal':
            # coastline_data = gpd.read_file(natural_earth_coastline_file)
            # coastline = gpd.GeoSeries(coastline_data.geometry.unary_union)
            #
            #
            # coast = np.zeros((len(customAlgoData["lon"]),2));
            # for i in range(len(customAlgoData["lon"])):
            #     coast[i,0] = calc_distance_to_coastline_km(customAlgoData.iloc[i]["lon"],customAlgoData.iloc[i]["lat"])
            #     g = np.where(np.abs(customAlgoData.iloc[i]["lon"] - lon) == np.min(np.abs(customAlgoData.iloc[i]["lon"] - lon)))[0]
            #     f = np.where(np.abs(customAlgoData.iloc[i]["lat"] - lat) == np.min(np.abs(customAlgoData.iloc[i]["lat"] - lat)))[0]
            #     coast[i,1] = depth[f[0],g[0]]
            # customAlgoData["distance_to_coast"] = coast[:,0]
            # customAlgoData["depth"] = coast[:,1]

            # customAlgoData = customAlgoData.loc[((customAlgoData["insitu_dist2coast_mean"] < coastal_distance_from_coast) & (customAlgoData["insitu_dist2coast_mean"] > 0)) | (customAlgoData["insitu_waterdepth_mean"] < coastal_depth_definition)]
            customAlgoData = customAlgoData.loc[((customAlgoData["distance_to_coast"] < coastal_distance_from_coast) & (customAlgoData["distance_to_coast"] > 0)) | (-customAlgoData["elevation"] < coastal_depth_definition)]

        elif customAlgo['coastal'] == 'open':
            # coastline_data = gpd.read_file(natural_earth_coastline_file)
            # coastline = gpd.GeoSeries(coastline_data.geometry.unary_union)
            #
            # coast = np.zeros((len(customAlgoData["lon"]),2));
            # for i in range(len(customAlgoData["lon"])):
            #     coast[i,0] = calc_distance_to_coastline_km(customAlgoData.iloc[i]["lon"],customAlgoData.iloc[i]["lat"])
            #     g = np.where(np.abs(customAlgoData.iloc[i]["lon"] - lon) == np.min(np.abs(customAlgoData.iloc[i]["lon"] - lon)))[0]
            #     f = np.where(np.abs(customAlgoData.iloc[i]["lat"] - lat) == np.min(np.abs(customAlgoData.iloc[i]["lat"] - lat)))[0]
            #     coast[i,1] = depth[f[0],g[0]]
            # customAlgoData["distance_to_coast"] = coast[:,0]
            # customAlgoData["depth"] = coast[:,1]
            #customAlgoData = customAlgoData.loc[(customAlgoData["insitu_dist2coast_mean"] > coastal_distance_from_coast) & (customAlgoData["insitu_waterdepth_mean"] > coastal_depth_definition)]
            customAlgoData = customAlgoData.loc[(customAlgoData["distance_to_coast"] > coastal_distance_from_coast) & (-customAlgoData["elevation"] > coastal_depth_definition)]

        # else:
        #     coast = np.zeros((len(customAlgoData["lon"]),2)); coast[:] = np.nan
        #     customAlgoData["distance_to_coast"] = coast[:,0]
        #     customAlgoData["depth"] = coast[:,1]
        ## DJF: The pH is corrected from 25 degC to in situ temp above, generating the QC'ed matchup database.
        # Then for some reason the data is loaded again at L691. So instead we copy the pandas table that we already generated,
        # and then continue the code. As the is an already QC'ed table, we dont need the QCing below.

        # # data are loaded in here again, delete rows that are removed by QC
        # A=customAlgoData.index#these are the rows in the subset
        # B=mdb_ind_rmv#these are the bad rows to removed
        # #this finds the indexes
        # B_unique_sorted, B_idx = np.unique(B, return_index=True)
        # B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
        # inthesubset=B_idx[B_in_A_bool]
        # C=B[inthesubset]
        # #this line then drops the bad rows
        #
        #
        # D=customAlgoData.index#these are the rows in the subset
        #
        # D_unique_sorted, D_idx = np.unique(D, return_index=True)
        # D_in_C_bool = np.in1d(D_unique_sorted, C, assume_unique=True)
        # inthesubset2=D_idx[D_in_C_bool]

        # customAlgoData = customAlgoData.drop(customAlgoData.index[inthesubset2]); #remove where there is no reference outputVar data

        # if customAlgo["name"] == "cmems_pco2": #unit conversion
        #     customAlgoData["cmems_pco2_mean"] = customAlgoData["cmems_pco2_mean"]*0.00000986923 * 1000000;

        algorithmOutput = {};
        algorithmOutput["instance"] = None;
        algorithmOutput["name"] = customAlgo["name"];
        algorithmOutput["propagatedInputUncertainty"] = None;
        algorithmOutput["rmsd"] = customAlgo["algoRMSD"]; #Goodness of fit from the original algorithm fit
        algorithmOutput["combinedUncertainty"] = customAlgo["combinedUncertainty"];
        algorithmOutput["outputVar"] = customAlgo["outputVar"]; #e.g. AT or DIC
        algorithmOutput["dataUsedIndices"] = customAlgoData.index;
        algorithmOutput["modelOutput"] = customAlgoData[customAlgo["matchupVariableName"]]; #The modelled output DIC or AT from matchup database
        algorithmOutput["outputVar_data"] = customAlgoData[customAlgo["outputVar"]]
        algorithmOutput["lon"] = customAlgoData["lon"]
        algorithmOutput["lat"] = customAlgoData["lat"]
        # algorithmOutput["distance_to_coast"] = customAlgoData["insitu_dist2coast_mean"]
        algorithmOutput["distance_to_coast"] = customAlgoData["distance_to_coast"]
        # algorithmOutput["depth"] = customAlgoData["insitu_waterdepth_mean"]
        algorithmOutput["depth"] = -customAlgoData["elevation"]
        print(algorithmOutput["modelOutput"])
        print(algorithmOutput["outputVar_data"])
        pickle.dump(algorithmOutput, open(path.join(outputRoot,'algorithm_outputs', customAlgo["name"]+".pickle"), 'wb'));
        #need to use new processed ph instead
        # if algorithmOutput["name"] == "ethz_ph" or algorithmOutput["name"] == "cmems_ph":
        #      #know the indexes that correspond to matchupdatabase, use these to access
        #      #the ph column added above with pyc02sys
        #      algorithmOutput["dataUsedIndices"]= matchupData["hydrogen_free"][customAlgoData.index];#was previously matchupData["ph_corr_insitu_temp"][customAlgoData.index];
        #      print("used correctly")
        #      #now remove entries from matchup where there isn't any temp corrected matchup
        #      Indexes_temp_adjusted_Ph=~(np.isnan(algorithmOutput["dataUsedIndices"]))#get the indexes where there is a pH
        #      algorithmOutput["dataUsedIndices"]=algorithmOutput["dataUsedIndices"][Indexes_temp_adjusted_Ph]
        #      algorithmOutput["modelOutput"]=10 ** (-1*algorithmOutput["modelOutput"][Indexes_temp_adjusted_Ph])# this is cmems/ethz ph- converted to H+
        #      algorithmOutput["outputVar"] = "hydrogen_free";
        #      algorithmOutput["dataUsedIndices"]=algorithmOutput["dataUsedIndices"].index
        # else:
        #     pass

        x=algorithmOutput["dataUsedIndices"]
        algorithmOutputs.append(algorithmOutput);
        dataUsedList.append((customAlgoData.reindex(x)));
        #dataUsedList.append((customAlgoData.iloc(customAlgoData.index)))
        if diagnosticPlots == True:
            basicMetrics = metrics.calc_basic_metrics(algorithmOutput, (customAlgoData.reindex(x)), settings);
            pickle.dump(basicMetrics, open(path.join(outputRoot,'algorithm_outputs', customAlgo["name"]+"_metrics.pickle"), 'wb'));
            savePath = path.join(outputRoot, "diagnostic_plots", customAlgo["name"]+".png");
            prediction_accuracy_plot(customAlgoData[customAlgo["outputVar"]], customAlgoData[customAlgo["matchupVariableName"]], customAlgo["name"], customAlgo["outputVar"],
                savePath=savePath,units = customAlgo["units"],variable=customAlgo['variable'],stats = basicMetrics,stats_p = True);


    basicMetricsMap = {};
    for i in range(0, len(algorithmOutputs)):
        basicMetrics = metrics.calc_basic_metrics(algorithmOutputs[i], dataUsedList[i], settings);
        basicMetricsMap[algorithmOutputs[i]["name"]] = basicMetrics;

#add two more dictionary entries then do the pH conversions
# import copy
# x=copy.copy(basicMetricsMap["cmems_ph"])
# y=copy.copy(basicMetricsMap["ethz_ph"])
# basicMetricsMap["cmems_ph_Hion_to_pH"] =x
# basicMetricsMap["ethz_ph_Hion_to_pH"] =y
#
# #pH conversions for these variables
# basicMetricsMap["ethz_ph_Hion_to_pH"]["model_output_mean"]=-np.log10(basicMetricsMap["ethz_ph"]["model_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["model_output_sd"]=-np.log10(basicMetricsMap["ethz_ph"]["model_output_mean"]+basicMetricsMap["ethz_ph"]["reference_output_sd"])+np.log10(basicMetricsMap["ethz_ph"]["model_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["reference_output_mean"]=-np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["reference_output_sd"]=-np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"]+basicMetricsMap["ethz_ph"]["reference_output_sd"])+np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["bias"]=-np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"]+basicMetricsMap["ethz_ph"]["bias"])+np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["mad"]=-np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"]+basicMetricsMap["ethz_ph"]["mad"])+np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"])
# basicMetricsMap["ethz_ph_Hion_to_pH"]["rmsd"]=-np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"]+basicMetricsMap["ethz_ph"]["rmsd"])+np.log10(basicMetricsMap["ethz_ph"]["reference_output_mean"])
#
#
# basicMetricsMap["cmems_ph_Hion_to_pH"]["model_output_mean"]=-np.log10(basicMetricsMap["cmems_ph"]["model_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["model_output_sd"]=-np.log10(basicMetricsMap["cmems_ph"]["model_output_mean"]+basicMetricsMap["cmems_ph"]["reference_output_sd"])+np.log10(basicMetricsMap["cmems_ph"]["model_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["reference_output_mean"]=-np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["reference_output_sd"]=-np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"]+basicMetricsMap["cmems_ph"]["reference_output_sd"])+np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["bias"]=-np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"]+basicMetricsMap["cmems_ph"]["bias"])+np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["mad"]=-np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"]+basicMetricsMap["cmems_ph"]["mad"])+np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"])
# basicMetricsMap["cmems_ph_Hion_to_pH"]["rmsd"]=-np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"]+basicMetricsMap["cmems_ph"]["rmsd"])+np.log10(basicMetricsMap["cmems_ph"]["reference_output_mean"])


#####Write basicMetrics object


#binary pickle dump
pickle.dump(basicMetricsMap, open(path.join(outputRoot, "basic_metrics_ETHZ.pickle"), 'wb'));
#basicMetricsRead = pickle.load(open(path.join(outputRoot, "basic_metrics_ETHZ.pickle"), 'rb'));

#For json file, remove long vectors
##For json files, convert pd.Series types to nunpy array for json serialisation
for key in basicMetricsMap.keys():
    del basicMetricsMap[key]["model_output"];
    del basicMetricsMap[key]["reference_output_uncertainty"];
    del basicMetricsMap[key]["weights"];
    del basicMetricsMap[key]["model_uncertainty"];


with open(path.join(outputRoot, "basic_metrics_ETHZ.json"), 'a') as file:
    json.dump(basicMetricsMap, file, indent=4);
