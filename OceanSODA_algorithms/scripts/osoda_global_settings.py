#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:12:50 2019

@author: tom holding
"""

from string import Template;
from os import path;

#encapsulates meta data about a data set, including file path info
class DatasetInfo:
    def __init__(self, commonName, datasetName, matchupVariableName, matchupDatabaseTemplate, matchupDatabaseError=None, predictionDatasetTemplate=None, predictionDatasetVariable=None, predictionDatasetError=None):
        self.commonName = commonName;
        self.datasetName = datasetName;
        self.matchupVariableName = matchupVariableName;
        self.matchupDatabaseTemplate = matchupDatabaseTemplate;
        self.predictionDatasetTemplate = predictionDatasetTemplate;
        self.predictionDatasetVariable = predictionDatasetVariable;
        self.predictionDatasetError = predictionDatasetError;
    
    def __str__(self):
        return "DatasetInfo: "+" ".join([self.commonName, self.dataset]);


#These are the names used within the scripts to refer to input and output variables.
#settings["columnMap"] maps these common names to dataset columns allowing different input datasets to be used
#   with different settings objects.
#Defined here for convenience. Your in-code values should match these, and you should change the settings["columnMap"] values
#   to refer to those in your input datasets.
COMMON_NAMES = ["date", #datetime object (converted from time in seconds since 1980-01-01)
                "lon", #longitude (degrees East)
                "lat", #latitude (degrees North)
                "SST", #sea surface temperature ()
                "SSS", #sea surface salinity (PSU?)
                "OC", #ocean colour ()
                "NO3", #Nitrate ()
                "DO", #Disolved oxygen ()
                "SiO4", #Silicate ()
                "PO4", #Phosphate ()
                "DIC", #dissolved inorganic carbon ()
                "AT", #total alkalinity ()
                "DIC_pred", #predicted DIC ()
                "AT_pred", #predicted AT ()
                ];

#Returns a dictionary containing the global settings
def get_default_settings():
    projectRoot = "../../"; #TODO: check these are correct
    repoRoot = "../";
    
    settings = {};
    settings["matchupDatasetTemplate"] = Template(path.join(projectRoot, "matchup_datasets/v1_with_DO", "${YYYY}_soda_mdb.nc")); #Location of the matchup dataset
    settings["regionMasksPath"] = path.join(repoRoot, "region_masks/osoda_region_masks_v2.nc"); #Where OSODA region masks can be found
    settings["outputPathRoot"] = path.join(repoRoot, "output"); #Where outputs are written
    settings["outputPathMetrics"] = path.join(settings["outputPathRoot"], "metric_outputs"); #Where algorithm metrics outputs are written
    
    #Mask file and variable for the depth mask and distance to coast mask. Set to None to turn off.
    settings["depthMaskPath"] = None;#path.join("../region_masks/depth_mask_500m.nc");
    settings["depthMaskVar"] = "depth_mask";
    settings["distToCoastMaskPath"] = None;#path.join("../region_masks/distance_to_land_mask_250km.nc");
    settings["distToCoastMaskVar"] = "dist_to_coast_mask";
    
    #Some algorithms use additional data sets (e.g. gridded masks). These can be specified here as key : value pairs,
    #   where the key is the class name of the algorithm and the value is the path to the data set
    settings["algorithmSpecificDataPaths"] = {"Hassoun2015_basins_at": "algorithms/algo_data/mediterranean_masks_tmh.nc",
                                              "Hassoun2015_basins_dic" : "algorithms/algo_data/mediterranean_masks_tmh.nc",
                                              "Lee2000_dic": "algorithms/algo_data/lee2000_masks_tmh.nc",
                                              "Lee2006_at": "algorithms/algo_data/lee2006_masks_tmh.nc",
                                              "Millero1998_at": "algorithms/algo_data/millero1998_masks_tmh.nc",
                                              "Rivaro2010_at": "algorithms/algo_data/mediterranean_masks_tmh.nc",
                                              "Sasse2013_at" : "algorithms/algo_data/sasse2013_masks_tmh.nc",
                                              "Sasse2013_dic" : "algorithms/algo_data/sasse2013_masks_tmh.nc",
                                              "Takahashi2013_at" : "algorithms/algo_data/takahashi2013_masks_tmh.nc",
                                              };
    
    #Which years to analysis for
    settings["years"] = range(1991, 2018); #[2010, 2011, 2012, 2013, 2014, 2016]; #range(2010, 2018); 
    
#    settings["usingTestMatchupDataset"] = False; #If true, this indicates the whole matchup dataset isn't being used
#                                                #and the driver script will skip subsetting for algorithm specific restrictions
    
    settings["algorithmInternalSpatialMasks"] = False; #Allow algorithms to use their own spatial masks to further subset matchup data (e.g. in addition to the OceanSODA region masks)
    settings["useErrorRatios"] = True; #rather than static error for in situ AT and DIC
    settings["assessUsingWeightedRMSDe"] = True; #When calculating the 'best' algorithm, should it use the weighted version of RMSDe? (True=yes, False=no, weighted RMSDe is not available for algorithms which do not report uncertainty)
    
    #nominal 'state-of-the-art' in situ measurement errors used in PATHFINDERS
    settings["insituError"] = {"DIC":2.5,
                               "AT":2.5};
    
    #Nominal 'state-of-the-art' in situ measurement errors as a percentage of the measured value. These are the same as used for PATHFINDERS
    settings["insituErrorRatio"] = {"DIC": 0.005, #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
                                    "AT": 0.005}; #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
    
    #Dictionary mapping common parameter names to dataset info
    #Each entry can be a single DatasetInfo object, or a list of them (where there is more than on data set for a single 'common name')
    settings["datasetInfoMap"] = {"date": DatasetInfo(commonName="date", datasetName="time", matchupVariableName="time", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "lon": DatasetInfo(commonName="lon", datasetName="lon", matchupVariableName="lon", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "lat": DatasetInfo(commonName="lat", datasetName="lat", matchupVariableName="lat", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "SST": [DatasetInfo(commonName="SST", datasetName="SST-ESACCI", matchupVariableName="OSTIA-ESACCI-L4-v02.1_mean", matchupDatabaseError="OSTIA-ESACCI-L4-v02.1_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/OISST_reynolds_SST/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_1.0x1.0.nc"), predictionDatasetVariable="sst_mean", predictionDatasetError="sst_stddev"),
                                           DatasetInfo(commonName="SST", datasetName="SST-insitu", matchupVariableName="SST_mean", matchupDatabaseError="SST_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/OISST_reynolds_SST/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_1.0x1.0.nc"), predictionDatasetVariable="sst_mean", predictionDatasetError="sst_stddev"),
                                           ],
                                   "SSS": [DatasetInfo(commonName="SSS", datasetName="SSS-ISAS", matchupVariableName="ISAS15PSAL_mean", matchupDatabaseError="ISAS15PSAL_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc"), predictionDatasetVariable="salinity_mean", predictionDatasetError="salinity_err"),
                                           DatasetInfo(commonName="SSS", datasetName="SSS-ESACCI", matchupVariableName="ESACCI-SSS-L4-SSS-MERGED-OI-7DAY-RUNNINGMEAN-DAILY-25km_mean", matchupDatabaseError="ESACCI-SSS-L4-SSS-MERGED-OI-7DAY-RUNNINGMEAN-DAILY-25km_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc"), predictionDatasetVariable="salinity_mean", predictionDatasetError="salinity_err"),
                                           DatasetInfo(commonName="SSS", datasetName="SSS-insitu", matchupVariableName="SSS_mean", matchupDatabaseError="SSS_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc"), predictionDatasetVariable="salinity_mean", predictionDatasetError="salinity_err"),
                                           ],
                                   "OC": DatasetInfo(commonName="OC", datasetName="OC-ESACCI", matchupVariableName="ESACCI-OC-L3S-OC_PRODUCTS-MERGED-1D_DAILY_4km_GEO_PML_OCx_QAA_mean", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]), #ocean colour
                                   "DO": DatasetInfo(commonName="DO", datasetName="DO-WOA", matchupVariableName="WOA18_Oxygen_o_an", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/WOA_dissolved_oxygen/woa18_all_o${MM}_01.nc"), predictionDatasetVariable="o_an", predictionDatasetError="o_sd"),
                                   "NO3": DatasetInfo(commonName="NO3", datasetName="NO3-WOA", matchupVariableName="WOA18_Nitrate_n_an", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/WOA_nitrate/woa18_all_n${MM}_01.nc"), predictionDatasetVariable="n_an", predictionDatasetError="n_sd"),
                                   "PO4": DatasetInfo(commonName="PO4", datasetName="PO4-WOA", matchupVariableName="WOA18_Phosphate_p_an", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/WOA_phosphate/woa18_all_p${MM}_01.nc"), predictionDatasetVariable="p_an", predictionDatasetError="p_sd"),
                                   "SiO4":DatasetInfo(commonName="SiO4", datasetName="SiO4-WOA", matchupVariableName="WOA18_Silicate_i_an", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template("../../prediction_datasets/WOA_silicate/woa18_all_i${MM}_01.nc"), predictionDatasetVariable="i_an", predictionDatasetError="i_sd"),
                                   "DIC": DatasetInfo(commonName="DIC", datasetName="DIC-matchup", matchupVariableName="DIC_mean", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "AT": DatasetInfo(commonName="AT", datasetName="AT-matchup", matchupVariableName="AT_mean", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   };
    
    ### Settings for prediction
    #Output location of gridded predicted timeseries
    settings["griddedPredictionOutputTemplate"] = Template(path.join(repoRoot, "output/gridded_predictions/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc"));
    
    #Input dataset locations for making predictions from. Dictionary containing inputParameterName:(netCDFVariableName, netCDFFileTemplate) where YYYY and MM are substituted for string representations of year and month
    settings["predictionDataPaths"] = {"SSS": ("salinity_mean", Template("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc")),
                                       "SSS_err": ("salinity_err", Template("../../prediction_datasets/smos_ifremer_salinity/processed/${YYYY}${MM}_smos_sss.nc")),
                                       "SST": ("sst_mean", Template("../../prediction_datasets/OISST_reynolds_SST/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_1.0x1.0.nc")),
                                       "DO": ("o_an", Template("../../prediction_datasets/WOA_dissolved_oxygen/woa18_all_o${MM}_01.nc")),
                                       "NO3": ("n_an", Template("../../prediction_datasets/WOA_nitrate/woa18_all_n${MM}_01.nc")),
                                       "PO4": ("p_an", Template("../../prediction_datasets/WOA_phosphate/woa18_all_p${MM}_01.nc")),
                                       "SiO4": ("i_an", Template("../../prediction_datasets/WOA_silicate/woa18_all_i${MM}_01.nc")),
                                       };
    
    #This defines both the region names and algorithms to use
    import algorithms.at_algorithms as at_algorithms;
    import algorithms.dic_algorithms as dic_algorithms;
    settings["algorithmRegionMapping"] = {"oceansoda_amazon_plume": [at_algorithms.Astor2017a_at,
                                                                     at_algorithms.Astor2017d_at,
                                                                     at_algorithms.Astor2017b_at,
                                                                     at_algorithms.Astor2017e_at,
                                                                     at_algorithms.Astor2017c_at,
                                                                     at_algorithms.Brewer1995_at,
                                                                     at_algorithms.Cai2010a_at,
                                                                     at_algorithms.Cooley2006a_at,
                                                                     at_algorithms.Goyet1998_at,
                                                                     at_algorithms.Lee2006_at,
                                                                     at_algorithms.Lefevre2010_at,
                                                                     at_algorithms.Millero1998_at,
                                                                     at_algorithms.Sasse2013_at,
                                                                     at_algorithms.Sasse2013_global_at,
                                                                     at_algorithms.Ternon2000_at,
                                                                     at_algorithms.Takahashi2013_at,
                                                                     dic_algorithms.Brewer1995_dic,
                                                                     dic_algorithms.Cooley2006a_dic,
                                                                     dic_algorithms.Lee2000_dic,
                                                                     dic_algorithms.Lefevre2010_dic,
                                                                     dic_algorithms.Lefevre2017_dic,
                                                                     dic_algorithms.Sasse2013_dic,
                                                                     dic_algorithms.Sasse2013_global_dic,
                                                                     dic_algorithms.Ternon2000_dic,
                                                                     ],
                                          "oceansoda_congo": [at_algorithms.Brewer1995_at,
                                                              at_algorithms.Goyet1998_at,
                                                              at_algorithms.Lee2006_at,
                                                              at_algorithms.Sasse2013_at,
                                                              at_algorithms.Sasse2013_global_at,
                                                              at_algorithms.Takahashi2013_at,
                                                              dic_algorithms.Bakker1999_lcr1_dic,
                                                              dic_algorithms.Bakker1999_lcr2_dic,
                                                              dic_algorithms.Bakker1999_outflow_dic,
                                                              dic_algorithms.Brewer1995_dic,
                                                              dic_algorithms.Lee2000_dic,
                                                              dic_algorithms.Sasse2013_dic,
                                                              dic_algorithms.Sasse2013_global_dic,
                                                              dic_algorithms.Vangriesheim2009_all_dic,
                                                              dic_algorithms.Vangriesheim2009_open_ocean_dic,
                                                              ],
                                          "oceansoda_mississippi": [at_algorithms.Astor2017a_at,
                                                                    at_algorithms.Astor2017d_at,
                                                                    at_algorithms.Astor2017b_at,
                                                                    at_algorithms.Astor2017e_at,
                                                                    at_algorithms.Astor2017c_at,
                                                                    at_algorithms.Brewer1995_at,
                                                                    at_algorithms.Cai2010b_at,
                                                                    at_algorithms.Huang2012_at,
                                                                    at_algorithms.Lee2006_at,
                                                                    at_algorithms.Millero1998_at,
                                                                    at_algorithms.Sasse2013_at,
                                                                    at_algorithms.Sasse2013_global_at,
                                                                    at_algorithms.Takahashi2013_at,
                                                                    dic_algorithms.Brewer1995_dic,
                                                                    dic_algorithms.Lee2000_dic,
                                                                    dic_algorithms.Sasse2013_dic,
                                                                    dic_algorithms.Sasse2013_global_dic,
                                                                    ],
                                          "oceansoda_st_lawrence": [at_algorithms.Brewer1995_at,
                                                                    at_algorithms.Cai2010c_at,
                                                                    at_algorithms.Cai2010d_at,
                                                                    at_algorithms.Cai2010e_at,
                                                                    at_algorithms.Corbiere2007_at,
                                                                    at_algorithms.Lee2006_at,
                                                                    at_algorithms.Millero1998_at,
                                                                    at_algorithms.Sasse2013_at,
                                                                    at_algorithms.Sasse2013_global_at,
                                                                    at_algorithms.Takahashi2013_at,
                                                                    dic_algorithms.Brewer1995_dic,
                                                                    dic_algorithms.Lee2000_dic,
                                                                    dic_algorithms.Sasse2013_dic,
                                                                    dic_algorithms.Sasse2013_global_dic,
                                                                    ],
                                          "oceansoda_mediterranean": [at_algorithms.Brewer1995_at,
                                                                      at_algorithms.CopinMontegut1993_at,
                                                                      at_algorithms.CopinMontegut2002_at,
                                                                      at_algorithms.Gemayel2015_at,
                                                                      at_algorithms.Hassoun2015_full_at,
                                                                      at_algorithms.Hassoun2015_basins_at,
                                                                      at_algorithms.Rivaro2010_at,
                                                                      at_algorithms.Schneider2007_at,
                                                                      at_algorithms.Touratier2009_at,
                                                                      at_algorithms.Touratier2011_at,
                                                                      dic_algorithms.Brewer1995_dic,
                                                                      dic_algorithms.CopinMontegut1993_dic,
                                                                      dic_algorithms.CopinMontegut2002a_dic,
                                                                      dic_algorithms.CopinMontegut2002b_dic,
                                                                      dic_algorithms.Gemayel2015_dic,
                                                                      dic_algorithms.Hassoun2015_full_dic,
                                                                      dic_algorithms.Hassoun2015_basins_dic,
                                                                      ],
                                          };
    
    ### Derived settings (do not change - these are lists/values derived from the settings above for convenience)
    settings["regions"] = settings["algorithmRegionMapping"].keys();
    
    
    return settings;


#TODO: Hannah - Paste your settings function here


#Helper function that searches algorithms.dic_algorithms and adds all child classes of BaseAlgorithm
def get_dic_algorithm_list():
    from algorithms.base_algorithm import BaseAlgorithm;
    import algorithms.dic_algorithms
    
    algoList = [];
    for attrName in dir(algorithms.dic_algorithms):
        attr = getattr(algorithms.dic_algorithms, attrName);
        if isinstance(attr, type):
            if issubclass(attr, BaseAlgorithm):
                if (attr == BaseAlgorithm) == False: #Don't add the base class itself
                    algoList.append(attr);
    return algoList;



#Helper function that searches algorithms.at_algorithms and adds all child classes of BaseAlgorithm
def get_at_algorithm_list():
    from algorithms.base_algorithm import BaseAlgorithm;
    import algorithms.at_algorithms
    
    algoList = [];
    for attrName in dir(algorithms.at_algorithms):
        attr = getattr(algorithms.at_algorithms, attrName);
        if isinstance(attr, type):
            if issubclass(attr, BaseAlgorithm):
                if (attr == BaseAlgorithm) == False: #Don't add the base class itself
                    algoList.append(attr);
    return algoList;






