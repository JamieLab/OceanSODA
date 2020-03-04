#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 15:12:50 2019

@author: tom holding
"""

from string import Template;
from os import path;

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
    settings["matchupDatasetTemplate"] = Template(path.join("../../matchup_datasets/v1_with_DO", "${YYYY}_soda_mdb.nc")); #Location of the matchup dataset
    settings["regionMasksPath"] = path.join("../region_masks/osoda_region_masks_v2.nc"); #Where OSODA region masks can be found
    settings["outputPathRoot"] = "../output"; #Where outputs are written
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
    settings["years"] = range(2002, 2016); #[2010, 2011, 2012, 2013, 2014, 2016]; #range(2010, 2018); 
    
#    settings["usingTestMatchupDataset"] = False; #If true, this indicates the whole matchup dataset isn't being used
#                                                #and the driver script will skip subsetting for algorithm specific restrictions
    
    settings["algorithmInternalSpatialMasks"] = False; #Allow algorithms to use their own spatial masks to further subset matchup data (e.g. in addition to the OceanSODA region masks)
    settings["useErrorRatios"] = True;
    
    #nominal 'state-of-the-art' in situ measurement errors used in PATHFINDERS
    settings["insituError"] = {"DIC":2.5,
                               "AT":2.5};
    
    #Nominal 'state-of-the-art' in situ measurement errors as a percentage of the measured value. These are the same as used for PATHFINDERS
    settings["insituErrorRatio"] = {"DIC": 0.005, #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
                                    "AT": 0.005}; #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
    
    #map between common name and netCDF variablename:
    settings["variableToDatabaseMap"] = {"date": "time",
                             "lon": "lon",
                             "lat": "lat",
                             "SST": ["OSTIA-ESACCI-L4-v02.1_mean", "SST_mean"],
                             "SSS": ["ISAS15PSAL_mean", "ESACCI-SSS-L4-SSS-MERGED-OI-7DAY-RUNNINGMEAN-DAILY-25km_mean", "SSS_mean"],
                             "OC": "ESACCI-OC-L3S-OC_PRODUCTS-MERGED-1D_DAILY_4km_GEO_PML_OCx_QAA_mean", #ocean colour
                             "DO": "WOA18_Oxygen_o_an",
                             "NO3": "WOA18_Nitrate_n_an",
                             "PO4": "WOA18_Phosphate_p_an",
                             "SiO4": "WOA18_Silicate_i_an",
                             "DIC": "DIC_mean",
                             "AT": "AT_mean",
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






