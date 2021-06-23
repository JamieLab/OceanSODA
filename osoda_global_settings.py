 #!/usr/bin/env python2
# -*- coding: utf-8 -*-


from string import Template;
from os import path;

#encapsulates meta data about a data set, including file path info
class DatasetInfo:
    def __init__(self, commonName, datasetName, matchupVariableName, matchupDatabaseTemplate, matchupDatabaseError=None, predictionDatasetTemplate=None, predictionDatasetVariable=None, predictionDatasetError=None):
        self.commonName = commonName;
        self.datasetName = datasetName;
        self.matchupVariableName = matchupVariableName;
        self.matchupDatabaseError = matchupDatabaseError;
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
                "chla", #ocean chlorophyll alpha
                "NO3", #Nitrate ()
                "DO", #Disolved oxygen ()
                "SiO4", #Silicate ()
                "PO4", #Phosphate ()
                "DIC", #dissolved inorganic carbon ()
                "AT", #total alkalinity ()
                #"DIC_err", #Uncertainty (std dev) in matchup dataset reference DIC
                #"AT_err", #Uncertainty (std dev) in matchup dataset reference AT
                "DIC_pred", #predicted DIC ()
                "AT_pred", #predicted AT ()
                ];

#Returns a dictionary containing the global settings
def get_default_settings():
    settings = {};
    
    ######################
    ### Set file paths ###
    projectRoot = path.dirname(__file__);
    settings["projectRoot"] = projectRoot;
    
    settings["dataPathRoot"] = path.join(projectRoot, "data"); #Where data is stored. This data is not included in the repository (e.g. matchup database and gridded prediction data sets)
    settings["matchupDatasetTemplate"] = Template(path.join(settings["dataPathRoot"], "matchup_datasets/v3.2/data/", "OCEANSODA-MMDB-${YYYY}-fv01.nc")); #Location of the matchup dataset
    settings["predictionDatasetsRoot"] = path.join(settings["dataPathRoot"], "prediction_datasets");
    settings["riverDischargeDataRoot"] = path.join(settings["dataPathRoot"], "gauging_station_discharge");
    settings["reefLocationsDataPath"] = path.join(settings["dataPathRoot"], "reefbase", "ReefLocations.csv");
    
    settings["auxDataPathRoot"] = path.join(projectRoot, "aux_data"); #Aux data includes small data files (e.g. masks) which are included with the repository
    settings["regionMasksPath"] = path.join(settings["auxDataPathRoot"], "osoda_region_masks_v2.nc"); #Where OSODA region masks can be found    
    settings["gridAreasPath"] = path.join(settings["auxDataPathRoot"], "grid_areas_1.0x1.0.csv"); #Where OSODA region masks can be found    
    #Mask files and variable for the depth mask and distance to coast mask. Set the boolean 'subsetWithDistToCoast' and 'subsetWithDepthMask' below to turn on/off turn off.
    settings["depthMaskPath"] = path.join(settings["auxDataPathRoot"], "depth_mask_500m.nc");
    settings["depthMaskVar"] = "depth_mask";
    settings["distToCoastMaskPath"] = path.join(settings["auxDataPathRoot"], "distance_to_land_mask_300km.nc");
    settings["distToCoastMaskVar"] = "dist_to_coast_mask";
    
    
    settings["outputPathRoot"] = path.join(projectRoot, "output"); #Root directory to write outputs
    settings["outputPathMetrics"] = path.join(settings["outputPathRoot"], "algo_metrics"); #Where algorithm metrics outputs are written
    settings["bestGriddedTimeSeriesPathTemplate"] = Template(path.join(settings["outputPathRoot"], "gridded_predictions/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc")); #path for predicted gridded time series using the 'best' version of the optimal algorithms
    settings["longGriddedTimeSeriesPathTemplate"] = Template(path.join(settings["outputPathRoot"], "gridded_predictions_min_year_range/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc")); #path for predicted gridded time series using the 'long' version of the optimal algorithms
    #other output directories here
    
    settings["logDirectoryRoot"] = path.join(settings["outputPathRoot"], "logs"); #Log files written here
    
    
    #Some algorithms use additional data sets (e.g. gridded masks). These can be specified here as key : value pairs,
    #   where the key is the class name of the algorithm and the value is the path to the data set
    algoDataPath = path.join(settings["auxDataPathRoot"], "algo_specific_masks");
    settings["algorithmSpecificDataPaths"] = {"Hassoun2015_basins_at": path.join(algoDataPath, "mediterranean_masks_tmh.nc"),
                                              "Hassoun2015_basins_dic" : path.join(algoDataPath, "mediterranean_masks_tmh.nc"),
                                              "Lee2000_dic": path.join(algoDataPath, "lee2000_masks_tmh.nc"),
                                              "Lee2006_at": path.join(algoDataPath, "lee2006_masks_tmh.nc"),
                                              "Millero1998_at": path.join(algoDataPath, "millero1998_masks_tmh.nc"),
                                              "Rivaro2010_at": path.join(algoDataPath, "mediterranean_masks_tmh.nc"),
                                              "Sasse2013_at" : path.join(algoDataPath, "sasse2013_masks_tmh.nc"),
                                              "Sasse2013_dic" : path.join(algoDataPath, "sasse2013_masks_tmh.nc"),
                                              "Takahashi2013_at" : path.join(algoDataPath, "takahashi2013_masks_tmh.nc"),
                                              };
    
    ##############################################
    ### Parameters to control script behaviour ###
    #Which years to analysis for
    settings["years"] = range(1957, 2021);
    
    settings["algorithmInternalSpatialMasks"] = False; #Allow algorithms to use their own spatial masks to further subset matchup data (e.g. in addition to the OceanSODA region masks)
    settings["useErrorRatios"] = True; #rather than static error for in situ AT and DIC
    settings["assessUsingWeightedRMSDe"] = True; #When calculating the 'best' algorithm, should it use the weighted version of RMSDe? (True=yes, False=no, weighted RMSDe is not available for algorithms which do not report uncertainty)
    
    #Use the ocean depth and distance to coast masks to remove shallow/coastal data?
    settings["subsetWithDepthMask"] = False
    settings["subsetWithDistToCoast"] = False;
    
    #nominal 'state-of-the-art' in situ measurement errors used in PATHFINDERS
    #These are substituted when no other uncertainty data is availablefor DIC and AT in the matchup database.
    settings["insituError"] = {"DIC": 2.5,
                               "AT": 2.5};
    
    #Nominal 'state-of-the-art' in situ measurement errors as a percentage of the measured value. These are the same as used for PATHFINDERS
    #These are substituted when no other uncertainty data is availablefor DIC and AT in the matchup database.
    settings["insituErrorRatio"] = {"DIC": 0.005, #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
                                    "AT": 0.005}; #0.5% nominal 'state-of-the-art' errors: Bockmon, E.E. and Dickson, A.G., 2015. An inter-laboratory comparison assessing the quality of seawater carbon dioxide measurements. Marine Chemistry, 171, pp.36-43.
    
    
    ############################
    ### data set definitions ###
    #'Common' names refer to a specific ocean parameter (e.g. SST), for which there may be multiple datasets available (multiple dataset names/specifications).
    #'datasetInfoMap' defines the mapping between common parameter names and available datasets in the matchup database and gridded prediction data.
    #Each entry can be a single DatasetInfo object, or a list of them (where there is more than on data set for a single 'common name')
    settings["datasetInfoMap"] = {"date": DatasetInfo(commonName="date", datasetName="time", matchupVariableName="time", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "lon": DatasetInfo(commonName="lon", datasetName="lon", matchupVariableName="lon", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "lat": DatasetInfo(commonName="lat", datasetName="lat", matchupVariableName="lat", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "SST": [DatasetInfo(commonName="SST", datasetName="SST-ESACCI-OSTIA", matchupVariableName="cci_sst_mean", matchupDatabaseError="cci_sst_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/ESACCI_SST_OSTIA/processed/${YYYY}/ESACCI-SST-OSTIA-LT-v02.0-fv01.1_${YYYY}_${MM}_processed.nc")), predictionDatasetVariable="sst", predictionDatasetError="sst_err"),
                                           DatasetInfo(commonName="SST", datasetName="SST-CORA", matchupVariableName="cora_temperature_mean", matchupDatabaseError="cora_temperature_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/CORA_SSS_SST/processed/${YYYY}/OA_CORA5.2_${YYYY}_${MM}_processed.nc")), predictionDatasetVariable="TEMP", predictionDatasetError="TEMP_err"),
                                           DatasetInfo(commonName="SST", datasetName="SST-OISST", matchupVariableName="noaa_sst_mean", matchupDatabaseError="noaa_sst_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/OISST_SST/processed/${YYYY}/${YYYY}${MM}01_OCF-SST-GLO-1M-100-REYNOLDS_1.0x1.0.nc")), predictionDatasetVariable="sst_mean", predictionDatasetError="sst_stddev"),
                                           ],
                                   "SSS": [DatasetInfo(commonName="SSS", datasetName="SSS-ESACCI-SMOS", matchupVariableName="cci_sss_mean", matchupDatabaseError="cci_sss_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/ESACCI_SSS_SMOS/processed/${YYYY}/ESACCI-SEASURFACESALINITY-L4-CENTRED15Day_${YYYY}_${MM}_processed.nc")), predictionDatasetVariable="sss", predictionDatasetError="sss_err"),
                                           DatasetInfo(commonName="SSS", datasetName="SSS-CORA", matchupVariableName="cora_salinity_mean", matchupDatabaseError="cora_salinity_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/CORA_SSS_SST/processed/${YYYY}/OA_CORA5.2_${YYYY}_${MM}_processed.nc")), predictionDatasetVariable="PSAL", predictionDatasetError="PSAL_err"),
                                           DatasetInfo(commonName="SSS", datasetName="SSS-RSS-SMAP", matchupVariableName="remss_smap_sss_mean", matchupDatabaseError="remss_smap_sss_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/RSS_SMAP_SSS/processed/RSS_smap_SSS_L3_monthly_${YYYY}_${MM}_FNL_v04.0_processed.nc")), predictionDatasetVariable="sss", predictionDatasetError="sss_err"),
                                           DatasetInfo(commonName="SSS", datasetName="SSS-ISAS", matchupVariableName="isas15_salinity_mean", matchupDatabaseError="isas15_salinity_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/ISAS_SSS_SST/processed/${YYYY}/ISAS15_DM_${YYYY}_${MM}_processed.nc")), predictionDatasetVariable="PSAL", predictionDatasetError="PSAL_err"),
                                           ],
                                   "OC": DatasetInfo(commonName="OC", datasetName="OC-ESACCI", matchupVariableName="cci_oc_chloro-a_mean", matchupDatabaseError="cci_oc_chloro-a_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]), #ocean colour
                                   "chla": DatasetInfo(commonName="chla", datasetName="OC-ESACCI", matchupVariableName="cci_oc_chloro-a_mean", matchupDatabaseError="cci_oc_chloro-a_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "DO": DatasetInfo(commonName="DO", datasetName="DO-WOA", matchupVariableName="woa18_oxygen_o_an", matchupDatabaseError="woa18_oxygen_o_sd", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/WOA_data_sets/WOA_dissolved_oxygen/woa18_all_o${MM}_processed.nc")), predictionDatasetVariable="o_an", predictionDatasetError="o_uncertainty"),
                                   "NO3": DatasetInfo(commonName="NO3", datasetName="NO3-WOA", matchupVariableName="woa18_nitrate_n_an", matchupDatabaseError="woa18_nitrate_n_sd", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/WOA_data_sets/WOA_nitrate/woa18_all_n${MM}_processed.nc")), predictionDatasetVariable="n_an", predictionDatasetError="n_uncertainty"),
                                   "PO4": DatasetInfo(commonName="PO4", datasetName="PO4-WOA", matchupVariableName="woa18_phosphate_p_an", matchupDatabaseError="woa18_phosphate_p_sd", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/WOA_data_sets/WOA_phosphate/woa18_all_p${MM}_processed.nc")), predictionDatasetVariable="p_an", predictionDatasetError="p_uncertainty"),
                                   "SiO4":DatasetInfo(commonName="SiO4", datasetName="SiO4-WOA", matchupVariableName="woa18_silicate_i_an", matchupDatabaseError="woa18_silicate_i_sd", matchupDatabaseTemplate=settings["matchupDatasetTemplate"], predictionDatasetTemplate=Template(path.join(projectRoot, "data", "prediction_datasets/WOA_data_sets/WOA_silicate/woa18_all_i${MM}_processed.nc")), predictionDatasetVariable="i_an", predictionDatasetError="i_uncertainty"),
                                   "DIC": DatasetInfo(commonName="DIC", datasetName="DIC-matchup", matchupVariableName="region_dic_mean", matchupDatabaseError="region_dic_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "AT": DatasetInfo(commonName="AT", datasetName="AT-matchup", matchupVariableName="region_at_mean", matchupDatabaseError="region_at_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "region_ph_mean": DatasetInfo(commonName="AT", datasetName="pH-matchup", matchupVariableName="region_ph_mean", matchupDatabaseError="region_ph_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   "region_pco2w_mean": DatasetInfo(commonName="AT", datasetName="pCO2-matchup", matchupVariableName="region_pco2w_mean", matchupDatabaseError="region_pco2w_stddev", matchupDatabaseTemplate=settings["matchupDatasetTemplate"]),
                                   };
    
    
    #######################
    # Defines both the region names and algorithms to use
    import os_algorithms.at_algorithms as at_algorithms;
    import os_algorithms.dic_algorithms as dic_algorithms;
    settings["algorithmRegionMapping"] = {
                                             "oceansoda_amazon_plume": [at_algorithms.Astor2017a_at,
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
                                          "oceansoda_mediterranean": [#at_algorithms.Brewer1995_at,
                                                                      at_algorithms.CopinMontegut1993_at,
                                                                      at_algorithms.CopinMontegut2002_at,
                                                                      at_algorithms.Gemayel2015_at,
                                                                      at_algorithms.Hassoun2015_full_at,
                                                                      at_algorithms.Hassoun2015_basins_at,
                                                                      at_algorithms.Rivaro2010_at,
                                                                      at_algorithms.Schneider2007_at,
                                                                      at_algorithms.Touratier2009_at,
                                                                      at_algorithms.Touratier2011_at,
                                                                      #dic_algorithms.Brewer1995_dic,
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
    
    
    ##########Refactoring ended here.
    ### Settings for prediction ###TODO: are these not needed now?
    #Output location of gridded predicted timeseries
    settings["griddedPredictionOutputTemplate"] = Template(path.join(projectRoot, "output/gridded_predictions/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc"));
    settings["griddedPredictionMinYearsOutputTemplate"] = Template(path.join(projectRoot, "output/gridded_predictions_min_year_range/gridded_${REGION}_${LATRES}x${LONRES}_${OUTPUTVAR}.nc"));
    
    return settings;



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






