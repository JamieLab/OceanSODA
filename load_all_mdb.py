# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 16:02:02 2022
@author: rps207
"""
#this is a simplified script for loading the entirety of the matchup database in to 
#perform statistics on it. This loops through all the combinations


#combination0__SST-ESACCI_SSS-ESACCI
# combination1__SST-CORA_SSS-ESACCI
# combination2__SST-OISST_SSS-ESACCI
# combination3__SST-ESACCI_SSS-CORA
# combination4__SST-CORA_SSS-CORA
# combination5__SST-OISST_SSS-CORA
# combination6__SST-ESACCI_SSS-RSS-SMAP
# combination7__SST-CORA_SSS-RSS-SMAP
# combination8__SST-OISST_SSS-RSS-SMAP
# combination9__SST-ESACCI_SSS-ISAS
# combination10__SST-CORA_SSS-ISAS
# combination11__SST-OISST_SSS-ISAS



import os_algorithms.utilities as utilities
import osoda_global_settings;
from netCDF4 import Dataset;
import numpy as np
settings = osoda_global_settings.get_default_settings(); #TODO: Hannah - Change this to call your settings function
specificVariableToDatabaseMaps, specificVariableToDatabaseMapNames = utilities.get_dataset_variable_map_combinations(settings);


for inputCombinationName, inputCombination in zip(specificVariableToDatabaseMapNames, specificVariableToDatabaseMaps):
    years = utilities.calculate_years_for_input_combination(settings, inputCombination);
    print(inputCombinationName)
    matchupData = utilities.load_matchup_to_dataframe(settings, inputCombination, years=years)
    
regionMaskPath=settings["regionMasksPath"]

## find all matchups in amazom region
region="oceansoda_amazon_plume"
if regionMaskPath != None:
    regionMaskNC = Dataset(regionMaskPath, 'r');
    subsetData_amazon = utilities.subset_from_mask(matchupData, regionMaskNC, region);
    regionMaskNC.close();    


np.nanstd(subsetData_amazon.AT)
np.nanstd(subsetData_amazon.DIC)
np.nanstd(subsetData_amazon.region_pco2w_mean)
np.nanstd(subsetData_amazon.region_ph_mean)
np.nanmean(subsetData_amazon.region_ph_mean)
# for ph in ha terms
(-np.log10(np.nanstd(10**-(subsetData_amazon.region_ph_mean))+np.nanmean(10**-(subsetData_amazon.region_ph_mean))))-(-np.log10(np.nanmean(10**-(subsetData_amazon.region_ph_mean))))
(-np.log10(np.nanmean(10**-(subsetData_amazon.region_ph_mean))))



## find all matchups in congo region
region="oceansoda_congo"
if regionMaskPath != None:
    regionMaskNC = Dataset(regionMaskPath, 'r');
    subsetData_congo = utilities.subset_from_mask(matchupData, regionMaskNC, region);
    regionMaskNC.close();    

np.nanstd(subsetData_congo.AT)
np.nanstd(subsetData_congo.DIC)
np.nanstd(subsetData_congo.region_pco2w_mean)
np.nanstd(subsetData_congo.region_ph_mean)
np.nanmean(subsetData_congo.region_ph_mean)

# for ph in ha terms
(-np.log10(np.nanstd(10**-(subsetData_congo.region_ph_mean))+np.nanmean(10**-(subsetData_congo.region_ph_mean))))-(-np.log10(np.nanmean(10**-(subsetData_congo.region_ph_mean))))
(-np.log10(np.nanmean(10**-(subsetData_congo.region_ph_mean))))

