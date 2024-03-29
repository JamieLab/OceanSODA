#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 07:55:00 2020

Example driver script to run all parts of the OceanSODA algorithm comparison,
gridded carbonate system time series, Amazon DIC outflow case study, and 
the coral reef vulnerability case study.

@author: tom holding
"""

from string import Template;
import osoda_global_settings;
from os import path;

import pandas as pd;
pd.set_option("mode.chained_assignment", None);

settings = osoda_global_settings.get_default_settings();

#########
# Compare algorithm performance using matchup data set 
# Compute all metrics and determine the 'best' and 'long' optimal algorithm for DIC and AT

import osoda_algorithm_comparison;
osoda_algorithm_comparison.main(settings);       

# # # ##########
# # # Download all prediction data sets and calculate gridded time series predictions

import osoda_calculate_gridded_predictions;
years = settings["years"];
regions = settings["regions"];
regionMaskPath= settings["regionMasksPath"];

#Run for the 'optimal algorithms' (aka min 8 year time series and n=30 matchups) - this is main dataset run. 
optAlgoTableLong = "output/algo_metrics/overall_best_algos_min_years=8.csv";
griddedTimeSeriesOutputPathLong = settings["longGriddedTimeSeriesPathTemplate"];
osoda_calculate_gridded_predictions.main(optAlgoTableLong, griddedTimeSeriesOutputPathLong, years, regions, regionMaskPath);


########
# Calculate DIC outflow for the Amazon (case study 1)
import osoda_dic_outflow;
osodaMasksPath = settings["regionMasksPath"];
precomputedGridAreaPath = settings["gridAreasPath"];

regions = ["oceansoda_amazon_plume"];

# #Run for the 'optimal algorithms' (aka min 8 year time series and n=30 matchups) - this is main dataset run. 
carbonateParametersTemplateBest = settings["longGriddedTimeSeriesPathTemplate"];
outputDirBest = path.join(settings["outputPathRoot"], "dic_outflow_amazon_best");
osoda_dic_outflow.main(carbonateParametersTemplateBest, outputDirBest, regions, osodaMasksPath, precomputedGridAreaPath);

##########
# # Assess and identify vulneral reefs (case study 2)
# import osoda_reef_vulnerability;
# osoda_reef_vulnerability.main(settings, useDistanceToLandMask=False);

import subprocess
subprocess.call("osoda_plot.py", shell=True)

