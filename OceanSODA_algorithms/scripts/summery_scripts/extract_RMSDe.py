#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:07:15 2020

@author: verwirrt
"""

import pandas as pd;
from os import path;
from string import Template;

inputPathTemplate = Template("/home/verwirrt/Projects/Work/20190816_OceanSODA/OceanSODA_algorithms/output/metric_outputs/$OOUTPUTVARP/${REGION}/final_scores.csv");

regions = ['oceansoda_amazon_plume', 'oceansoda_congo', 'oceansoda_st_lawrence', 'oceansoda_mississippi', 'oceansoda_mediterranean'];
outputVars = ["AT", "DIC"];

for outputVar in outputVars:
    outdf = pd.DataFrame();
    for region in regions:
        inputPath = inputPathTemplate.safe_substitute(OUTPUTVAR=outputVar, REGION=region);
        inputdf = pd.read_csv(inputPath);
        
        