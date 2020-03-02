#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:29:41 2020

@author: tom holding
"""

from string import Template;
import numpy as np;
import pandas as pd;
from utilities import run_sea_carb;
from os import path;

from rpy2.robjects import pandas2ri;
from rpy2.robjects.conversion import localconverter;
import rpy2;
import rpy2.robjects as ro
from rpy2.robjects.packages import importr

try:
    seacarb = importr("seacarb");
except rpy2.RRuntimeError:
    #install seacarb in the rpy2 version of r
    print("Installing R package: seacarb");
    utils = importr('utils');
    utils.install_packages('seacarb', repos='https://cloud.r-project.org');
    seacarb = importr("seacarb");

reefDataPathTemplate = Template("../output/reef_time_series/individual/reef_${ID}_${REGION}.csv");

reefDataPath = reefDataPathTemplate.safe_substitute(ID="510", REGION="oceansoda_amazon_plume");
inputData = pd.read_csv(reefDataPath);

#convert to rpy objects
with localconverter(ro.default_converter + pandas2ri.converter):
    r_at = ro.conversion.py2rpy(inputData["AT_pred"]);
    r_dic = ro.conversion.py2rpy(inputData["DIC_pred"]);
    r_sss = ro.conversion.py2rpy(inputData["SSS"]);
    r_sst = ro.conversion.py2rpy(inputData["SST"]);
    r_time = ro.conversion.py2rpy(inputData["time"]);
    #base = importr('base');
    #print(base.summary(r_dic));


tflag=15; #ALK and DIC

TODO: run seacarb directly using the 'seacarb' importr object
seacarb.carb(......);
#output = run_sea_carb(tflag, r_at, r_dic, r_sss, r_sst, verbose = False);



