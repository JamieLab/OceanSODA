#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 12:38:19 2020

@author: tom holding
"""

import pickle;
import numpy as np;

#For amazon find lowest RMSDe algorithm and report the RMSDe and bias

f = open("/home/verwirrt/Projects/Work/20190816_OceanSODA/metric_outputs/oceansoda_amazon_plume/algorithms_used.csv", 'r');
algoNames = f.readline().split(',')[:-1];

metrics = pickle.load(open("/home/verwirrt/Projects/Work/20190816_OceanSODA/metric_outputs/oceansoda_amazon_plume/basic_metrics.json", 'rb'));

wrmsds = np.genfromtxt("/home/verwirrt/Projects/Work/20190816_OceanSODA/metric_outputs/oceansoda_amazon_plume/paired_wrmsd_matrix.csv", delimiter=',');
wscores = np.genfromtxt("/home/verwirrt/Projects/Work/20190816_OceanSODA/metric_outputs/oceansoda_amazon_plume/paired_wscore_matrix.csv", delimiter=',');
nintersect = np.genfromtxt("/home/verwirrt/Projects/Work/20190816_OceanSODA/metric_outputs/oceansoda_amazon_plume/n_intersect_matrix.csv", delimiter=',');

finalScores = np.nanmean(wscores, axis=1);
rmsd_rep_loc = np.array([finalScores[i] / metrics[i]["n"] for i in range(len(algoNames))]).argmin();

wrmsdIndiv = np.array([metrics[i]["wrmsd"] for i in range(len(algoNames))]);

rmsd_rep = wrmsdIndiv[rmsd_rep_loc];

rmsdes = rmsd_rep*finalScores/finalScores[rmsd_rep_loc];
rmsdes[rmsd_rep_loc] = rmsd_rep;


print("RMSDe values for Amazon plume region using SMOS salinity:");
for i in range(len(algoNames)):
    print("\t"+algoNames[i]+":", rmsdes[i], "(bias: "+str(metrics[i]["bias"])+")");