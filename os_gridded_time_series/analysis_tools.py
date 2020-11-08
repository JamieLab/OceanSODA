#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:49:57 2020

@author: tom holding
"""

from os import path;
import pandas as pd;
import matplotlib.pyplot as plt;



def prediction_accuracy_plot(insituOutput, modelOutput, title, outputVariable):
    guidelineX1 = min(insituOutput.min(), modelOutput.min())#*0.9;
    guidelineX2 = max(insituOutput.max(), modelOutput.max())#*1.1;
    
    plt.figure();
    plt.plot([guidelineX1, guidelineX2], [guidelineX1, guidelineX2]);
    plt.scatter(insituOutput, modelOutput);
    plt.xlabel("in situ %s (umol/kg)" % outputVariable);
    plt.ylabel("predicted %s (umol/kg)" % outputVariable);
    plt.title(title);


if __name__ == "__main__":
######tmp scratchpad
    import osoda_global_settings;
    
    
    settings = osoda_global_settings.get_default_settings();
    tmp = osoda_global_settings.get_dic_algorithm_list();
    
    algoNames = [algo.__name__ for algo in settings["dic_algorithms"]];
    algoNames = ["Ternon2000_dic"];
    
    for algoName in algoNames:
        #read algo data
        readPath = path.join(settings["outputPath"], "global", "matchup_appended_"+algoName+".csv");
        try:
            df = pd.read_csv(readPath, sep=",", index_col=0);
        except FileNotFoundError:
            continue;
        
        prediction_accuracy_plot(df["DIC"], df["DIC_pred"], algoName, "DIC");
    
