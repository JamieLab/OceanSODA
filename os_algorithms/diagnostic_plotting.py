#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 13:49:57 2020

@author: tom holding
"""

from os import path;
import pandas as pd;
import matplotlib.pyplot as plt;
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from math import log10, floor

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)
# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-200, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=2048)


def prediction_accuracy_plot(insituOutput, modelOutput, title, outputVariable, savePath=None,units = 'umol/kg', variable = 'pCO$_{2 (sw)}$',stats_p=False,stats = []):

    mean = np.nanmean([np.nanmean(insituOutput),np.nanmean(modelOutput)])
    std = np.nanmean([np.nanstd(insituOutput),np.nanstd(modelOutput)])
    guidelineX1 = mean-5*std
    guidelineX2 = mean+5*std

    guidelineX1 = np.nanmin([np.nanmin(insituOutput),np.nanmin(modelOutput)])-std*2
    guidelineX2 = np.nanmax([np.nanmax(insituOutput),np.nanmax(modelOutput)]) + std*2
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1, 1, 1,projection='scatter_density')
    ax.plot([guidelineX1, guidelineX2], [guidelineX1, guidelineX2],'k-',linewidth=0.5);
    ax.plot([guidelineX1, guidelineX2],np.array([guidelineX1, guidelineX2])*stats['slope'] + stats['intercept'],'k--')
    density = ax.scatter_density(insituOutput, modelOutput,cmap = white_viridis,dpi=48);
    fig.colorbar(density, label='Number of points per pixel')
    # ax.scatter(insituOutput, modelOutput)
    ax.set_xlabel(f"in situ %s ({units})" % variable);
    ax.set_ylabel(f"predicted %s ({units})" % variable);
    ax.set_title(title);
    ax.set_xlim([guidelineX1, guidelineX2])
    ax.set_ylim([guidelineX1, guidelineX2])
    if stats_p:
        rmsd = str(round_sig(stats['rmsd'],2))
        bias = str(round_sig(stats['bias'],2));
        mad = str(round_sig(stats['mad'],2))
        slope = str(round_sig(stats['slope'],3));
        n_val = stats['n']
        ax.text(0.57,0.30,f'Unweighted Stats\nRMSD = {rmsd} {units}\nBias = {bias} {units}\nMAD = {mad} {units}\nSlope = {slope}\nN = {n_val}',transform=ax.transAxes,va='top',fontsize=12)

    if savePath is not None:
        fig.savefig(savePath);
        #fig.close();
        plt.close()

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
