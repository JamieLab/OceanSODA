# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 17:22:37 2021

@author: rps207
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 14:49:35 2020

@author: tom holding
"""

from string import Template;
from netCDF4 import Dataset;
import numpy as np;
import matplotlib.pyplot as plt;
from datetime import datetime, timedelta;
from matplotlib import colorbar, colors;
import os
from os import path;
import cartopy.crs as ccrs;
import cartopy.feature as cfeature
from cartopy import config;
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER;
import cv2
import csv
import json


inputTemplate = Template("C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\output\\gridded_predictions_min_year_range\\gridded_${REGION}_1.0x1.0_${VAR}.nc");
figsize=(7,6.5);
figsize2=(figsize[0],figsize[1]*0.75);
labelsize = 15;

def extract_var_means(nc, var):
    def calc_annual_mean(vals, time):
        groupLabels = np.array([dt.year for dt in time]);
        years = np.unique(groupLabels);
        
        annualMeanVals = [];
        for year in years:
            w = np.where(groupLabels == year);
            if np.any(np.isfinite(vals[w])):
                annualMeanVals.append(np.nanmean(vals[w]));
            else:
                annualMeanVals.append(np.nan);
        
        annualMeanVals = np.array(annualMeanVals);
        yearDTs = np.array([datetime(year, 6, 1) for year in years]);
        return annualMeanVals, yearDTs;
    
    var = nc.variables[var][:];
    plumeMask = np.zeros(var.shape, dtype=int);
    plumeMask[nc.variables["SSS"][:] < 35] = 1;
    
    plumeVar = var.copy();
    plumeVar[plumeMask!=1] = np.nan;
    meanPlumeVar = np.nanmean(plumeVar, axis=(1,2));
    del plumeVar;
    
    notPlumeVar = var.copy();
    notPlumeVar[plumeMask==1] = np.nan;
    meanNotPlumeVar = np.nanmean(notPlumeVar, axis=(1,2));
    del notPlumeVar;
    
    meanAllVar = np.nanmean(var, axis=(1,2));
    
    #annual means
    time = get_datetimes(nc["time"][:]);
    
    annualMeanPlumeVar, years = calc_annual_mean(meanPlumeVar, time);
    annualMeanNotPlumeVar, years = calc_annual_mean(meanNotPlumeVar, time);
    annualMeanAllVar, years = calc_annual_mean(meanAllVar, time);
    
    return meanPlumeVar, meanNotPlumeVar, meanAllVar, years, annualMeanPlumeVar, annualMeanNotPlumeVar, annualMeanAllVar;

def extract_var_stds(nc, var):
    def calc_annual_std(vals, time):
        groupLabels = np.array([dt.year for dt in time]);
        years = np.unique(groupLabels);
        
        annualstdVals = [];
        for year in years:
            w = np.where(groupLabels == year);
            if np.any(np.isfinite(vals[w])):
                annualstdVals.append(np.nanstd(vals[w]));
            else:
                annualstdVals.append(np.nan);
        
        annualstdVals = np.array(annualstdVals);
        yearDTs = np.array([datetime(year, 6, 1) for year in years]);
        return annualstdVals, yearDTs;
    
    var = nc.variables[var][:];
    plumeMask = np.zeros(var.shape, dtype=int);
    plumeMask[nc.variables["SSS"][:] < 35] = 1;
    
    plumeVar = var.copy();
    plumeVar[plumeMask!=1] = np.nan;
    stdPlumeVar = np.nanstd(plumeVar, axis=(1,2));
    del plumeVar;
    
    notPlumeVar = var.copy();
    notPlumeVar[plumeMask==1] = np.nan;
    stdNotPlumeVar = np.nanstd(notPlumeVar, axis=(1,2));
    del notPlumeVar;
    
    stdAllVar = np.nanstd(var, axis=(1,2));
    
    #annual stds
    time = get_datetimes(nc["time"][:]);
    
    annualstdPlumeVar, years = calc_annual_std(stdPlumeVar, time);
    annualstdNotPlumeVar, years = calc_annual_std(stdNotPlumeVar, time);
    annualstdAllVar, years = calc_annual_std(stdAllVar, time);
    
    return stdPlumeVar, stdNotPlumeVar, stdAllVar, years, annualstdPlumeVar, annualstdNotPlumeVar, annualstdAllVar;

def extract_var_mins(nc, var):
    def calc_annual_min(vals, time):
        groupLabels = np.array([dt.year for dt in time]);
        years = np.unique(groupLabels);
        
        annualminVals = [];
        for year in years:
            w = np.where(groupLabels == year);
            if np.any(np.isfinite(vals[w])):
                annualminVals.append(np.nanmin(vals[w]));
            else:
                annualminVals.append(np.nan);
        
        annualminVals = np.array(annualminVals);
        yearDTs = np.array([datetime(year, 6, 1) for year in years]);
        return annualminVals, yearDTs;
    
    var = nc.variables[var][:];
    plumeMask = np.zeros(var.shape, dtype=int);
    plumeMask[nc.variables["SSS"][:] < 35] = 1;
    
    plumeVar = var.copy();
    plumeVar[plumeMask!=1] = np.nan;
    minPlumeVar = np.nanmin(plumeVar, axis=(1,2));
    del plumeVar;
    
    notPlumeVar = var.copy();
    notPlumeVar[plumeMask==1] = np.nan;
    minNotPlumeVar = np.nanmin(notPlumeVar, axis=(1,2));
    del notPlumeVar;
    
    minAllVar = np.nanmin(var, axis=(1,2));
    
    #annual mins
    time = get_datetimes(nc["time"][:]);
    
    annualminPlumeVar, years = calc_annual_min(minPlumeVar, time);
    annualminNotPlumeVar, years = calc_annual_min(minNotPlumeVar, time);
    annualminAllVar, years = calc_annual_min(minAllVar, time);
    
    return minPlumeVar, minNotPlumeVar, minAllVar, years, annualminPlumeVar, annualminNotPlumeVar, annualminAllVar;

def extract_var_maxs(nc, var):
    def calc_annual_max(vals, time):
        groupLabels = np.array([dt.year for dt in time]);
        years = np.unique(groupLabels);
        
        annualmaxVals = [];
        for year in years:
            w = np.where(groupLabels == year);
            if np.any(np.isfinite(vals[w])):
                annualmaxVals.append(np.nanmax(vals[w]));
            else:
                annualmaxVals.append(np.nan);
        
        annualmaxVals = np.array(annualmaxVals);
        yearDTs = np.array([datetime(year, 6, 1) for year in years]);
        return annualmaxVals, yearDTs;
    
    var = nc.variables[var][:];
    plumeMask = np.zeros(var.shape, dtype=int);
    plumeMask[nc.variables["SSS"][:] < 35] = 1;
    
    plumeVar = var.copy();
    plumeVar[plumeMask!=1] = np.nan;
    maxPlumeVar = np.nanmax(plumeVar, axis=(1,2));
    del plumeVar;
    
    notPlumeVar = var.copy();
    notPlumeVar[plumeMask==1] = np.nan;
    maxNotPlumeVar = np.nanmax(notPlumeVar, axis=(1,2));
    del notPlumeVar;
    
    maxAllVar = np.nanmax(var, axis=(1,2));
    
    #annual maxs
    time = get_datetimes(nc["time"][:]);
    
    annualmaxPlumeVar, years = calc_annual_max(maxPlumeVar, time);
    annualmaxNotPlumeVar, years = calc_annual_max(maxNotPlumeVar, time);
    annualmaxAllVar, years = calc_annual_max(maxAllVar, time);
    
    return maxPlumeVar, maxNotPlumeVar, maxAllVar, years, annualmaxPlumeVar, annualmaxNotPlumeVar, annualmaxAllVar;

def get_time_index(target, firstIndexSecs):
    baseDatetime = datetime(1980, 1, 1);
    firstDatetime = baseDatetime + timedelta(seconds=firstIndexSecs);
    
    #assumes monthly temporal resolution
    targetIndex = (target.year - firstDatetime.year)*12 + (target.month-firstDatetime.month);
    return targetIndex;


def get_datetimes(secondsSince1980):
    base = datetime(1980,1,1);
    return np.array([base+timedelta(seconds=int(t)) for t in secondsSince1980]);

def get_datetime(secondsSince1980):
    baseDate = datetime(1980,1,1);
    return baseDate+timedelta(seconds=secondsSince1980);

def get_extents(data, lats, lons, pad=2):
    w = np.where(np.isfinite(data));
    latMin = lats[min(w[0])];
    latMax = lats[max(w[0])];
    lonMin = lons[min(w[1])];
    lonMax = lons[max(w[1])];
    
    return (lonMin-pad, lonMax+pad, latMin-pad, latMax+pad);







#### look at regions individually

#at
atAmazonNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_amazon_plume", VAR="AT"));
meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, years, annualMeanPlumeATAmazon, annualMeanNotPlumeATAmazon, annualMeanAllATAmazon = extract_var_means(atAmazonNC, "AT");
dates = get_datetimes(atAmazonNC.variables["time"][:]);

# atMississippiNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mississippi", VAR="AT"));
# meanPlumeATMississippi, meanNotPlumeATMississippi, meanAllATMississippi, years, annualMeanPlumeATMississippi, annualMeanNotPlumeATMississippi, annualMeanAllATMississippi = extract_var_means(atMississippiNC, "AT");

atCongoNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_congo", VAR="AT"));
meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo, years, annualMeanPlumeATCongo, annualMeanNotPlumeATCongo, annualMeanAllATCongo = extract_var_means(atCongoNC,"AT");

# atStLawrenceNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_st_lawrence", VAR="AT"));
# meanPlumeATStLawrence, meanNotPlumeATStLawrence, meanAllATStLawrence, years, annualMeanPlumeATStLawrence, annualMeanNotPlumeATStLawrence, annualMeanAllATStLawrence = extract_var_means(atStLawrenceNC, "AT");

atMediterraneanNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mediterranean", VAR="AT"));
meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean, years, annualMeanPlumeATMediterranean, annualMeanNotPlumeATMediterranean, annualMeanAllATMediterranean = extract_var_means(atMediterraneanNC, "AT");

#min and max values for y scale:
#listAll = [meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, meanPlumeATMississippi, meanNotPlumeATMississippi, meanAllATMississippi, meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo, meanPlumeATStLawrence, meanNotPlumeATStLawrence, meanAllATStLawrence,meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean];
listAll = [meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean, meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo];

minY = np.nanmin(listAll)*0.95;
maxY = np.nanmax(listAll)*1.05;

#dic
dicAmazonNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_amazon_plume", VAR="DIC"));
meanPlumeDICAmazon, meanNotPlumeDICAmazon, meanAllDICAmazon, years, annualMeanPlumeDICAmazon, annualMeanNotPlumeDICAmazon, annualMeanAllDICAmazon= extract_var_means(dicAmazonNC, "DIC");
dates = get_datetimes(dicAmazonNC.variables["time"][:]);

# dicMississippiNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mississippi", VAR="DIC"));
# meanPlumeDICMississippi, meanNotPlumeDICMississippi, meanAllDICMississippi, years, annualMeanPlumeDICMississippi, annualMeanNotPlumeDICMississippi, annualMeanAllDICMississippi = extract_var_means(dicMississippiNC, "DIC");

dicCongoNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_congo", VAR="DIC"));
meanPlumeDICCongo, meanNotPlumeDICCongo, meanAllDICCongo,years, annualMeanPlumeDICCongo, annualMeanNotPlumeDICCongo, annualMeanAllDICCongo = extract_var_means(dicCongoNC, "DIC");

# dicStLawrenceNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_st_lawrence", VAR="DIC"));
# meanPlumeDICStLawrence, meanNotPlumeDICStLawrence, meanAllDICStLawrence, years, annualMeanPlumeDICStLawrence, annualMeanNotPlumeDICStLawrence, annualMeanAllDICStLawrence = extract_var_means(dicStLawrenceNC, "DIC");

dicMediterraneanNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mediterranean", VAR="DIC"));
meanPlumeDICMediterranean, meanNotPlumeDICMediterranean, meanAllDICMediterranean, years, annualMeanPlumeDICMediterranean, annualMeanNotPlumeDICMediterranean, annualMeanAllDICMediterranean = extract_var_means(dicMediterraneanNC, "DIC");

#min and max values for y scale:
#listAll = [meanPlumeDICAmazon, meanNotPlumeDICAmazon, meanAllDICAmazon, meanPlumeDICMississippi, meanNotPlumeDICMississippi, meanAllDICMississippi, meanPlumeDICCongo, meanNotPlumeDICCongo, meanAllDICCongo];
listAll = [meanPlumeDICAmazon, meanNotPlumeDICAmazon, meanAllDICAmazon, meanPlumeDICMediterranean, meanNotPlumeDICMediterranean, meanAllDICMediterranean, meanPlumeDICCongo, meanNotPlumeDICCongo, meanAllDICCongo];

minY = np.nanmin(listAll)*0.95;
maxY = np.nanmax(listAll)*1.05;



#61 years 12 months =  768 months 
#probably shouldn't hard code this in!
number_list_jfm=[]
number_list_amj=[]
number_list_jas=[]
number_list_ond=[]
A = np.arange(768)
for i in range(64):
    number_list_jfm=np.append(number_list_jfm,[i*12 + A[0:3:1]])
    number_list_amj=np.append(number_list_amj,[i*12 + A[3:6:1]])
    number_list_jas=np.append(number_list_jas,[i*12 + A[6:9:1]])
    number_list_ond=np.append(number_list_ond,[i*12 + A[9:12:1]])

number_list_jfm=number_list_jfm.astype(int)
number_list_amj=number_list_amj.astype(int)
number_list_jas=number_list_jas.astype(int)
number_list_ond=number_list_ond.astype(int)

    #### plot of AT in all three regions
fig1=plt.figure(figsize=figsize);
plt.subplot(3,1,1);
# plt.plot(years, annualMeanPlumeATAmazon, 'k--');
# plt.plot(years, annualMeanNotPlumeATAmazon, 'k:');
# plt.plot(years, annualMeanAllATAmazon, 'k');
plt.plot(dates, meanPlumeATAmazon, 'k--');
plt.plot(dates, meanNotPlumeATAmazon, 'k:');
plt.plot(dates, meanAllATAmazon, 'k');
plt.locator_params(axis='y', nbins=4);
plt.ylabel("TA ($\mu mol \ kg^{-1}$)", fontsize=labelsize);
plt.title('Amazon');

plt.subplot(3,1,2);
# plt.plot(years, annualMeanPlumeATCongo, 'k--');
# plt.plot(years, annualMeanNotPlumeATCongo, 'k:');
# plt.plot(years, annualMeanAllATCongo, 'k');
plt.plot(dates, meanPlumeATCongo, 'k--');
plt.plot(dates, meanNotPlumeATCongo, 'k:');
plt.plot(dates, meanAllATCongo, 'k');
plt.locator_params(axis='y', nbins=4);
plt.ylabel("TA ($\mu mol \ kg^{-1}$)", fontsize=labelsize);
plt.title('Congo');

plt.subplot(3,1,3);
# plt.plot(years, annualMeanPlumeATMediterranean, 'k--');
# plt.plot(years, annualMeanNotPlumeATMediterranean, 'k:');
# plt.plot(years, annualMeanAllATMediterranean, 'k');
plt.plot(dates, meanPlumeATMediterranean, 'k--');
plt.plot(dates, meanNotPlumeATMediterranean, 'k:');
plt.plot(dates, meanAllATMediterranean, 'k');
plt.locator_params(axis='y', nbins=4);
plt.xlabel("time (year)", fontsize=labelsize);
plt.ylabel("TA ($\mu mol \ kg^{-1}$)", fontsize=labelsize);
plt.title('Mediterranean');

plt.tight_layout();
plt.savefig("os_plots\\at_timeseries.pdf");
plt.savefig("os_plots\\at_timeseries.png");


    #### plot of DIC in all three regions
fig2=plt.figure(figsize=figsize2);
plt.subplot(3,1,1);
plt.plot(dates, meanPlumeDICAmazon, 'k--'); #plt.ylim(minY, maxY);
plt.plot(dates, meanNotPlumeDICAmazon, 'k:'); #plt.ylim(minY, maxY);
plt.plot(dates, meanAllDICAmazon, 'k'); #plt.ylim(minY, maxY);
plt.ylabel("DIC ($\mu mol \ kg^{-1}$)", fontsize=labelsize);
plt.title('Amazon');

plt.subplot(3,1,2);
plt.plot(dates, meanPlumeDICCongo, 'k--'); #plt.ylim(minY, maxY);
plt.plot(dates, meanNotPlumeDICCongo, 'k:'); #plt.ylim(minY, maxY);
plt.plot(dates, meanAllDICCongo, 'k'); #plt.ylim(minY, maxY);
plt.ylabel("DIC($\mu mol \ kg^{-1}$)", fontsize=labelsize);
plt.title('Congo');

plt.subplot(3,1,3);
plt.plot(dates, meanPlumeDICMediterranean, 'k--'); #plt.ylim(minY, maxY);
plt.plot(dates, meanNotPlumeDICMediterranean, 'k:'); #plt.ylim(minY, maxY);
plt.plot(dates, meanAllDICMediterranean, 'k'); #plt.ylim(minY, maxY);
plt.title('Mediterranean');
plt.ylabel("DIC ($\mu mol \ kg^{-1}$)", fontsize=labelsize);

plt.xlabel("time (year)", fontsize=labelsize);

plt.tight_layout();
plt.savefig("os_plots\\dic_timeseries.pdf");
plt.savefig("os_plots\\dic_timeseries.png");

del atAmazonNC,atCongoNC,atMediterraneanNC,dicAmazonNC,dicCongoNC,dicMediterraneanNC;


#### loop through regions

    #### load data, calculate stats and make tables

regions = ["oceansoda_amazon_plume", "oceansoda_congo", "oceansoda_mediterranean"];#, "oceansoda_st_lawrence"];
video_vars = ["DIC", "SSS","SST", "pH_free","pCO2", "saturation_aragonite", "saturation_calcite"];
create_animations="False";

for region in regions:
        #### Define the two netCDF files for each region
    dicNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="DIC"), 'r');
    atNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="AT"), 'r');
    time = get_datetimes(dicNC["time"][:]);
    
        #### Calculate stats for all key variables
    #DIC
    meanPlumeDIC, meanNotPlumeDIC, meanAllDIC, _, _, _, _ = extract_var_means(dicNC, "DIC");
    stdPlumeDIC, stdNotPlumeDIC, stdAllDIC, _, _, _, _ = extract_var_stds(dicNC, "DIC");
    minPlumeDIC, minNotPlumeDIC, minAllDIC, _, _, _, _ = extract_var_mins(dicNC, "DIC");
    maxPlumeDIC, maxNotPlumeDIC, maxAllDIC, _, _, _, _ = extract_var_maxs(dicNC, "DIC");
    meanPlumeDICuncertainty, meanNotPlumeDICuncertainty, meanAllDICuncertainty, _, _, _, _ = extract_var_means(dicNC, "DIC_Combined_uncertainty_dueto_RMSD_and_Bias");
    meanPlumeDICuncertainty_botup, meanNotPlumeDICuncertainty_botup, meanAllDICuncertainty_botup, _, _, _, _ = extract_var_means(dicNC, "DIC_pred_combined_uncertainty");

    #AT
    meanPlumeAT, meanNotPlumeAT, meanAllAT, _, _, _, _ = extract_var_means(atNC, "AT");
    stdPlumeAT, stdNotPlumeAT, stdAllAT, _, _, _, _ = extract_var_stds(atNC, "AT");
    minPlumeAT, minNotPlumeAT, minAllAT, _, _, _, _ = extract_var_mins(atNC, "AT");
    maxPlumeAT, maxNotPlumeAT, maxAllAT, _, _, _, _ = extract_var_maxs(atNC, "AT");
    meanPlumeATuncertainty, meanNotPlumeATuncertainty, meanAllATuncertainty, _, _, _, _ = extract_var_means(atNC, "AT_Combined_uncertainty_dueto_RMSD_and_Bias");
    meanPlumeATuncertainty_botup, meanNotPlumeATuncertainty_botup, meanAllATuncertainty_botup, _, _, _, _ = extract_var_means(atNC, "AT_pred_combined_uncertainty");

    #pH
    meanPlumepH_dic, meanNotPlumepH_dic, meanAllpH_dic, _, _, _, _ = extract_var_means(dicNC, "pH");
    stdPlumepH_dic, stdNotPlumepH_dic, stdAllpH_dic, _, _, _, _ = extract_var_stds(dicNC, "pH");
    minPlumepH_dic, minNotPlumepH_dic, minAllpH_dic, _, _, _, _ = extract_var_mins(dicNC, "pH");
    maxPlumepH_dic, maxNotPlumepH_dic, maxAllpH_dic, _, _, _, _ = extract_var_maxs(dicNC, "pH");
    meanPlumepHuncertainty_dic, meanNotPlumepHuncertainty_dic, meanAllpHuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "pH_uncertainty");
    
    meanPlumepH_at, meanNotPlumepH_at, meanAllpH_at, _, _, _, _ = extract_var_means(atNC, "pH");
    stdPlumepH_at, stdNotPlumepH_at, stdAllpH_at, _, _, _, _ = extract_var_stds(atNC, "pH");
    minPlumepH_at, minNotPlumepH_at, minAllpH_at, _, _, _, _ = extract_var_mins(atNC, "pH");
    maxPlumepH_at, maxNotPlumepH_at, maxAllpH_at, _, _, _, _ = extract_var_maxs(atNC, "pH");
    meanPlumepHuncertainty_at, meanNotPlumepHuncertainty_at, meanAllpHuncertainty_at, _, _, _, _ = extract_var_means(atNC, "pH_uncertainty");

    #hydrogen_free
    meanPlumehydrogen_free_dic, meanNotPlumehydrogen_free_dic, meanAllhydrogen_free_dic, _, _, _, _ = extract_var_means(dicNC, "hydrogen_free");
    stdPlumehydrogen_free_dic, stdNotPlumehydrogen_free_dic, stdAllhydrogen_free_dic, _, _, _, _ = extract_var_stds(dicNC, "hydrogen_free");
    minPlumehydrogen_free_dic, minNotPlumehydrogen_free_dic, minAllhydrogen_free_dic, _, _, _, _ = extract_var_mins(dicNC, "hydrogen_free");
    maxPlumehydrogen_free_dic, maxNotPlumehydrogen_free_dic, maxAllhydrogen_free_dic, _, _, _, _ = extract_var_maxs(dicNC, "hydrogen_free");
    meanPlumehydrogen_freeuncertainty_dic, meanNotPlumehydrogen_freeuncertainty_dic, meanAllhydrogen_freeuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "hydrogen_free_uncertainty");
    
    meanPlumehydrogen_free_at, meanNotPlumehydrogen_free_at, meanAllhydrogen_free_at, _, _, _, _ = extract_var_means(atNC, "hydrogen_free");
    stdPlumehydrogen_free_at, stdNotPlumehydrogen_free_at, stdAllhydrogen_free_at, _, _, _, _ = extract_var_stds(atNC, "hydrogen_free");
    minPlumehydrogen_free_at, minNotPlumehydrogen_free_at, minAllhydrogen_free_at, _, _, _, _ = extract_var_mins(atNC, "hydrogen_free");
    maxPlumehydrogen_free_at, maxNotPlumehydrogen_free_at, maxAllhydrogen_free_at, _, _, _, _ = extract_var_maxs(atNC, "hydrogen_free");
    meanPlumehydrogen_freeuncertainty_at, meanNotPlumehydrogen_freeuncertainty_at, meanAllhydrogen_freeuncertainty_at, _, _, _, _ = extract_var_means(atNC, "hydrogen_free_uncertainty");

    #CO3
    meanPlumeCO3_dic, meanNotPlumeCO3_dic, meanAllCO3_dic, _, _, _, _ = extract_var_means(dicNC, "CO3");
    stdPlumeCO3_dic, stdNotPlumeCO3_dic, stdAllCO3_dic, _, _, _, _ = extract_var_stds(dicNC, "CO3");
    
    meanPlumeCO3_at, meanNotPlumeCO3_at, meanAllCO3_at, _, _, _, _ = extract_var_means(atNC, "CO3");
    stdPlumeCO3_at, stdNotPlumeCO3_at, stdAllCO3_at, _, _, _, _ = extract_var_stds(atNC, "CO3");

    # pco2
    meanPlumepco2_dic, meanNotPlumepco2_dic, meanAllpco2_dic, _, _, _, _ = extract_var_means(dicNC, "pCO2");
    stdPlumepco2_dic, stdNotPlumepco2_dic, stdAllpco2_dic, _, _, _, _ = extract_var_stds(dicNC, "pCO2");
    minPlumepco2_dic, minNotPlumepco2_dic, minAllpco2_dic, _, _, _, _ = extract_var_mins(dicNC, "pCO2");
    maxPlumepco2_dic, maxNotPlumepco2_dic, maxAllpco2_dic, _, _, _, _ = extract_var_maxs(dicNC, "pCO2");
    meanPlumepco2uncertainty_dic, meanNotPlumepco2uncertainty_dic, meanAllpco2uncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "pCO2_uncertainty");

    meanPlumepco2_at, meanNotPlumepco2_at, meanAllpco2_at, _, _, _, _ = extract_var_means(atNC, "pCO2");
    stdPlumepco2_at, stdNotPlumepco2_at, stdAllpco2_at, _, _, _, _ = extract_var_stds(atNC, "pCO2");
    minPlumepco2_at, minNotPlumepco2_at, minAllpco2_at, _, _, _, _ = extract_var_mins(atNC, "pCO2");
    maxPlumepco2_at, maxNotPlumepco2_at, maxAllpco2_at, _, _, _, _ = extract_var_maxs(atNC, "pCO2");
    meanPlumepco2uncertainty_at, meanNotPlumepco2uncertainty_at, meanAllpco2uncertainty_at, _, _, _, _ = extract_var_means(atNC, "pCO2_uncertainty");
    
    #Aragonite saturation state
    meanPlumeOmegaAragonite_dic, meanNotPlumeOmegaAragonite_dic, meanAllOmegaAragonite_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_aragonite");
    stdPlumeOmegaAragonite_dic, stdNotPlumeOmegaAragonite_dic, stdAllOmegaAragonite_dic, _, _, _, _ = extract_var_stds(dicNC, "saturation_aragonite");
    minPlumeOmegaAragonite_dic, minNotPlumeOmegaAragonite_dic, minAllOmegaAragonite_dic, _, _, _, _ = extract_var_mins(dicNC, "saturation_aragonite");
    maxPlumeOmegaAragonite_dic, maxNotPlumeOmegaAragonite_dic, maxAllOmegaAragonite_dic, _, _, _, _ = extract_var_maxs(dicNC, "saturation_aragonite");
    meanPlumeOmegaAragoniteuncertainty_dic, meanNotPlumeOmegaAragoniteuncertainty_dic, meanAllOmegaAragoniteuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_aragonite_uncertainty");

    meanPlumeOmegaAragonite_at, meanNotPlumeOmegaAragonite_at, meanAllOmegaAragonite_at, _, _, _, _ = extract_var_means(atNC, "saturation_aragonite");
    stdPlumeOmegaAragonite_at, stdNotPlumeOmegaAragonite_at, stdAllOmegaAragonite_at, _, _, _, _ = extract_var_stds(atNC, "saturation_aragonite");
    minPlumeOmegaAragonite_at, minNotPlumeOmegaAragonite_at, minAllOmegaAragonite_at, _, _, _, _ = extract_var_mins(atNC, "saturation_aragonite");
    maxPlumeOmegaAragonite_at, maxNotPlumeOmegaAragonite_at, maxAllOmegaAragonite_at, _, _, _, _ = extract_var_maxs(atNC, "saturation_aragonite");
    meanPlumeOmegaAragoniteuncertainty_at, meanNotPlumeOmegaAragoniteuncertainty_at, meanAllOmegaAragoniteuncertainty_at, _, _, _, _ = extract_var_means(atNC, "saturation_aragonite_uncertainty");

    #Calcite saturation state
    meanPlumeOmegaCalcite_dic, meanNotPlumeOmegaCalcite_dic, meanAllOmegaCalcite_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_calcite");
    stdPlumeOmegaCalcite_dic, stdNotPlumeOmegaCalcite_dic, stdAllOmegaCalcite_dic, _, _, _, _ = extract_var_stds(dicNC, "saturation_calcite");
    minPlumeOmegaCalcite_dic, minNotPlumeOmegaCalcite_dic, minAllOmegaCalcite_dic, _, _, _, _ = extract_var_mins(dicNC, "saturation_calcite");
    maxPlumeOmegaCalcite_dic, maxNotPlumeOmegaCalcite_dic, maxAllOmegaCalcite_dic, _, _, _, _ = extract_var_maxs(dicNC, "saturation_calcite");
    meanPlumeOmegaCalciteuncertainty_dic, meanNotPlumeOmegaCalciteuncertainty_dic, meanAllOmegaCalciteuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_calcite_uncertainty");

    meanPlumeOmegaCalcite_at, meanNotPlumeOmegaCalcite_at, meanAllOmegaCalcite_at, _, _, _, _ = extract_var_means(atNC, "saturation_calcite");
    stdPlumeOmegaCalcite_at, stdNotPlumeOmegaCalcite_at, stdAllOmegaCalcite_at, _, _, _, _ = extract_var_stds(atNC, "saturation_calcite");
    minPlumeOmegaCalcite_at, minNotPlumeOmegaCalcite_at, minAllOmegaCalcite_at, _, _, _, _ = extract_var_mins(atNC, "saturation_calcite");
    maxPlumeOmegaCalcite_at, maxNotPlumeOmegaCalcite_at, maxAllOmegaCalcite_at, _, _, _, _ = extract_var_maxs(atNC, "saturation_calcite");
    meanPlumeOmegaCalciteuncertainty_at, meanNotPlumeOmegaCalciteuncertainty_at, meanAllOmegaCalciteuncertainty_at, _, _, _, _ = extract_var_means(atNC, "saturation_calcite_uncertainty");

        #### Create summary dictionary for key stats
    summary_dict={
                    #### - TA
                    "TA_mean": np.mean(meanAllAT),
                    "TA_mean_JFM": np.mean(meanAllAT[number_list_jfm]),
                    "TA_mean_AMJ": np.mean(meanAllAT[number_list_amj]),
                    "TA_mean_JAS": np.mean(meanAllAT[number_list_jas]),
                    "TA_mean_OND": np.mean(meanAllAT[number_list_ond]),
                    
                    "TA_std": np.std(stdAllAT),
                    "TA_std_JFM": np.std(stdAllAT[number_list_jfm]),
                    "TA_std_AMJ": np.std(stdAllAT[number_list_amj]),
                    "TA_std_JAS": np.std(stdAllAT[number_list_jas]),
                    "TA_std_OND": np.std(stdAllAT[number_list_ond]),
                    
                    "TA_min": np.min(minAllAT),
                    "TA_min_JFM": np.min(minAllAT[number_list_jfm]),
                    "TA_min_AMJ": np.min(minAllAT[number_list_amj]),
                    "TA_min_JAS": np.min(minAllAT[number_list_jas]),
                    "TA_min_OND": np.min(minAllAT[number_list_ond]),
                    
                    "TA_max": np.max(maxAllAT),
                    "TA_max_JFM": np.max(maxAllAT[number_list_jfm]),
                    "TA_max_AMJ": np.max(maxAllAT[number_list_amj]),
                    "TA_max_JAS": np.max(maxAllAT[number_list_jas]),
                    "TA_max_OND": np.max(maxAllAT[number_list_ond]),
                    
                    "TA_uncert": np.mean(meanAllATuncertainty),
                    "TA_uncert_JFM": np.mean(meanAllATuncertainty[number_list_jfm]),
                    "TA_uncert_AMJ": np.mean(meanAllATuncertainty[number_list_amj]),
                    "TA_uncert_JAS": np.mean(meanAllATuncertainty[number_list_jas]),
                    "TA_uncert_OND": np.mean(meanAllATuncertainty[number_list_ond]),
                    
                    "TA_uncert_botup": np.mean(meanAllATuncertainty_botup),
                    "TA_uncert_JFM_botup": np.mean(meanAllATuncertainty_botup[number_list_jfm]),
                    "TA_uncert_AMJ_botup": np.mean(meanAllATuncertainty_botup[number_list_amj]),
                    "TA_uncert_JAS_botup": np.mean(meanAllATuncertainty_botup[number_list_jas]),
                    "TA_uncert_OND_botup": np.mean(meanAllATuncertainty_botup[number_list_ond]),
                    
                    #### - DIC                 
                    "DIC_mean": np.mean(meanAllDIC),
                    "DIC_mean_JFM": np.mean(meanAllDIC[number_list_jfm]),
                    "DIC_mean_AMJ": np.mean(meanAllDIC[number_list_amj]),
                    "DIC_mean_JAS": np.mean(meanAllDIC[number_list_jas]),
                    "DIC_mean_OND": np.mean(meanAllDIC[number_list_ond]),
                    
                    "DIC_std": np.std(stdAllDIC),
                    "DIC_std_JFM": np.std(stdAllDIC[number_list_jfm]),
                    "DIC_std_AMJ": np.std(stdAllDIC[number_list_amj]),
                    "DIC_std_JAS": np.std(stdAllDIC[number_list_jas]),
                    "DIC_std_OND": np.std(stdAllDIC[number_list_ond]),
                    
                    "DIC_min": np.min(minAllDIC),
                    "DIC_min_JFM": np.min(minAllDIC[number_list_jfm]),
                    "DIC_min_AMJ": np.min(minAllDIC[number_list_amj]),
                    "DIC_min_JAS": np.min(minAllDIC[number_list_jas]),
                    "DIC_min_OND": np.min(minAllDIC[number_list_ond]),
                    
                    "DIC_max": np.max(maxAllDIC),
                    "DIC_max_JFM": np.max(maxAllDIC[number_list_jfm]),
                    "DIC_max_AMJ": np.max(maxAllDIC[number_list_amj]),
                    "DIC_max_JAS": np.max(maxAllDIC[number_list_jas]),
                    "DIC_max_OND": np.max(maxAllDIC[number_list_ond]),
                    
                    "DIC_uncert": np.mean(meanAllDICuncertainty),
                    "DIC_uncert_JFM": np.mean(meanAllDICuncertainty[number_list_jfm]),
                    "DIC_uncert_AMJ": np.mean(meanAllDICuncertainty[number_list_amj]),
                    "DIC_uncert_JAS": np.mean(meanAllDICuncertainty[number_list_jas]),
                    "DIC_uncert_OND": np.mean(meanAllDICuncertainty[number_list_ond]),             
                    
                    "DIC_uncert_botup": np.mean(meanAllDICuncertainty_botup),
                    "DIC_uncert_JFM_botup": np.mean(meanAllDICuncertainty_botup[number_list_jfm]),
                    "DIC_uncert_AMJ_botup": np.mean(meanAllDICuncertainty_botup[number_list_amj]),
                    "DIC_uncert_JAS_botup": np.mean(meanAllDICuncertainty_botup[number_list_jas]),
                    "DIC_uncert_OND_botup": np.mean(meanAllDICuncertainty_botup[number_list_ond]),  
                    
                    #### - pH from TA      
                    "pH_mean_SSSSST_from_ta": np.mean(meanAllpH_at),
                    "pH_mean_SSSSST_from_ta_JFM": np.mean(meanAllpH_at[number_list_jfm]),
                    "pH_mean_SSSSST_from_ta_AMJ": np.mean(meanAllpH_at[number_list_amj]),
                    "pH_mean_SSSSST_from_ta_JAS": np.mean(meanAllpH_at[number_list_jas]),
                    "pH_mean_SSSSST_from_ta_OND": np.mean(meanAllpH_at[number_list_ond]),
                    
                    "pH_std_SSSSST_from_ta": np.std(stdAllpH_at),
                    "pH_std_SSSSST_from_ta_JFM": np.std(stdAllpH_at[number_list_jfm]),
                    "pH_std_SSSSST_from_ta_AMJ": np.std(stdAllpH_at[number_list_amj]),
                    "pH_std_SSSSST_from_ta_JAS": np.std(stdAllpH_at[number_list_jas]),
                    "pH_std_SSSSST_from_ta_OND": np.std(stdAllpH_at[number_list_ond]),
                    
                    "pH_min_SSSSST_from_ta": np.min(minAllpH_at),
                    "pH_min_SSSSST_from_ta_JFM": np.min(minAllpH_at[number_list_jfm]),
                    "pH_min_SSSSST_from_ta_AMJ": np.min(minAllpH_at[number_list_amj]),
                    "pH_min_SSSSST_from_ta_JAS": np.min(minAllpH_at[number_list_jas]),
                    "pH_min_SSSSST_from_ta_OND": np.min(minAllpH_at[number_list_ond]),
                    
                    "pH_max_SSSSST_from_ta": np.max(maxAllpH_at),
                    "pH_max_JFM_SSSSST_from_ta": np.max(maxAllpH_at[number_list_jfm]),
                    "pH_max_AMJ_SSSSST_from_ta": np.max(maxAllpH_at[number_list_amj]),
                    "pH_max_JAS_SSSSST_from_ta": np.max(maxAllpH_at[number_list_jas]),
                    "pH_max_OND_SSSSST_from_ta": np.max(maxAllpH_at[number_list_ond]),
                    
                    "pH_uncert_SSSSST_from_ta": np.mean(meanAllpHuncertainty_at),
                    "pH_uncert_SSSSST_from_ta_JFM": np.mean(meanAllpHuncertainty_at[number_list_jfm]),
                    "pH_uncert_SSSSST_from_ta_AMJ": np.mean(meanAllpHuncertainty_at[number_list_amj]),
                    "pH_uncert_SSSSST_from_ta_JAS": np.mean(meanAllpHuncertainty_at[number_list_jas]),
                    "pH_uncert_SSSSST_from_ta_OND": np.mean(meanAllpHuncertainty_at[number_list_ond]),                      
                    
                    #### - pH from DIC
                    "pH_mean_SSSSST_from_dic": np.mean(meanAllpH_dic),
                    "pH_mean_SSSSST_from_dic_JFM": np.mean(meanAllpH_dic[number_list_jfm]),
                    "pH_mean_SSSSST_from_dic_AMJ": np.mean(meanAllpH_dic[number_list_amj]),
                    "pH_mean_SSSSST_from_dic_JAS": np.mean(meanAllpH_dic[number_list_jas]),
                    "pH_mean_SSSSST_from_dic_OND": np.mean(meanAllpH_dic[number_list_ond]),
                    
                    "pH_std_SSSSST_from_dic": np.std(stdAllpH_dic),
                    "pH_std_SSSSST_from_dic_JFM": np.std(stdAllpH_dic[number_list_jfm]),
                    "pH_std_SSSSST_from_dic_AMJ": np.std(stdAllpH_dic[number_list_amj]),
                    "pH_std_SSSSST_from_dic_JAS": np.std(stdAllpH_dic[number_list_jas]),
                    "pH_std_SSSSST_from_dic_OND": np.std(stdAllpH_dic[number_list_ond]),
                    
                    "pH_min_SSSSST_from_dic": np.min(minAllpH_dic),
                    "pH_min_SSSSST_from_dic_JFM": np.min(minAllpH_dic[number_list_jfm]),
                    "pH_min_SSSSST_from_dic_AMJ": np.min(minAllpH_dic[number_list_amj]),
                    "pH_min_SSSSST_from_dic_JAS": np.min(minAllpH_dic[number_list_jas]),
                    "pH_min_SSSSST_from_dic_OND": np.min(minAllpH_dic[number_list_ond]),
                    
                    "pH_max_SSSSST_from_dic": np.max(maxAllpH_dic),
                    "pH_max_JFM_SSSSST_from_dic": np.max(maxAllpH_dic[number_list_jfm]),
                    "pH_max_AMJ_SSSSST_from_dic": np.max(maxAllpH_dic[number_list_amj]),
                    "pH_max_JAS_SSSSST_from_dic": np.max(maxAllpH_dic[number_list_jas]),
                    "pH_max_OND_SSSSST_from_dic": np.max(maxAllpH_dic[number_list_ond]),
                    
                    "pH_uncert_SSSSST_from_dic": np.mean(meanAllpHuncertainty_dic),
                    "pH_uncert_SSSSST_from_dic_JFM": np.mean(meanAllpHuncertainty_dic[number_list_jfm]),
                    "pH_uncert_SSSSST_from_dic_AMJ": np.mean(meanAllpHuncertainty_dic[number_list_amj]),
                    "pH_uncert_SSSSST_from_dic_JAS": np.mean(meanAllpHuncertainty_dic[number_list_jas]),
                    "pH_uncert_SSSSST_from_dic_OND": np.mean(meanAllpHuncertainty_dic[number_list_ond]),  
                    
                    #### - hydrogen_free from TA      
                    "hydrogen_free_mean_SSSSST_from_ta": np.mean(meanAllhydrogen_free_at),
                    "hydrogen_free_mean_SSSSST_from_ta_JFM": np.mean(meanAllhydrogen_free_at[number_list_jfm]),
                    "hydrogen_free_mean_SSSSST_from_ta_AMJ": np.mean(meanAllhydrogen_free_at[number_list_amj]),
                    "hydrogen_free_mean_SSSSST_from_ta_JAS": np.mean(meanAllhydrogen_free_at[number_list_jas]),
                    "hydrogen_free_mean_SSSSST_from_ta_OND": np.mean(meanAllhydrogen_free_at[number_list_ond]),
                    
                    "hydrogen_free_std_SSSSST_from_ta": np.std(stdAllhydrogen_free_at),
                    "hydrogen_free_std_SSSSST_from_ta_JFM": np.std(stdAllhydrogen_free_at[number_list_jfm]),
                    "hydrogen_free_std_SSSSST_from_ta_AMJ": np.std(stdAllhydrogen_free_at[number_list_amj]),
                    "hydrogen_free_std_SSSSST_from_ta_JAS": np.std(stdAllhydrogen_free_at[number_list_jas]),
                    "hydrogen_free_std_SSSSST_from_ta_OND": np.std(stdAllhydrogen_free_at[number_list_ond]),
                    
                    "hydrogen_free_min_SSSSST_from_ta": np.min(minAllhydrogen_free_at),
                    "hydrogen_free_min_SSSSST_from_ta_JFM": np.min(minAllhydrogen_free_at[number_list_jfm]),
                    "hydrogen_free_min_SSSSST_from_ta_AMJ": np.min(minAllhydrogen_free_at[number_list_amj]),
                    "hydrogen_free_min_SSSSST_from_ta_JAS": np.min(minAllhydrogen_free_at[number_list_jas]),
                    "hydrogen_free_min_SSSSST_from_ta_OND": np.min(minAllhydrogen_free_at[number_list_ond]),
                    
                    "hydrogen_free_max_SSSSST_from_ta": np.max(maxAllhydrogen_free_at),
                    "hydrogen_free_max_JFM_SSSSST_from_ta": np.max(maxAllhydrogen_free_at[number_list_jfm]),
                    "hydrogen_free_max_AMJ_SSSSST_from_ta": np.max(maxAllhydrogen_free_at[number_list_amj]),
                    "hydrogen_free_max_JAS_SSSSST_from_ta": np.max(maxAllhydrogen_free_at[number_list_jas]),
                    "hydrogen_free_max_OND_SSSSST_from_ta": np.max(maxAllhydrogen_free_at[number_list_ond]),
                    
                    "hydrogen_free_uncert_SSSSST_from_ta": np.mean(meanAllhydrogen_freeuncertainty_at),
                    "hydrogen_free_uncert_SSSSST_from_ta_JFM": np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jfm]),
                    "hydrogen_free_uncert_SSSSST_from_ta_AMJ": np.mean(meanAllhydrogen_freeuncertainty_at[number_list_amj]),
                    "hydrogen_free_uncert_SSSSST_from_ta_JAS": np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jas]),
                    "hydrogen_free_uncert_SSSSST_from_ta_OND": np.mean(meanAllhydrogen_freeuncertainty_at[number_list_ond]),                      
                    
                    #### - hydrogen_free from DIC
                    "hydrogen_free_mean_SSSSST_from_dic": np.mean(meanAllhydrogen_free_dic),
                    "hydrogen_free_mean_SSSSST_from_dic_JFM": np.mean(meanAllhydrogen_free_dic[number_list_jfm]),
                    "hydrogen_free_mean_SSSSST_from_dic_AMJ": np.mean(meanAllhydrogen_free_dic[number_list_amj]),
                    "hydrogen_free_mean_SSSSST_from_dic_JAS": np.mean(meanAllhydrogen_free_dic[number_list_jas]),
                    "hydrogen_free_mean_SSSSST_from_dic_OND": np.mean(meanAllhydrogen_free_dic[number_list_ond]),
                    
                    "hydrogen_free_std_SSSSST_from_dic": np.std(stdAllhydrogen_free_dic),
                    "hydrogen_free_std_SSSSST_from_dic_JFM": np.std(stdAllhydrogen_free_dic[number_list_jfm]),
                    "hydrogen_free_std_SSSSST_from_dic_AMJ": np.std(stdAllhydrogen_free_dic[number_list_amj]),
                    "hydrogen_free_std_SSSSST_from_dic_JAS": np.std(stdAllhydrogen_free_dic[number_list_jas]),
                    "hydrogen_free_std_SSSSST_from_dic_OND": np.std(stdAllhydrogen_free_dic[number_list_ond]),
                    
                    "hydrogen_free_min_SSSSST_from_dic": np.min(minAllhydrogen_free_dic),
                    "hydrogen_free_min_SSSSST_from_dic_JFM": np.min(minAllhydrogen_free_dic[number_list_jfm]),
                    "hydrogen_free_min_SSSSST_from_dic_AMJ": np.min(minAllhydrogen_free_dic[number_list_amj]),
                    "hydrogen_free_min_SSSSST_from_dic_JAS": np.min(minAllhydrogen_free_dic[number_list_jas]),
                    "hydrogen_free_min_SSSSST_from_dic_OND": np.min(minAllhydrogen_free_dic[number_list_ond]),
                    
                    "hydrogen_free_max_SSSSST_from_dic": np.max(maxAllhydrogen_free_dic),
                    "hydrogen_free_max_JFM_SSSSST_from_dic": np.max(maxAllhydrogen_free_dic[number_list_jfm]),
                    "hydrogen_free_max_AMJ_SSSSST_from_dic": np.max(maxAllhydrogen_free_dic[number_list_amj]),
                    "hydrogen_free_max_JAS_SSSSST_from_dic": np.max(maxAllhydrogen_free_dic[number_list_jas]),
                    "hydrogen_free_max_OND_SSSSST_from_dic": np.max(maxAllhydrogen_free_dic[number_list_ond]),
                    
                    "hydrogen_free_uncert_SSSSST_from_dic": np.mean(meanAllhydrogen_freeuncertainty_dic),
                    "hydrogen_free_uncert_SSSSST_from_dic_JFM": np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jfm]),
                    "hydrogen_free_uncert_SSSSST_from_dic_AMJ": np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_amj]),
                    "hydrogen_free_uncert_SSSSST_from_dic_JAS": np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jas]),
                    "hydrogen_free_uncert_SSSSST_from_dic_OND": np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_ond]),  
                                        
                    #### - pCO2 from TA      
                    "pco2_mean_SSSSST_from_ta": np.mean(meanAllpco2_at),
                    "pco2_mean_SSSSST_from_ta_JFM": np.mean(meanAllpco2_at[number_list_jfm]),
                    "pco2_mean_SSSSST_from_ta_AMJ": np.mean(meanAllpco2_at[number_list_amj]),
                    "pco2_mean_SSSSST_from_ta_JAS": np.mean(meanAllpco2_at[number_list_jas]),
                    "pco2_mean_SSSSST_from_ta_OND": np.mean(meanAllpco2_at[number_list_ond]),
                    
                    "pco2_std_SSSSST_from_ta": np.std(stdAllpco2_at),
                    "pco2_std_SSSSST_from_ta_JFM": np.std(stdAllpco2_at[number_list_jfm]),
                    "pco2_std_SSSSST_from_ta_AMJ": np.std(stdAllpco2_at[number_list_amj]),
                    "pco2_std_SSSSST_from_ta_JAS": np.std(stdAllpco2_at[number_list_jas]),
                    "pco2_std_SSSSST_from_ta_OND": np.std(stdAllpco2_at[number_list_ond]),
                    
                    "pco2_min_SSSSST_from_ta": np.min(minAllpco2_at),
                    "pco2_min_SSSSST_from_ta_JFM": np.min(minAllpco2_at[number_list_jfm]),
                    "pco2_min_SSSSST_from_ta_AMJ": np.min(minAllpco2_at[number_list_amj]),
                    "pco2_min_SSSSST_from_ta_JAS": np.min(minAllpco2_at[number_list_jas]),
                    "pco2_min_SSSSST_from_ta_OND": np.min(minAllpco2_at[number_list_ond]),
                    
                    "pco2_max_SSSSST_from_ta": np.max(maxAllpco2_at),
                    "pco2_max_JFM_SSSSST_from_ta": np.max(maxAllpco2_at[number_list_jfm]),
                    "pco2_max_AMJ_SSSSST_from_ta": np.max(maxAllpco2_at[number_list_amj]),
                    "pco2_max_JAS_SSSSST_from_ta": np.max(maxAllpco2_at[number_list_jas]),
                    "pco2_max_OND_SSSSST_from_ta": np.max(maxAllpco2_at[number_list_ond]),
                    
                    "pco2_uncert_SSSSST_from_ta": np.mean(meanAllpco2uncertainty_at),
                    "pco2_uncert_SSSSST_from_ta_JFM": np.mean(meanAllpco2uncertainty_at[number_list_jfm]),
                    "pco2_uncert_SSSSST_from_ta_AMJ": np.mean(meanAllpco2uncertainty_at[number_list_amj]),
                    "pco2_uncert_SSSSST_from_ta_JAS": np.mean(meanAllpco2uncertainty_at[number_list_jas]),
                    "pco2_uncert_SSSSST_from_ta_OND": np.mean(meanAllpco2uncertainty_at[number_list_ond]),                      
                    
                    #### - pCO2 from DIC
                    "pco2_mean_SSSSST_from_dic": np.mean(meanAllpco2_dic),
                    "pco2_mean_SSSSST_from_dic_JFM": np.mean(meanAllpco2_dic[number_list_jfm]),
                    "pco2_mean_SSSSST_from_dic_AMJ": np.mean(meanAllpco2_dic[number_list_amj]),
                    "pco2_mean_SSSSST_from_dic_JAS": np.mean(meanAllpco2_dic[number_list_jas]),
                    "pco2_mean_SSSSST_from_dic_OND": np.mean(meanAllpco2_dic[number_list_ond]),
                    
                    "pco2_std_SSSSST_from_dic": np.std(stdAllpco2_dic),
                    "pco2_std_SSSSST_from_dic_JFM": np.std(stdAllpco2_dic[number_list_jfm]),
                    "pco2_std_SSSSST_from_dic_AMJ": np.std(stdAllpco2_dic[number_list_amj]),
                    "pco2_std_SSSSST_from_dic_JAS": np.std(stdAllpco2_dic[number_list_jas]),
                    "pco2_std_SSSSST_from_dic_OND": np.std(stdAllpco2_dic[number_list_ond]),
                    
                    "pco2_min_SSSSST_from_dic": np.min(minAllpco2_dic),
                    "pco2_min_SSSSST_from_dic_JFM": np.min(minAllpco2_dic[number_list_jfm]),
                    "pco2_min_SSSSST_from_dic_AMJ": np.min(minAllpco2_dic[number_list_amj]),
                    "pco2_min_SSSSST_from_dic_JAS": np.min(minAllpco2_dic[number_list_jas]),
                    "pco2_min_SSSSST_from_dic_OND": np.min(minAllpco2_dic[number_list_ond]),
                    
                    "pco2_max_SSSSST_from_dic": np.max(maxAllpco2_dic),
                    "pco2_max_JFM_SSSSST_from_dic": np.max(maxAllpco2_dic[number_list_jfm]),
                    "pco2_max_AMJ_SSSSST_from_dic": np.max(maxAllpco2_dic[number_list_amj]),
                    "pco2_max_JAS_SSSSST_from_dic": np.max(maxAllpco2_dic[number_list_jas]),
                    "pco2_max_OND_SSSSST_from_dic": np.max(maxAllpco2_dic[number_list_ond]),
                    
                    "pco2_uncert_SSSSST_from_dic": np.mean(meanAllpco2uncertainty_dic),
                    "pco2_uncert_SSSSST_from_dic_JFM": np.mean(meanAllpco2uncertainty_dic[number_list_jfm]),
                    "pco2_uncert_SSSSST_from_dic_AMJ": np.mean(meanAllpco2uncertainty_dic[number_list_amj]),
                    "pco2_uncert_SSSSST_from_dic_JAS": np.mean(meanAllpco2uncertainty_dic[number_list_jas]),
                    "pco2_uncert_SSSSST_from_dic_OND": np.mean(meanAllpco2uncertainty_dic[number_list_ond]),  
                    
                      
                    #### - Omega Aragonite from TA      
                    "OmegaAragonite_mean_SSSSST_from_ta": np.mean(meanAllOmegaAragonite_at),
                    "OmegaAragonite_mean_SSSSST_from_ta_JFM": np.mean(meanAllOmegaAragonite_at[number_list_jfm]),
                    "OmegaAragonite_mean_SSSSST_from_ta_AMJ": np.mean(meanAllOmegaAragonite_at[number_list_amj]),
                    "OmegaAragonite_mean_SSSSST_from_ta_JAS": np.mean(meanAllOmegaAragonite_at[number_list_jas]),
                    "OmegaAragonite_mean_SSSSST_from_ta_OND": np.mean(meanAllOmegaAragonite_at[number_list_ond]),
                    
                    "OmegaAragonite_std_SSSSST_from_ta": np.std(stdAllOmegaAragonite_at),
                    "OmegaAragonite_std_SSSSST_from_ta_JFM": np.std(stdAllOmegaAragonite_at[number_list_jfm]),
                    "OmegaAragonite_std_SSSSST_from_ta_AMJ": np.std(stdAllOmegaAragonite_at[number_list_amj]),
                    "OmegaAragonite_std_SSSSST_from_ta_JAS": np.std(stdAllOmegaAragonite_at[number_list_jas]),
                    "OmegaAragonite_std_SSSSST_from_ta_OND": np.std(stdAllOmegaAragonite_at[number_list_ond]),
                    
                    "OmegaAragonite_min_SSSSST_from_ta": np.min(minAllOmegaAragonite_at),
                    "OmegaAragonite_min_SSSSST_from_ta_JFM": np.min(minAllOmegaAragonite_at[number_list_jfm]),
                    "OmegaAragonite_min_SSSSST_from_ta_AMJ": np.min(minAllOmegaAragonite_at[number_list_amj]),
                    "OmegaAragonite_min_SSSSST_from_ta_JAS": np.min(minAllOmegaAragonite_at[number_list_jas]),
                    "OmegaAragonite_min_SSSSST_from_ta_OND": np.min(minAllOmegaAragonite_at[number_list_ond]),
                    
                    "OmegaAragonite_max_SSSSST_from_ta": np.max(maxAllOmegaAragonite_at),
                    "OmegaAragonite_max_JFM_SSSSST_from_ta": np.max(maxAllOmegaAragonite_at[number_list_jfm]),
                    "OmegaAragonite_max_AMJ_SSSSST_from_ta": np.max(maxAllOmegaAragonite_at[number_list_amj]),
                    "OmegaAragonite_max_JAS_SSSSST_from_ta": np.max(maxAllOmegaAragonite_at[number_list_jas]),
                    "OmegaAragonite_max_OND_SSSSST_from_ta": np.max(maxAllOmegaAragonite_at[number_list_ond]),
                    
                    "OmegaAragonite_uncert_SSSSST_from_ta": np.mean(meanAllOmegaAragoniteuncertainty_at),
                    "OmegaAragonite_uncert_SSSSST_from_ta_JFM": np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_jfm]),
                    "OmegaAragonite_uncert_SSSSST_from_ta_AMJ": np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_amj]),
                    "OmegaAragonite_uncert_SSSSST_from_ta_JAS": np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_jas]),
                    "OmegaAragonite_uncert_SSSSST_from_ta_OND": np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_ond]),                      
                    
                    #### - Omega Aragonite from DIC
                    "OmegaAragonite_mean_SSSSST_from_dic": np.mean(meanAllOmegaAragonite_dic),
                    "OmegaAragonite_mean_SSSSST_from_dic_JFM": np.mean(meanAllOmegaAragonite_dic[number_list_jfm]),
                    "OmegaAragonite_mean_SSSSST_from_dic_AMJ": np.mean(meanAllOmegaAragonite_dic[number_list_amj]),
                    "OmegaAragonite_mean_SSSSST_from_dic_JAS": np.mean(meanAllOmegaAragonite_dic[number_list_jas]),
                    "OmegaAragonite_mean_SSSSST_from_dic_OND": np.mean(meanAllOmegaAragonite_dic[number_list_ond]),
                    
                    "OmegaAragonite_std_SSSSST_from_dic": np.std(stdAllOmegaAragonite_dic),
                    "OmegaAragonite_std_SSSSST_from_dic_JFM": np.std(stdAllOmegaAragonite_dic[number_list_jfm]),
                    "OmegaAragonite_std_SSSSST_from_dic_AMJ": np.std(stdAllOmegaAragonite_dic[number_list_amj]),
                    "OmegaAragonite_std_SSSSST_from_dic_JAS": np.std(stdAllOmegaAragonite_dic[number_list_jas]),
                    "OmegaAragonite_std_SSSSST_from_dic_OND": np.std(stdAllOmegaAragonite_dic[number_list_ond]),
                    
                    "OmegaAragonite_min_SSSSST_from_dic": np.min(minAllOmegaAragonite_dic),
                    "OmegaAragonite_min_SSSSST_from_dic_JFM": np.min(minAllOmegaAragonite_dic[number_list_jfm]),
                    "OmegaAragonite_min_SSSSST_from_dic_AMJ": np.min(minAllOmegaAragonite_dic[number_list_amj]),
                    "OmegaAragonite_min_SSSSST_from_dic_JAS": np.min(minAllOmegaAragonite_dic[number_list_jas]),
                    "OmegaAragonite_min_SSSSST_from_dic_OND": np.min(minAllOmegaAragonite_dic[number_list_ond]),
                    
                    "OmegaAragonite_max_SSSSST_from_dic": np.max(maxAllOmegaAragonite_dic),
                    "OmegaAragonite_max_JFM_SSSSST_from_dic": np.max(maxAllOmegaAragonite_dic[number_list_jfm]),
                    "OmegaAragonite_max_AMJ_SSSSST_from_dic": np.max(maxAllOmegaAragonite_dic[number_list_amj]),
                    "OmegaAragonite_max_JAS_SSSSST_from_dic": np.max(maxAllOmegaAragonite_dic[number_list_jas]),
                    "OmegaAragonite_max_OND_SSSSST_from_dic": np.max(maxAllOmegaAragonite_dic[number_list_ond]),
                    
                    "OmegaAragonite_uncert_SSSSST_from_dic": np.mean(meanAllOmegaAragoniteuncertainty_dic),
                    "OmegaAragonite_uncert_SSSSST_from_dic_JFM": np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_jfm]),
                    "OmegaAragonite_uncert_SSSSST_from_dic_AMJ": np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_amj]),
                    "OmegaAragonite_uncert_SSSSST_from_dic_JAS": np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_jas]),
                    "OmegaAragonite_uncert_SSSSST_from_dic_OND": np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_ond]),  
                    
                
                    #### - Omega Calcilte from TA      
                    "OmegaCalcite_mean_SSSSST_from_ta": np.mean(meanAllOmegaCalcite_at),
                    "OmegaCalcite_mean_SSSSST_from_ta_JFM": np.mean(meanAllOmegaCalcite_at[number_list_jfm]),
                    "OmegaCalcite_mean_SSSSST_from_ta_AMJ": np.mean(meanAllOmegaCalcite_at[number_list_amj]),
                    "OmegaCalcite_mean_SSSSST_from_ta_JAS": np.mean(meanAllOmegaCalcite_at[number_list_jas]),
                    "OmegaCalcite_mean_SSSSST_from_ta_OND": np.mean(meanAllOmegaCalcite_at[number_list_ond]),
                    
                    "OmegaCalcite_std_SSSSST_from_ta": np.std(stdAllOmegaCalcite_at),
                    "OmegaCalcite_std_SSSSST_from_ta_JFM": np.std(stdAllOmegaCalcite_at[number_list_jfm]),
                    "OmegaCalcite_std_SSSSST_from_ta_AMJ": np.std(stdAllOmegaCalcite_at[number_list_amj]),
                    "OmegaCalcite_std_SSSSST_from_ta_JAS": np.std(stdAllOmegaCalcite_at[number_list_jas]),
                    "OmegaCalcite_std_SSSSST_from_ta_OND": np.std(stdAllOmegaCalcite_at[number_list_ond]),
                    
                    "OmegaCalcite_min_SSSSST_from_ta": np.min(minAllOmegaCalcite_at),
                    "OmegaCalcite_min_SSSSST_from_ta_JFM": np.min(minAllOmegaCalcite_at[number_list_jfm]),
                    "OmegaCalcite_min_SSSSST_from_ta_AMJ": np.min(minAllOmegaCalcite_at[number_list_amj]),
                    "OmegaCalcite_min_SSSSST_from_ta_JAS": np.min(minAllOmegaCalcite_at[number_list_jas]),
                    "OmegaCalcite_min_SSSSST_from_ta_OND": np.min(minAllOmegaCalcite_at[number_list_ond]),
                    
                    "OmegaCalcite_max_SSSSST_from_ta": np.max(maxAllOmegaCalcite_at),
                    "OmegaCalcite_max_JFM_SSSSST_from_ta": np.max(maxAllOmegaCalcite_at[number_list_jfm]),
                    "OmegaCalcite_max_AMJ_SSSSST_from_ta": np.max(maxAllOmegaCalcite_at[number_list_amj]),
                    "OmegaCalcite_max_JAS_SSSSST_from_ta": np.max(maxAllOmegaCalcite_at[number_list_jas]),
                    "OmegaCalcite_max_OND_SSSSST_from_ta": np.max(maxAllOmegaCalcite_at[number_list_ond]),
                    
                    "OmegaCalcite_uncert_SSSSST_from_ta": np.mean(meanAllOmegaCalciteuncertainty_at),
                    "OmegaCalcite_uncert_SSSSST_from_ta_JFM": np.mean(meanAllOmegaCalciteuncertainty_at[number_list_jfm]),
                    "OmegaCalcite_uncert_SSSSST_from_ta_AMJ": np.mean(meanAllOmegaCalciteuncertainty_at[number_list_amj]),
                    "OmegaCalcite_uncert_SSSSST_from_ta_JAS": np.mean(meanAllOmegaCalciteuncertainty_at[number_list_jas]),
                    "OmegaCalcite_uncert_SSSSST_from_ta_OND": np.mean(meanAllOmegaCalciteuncertainty_at[number_list_ond]),                      
                    
                    #### - Omega Calcite from DIC
                    "OmegaCalcite_mean_SSSSST_from_dic": np.mean(meanAllOmegaCalcite_dic),
                    "OmegaCalcite_mean_SSSSST_from_dic_JFM": np.mean(meanAllOmegaCalcite_dic[number_list_jfm]),
                    "OmegaCalcite_mean_SSSSST_from_dic_AMJ": np.mean(meanAllOmegaCalcite_dic[number_list_amj]),
                    "OmegaCalcite_mean_SSSSST_from_dic_JAS": np.mean(meanAllOmegaCalcite_dic[number_list_jas]),
                    "OmegaCalcite_mean_SSSSST_from_dic_OND": np.mean(meanAllOmegaCalcite_dic[number_list_ond]),
                    
                    "OmegaCalcite_std_SSSSST_from_dic": np.std(stdAllOmegaCalcite_dic),
                    "OmegaCalcite_std_SSSSST_from_dic_JFM": np.std(stdAllOmegaCalcite_dic[number_list_jfm]),
                    "OmegaCalcite_std_SSSSST_from_dic_AMJ": np.std(stdAllOmegaCalcite_dic[number_list_amj]),
                    "OmegaCalcite_std_SSSSST_from_dic_JAS": np.std(stdAllOmegaCalcite_dic[number_list_jas]),
                    "OmegaCalcite_std_SSSSST_from_dic_OND": np.std(stdAllOmegaCalcite_dic[number_list_ond]),
                    
                    "OmegaCalcite_min_SSSSST_from_dic": np.min(minAllOmegaCalcite_dic),
                    "OmegaCalcite_min_SSSSST_from_dic_JFM": np.min(minAllOmegaCalcite_dic[number_list_jfm]),
                    "OmegaCalcite_min_SSSSST_from_dic_AMJ": np.min(minAllOmegaCalcite_dic[number_list_amj]),
                    "OmegaCalcite_min_SSSSST_from_dic_JAS": np.min(minAllOmegaCalcite_dic[number_list_jas]),
                    "OmegaCalcite_min_SSSSST_from_dic_OND": np.min(minAllOmegaCalcite_dic[number_list_ond]),
                    
                    "OmegaCalcite_max_SSSSST_from_dic": np.max(maxAllOmegaCalcite_dic),
                    "OmegaCalcite_max_JFM_SSSSST_from_dic": np.max(maxAllOmegaCalcite_dic[number_list_jfm]),
                    "OmegaCalcite_max_AMJ_SSSSST_from_dic": np.max(maxAllOmegaCalcite_dic[number_list_amj]),
                    "OmegaCalcite_max_JAS_SSSSST_from_dic": np.max(maxAllOmegaCalcite_dic[number_list_jas]),
                    "OmegaCalcite_max_OND_SSSSST_from_dic": np.max(maxAllOmegaCalcite_dic[number_list_ond]),
                    
                    "OmegaCalcite_uncert_SSSSST_from_dic": np.mean(meanAllOmegaCalciteuncertainty_dic),
                    "OmegaCalcite_uncert_SSSSST_from_dic_JFM": np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_jfm]),
                    "OmegaCalcite_uncert_SSSSST_from_dic_AMJ": np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_amj]),
                    "OmegaCalcite_uncert_SSSSST_from_dic_JAS": np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_jas]),
                    "OmegaCalcite_uncert_SSSSST_from_dic_OND": np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_ond]),  
                    
                    };
    
    
    # alternatively save as .json files
    # json = json.dumps(summary_dict, indent=4)    # create json object from dictionary
    # f = open("uncert_summary_{0}.json".format(region),"w")# open file for writing, "w" 
    # f.write(json)# write json object to file
    # f.close()# close file

        #### Save summary dictionary for key stats    
    with open('output\\gridded_predictions_min_year_range\\uncert_summary_{0}.csv'.format(region), 'w') as f:
        for key in summary_dict.keys():
            f.write("%s,%s\n"%(key,summary_dict[key]))



        #### Formatted table of summary stats
    field_names = ['Timeperiod', 'mean', 'std','min','max','uncertainty','uncertainty_botup']
    dp_to_round=2;
    TA = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllAT),dp_to_round), 'std': np.round(np.std(stdAllAT),dp_to_round), 'min': np.round(np.min(minAllAT),dp_to_round), 'max':np.round( np.max(maxAllAT),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllATuncertainty_botup),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllAT[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_jfm]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllATuncertainty_botup[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllAT[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_amj]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllATuncertainty_botup[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllAT[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_jas]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllATuncertainty_botup[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllAT[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_ond]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllATuncertainty_botup[number_list_ond]),dp_to_round)},
    ]
    
    DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllDIC),dp_to_round), 'std': np.round(np.std(stdAllDIC),dp_to_round), 'min': np.round(np.min(minAllDIC),dp_to_round), 'max': np.round(np.max(maxAllDIC),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllDICuncertainty_botup),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllDIC[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_jfm]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllDICuncertainty_botup[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllDIC[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_amj]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllDICuncertainty_botup[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllDIC[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_jas]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllDICuncertainty_botup[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllDIC[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_ond]),dp_to_round), 'uncertainty_botup': np.round(np.mean(meanAllDICuncertainty_botup[number_list_ond]),dp_to_round)},
    ]
    
    pH_AT = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllpH_at),dp_to_round), 'std': np.round(np.std(stdAllpH_at),dp_to_round), 'min': np.round(np.min(minAllpH_at),dp_to_round), 'max': np.round(np.max(maxAllpH_at),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_at),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllpH_at[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllpH_at[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllpH_at[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllpH_at[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_at[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllpH_at[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllpH_at[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllpH_at[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllpH_at[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_at[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllpH_at[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllpH_at[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllpH_at[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllpH_at[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_at[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllpH_at[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllpH_at[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllpH_at[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllpH_at[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_at[number_list_ond]),dp_to_round)},
    ]
    
    pH_DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllpH_dic),dp_to_round), 'std': np.round(np.std(stdAllpH_dic),dp_to_round), 'min': np.round(np.min(minAllpH_dic),dp_to_round), 'max': np.round(np.max(maxAllpH_dic),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_dic),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllpH_dic[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllpH_dic[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllpH_dic[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllpH_dic[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_dic[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllpH_dic[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllpH_dic[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllpH_dic[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllpH_dic[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_dic[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllpH_dic[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllpH_dic[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllpH_dic[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllpH_dic[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_dic[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllpH_dic[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllpH_dic[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllpH_dic[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllpH_dic[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpHuncertainty_dic[number_list_ond]),dp_to_round)},
    ]
    
    hydrogen_free_AT = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllhydrogen_free_at),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_at),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_at),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_at),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_at),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllhydrogen_free_at[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_at[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_at[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_at[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllhydrogen_free_at[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_at[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_at[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_at[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_at[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllhydrogen_free_at[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_at[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_at[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_at[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllhydrogen_free_at[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_at[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_at[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_at[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_at[number_list_ond]),dp_to_round)},
    ]
    
    hydrogen_free_DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllhydrogen_free_dic),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_dic),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_dic),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_dic),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_dic),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllhydrogen_free_dic[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_dic[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_dic[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_dic[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllhydrogen_free_dic[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_dic[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_dic[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_dic[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllhydrogen_free_dic[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_dic[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_dic[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_dic[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllhydrogen_free_dic[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllhydrogen_free_dic[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllhydrogen_free_dic[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllhydrogen_free_dic[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_ond]),dp_to_round)},
    ]
    
    
    ph_from_hydrogen_free_AT = [
    {'Timeperiod':'Annual' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_at))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_at)+np.std(stdAllhydrogen_free_at))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_at))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_at))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_at)+np.mean(meanAllhydrogen_freeuncertainty_at))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at))*1e-6)},
    {'Timeperiod':'JFM', 'mean': -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jfm]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jfm])+np.std(stdAllhydrogen_free_at[number_list_jfm]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_jfm]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_at[number_list_jfm]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_at[number_list_jfm]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jfm])+np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jfm]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_jfm]))*1e-6)},
    {'Timeperiod':'AMJ' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_at[number_list_amj]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_amj])+np.std(stdAllhydrogen_free_at[number_list_amj]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_amj]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_at[number_list_amj]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_at[number_list_amj]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_amj])+np.mean(meanAllhydrogen_freeuncertainty_at[number_list_amj]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_amj]))*1e-6)},
    {'Timeperiod':'JAS', 'mean': -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jas]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jas])+np.std(stdAllhydrogen_free_at[number_list_jas]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_jas]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_at[number_list_jas]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_at[number_list_jas]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_jas])+np.mean(meanAllhydrogen_freeuncertainty_at[number_list_jas]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_jas]))*1e-6)},
    {'Timeperiod':'OND' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_at[number_list_ond]))*1e-6), 'std':    - np.log10((np.mean(meanAllhydrogen_free_at[number_list_ond])+np.std(stdAllhydrogen_free_at[number_list_ond]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_ond]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_at[number_list_ond]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_at[number_list_ond]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_at[number_list_ond])+np.mean(meanAllhydrogen_freeuncertainty_at[number_list_ond]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_at[number_list_ond]))*1e-6)},
    ]
    
    ph_from_hydrogen_free_DIC = [
    {'Timeperiod':'Annual' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_dic))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_dic)+np.std(stdAllhydrogen_free_dic))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_dic))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_dic))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_dic)+np.mean(meanAllhydrogen_freeuncertainty_dic))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic))*1e-6)},
    {'Timeperiod':'JFM', 'mean': -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jfm]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jfm])+np.std(stdAllhydrogen_free_dic[number_list_jfm]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jfm]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_dic[number_list_jfm]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_dic[number_list_jfm]))*1e-6), 'uncertainty':    - np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jfm])+np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jfm]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jfm]))*1e-6)},
    {'Timeperiod':'AMJ' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_amj]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_amj])+np.std(stdAllhydrogen_free_dic[number_list_amj]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_amj]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_dic[number_list_amj]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_dic[number_list_amj]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_amj])+np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_amj]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_amj]))*1e-6)},
    {'Timeperiod':'JAS', 'mean': -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jas]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jas])+np.std(stdAllhydrogen_free_dic[number_list_jas]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jas]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_dic[number_list_jas]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_dic[number_list_jas]))*1e-6), 'uncertainty':    -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jas])+np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_jas]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_jas]))*1e-6)},
    {'Timeperiod':'OND' , 'mean': -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_ond]))*1e-6), 'std':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_ond])+np.std(stdAllhydrogen_free_dic[number_list_ond]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_ond]))*1e-6), 'max': -np.log10((np.min(minAllhydrogen_free_dic[number_list_ond]))*1e-6), 'min': -np.log10((np.max(maxAllhydrogen_free_dic[number_list_ond]))*1e-6), 'uncertainty':     -np.log10((np.mean(meanAllhydrogen_free_dic[number_list_ond])+np.mean(meanAllhydrogen_freeuncertainty_dic[number_list_ond]))*1e-6)+np.log10((np.mean(meanAllhydrogen_free_dic[number_list_ond]))*1e-6)},
    ]
    
    pco2_AT = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllpco2_at),dp_to_round), 'std': np.round(np.std(stdAllpco2_at),dp_to_round), 'min': np.round(np.min(minAllpco2_at),dp_to_round), 'max': np.round(np.max(maxAllpco2_at),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_at),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllpco2_at[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllpco2_at[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllpco2_at[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllpco2_at[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_at[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllpco2_at[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllpco2_at[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllpco2_at[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllpco2_at[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_at[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllpco2_at[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllpco2_at[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllpco2_at[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllpco2_at[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_at[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllpco2_at[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllpco2_at[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllpco2_at[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllpco2_at[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_at[number_list_ond]),dp_to_round)},
    ]
    
    pco2_DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllpco2_dic),dp_to_round), 'std': np.round(np.std(stdAllpco2_dic),dp_to_round), 'min': np.round(np.min(minAllpco2_dic),dp_to_round), 'max': np.round(np.max(maxAllpco2_dic),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_dic),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllpco2_dic[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllpco2_dic[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllpco2_dic[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllpco2_dic[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_dic[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllpco2_dic[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllpco2_dic[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllpco2_dic[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllpco2_dic[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_dic[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllpco2_dic[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllpco2_dic[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllpco2_dic[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllpco2_dic[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_dic[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllpco2_dic[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllpco2_dic[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllpco2_dic[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllpco2_dic[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllpco2uncertainty_dic[number_list_ond]),dp_to_round)},
    ]
    
    OmegaAragonite_AT = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllOmegaAragonite_at),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_at),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_at),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_at),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_at),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllOmegaAragonite_at[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_at[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_at[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_at[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllOmegaAragonite_at[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_at[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_at[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_at[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllOmegaAragonite_at[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_at[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_at[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_at[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllOmegaAragonite_at[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_at[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_at[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_at[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_at[number_list_ond]),dp_to_round)},
    ]
    
    OmegaAragonite_DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllOmegaAragonite_dic),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_dic),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_dic),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_dic),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_dic),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllOmegaAragonite_dic[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_dic[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_dic[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_dic[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllOmegaAragonite_dic[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_dic[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_dic[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_dic[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllOmegaAragonite_dic[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_dic[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_dic[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_dic[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllOmegaAragonite_dic[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllOmegaAragonite_dic[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllOmegaAragonite_dic[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllOmegaAragonite_dic[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaAragoniteuncertainty_dic[number_list_ond]),dp_to_round)},
    ]
    
    OmegaCalcite_AT = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllOmegaCalcite_at),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_at),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_at),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_at),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_at),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllOmegaCalcite_at[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_at[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_at[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_at[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_at[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllOmegaCalcite_at[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_at[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_at[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_at[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_at[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllOmegaCalcite_at[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_at[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_at[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_at[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_at[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllOmegaCalcite_at[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_at[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_at[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_at[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_at[number_list_ond]),dp_to_round)},
    ]
    
    OmegaCalcite_DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllOmegaCalcite_dic),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_dic),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_dic),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_dic),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_dic),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllOmegaCalcite_dic[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_dic[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_dic[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_dic[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllOmegaCalcite_dic[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_dic[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_dic[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_dic[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllOmegaCalcite_dic[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_dic[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_dic[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_dic[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllOmegaCalcite_dic[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllOmegaCalcite_dic[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllOmegaCalcite_dic[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllOmegaCalcite_dic[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllOmegaCalciteuncertainty_dic[number_list_ond]),dp_to_round)},
    ]
    
        #### Save formatted table of sumnmary stats
    with open('output\\gridded_predictions_min_year_range\\Gridded_data_summary_stats_{0}.csv'.format(region), 'w',newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames = field_names)
        
        csvfile.write('TA ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(TA)
        csvfile.write('\n') #Give your csv text here.
        
        csvfile.write('DIC ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(DIC)
        csvfile.write('\n') #Give your csv text here.
        
        csvfile.write('pH using SSSSST from TA algorithm (note: averages/std etc calculated in logs!!!) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(pH_AT)
        csvfile.write('\n') #Give your csv text here
        
        csvfile.write('pH using SSSSST from DIC algorithm (note: averages/std/uncertainty etc calculated in logs!!!) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(pH_DIC)
        csvfile.write('\n') #Give your csv text here.
    
        csvfile.write('Hydrogen free aka H+ using SSSSST from TA algorithm ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(hydrogen_free_AT)
        csvfile.write('\n') #Give your csv text here
        
        csvfile.write('Hydrogen free aka H+ using SSSSST from DIC algorithm ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(hydrogen_free_DIC)
        csvfile.write('\n') #Give your csv text here.
        
        csvfile.write('ph calc from Hydrogen free aka H+ using SSSSST from DIC algorithm ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        csvfile.write('umol converted to mol and logged to give mean, min and max \n') #Give your csv text here.
        csvfile.write('uncertainty and std are meaningless without refernce, so ro calc unc and std did the following log10(H+ mean and H+std) - log10 (H+)  \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(ph_from_hydrogen_free_AT)
        csvfile.write('\n') #Give your csv text here.
    
        csvfile.write('pH calc from Hydrogen free aka H+ using SSSSST from DIC algorithm ($\mu mol \ kg^{-1}$) \n') #Give your csv text here.
        csvfile.write('umol converted to mol and logged to give mean, min and max \n') #Give your csv text here.
        csvfile.write('uncertainty and std are meaningless without refernce, so ro calc unc and std did the following log10(H+ mean and H+std) - log10 (H+)  \n') #Give your csv text here.writer.writeheader()
        writer.writerows(ph_from_hydrogen_free_DIC)
        csvfile.write('\n') #Give your csv text here.
    
        csvfile.write('pco2 using SSSSST from TA algorithm (ppm) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(pco2_AT)
        csvfile.write('\n') #Give your csv text here
        
        csvfile.write('pco2 using SSSSST from DIC algorithm  (ppm) \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(pco2_DIC)
        csvfile.write('\n') #Give your csv text here.
    
        csvfile.write('OmegaAragonite using SSSSST from TA algorithm \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(OmegaAragonite_AT)
        csvfile.write('\n') #Give your csv text here
        
        csvfile.write('Omegaaragonite using SSSSST from DIC algorithm \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(OmegaAragonite_DIC)
        csvfile.write('\n') #Give your csv text here.
    
        csvfile.write('OmegaCalcite using SSSSST from TA algorithm \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(OmegaCalcite_AT)
        csvfile.write('\n') #Give your csv text here
        
        csvfile.write('OmegaCalcite using SSSSST from DIC algorithm \n') #Give your csv text here.
        writer.writeheader()
        writer.writerows(OmegaCalcite_DIC)
        csvfile.write('\n') #Give your csv text here.

    #### Videos
        #### Create png images of the region every month
    if create_animations=="True": #this can be slow so make it an on/off switch 
        for video_var in video_vars:
            outputPath = "plots/";
            maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
            maskNC = Dataset(maskPath, 'r');
        
            #### For DIC
            mask = maskNC.variables[region][:];
            xmin = min(np.where(np.flipud(mask)==1)[1]);
            xmax = max(np.where(np.flipud(mask)==1)[1]);
            ymin = min(np.where(np.flipud(mask)==1)[0]);
            ymax = max(np.where(np.flipud(mask)==1)[0]);
            pad = 4; #size of xlim/ylim extents
            
            figsize = (8,4);
            ticksize = 10;
            dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
            varALL = dicNC.variables["{0}".format(video_var)][:];
            hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
            
            # these cover the full range in time
            # maxvar = np.nanmax(varALL);
            # minvar = np.nanmin(varALL);
            
            # arbitary limits
            if region == "oceansoda_amazon_plume":
                if video_var == "DIC":
                    maxvar = 2000;
                    minvar = 1400; 
                elif video_var == "SSS":
                    maxvar = 37;
                    minvar = 32; 
                elif video_var == "SST":
                    maxvar = 30;
                    minvar = 22; 
                elif video_var == "pH_free": 
                    maxvar = 8.5;
                    minvar = 7.5; 
                elif video_var == "pCO2":
                    maxvar = 3000;
                    minvar = 200; 
                elif video_var == "saturation_aragonite":
                    maxvar = 10;
                    minvar = 0; 
                elif video_var == "saturation_calcite":                 
                    maxvar = 10;
                    minvar = 0; 
                else:
                    pass
                    
            elif region == "oceansoda_congo":
                if video_var == "DIC":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "SSS":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "SST":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "pH_free": 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "pCO2":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "saturation_aragonite":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif video_var == "saturation_calcite":                 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                else:
                    pass
            elif region == "oceansoda_mediterranean":
                if video_var == "DIC":
                    maxvar = 2300;
                    minvar = 1950; 
                elif video_var == "SSS":
                    maxvar = 37;
                    minvar = 32; 
                elif video_var == "SST":
                    maxvar = 25;
                    minvar = 10; 
                elif video_var == "pH_free": 
                    maxvar = 8.2;
                    minvar = 8; 
                elif video_var == "pCO2":
                    maxvar = 600;
                    minvar = 200; 
                elif video_var == "saturation_aragonite":
                    maxvar = 10;
                    minvar = 0; 
                elif video_var == "saturation_calcite":                 
                    maxvar = 10;
                    minvar = 0; 
                else:
                    pass
            else:
                pass
            
    
            lats = dicNC.variables["lat"][:];
            lons = dicNC.variables["lon"][:];
            
            #check that folder exists to save images, otherwise create it
            from pathlib import Path
            Path("os_plots\{0}_{1}_animation_images".format(region,video_var)).mkdir(parents=True, exist_ok=True)
            
            iframe = 0;
            for t in range(0, varALL.shape[0]):
                if hasData[t] == False:
                    continue;
                
                date = dates[t];
                plotvar = varALL[t,:,:];
                if video_var == "SST":
                    plotvar=plotvar-273.15;
                plotvar.mask[mask!=1] = True;
                
                #print("Plotting frame for {0} {1}".format(date.year, format(date.month, "02d")));
                
                plt.figure(figsize=figsize);
                ax = plt.axes(projection=ccrs.PlateCarree());
                
                if video_var == "pCO2":
                    cmap_reversed = plt.cm.get_cmap('inferno_r')
    
                    #this plots the data as 1 x 1 grids 
                    data_proj = ccrs.PlateCarree()
                    contPlot = plt.pcolormesh(lons, lats, plotvar, transform=data_proj , vmin=minvar, vmax=maxvar, cmap=cmap_reversed);
                    
                    #this contours the data and smooths it
                    #contPlot = plt.contourf(lons, lats, plotvar,100, transform=ccrs.PlateCarree() , vmin=minvar, vmax=maxvar, cmap=cmap_reversed);
                else:
                    
                    #this plots the data as 1 x 1 grids 
                    data_proj = ccrs.PlateCarree()
                    contPlot = plt.pcolormesh(lons, lats, plotvar, transform=data_proj , vmin=minvar, vmax=maxvar, cmap=plt.cm.inferno);
    
                    #this contours the data and smooths it
                    #contPlot = plt.contourf(lons, lats, plotvar, 100, transform=ccrs.PlateCarree(), vmin=minvar, vmax=maxvar, cmap=plt.cm.inferno);
                
                gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
                gl.xlabels_top = False;
                gl.ylabels_right = False;
                gl.xformatter = LONGITUDE_FORMATTER;
                gl.yformatter = LATITUDE_FORMATTER;
                gl.xlabel_style = {'size': ticksize};
                gl.ylabel_style = {'size': ticksize};
                            
                plt.xlabel("Longitude($^\circ$)", fontsize=labelsize);
                plt.ylabel("Latitude($^\circ$)", fontsize=labelsize);
        
                if region == "oceansoda_amazon_plume":
                    plt.xlim(-74, -30);
                    plt.ylim(-4, 26);
                elif region == "oceansoda_congo":
                    plt.xlim(-4, 16);
                    plt.ylim(-12, 6);
                elif region == "oceansoda_mediterranean":
                    plt.xlim(-6, 40);
                    plt.ylim(28, 48);
                else:
                    pass
                
                ax.coastlines();
                resol = '50m'  # use data at this scale
                land = cfeature.NaturalEarthFeature('physical', 'land', \
                scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
                ax.add_feature(land, facecolor='beige')
            
                        
                if video_var == "pCO2":
                    m = plt.cm.ScalarMappable(cmap=cmap_reversed)
                else:
                    m = plt.cm.ScalarMappable(cmap=plt.cm.inferno)
    
                m.set_array(plotvar)
                m.set_clim(minvar, maxvar)
                cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))
                
                if video_var == "DIC":
                    cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
                    plt.title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
                elif video_var == "SSS":
                    cb.set_label("Salinity (PSU)");
                    plt.title("Salinity for {0} {1}".format(date.year, format(date.month, "02d")));
                elif video_var == "SST":
                    cb.set_label("SST (($^\circ$ C))");
                    plt.title("SST for {0} {1}".format(date.year, format(date.month, "02d")));
                elif video_var == "pH_free":
                    cb.set_label("pH");
                    plt.title("pH for {0} {1}".format(date.year, format(date.month, "02d")));
                elif video_var == "pCO2":
                    cb.set_label("pCO$_{2}$ (ppm)");
                    plt.title("pCO2 (ppm) for {0} {1}".format(date.year, format(date.month, "02d")));
                elif video_var == "saturation_aragonite":
                    cb.set_label("$\Omega$ Aragonite");
                    plt.title("$\Omega$ Aragonite for {0} {1}".format(date.year, format(date.month, "02d")));
                elif video_var == "saturation_calcite":
                    cb.set_label("$\Omega$ Calcite");
                    plt.title("$\Omega$ Calcite for {0} {1}".format(date.year, format(date.month, "02d")));
                else:
                    pass
    
                #### Save plots for animations
                plt.savefig(path.join("os_plots\{0}_monthly_plots_for_animation\{0}_{1}_animation_images".format(region,video_var), "{0}_frame.png".format(format(iframe, "03d"))));
                plt.close();
                iframe=iframe+1;
            del varALL;
            
            #### Create video from images    
            image_folder = "os_plots\{0}_{1}_animation_images".format(region,video_var)
            video_name = "{0}_{1}.avi".format(region,video_var)
            
            images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
            frame = cv2.imread(os.path.join(image_folder, images[0]))
            height, width, layers = frame.shape
            
            video = cv2.VideoWriter(video_name, 0, 2, (width,height))
            
            for image in images:
                video.write(cv2.imread(os.path.join(image_folder, image)))
            
            cv2.destroyAllWindows()
            video.release()
    else:
        pass #don't make animations
    
    

    #### Maps of the variables
    for plot_var in video_vars:
        outputPath = "plots/";
        maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
        maskNC = Dataset(maskPath, 'r');
        
        dateToPlot1 = datetime(2015, 1, 1); #Which time slice to plot?
        dateToPlot2 = datetime(2015, 4, 1); #Which time slice to plot?
        dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
        dateToPlot4 = datetime(2015, 10, 1); #Which time slice to plot?
        
        datetoplotall=[dateToPlot1,dateToPlot2,dateToPlot3,dateToPlot4]
        
        mask = maskNC.variables[region][:];
        xmin = min(np.where(np.flipud(mask)==1)[1]);
        xmax = max(np.where(np.flipud(mask)==1)[1]);
        ymin = min(np.where(np.flipud(mask)==1)[0]);
        ymax = max(np.where(np.flipud(mask)==1)[0]);
        pad = 4; #size of xlim/ylim extents
        
        fig3 = plt.figure(figsize=(16,24))
        gs = fig3.add_gridspec(4, 1)
        ticksize = 10;

        dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
        varALL = dicNC.variables["{0}".format(plot_var)][:];
        hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
        
        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        lats = dicNC.variables["lat"];
        lons = dicNC.variables["lon"];
        
        counter=0
        for subplot_loop_no in datetoplotall:
        
            timeIndex = get_time_index(datetoplotall[counter], firstIndexSecs = int(dicNC.variables["time"][0]));
            data = dicNC.variables[plot_var][timeIndex, :, :];
            data[mask != 1] = np.nan;
            date = dates[timeIndex]

            if plot_var == "pCO2":
                cmap_touse = plt.cm.get_cmap('inferno_r');
            else:
                cmap_touse = plt.cm.inferno;
                
            if plot_var == "SST":
                data=data-273.15;    
                
            # these cover the full range in time
            # maxvar = np.nanmax(varALL);
            # minvar = np.nanmin(varALL);
            
            # arbitary limits
            if region == "oceansoda_amazon_plume":
                if plot_var == "DIC":
                    maxvar = 2000;
                    minvar = 1400; 
                elif plot_var == "SSS":
                    maxvar = 37;
                    minvar = 32; 
                elif plot_var == "SST":
                    maxvar = 30;
                    minvar = 22; 
                elif plot_var == "pH_free": 
                    maxvar = 8.5;
                    minvar = 7.5; 
                elif plot_var == "pCO2":
                    maxvar = 3000;
                    minvar = 200; 
                elif plot_var == "saturation_aragonite":
                    maxvar = 10;
                    minvar = 0; 
                elif plot_var == "saturation_calcite":                 
                    maxvar = 10;
                    minvar = 0; 
                else:
                    pass
                    
            elif region == "oceansoda_congo":
                if plot_var == "DIC":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SSS":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SST":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pH_free": 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pCO2":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_aragonite":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_calcite":                 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                else:
                    pass
            elif region == "oceansoda_mediterranean":
                if plot_var == "DIC":
                    maxvar = 2300;
                    minvar = 1950; 
                elif plot_var == "SSS":
                    maxvar = 37;
                    minvar = 32; 
                elif plot_var == "SST":
                    maxvar = 25;
                    minvar = 10; 
                elif plot_var == "pH_free": 
                    maxvar = 8.2;
                    minvar = 8; 
                elif plot_var == "pCO2":
                    maxvar = 600;
                    minvar = 200; 
                elif plot_var == "saturation_aragonite":
                    maxvar = 10;
                    minvar = 0; 
                elif plot_var == "saturation_calcite":                 
                    maxvar = 10;
                    minvar = 0; 
                else:
                    pass
            else:
                pass
            
            f3_ax = fig3.add_subplot(gs[counter:counter+1,0:1],projection=ccrs.PlateCarree())
            gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
            contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            gl.xlabels_top = False;
            gl.ylabels_right = False;
            gl.xformatter = LONGITUDE_FORMATTER;
            gl.yformatter = LATITUDE_FORMATTER;
            gl.xlabel_style = {'size': ticksize};
            gl.ylabel_style = {'size': ticksize};
            f3_ax.set_xlabel("Longitude($^\circ$)", fontsize=labelsize);
            f3_ax.set_ylabel("Latitude($^\circ$)", fontsize=labelsize);
            
            if region == "oceansoda_amazon_plume":
                f3_ax.set_xlim(-74, -30);
                f3_ax.set_ylim(-4, 26);
            elif region == "oceansoda_congo":
                f3_ax.set_xlim(-4, 16);
                f3_ax.set_ylim(-12, 6);
            elif region == "oceansoda_mediterranean":
                f3_ax.set_xlim(-6, 40);
                f3_ax.set_ylim(28, 48);
            else:
                pass
            
            f3_ax.coastlines();
            resol = '50m'  # use data at this scale
            land = cfeature.NaturalEarthFeature('physical', 'land', \
            scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
            f3_ax.add_feature(land, facecolor='beige')

            

            m = plt.cm.ScalarMappable(cmap=cmap_touse)
            m.set_array(plot_var)
            m.set_clim(minvar, maxvar)
            cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))   
        
            if plot_var == "DIC":
                cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
                f3_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
            elif plot_var == "SSS":
                cb.set_label("Salinity (PSU)");
                f3_ax.set_title("Salinity for {0} {1}".format(date.year, format(date.month, "02d")));
            elif plot_var == "SST":
                cb.set_label("SST (($^\circ$ C))");
                f3_ax.set_title("SST for {0} {1}".format(date.year, format(date.month, "02d")));
            elif plot_var == "pH_free":
                cb.set_label("pH");
                f3_ax.set_title("pH for {0} {1}".format(date.year, format(date.month, "02d")));
            elif plot_var == "pCO2":
                cb.set_label("pCO$_{2}$ (ppm)");
                f3_ax.set_title("pCO2 (ppm) for {0} {1}".format(date.year, format(date.month, "02d")));
            elif plot_var == "saturation_aragonite":
                cb.set_label("$\Omega$ Aragonite");
                f3_ax.set_title("$\Omega$ Aragonite for {0} {1}".format(date.year, format(date.month, "02d")));
            elif plot_var == "saturation_calcite":
                cb.set_label("$\Omega$ Calcite");
                f3_ax.set_title("$\Omega$ Calcite for {0} {1}".format(date.year, format(date.month, "02d")));
            else:
                pass
 
            counter=counter+1;
            
        #check that folder exists to save images, otherwise create it
        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        fig3.savefig(path.join("os_plots\{0}_seasonal".format(region,plot_var), "{0}.png".format(plot_var)));
        plt.close();
        del varALL;
        
    #### Hovmuller plots
    for plot_var in video_vars:
        outputPath = "plots/";
        maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
        maskNC = Dataset(maskPath, 'r');
        
        if region == "oceansoda_amazon_plume":
            long_to_slice1 = 55; #Which lat slice to plot?
            long_to_slice2 = 52; #Which lat slice to plot?
            long_to_slice3 = 50; #Which lat slice to plot?
            long_to_slice4 = 45; #Which lat slice to plot?
        elif region == "oceansoda_congo":
            long_to_slice1 = -10; #Which lat slice to plot?
            long_to_slice2 = -7.5; #Which lat slice to plot?
            long_to_slice3 = -5; #Which lat slice to plot?
            long_to_slice4 = -2.5; #Which lat slice to plot?
        elif region == "oceansoda_mediterranean":
            long_to_slice1 = 0; #Which lat slice to plot?
            long_to_slice2 = -10; #Which lat slice to plot?
            long_to_slice3 = -20; #Which lat slice to plot?
            long_to_slice4 = -30; #Which lat slice to plot? 
            
        long_to_slice=[long_to_slice1,long_to_slice2,long_to_slice3,long_to_slice4]
        
        mask = maskNC.variables[region][:];
        xmin = min(np.where(np.flipud(mask)==1)[1]);
        xmax = max(np.where(np.flipud(mask)==1)[1]);
        ymin = min(np.where(np.flipud(mask)==1)[0]);
        ymax = max(np.where(np.flipud(mask)==1)[0]);
        pad = 4; #size of xlim/ylim extents
        
        fig41 = plt.figure(figsize=(16,24))
        gs = fig41.add_gridspec(4, 1)
        ticksize = 10;

        dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
        varALL = dicNC.variables["{0}".format(plot_var)][:];
        #Create second variable to get the shorter length where TA and DIC overlap
        #which we need for plotting
        varALL_TADIC_overlap = atNC.variables["{0}".format("pCO2")][:];
        hasData = np.array([np.any(varALL_TADIC_overlap[t,:,:].mask==False) for t in range(0, varALL_TADIC_overlap.shape[0])]);

        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        lats = dicNC.variables["lat"];
        lons = dicNC.variables["lon"];
        
        counter=0
        for subplot_loop_no in long_to_slice:
        
            data = dicNC.variables[plot_var][:, :,180-subplot_loop_no ];
            data=data.T
            data.filled(np.nan) 
            
            if plot_var == "pCO2":
                cmap_touse = plt.cm.get_cmap('inferno_r');
            else:
                cmap_touse = plt.cm.inferno;
                
            if plot_var == "SST":
                data=data-273.15;    
                
            # these cover the full range in time
            # maxvar = np.nanmax(varALL);
            # minvar = np.nanmin(varALL);
            
            # arbitary limits
            if region == "oceansoda_amazon_plume":
                if plot_var == "DIC":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SSS":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SST":
                    maxvar = np.nanmax(varALL)-273.15;
                    minvar = np.nanmin(varALL)-273.15; 
                elif plot_var == "pH_free": 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pCO2":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_aragonite":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_calcite":                 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                else:
                    pass
                    
            elif region == "oceansoda_congo":
                if plot_var == "DIC":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SSS":
                    maxvar = np.nanmax(varALL)-273.15;
                    minvar = np.nanmin(varALL)-273.15; 
                elif plot_var == "SST":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pH_free": 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pCO2":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_aragonite":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_calcite":                 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                else:
                    pass
            elif region == "oceansoda_mediterranean":
                if plot_var == "DIC":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SSS":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "SST":
                    maxvar = np.nanmax(varALL)-273.15;
                    minvar = np.nanmin(varALL)-273.15; 
                elif plot_var == "pH_free": 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "pCO2":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_aragonite":
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                elif plot_var == "saturation_calcite":                 
                    maxvar = np.nanmax(varALL);
                    minvar = np.nanmin(varALL); 
                else:
                    pass
            else:
                pass
            
            
            f3_ax = fig41.add_subplot(gs[counter:counter+1,0:1])
            contPlot1 = f3_ax.pcolor(dates,lats,  data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            f3_ax.set_xlabel("Time", fontsize=labelsize);
            f3_ax.set_ylabel("Latitude($^\circ$N)", fontsize=labelsize);
            
            
            if region == "oceansoda_mediterranean":
                f3_ax.set_title("Longitude slice ({0}$^\circ$E)".format(subplot_loop_no*-1), fontsize=labelsize);
            elif region == "oceansoda_congo":
                f3_ax.set_title("Longitude slice ({0}$^\circ$E)".format(subplot_loop_no*-1), fontsize=labelsize);
            else:
                f3_ax.set_title("Longitude slice ({0}$^\circ$W)".format(subplot_loop_no), fontsize=labelsize);

            
            x=np.where(hasData)[0];
            firsttimeind=x[0];
            lastimeind=x[-1];

            if region == "oceansoda_amazon_plume":
                f3_ax.set_ylim(0, 24);
                f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
            elif region == "oceansoda_congo":
                f3_ax.set_ylim(-10, 3);
                f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
            elif region == "oceansoda_mediterranean":
                f3_ax.set_ylim(28, 48);
                f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
            else:
                pass
            cb = plt.colorbar(contPlot1)
        
            if plot_var == "DIC":
                cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
            elif plot_var == "SSS":
                cb.set_label("Salinity (PSU)");
            elif plot_var == "SST":
                cb.set_label("SST (($^\circ$ C))");
            elif plot_var == "pH_free":
                cb.set_label("pH");
            elif plot_var == "pCO2":
                cb.set_label("pCO$_{2}$ (ppm)");
            elif plot_var == "saturation_aragonite":
                cb.set_label("$\Omega$ Aragonite");
            elif plot_var == "saturation_calcite":
                cb.set_label("$\Omega$ Calcite");
            else:
                pass
 
            counter=counter+1;
            
            #check that folder exists to save images, otherwise create it
            from pathlib import Path
            Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
            
        fig41.savefig(path.join("os_plots\{0}_seasonal".format(region,plot_var), "{0}_hov.png".format(plot_var)));
        plt.close();
        del varALL;        
    
    
    hovmoll_vars = ["DIC", "AT", "pH_free","pCO2", "saturation_aragonite", "saturation_calcite"];
    
    #### Figure 3 Hovmuller plots combined
    counter_hov=0;
    
    fig42 = plt.figure(figsize=(16,24))
    gs = fig42.add_gridspec(6, 1)
    ticksize = 10;

    for plot_var in hovmoll_vars:
        outputPath = "plots/";
        maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
        maskNC = Dataset(maskPath, 'r');
        
        if region == "oceansoda_amazon_plume":
            long_to_slice = 52; #Which lat slice to plot?
        elif region == "oceansoda_congo":
            lat_to_slice = -6; #Which lat slice to plot?
        elif region == "oceansoda_mediterranean":
            long_to_slice = -20; #Which lat slice to plot?
                  
        mask = maskNC.variables[region][:];
        xmin = min(np.where(np.flipud(mask)==1)[1]);
        xmax = max(np.where(np.flipud(mask)==1)[1]);
        ymin = min(np.where(np.flipud(mask)==1)[0]);
        ymax = max(np.where(np.flipud(mask)==1)[0]);
        pad = 4; #size of xlim/ylim extents
        
        if plot_var == "AT":
            dates = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
            varALL = atNC.variables["{0}".format(plot_var)][:];
            #Create second variable to get the shorter length where TA and DIC overlap
            #which we need for plotting
            varALL_TADIC_overlap = atNC.variables["{0}".format("pCO2")][:];
            hasData = np.array([np.any(varALL_TADIC_overlap[t,:,:].mask==False) for t in range(0, varALL_TADIC_overlap.shape[0])]);
        else:
            dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
            varALL = dicNC.variables["{0}".format(plot_var)][:];
            #Create second variable to get the shorter length where TA and DIC overlap
            #which we need for plotting
            varALL_TADIC_overlap = atNC.variables["{0}".format("pCO2")][:];
            hasData = np.array([np.any(varALL_TADIC_overlap[t,:,:].mask==False) for t in range(0, varALL_TADIC_overlap.shape[0])]);


        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        lats = dicNC.variables["lat"];
        lons = dicNC.variables["lon"];

        

        if region == "oceansoda_congo": #latitude slice
            if plot_var == "AT":
                data = atNC.variables[plot_var][:, 90+lat_to_slice,: ]; 
                data=data.T
                data.filled(np.nan) 
            else:
                data = dicNC.variables[plot_var][:,90+lat_to_slice,: ];
                data=data.T
                data.filled(np.nan) 
        else: #longitude slice
            if plot_var == "AT":
                data = atNC.variables[plot_var][:, :,180-long_to_slice ];
                data=data.T
                data.filled(np.nan) 
            else:
                data = dicNC.variables[plot_var][:, :,180-long_to_slice ];
                data=data.T
                data.filled(np.nan) 
        
        if plot_var == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var == "DIC":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            if plot_var == "AT":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "SST":
                maxvar = np.nanmax(data)-273.15;
                minvar = np.nanmin(data)-273.15; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var == "DIC":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);
            elif plot_var == "AT":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);   
            elif plot_var == "SSS":
                maxvar = np.nanmax(data)-273.15;
                minvar = np.nanmin(data)-273.15; 
            elif plot_var == "SST":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);  
            elif plot_var == "pCO2":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);  
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var == "DIC":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "AT":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);  
            elif plot_var == "SST":
                maxvar = np.nanmax(data)-273.15;
                minvar = np.nanmin(data)-273.15; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data);  
            elif plot_var == "pCO2":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(data);
                minvar = np.nanmin(data); 
            else:
                pass
        else:
            pass
        
        
        f3_ax = fig42.add_subplot(gs[counter_hov:counter_hov+1,0:1])
        plt.xticks(fontsize=24)
        plt.yticks(fontsize=24)

        if region == "oceansoda_congo":
            contPlot1 = f3_ax.pcolor(dates,lons,  data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            f3_ax.set_xlabel("Time (years)", fontsize=24);
            f3_ax.set_ylabel("Longitude ($^\circ$E)", fontsize=24);
        else:
            contPlot1 = f3_ax.pcolor(dates,lats,  data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            f3_ax.set_xlabel("Time (years)", fontsize=24);
            f3_ax.set_ylabel("Latitude ($^\circ$N)", fontsize=24);
        
        # if region == "oceansoda_mediterranean":
        #     f3_ax.set_title("Longitude slice ({0}$^\circ$E)".format(long_to_slice*-1), fontsize=labelsize);
        # elif region == "oceansoda_congo":
        #     f3_ax.set_title("Latitude slice ({0}$^\circ$S)".format(lat_to_slice), fontsize=labelsize);
        # else:
        #     f3_ax.set_title("Longitude slice ({0}$^\circ$W)".format(long_to_slice), fontsize=labelsize);

        
        x=np.where(hasData)[0];
        firsttimeind=x[0];
        lastimeind=x[-1];

        if region == "oceansoda_amazon_plume":
            f3_ax.set_ylim(4, 24);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        elif region == "oceansoda_congo":
            f3_ax.set_ylim(-2, 12);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_ylim(32, 36);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        else:
            pass
        cb = plt.colorbar(contPlot1)
        cb.ax.tick_params(labelsize=24)

        if plot_var == "DIC":
            cb.set_label("DIC ($\mu mol \ kg^{-1}$)", fontsize=24);
        elif plot_var == "AT":
            cb.set_label("TA ($\mu mol \ kg^{-1}$)", fontsize=24);
        elif plot_var == "SSS":
            cb.set_label("Salinity (PSU)", fontsize=24);
        elif plot_var == "SST":
            cb.set_label("SST (($^\circ$ C))", fontsize=24);
        elif plot_var == "pH_free":
            cb.set_label("pH", fontsize=24);
        elif plot_var == "pCO2":
            cb.set_label("pCO$_{2}$ (ppm)", fontsize=24);
        elif plot_var == "saturation_aragonite":
            cb.set_label("$\Omega$ Aragonite", fontsize=24);
        elif plot_var == "saturation_calcite":
            cb.set_label("$\Omega$ Calcite", fontsize=24);
        else:
            pass
 
        if counter_hov==0:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(a)',transform=f3_ax.transAxes, fontsize=24, weight="bold");
        elif counter_hov==1:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(b)', transform=f3_ax.transAxes,fontsize=24, weight="bold");
        elif counter_hov==2:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(c)',transform=f3_ax.transAxes,  fontsize=24, weight="bold");
        elif counter_hov==3:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(d)',transform=f3_ax.transAxes,  fontsize=24, weight="bold");
        elif counter_hov==4:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(e)',transform=f3_ax.transAxes, fontsize=24, weight="bold");
        elif counter_hov==5:
            f3_ax_txt = f3_ax.text(-0.1, 1.05, '(f)',transform=f3_ax.transAxes,  fontsize=24, weight="bold");
        
        counter_hov=counter_hov+1;
        
        #check that folder exists to save images, otherwise create it
        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    #fig42.tight_layout()

    fig42.savefig(path.join("os_plots\{0}_seasonal".format(region), "_hov_multivar_all_vars.png"));
    plt.close();
    del varALL;        
    #### Figure 4 - Timeseries plots
    
    #plotting parameters
    fontsize = 12;
    figsizex = 9.0;
    figsizey = 2.0*5;
    figsize = (16,24);
    
    plt.figure(figsize=figsize);
    plt.title("carbonate_timeseries_{0}".format(region))

    ax11 = plt.subplot(6,1,1);
    ax11.plot(time, meanPlumeDIC, 'r--');#, label="plume");
    ax11.plot(time, meanNotPlumeDIC, 'r:');#, label="non-plume");
    ax11.plot(time, meanAllDIC, 'r');#, label="whole region");
    ax11.set_ylabel("DIC\n$(\mu mol \ kg^{-1})$", fontsize=24);
    ax11.text(-0.12, 1.06, '(a)', transform=ax11.transAxes, fontsize=24, weight="bold");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax11.set_xlim(dates[firsttimeind], dates[lastimeind]);
            
    ax21 = plt.subplot(6,1,2);
    ax21.plot(time, meanPlumeAT, 'b--', label="plume");
    ax21.plot(time, meanNotPlumeAT, 'b:', label="non-plume");
    ax21.plot(time, meanAllAT, 'b', label="whole region");
    ax21.set_ylabel("TA\n$(\mu mol \ kg^{-1})$", fontsize=24);
    ax21.text(-0.12, 1.06, '(b)', transform=ax21.transAxes, fontsize=24, weight="bold");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax21.set_xlim(dates[firsttimeind], dates[lastimeind]);

    ax31 = plt.subplot(6,1,3);
    ax31.plot(time, meanPlumepH_dic, 'r--', label="plume");
    ax31.plot(time, meanNotPlumepH_dic, 'r:', label="non-plume");
    ax31.plot(time, meanAllpH_dic, 'r', label="whole region");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)

    ax31.plot(time,meanPlumepH_at , 'b--', label="plume");
    ax31.plot(time,meanNotPlumepH_at, 'b:', label="non-plume");
    ax31.plot(time,meanAllpH_at , 'b', label="whole region");
    ax31.set_ylabel("pH", fontsize=24);
    ax31.text(-0.12, 1.06, '(c)', transform=ax31.transAxes, fontsize=24, weight="bold");
    ax31.set_xlim(dates[firsttimeind], dates[lastimeind]);

    ax41 = plt.subplot(6,1,4);
    ax41.plot(time, meanPlumepco2_dic, 'r--', label="plume");
    ax41.plot(time, meanNotPlumepco2_dic, 'r:', label="non-plume");
    ax41.plot(time,meanAllpco2_dic , 'r', label="whole region");
    
    ax41.plot(time, meanPlumepco2_at , 'b--', label="plume");
    ax41.plot(time, meanNotPlumepco2_at, 'b:', label="non-plume");
    ax41.plot(time, meanAllpco2_at, 'b', label="whole region");
    ax41.set_ylabel("pCO$_{2}$ \n (ppm)", fontsize=24);
    ax41.text(-0.12, 1.06, '(d)', transform=ax41.transAxes, fontsize=24, weight="bold");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax41.set_xlim(dates[firsttimeind], dates[lastimeind]);

    ax51 = plt.subplot(6,1,5);
    ax51.plot(time,meanPlumeOmegaAragonite_dic , 'r--', label="plume");
    ax51.plot(time, meanNotPlumeOmegaAragonite_dic, 'r:', label="non-plume");
    ax51.plot(time, meanAllOmegaAragonite_dic, 'r', label="whole region");
    
    ax51.plot(time,  meanPlumeOmegaAragonite_at, 'b--', label="plume");
    ax51.plot(time, meanNotPlumeOmegaAragonite_at, 'b:', label="non-plume");
    ax51.plot(time, meanAllOmegaAragonite_at, 'b', label="whole region");
    ax51.set_ylabel("Aragonite \n saturation \n state", fontsize=24);
    ax51.text(-0.12, 1.06, '(e)', transform=ax51.transAxes, fontsize=24, weight="bold");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax51.set_xlim(dates[firsttimeind], dates[lastimeind]);

    ax61 = plt.subplot(6,1,6);
    ax61.plot(time,meanPlumeOmegaCalcite_dic , 'r--', label="plume");
    ax61.plot(time, meanNotPlumeOmegaCalcite_dic, 'r:', label="non-plume");
    ax61.plot(time, meanAllOmegaCalcite_dic, 'r', label="whole region");
    
    ax61.plot(time,  meanPlumeOmegaCalcite_at, 'b--', label="plume");
    ax61.plot(time, meanNotPlumeOmegaCalcite_at, 'b:', label="non-plume");
    ax61.plot(time, meanAllOmegaCalcite_at, 'b', label="whole region");
    ax61.set_ylabel("Calcite\n saturation \n state", fontsize=24);
    ax61.set_xlabel("Time (years)", fontsize=24);
    ax61.text(-0.12, 1.06, '(f)', transform=ax61.transAxes, fontsize=24, weight="bold");
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    ax61.set_xlim(dates[firsttimeind], dates[lastimeind]);

    #legend plotting
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k--', label="Plume");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k:', label="Non-plume");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k', label="Whole region");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'rs', label="Using best SST & SSS inputs from DIC algorithm");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'bs', label="Using best SST & SSS inputs from TA algorithm");
    ax11.legend(loc="center",bbox_to_anchor=(0.5, 1.5), fontsize=24,ncol=2);
    
    # plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    #plt.tight_layout()

    plt.savefig("os_plots\\fig4_carbonate_timeseries_{0}.png".format(region));
    plt.savefig("os_plots\\fig4_carbonate_timeseries_{0}.pdf".format(region));    
    
    

    #### DIC and TA in 2015
    plot_var="DIC"
    plot_var2="AT"
    outputPath = "plots/";
    maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
    maskNC = Dataset(maskPath, 'r');
    
    dateToPlot1 = datetime(2015, 1, 1); #Which time slice to plot?
    dateToPlot2 = datetime(2015, 4, 1); #Which time slice to plot?
    dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
    dateToPlot4 = datetime(2015, 10, 1); #Which time slice to plot?
    
    datetoplotall=[dateToPlot1,dateToPlot2,dateToPlot3,dateToPlot4]
    
    mask = maskNC.variables[region][:];
    xmin = min(np.where(np.flipud(mask)==1)[1]);
    xmax = max(np.where(np.flipud(mask)==1)[1]);
    ymin = min(np.where(np.flipud(mask)==1)[0]);
    ymax = max(np.where(np.flipud(mask)==1)[0]);
    pad = 4; #size of xlim/ylim extents
    
    fig5 = plt.figure(figsize=(12,15))
    gs = fig5.add_gridspec(4, 2)
    ticksize = 10;

    dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
    varALL = dicNC.variables["{0}".format(plot_var)][:];
    hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
    
    dates2 = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
    varALL2 = atNC.variables["{0}".format(plot_var2)][:];
    hasData2 = np.array([np.any(varALL2[t,:,:].mask==False) for t in range(0, varALL2.shape[0])]);
    
    
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    lats = dicNC.variables["lat"];
    lons = dicNC.variables["lon"];
    
    counter=0
    for subplot_loop_no in datetoplotall:
        
        #### DIC half of plot first
        timeIndex = get_time_index(datetoplotall[counter], firstIndexSecs = int(dicNC.variables["time"][0]));
        data = dicNC.variables[plot_var][timeIndex, :, :];
        data[mask != 1] = np.nan;
        date = dates[timeIndex]

        if plot_var == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var == "DIC":
                maxvar = 2000;
                minvar = 1400; 
            elif plot_var == "AT":
                maxvar = 2200;
                minvar = 1400; 
            elif plot_var == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var == "SST":
                maxvar = 30;
                minvar = 22; 
            elif plot_var == "pH_free": 
                maxvar = 8.5;
                minvar = 7.5; 
            elif plot_var == "pCO2":
                maxvar = 3000;
                minvar = 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var == "DIC":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "SST":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var == "DIC":
                maxvar = 2300;
                minvar = 1950; 
            elif plot_var == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var == "SST":
                maxvar = 25;
                minvar = 10; 
            elif plot_var == "pH_free": 
                maxvar = 8.2;
                minvar = 8; 
            elif plot_var == "pCO2":
                maxvar = 600;
                minvar = 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig5.add_subplot(gs[counter:counter+1,0:1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-4, 16);
            f3_ax.set_ylim(-12, 6);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')

        
        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))  

    
        if plot_var == "DIC":
            cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
            f3_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
        elif plot_var == "AT":
            cb.set_label("TA ($\mu mol \ kg^{-1}$)");
            f3_ax.set_title("TA for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "SSS":
            cb.set_label("Salinity (PSU)");
            f3_ax.set_title("Salinity for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "SST":
            cb.set_label("SST (($^\circ$ C))");
            f3_ax.set_title("SST for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "pH_free":
            cb.set_label("pH");
            f3_ax.set_title("pH for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "pCO2":
            cb.set_label("pCO$_{2}$ (ppm)");
            f3_ax.set_title("pCO2 (ppm) for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "saturation_aragonite":
            cb.set_label("$\Omega$ Aragonite");
            f3_ax.set_title("$\Omega$ Aragonite for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var == "saturation_calcite":
            cb.set_label("$\Omega$ Calcite");
            f3_ax.set_title("$\Omega$ Calcite for {0} {1}".format(date.year, format(date.month, "02d")));
        else:
            pass
  
        counter=counter+1;
    
        #### TA half of plot first
    counter=0
    for subplot_loop_no in datetoplotall:
    
        timeIndex = get_time_index(datetoplotall[counter], firstIndexSecs = int(atNC.variables["time"][0]));
        data = atNC.variables[plot_var2][timeIndex, :, :];
        data[mask != 1] = np.nan;
        date = dates[timeIndex]

        if plot_var2 == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var2 == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var2 == "DIC":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT":
                maxvar = 2400;
                minvar = 600; 
            elif plot_var2 == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var2 == "SST":
                maxvar = 30;
                minvar = 22; 
            elif plot_var2 == "pH_free": 
                maxvar = 8.5;
                minvar = 7.5; 
            elif plot_var2 == "pCO2":
                maxvar = 3000;
                minvar = 200; 
            elif plot_var2 == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var2 == "DIC":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "SSS":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "SST":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "pH_free": 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "pCO2":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "saturation_aragonite":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var2 == "DIC":
                maxvar = 2300;
                minvar = 1950; 
            elif plot_var2 == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var2 == "SST":
                maxvar = 25;
                minvar = 10; 
            elif plot_var2 == "pH_free": 
                maxvar = 8.2;
                minvar = 8; 
            elif plot_var2 == "pCO2":
                maxvar = 600;
                minvar = 200; 
            elif plot_var2 == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig5.add_subplot(gs[counter:counter+1,1:2],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-4, 16);
            f3_ax.set_ylim(-12, 6);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')
        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var2)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))    
        if plot_var2 == "DIC":
            cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
            f3_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
        elif plot_var2 == "AT":
            cb.set_label("TA ($\mu mol \ kg^{-1}$)");
            f3_ax.set_title("TA for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "SSS":
            cb.set_label("Salinity (PSU)");
            f3_ax.set_title("Salinity for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "SST":
            cb.set_label("SST (($^\circ$ C))");
            f3_ax.set_title("SST for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "pH_free":
            cb.set_label("pH");
            f3_ax.set_title("pH for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "pCO2":
            cb.set_label("pCO$_{2}$ (ppm)");
            f3_ax.set_title("pCO2 (ppm) for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "saturation_aragonite":
            cb.set_label("$\Omega$ Aragonite");
            f3_ax.set_title("$\Omega$ Aragonite for {0} {1}".format(date.year, format(date.month, "02d")));
        elif plot_var2 == "saturation_calcite":
            cb.set_label("$\Omega$ Calcite");
            f3_ax.set_title("$\Omega$ Calcite for {0} {1}".format(date.year, format(date.month, "02d")));
        else:
            pass

        counter=counter+1;
            
    #check that folder exists to save images, otherwise create it
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    fig5.savefig(path.join("os_plots\{0}_seasonal".format(region), "{0}_{1}.png".format(plot_var,plot_var2)));
    plt.close();
    del varALL;        
    
    
    #### Figure 2 - DIC and TA mmm
    plot_var="DIC"
    plot_var2="AT"
    outputPath = "plots/";
    maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
    maskNC = Dataset(maskPath, 'r');
    
    dateToPlot1 = datetime(2015, 1, 1); #Which time slice to plot?
    dateToPlot2 = datetime(2015, 4, 1); #Which time slice to plot?
    dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
    dateToPlot4 = datetime(2015, 10, 1); #Which time slice to plot?
    
    datetoplotall=[dateToPlot1,dateToPlot2,dateToPlot3,dateToPlot4]
    
    mask = maskNC.variables[region][:];
    xmin = min(np.where(np.flipud(mask)==1)[1]);
    xmax = max(np.where(np.flipud(mask)==1)[1]);
    ymin = min(np.where(np.flipud(mask)==1)[0]);
    ymax = max(np.where(np.flipud(mask)==1)[0]);
    pad = 4; #size of xlim/ylim extents
    
    fig6 = plt.figure(figsize=(12,15))
    gs = fig6.add_gridspec(4, 2)
    ticksize = 10;

    dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
    varALL = dicNC.variables["{0}".format(plot_var)][:];
    hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
    
    dates2 = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
    varALL2 = atNC.variables["{0}".format(plot_var2)][:];
    hasData2 = np.array([np.any(varALL2[t,:,:].mask==False) for t in range(0, varALL2.shape[0])]);
    
    
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    lats = dicNC.variables["lat"];
    lons = dicNC.variables["lon"];
    
    counter=0
    for subplot_loop_no in datetoplotall:
        
        #### DIC half of plot first

        if counter==0:
            months_average=number_list_jfm
        elif counter==1:
            months_average=number_list_amj
        elif counter==2:
            months_average=number_list_jas
        elif counter==3:
            months_average=number_list_ond

        data = dicNC.variables[plot_var][months_average, :, :];
        data[:,mask != 1] = np.nan;
        data= np.nanmean(data, axis=0);

        
        if plot_var == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var == "DIC":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT":
                maxvar = 2400;
                minvar = 600; 
            elif plot_var == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var == "SST":
                maxvar = 30;
                minvar = 22; 
            elif plot_var == "pH_free": 
                maxvar = 8.5;
                minvar = 7.5; 
            elif plot_var == "pCO2":
                maxvar = 3000;
                minvar = 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var == "DIC":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "AT":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "SST":
                maxvar = np.nanmax(varALL)-273.15;
                minvar = np.nanmin(varALL)-273.15; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var == "DIC":
                maxvar = 2300;
                minvar = 1950; 
            elif plot_var == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var == "SST":
                maxvar = 25;
                minvar = 10; 
            elif plot_var == "pH_free": 
                maxvar = 8.2;
                minvar = 8; 
            elif plot_var == "pCO2":
                maxvar = 600;
                minvar = 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1,0:1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-4, 16);
            f3_ax.set_ylim(-12, 6);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')

        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))  

        if plot_var == "DIC":
            cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
        elif plot_var == "AT":
            cb.set_label("TA ($\mu mol \ kg^{-1}$)");
        elif plot_var == "SSS":
            cb.set_label("Salinity (PSU)");
        elif plot_var == "SST":
            cb.set_label("SST (($^\circ$ C))");
        elif plot_var == "pH_free":
            cb.set_label("pH");
        elif plot_var == "pCO2":
            cb.set_label("pCO$_{2}$ (ppm)");
        elif plot_var == "saturation_aragonite":
            cb.set_label("$\Omega$ Aragonite");
        elif plot_var == "saturation_calcite":
            cb.set_label("$\Omega$ Calcite");
        else:
            pass

        if counter==0:
            f3_ax.set_title("JFM");
        elif counter==1:
            f3_ax.set_title("AMJ");
        elif counter==2:
            f3_ax.set_title("JAS");
        elif counter==3:
            f3_ax.set_title("OND");

        counter=counter+1;
    
        #### TA half of plot first
    counter=0
    for subplot_loop_no in datetoplotall:
    
        if counter==0:
            months_average=number_list_jfm
        elif counter==1:
            months_average=number_list_amj
        elif counter==2:
            months_average=number_list_jas
        elif counter==3:
            months_average=number_list_ond

        data = atNC.variables[plot_var2][months_average, :, :];
        data[:,mask != 1] = np.nan;
        data= np.nanmean(data, axis=0);

        if plot_var2 == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var2 == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var2 == "DIC":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var2 == "AT":
                maxvar = 2400;
                minvar = 600; 
            elif plot_var2 == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var2 == "SST":
                maxvar = 30;
                minvar = 22; 
            elif plot_var2 == "pH_free": 
                maxvar = 8.5;
                minvar = 7.5; 
            elif plot_var2 == "pCO2":
                maxvar = 3000;
                minvar = 200; 
            elif plot_var2 == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var2 == "DIC":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "SSS":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "SST":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "pH_free": 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "pCO2":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "saturation_aragonite":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var2 == "DIC":
                maxvar = 2300;
                minvar = 1950; 
            elif plot_var2 == "SSS":
                maxvar = 37;
                minvar = 32; 
            elif plot_var2 == "SST":
                maxvar = 25;
                minvar = 10; 
            elif plot_var2 == "pH_free": 
                maxvar = 8.2;
                minvar = 8; 
            elif plot_var2 == "pCO2":
                maxvar = 600;
                minvar = 200; 
            elif plot_var2 == "saturation_aragonite":
                maxvar = 10;
                minvar = 0; 
            elif plot_var2 == "saturation_calcite":                 
                maxvar = 10;
                minvar = 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1,1:2],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'latitude', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'longitude', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-4, 16);
            f3_ax.set_ylim(-12, 6);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')
        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var2)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))    
        if plot_var2 == "DIC":
            cb.set_label("DIC ($\mu mol \ kg^{-1}$)");
        elif plot_var2 == "AT":
            cb.set_label("TA ($\mu mol \ kg^{-1}$)");
        elif plot_var2 == "SSS":
            cb.set_label("Salinity (PSU)");
        elif plot_var2 == "SST":
            cb.set_label("SST (($^\circ$ C))");
        elif plot_var2 == "pH_free":
            cb.set_label("pH");
        elif plot_var2 == "pCO2":
            cb.set_label("pCO$_{2}$ (ppm)");
        elif plot_var2 == "saturation_aragonite":
            cb.set_label("$\Omega$ Aragonite");
        elif plot_var2 == "saturation_calcite":
            cb.set_label("$\Omega$ Calcite");
        else:
            pass
        
        if counter==0:
            f3_ax.set_title("JFM");
        elif counter==1:
            f3_ax.set_title("AMJ");
        elif counter==2:
            f3_ax.set_title("JAS");
        elif counter==3:
            f3_ax.set_title("OND");
        counter=counter+1;
            
    #check that folder exists to save images, otherwise create it
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    fig6.savefig(path.join("os_plots\{0}_seasonal".format(region), "{0}_{1}_monthly_averages.png".format(plot_var,plot_var2)));
    plt.close();
    del varALL;           
    
    #### Figure 2a - DIC mmm
    plot_var="DIC"
    outputPath = "plots/";
    maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
    maskNC = Dataset(maskPath, 'r');
    
    dateToPlot1 = datetime(2015, 1, 1); #Which time slice to plot?
    dateToPlot2 = datetime(2015, 4, 1); #Which time slice to plot?
    dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
    dateToPlot4 = datetime(2015, 10, 1); #Which time slice to plot?
    
    datetoplotall=[dateToPlot1,dateToPlot2,dateToPlot3,dateToPlot4]
    
    mask = maskNC.variables[region][:];
    xmin = min(np.where(np.flipud(mask)==1)[1]);
    xmax = max(np.where(np.flipud(mask)==1)[1]);
    ymin = min(np.where(np.flipud(mask)==1)[0]);
    ymax = max(np.where(np.flipud(mask)==1)[0]);
    pad = 4; #size of xlim/ylim extents
    
                
    if region == "oceansoda_amazon_plume":
        fig6 = plt.figure(figsize=(20,12))
    elif region == "oceansoda_congo":
        fig6 = plt.figure(figsize=(20,12))#15,10 IS ALMOST SQUARE 
    elif region == "oceansoda_mediterranean":
        fig6 = plt.figure(figsize=(15,10))
    else:
        pass
        
    gs = fig6.add_gridspec(2, 2)
    ticksize = 10;

    dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
    varALL = dicNC.variables["{0}".format(plot_var)][:];
    hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
    
    dates2 = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
    varALL2 = atNC.variables["{0}".format(plot_var2)][:];
    hasData2 = np.array([np.any(varALL2[t,:,:].mask==False) for t in range(0, varALL2.shape[0])]);
    
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    lats = dicNC.variables["lat"];
    lons = dicNC.variables["lon"];
    
    counter=0
    for subplot_loop_no in datetoplotall:
        
        #### DIC half of plot first

        if counter==0:
            months_average=number_list_jfm
        elif counter==1:
            months_average=number_list_amj
        elif counter==2:
            months_average=number_list_jas
        elif counter==3:
            months_average=number_list_ond

        data = dicNC.variables[plot_var][months_average, :, :];
        data[:,mask != 1] = np.nan;
        data= np.nanmean(data, axis=0);

        #find the maximum and minimum of data in the subplots so they are all
        #standardised and         
        data1 = dicNC.variables[plot_var][number_list_jfm, :, :];
        data1[:,mask != 1] = np.nan;
        data1= np.nanmean(data1, axis=0);
        data2 = dicNC.variables[plot_var][number_list_amj, :, :];
        data2[:,mask != 1] = np.nan;
        data2= np.nanmean(data2, axis=0);
        data3= dicNC.variables[plot_var][number_list_jas, :, :];
        data3[:,mask != 1] = np.nan;
        data3= np.nanmean(data3, axis=0);        
        data4 = dicNC.variables[plot_var][number_list_ond, :, :];
        data4[:,mask != 1] = np.nan;
        data4= np.nanmean(data4, axis=0);

        dataall_subplots=([data1,data2,data3,data4]);

        if plot_var == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);#2100;
                minvar = np.nanmin(dataall_subplots);#1000; 
            elif plot_var == "AT":
                maxvar = np.nanmax(dataall_subplots);#2400;
                minvar = np.nanmin(dataall_subplots);#600; 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);#37;
                minvar = np.nanmin(dataall_subplots);#32; 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots);#30;
                minvar = np.nanmin(dataall_subplots);#22; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);#8.5;
                minvar = np.nanmin(dataall_subplots);#7.5; 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);#3000;
                minvar = np.nanmin(dataall_subplots);#200; 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);#0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);#0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "AT":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots)-273.15;
                minvar = np.nanmin(dataall_subplots)-273.15; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);#
                minvar = np.nanmin(dataall_subplots);# 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);#2300;
                minvar = np.nanmin(dataall_subplots);# 1950; 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);#37;
                minvar = np.nanmin(dataall_subplots);# 32; 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots);#25;
                minvar = np.nanmin(dataall_subplots);# 10; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);#8.2;
                minvar = np.nanmin(dataall_subplots);# 8; 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);#600;
                minvar = np.nanmin(dataall_subplots);# 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);# 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);# 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': 16};
        gl.ylabel_style = {'size': 16};
        f3_ax.text(-0.12, 0.55, 'Latitude $(^{\circ})$', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',fontsize=16,
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.13, 'Longitude $(^{\circ})$', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',fontsize=16,
            transform=f3_ax.transAxes)
        
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-3, 18);
            f3_ax.set_ylim(-11, 5);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')

        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))  
        cb.ax.tick_params(labelsize=16)
        
        if plot_var == "DIC":
            cb.ax.set_title("DIC ($\mu mol \ kg^{-1}$)",fontsize=16);
        elif plot_var == "AT":
            cb.ax.set_title("TA ($\mu mol \ kg^{-1}$)",fontsize=16);
        elif plot_var == "SSS":
            cb.ax.set_title("Salinity (PSU)",fontsize=16);
        elif plot_var == "SST":
            cb.ax.set_title("SST (($^\circ$ C))",fontsize=16);
        elif plot_var == "pH_free":
            cb.ax.set_title("pH",fontsize=16);
        elif plot_var == "pCO2":
            cb.ax.set_title("pCO$_{2}$ (ppm)",fontsize=16);
        elif plot_var == "saturation_aragonite":
            cb.ax.set_title("$\Omega$ Aragonite",fontsize=16);
        elif plot_var == "saturation_calcite":
            cb.ax.set_title("$\Omega$ Calcite",fontsize=16);
        else:
            pass

        # #add a marker for the mouth of the Amazon
        if region == "oceansoda_amazon_plume":
            # Mark some particular places with a small circle and a name label...
            # Define some test points with latitude and longitude coordinates.
            city_data = [('Amazon River mouth', -1.455833,-48.503889 ), 
                          ('Orinoco River mouth', 8.616667, -62.250000),
                          ('Maroni River mouth', 5.745833, -53.968333)]
                            #('Essequibo River mouth', 7.033333, -58.450000)
            crs_latlon = ccrs.PlateCarree()
            # Place a single marker point and a text annotation at each place.
            for name, lat, lon in city_data:
            
                plt.plot(lon, lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                          markerfacecolor='black', markeredgecolor='white',transform=crs_latlon)
                # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
                # so for this one we transform the coordinates with a Cartopy call.
                at_x, at_y = f3_ax.projection.transform_point(lon, lat,src_crs=crs_latlon)
                if name =="Maroni River mouth":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-150, -40), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Orinoco River mouth":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-100, -30), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                else:
                # place_labels=plt.text(lon, lat,name,fontsize=15,weight='bold')
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-200, -10), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                    
        elif region == "oceansoda_congo":
            # Mark some particular places with a small circle and a name label...
            # Define some test points with latitude and longitude coordinates.
            city_data = [('Congo River', -6.060116, 12.494964 ), 
                          ('Niger Delta', 4.30996655147119, 6.22693906885239),
                          ('Ogooué River', -1.027200, 8.884800),
                          ('Sanaga River', 3.559338, 9.652175)]
            crs_latlon = ccrs.PlateCarree()
            # Place a single marker point and a text annotation at each place.
            for name, lat, lon in city_data:
            
                plt.plot(lon, lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                          markerfacecolor='black', markeredgecolor='white',transform=crs_latlon)
                # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
                # so for this one we transform the coordinates with a Cartopy call.
                at_x, at_y = f3_ax.projection.transform_point(lon, lat,src_crs=crs_latlon)
                if name =="Congo River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-20, 50), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Niger Delta":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(130, -5), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Ogooué River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(50, 0), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))   
                elif name =="Sanaga River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(50, -40), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5)) 
                    
        # # add the slice as a line
        if region == "oceansoda_amazon_plume":
            ad_lat, ad_lon = 4, -52
            liv_lat, liv_lon = 23, -52
            x1, y1 = f3_ax.projection.transform_point(ad_lon, ad_lat,src_crs=crs_latlon)
            x2, y2 = f3_ax.projection.transform_point(liv_lon, liv_lat,src_crs=crs_latlon)
    
            f3_ax.plot([x1, x2], [y1, y2],
                      color='black', linewidth=1, marker='s', markersize=3, markerfacecolor='black',linestyle='--',transform=crs_latlon)
        elif region == "oceansoda_congo":
            ad_lat, ad_lon = -6, -1.5
            liv_lat, liv_lon = -6, 11.5
            x1, y1 = f3_ax.projection.transform_point(ad_lon, ad_lat,src_crs=crs_latlon)
            x2, y2 = f3_ax.projection.transform_point(liv_lon, liv_lat,src_crs=crs_latlon)
    
            f3_ax.plot([x1, x2], [y1, y2],
                      color='black', linewidth=1, marker='s', markersize=3, markerfacecolor='black',linestyle='--',transform=crs_latlon)
           

        if counter==0:
            f3_ax.set_title("JFM",fontsize=16);
        elif counter==1:
            f3_ax.set_title("AMJ",fontsize=16);
        elif counter==2:
            f3_ax.set_title("JAS",fontsize=16);
        elif counter==3:
            f3_ax.set_title("OND",fontsize=16);
            
        if counter==0:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(a)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==1:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(b)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==2:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(c)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==3:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(d)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
            
        counter=counter+1;
            
    #check that folder exists to save images, otherwise create it
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    fig6.savefig(path.join("os_plots\{0}_seasonal".format(region), "{0}_monthly_averages.png".format(plot_var)));
    fig6.tight_layout()
    plt.close();
    del varALL;     
    
    #### Figure 2b - TA mmm
    plot_var="AT"
    outputPath = "plots/";
    maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
    maskNC = Dataset(maskPath, 'r');
    
    dateToPlot1 = datetime(2015, 1, 1); #Which time slice to plot?
    dateToPlot2 = datetime(2015, 4, 1); #Which time slice to plot?
    dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
    dateToPlot4 = datetime(2015, 10, 1); #Which time slice to plot?
    
    datetoplotall=[dateToPlot1,dateToPlot2,dateToPlot3,dateToPlot4]
    
    mask = maskNC.variables[region][:];
    xmin = min(np.where(np.flipud(mask)==1)[1]);
    xmax = max(np.where(np.flipud(mask)==1)[1]);
    ymin = min(np.where(np.flipud(mask)==1)[0]);
    ymax = max(np.where(np.flipud(mask)==1)[0]);
    pad = 4; #size of xlim/ylim extents
    
    if region == "oceansoda_amazon_plume":
        fig6 = plt.figure(figsize=(20,12))
    elif region == "oceansoda_congo":
        fig6 = plt.figure(figsize=(20,12))
    elif region == "oceansoda_mediterranean":
        fig6 = plt.figure(figsize=(15,10))
    else:
        pass
    
    gs = fig6.add_gridspec(2, 2)
    ticksize = 10;

    dates = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
    varALL = atNC.variables["{0}".format(plot_var)][:];
    hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
    
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    lats = atNC.variables["lat"];
    lons = atNC.variables["lon"];
    
    counter=0
    for subplot_loop_no in datetoplotall:
        
        #### DIC half of plot first

        if counter==0:
            months_average=number_list_jfm
        elif counter==1:
            months_average=number_list_amj
        elif counter==2:
            months_average=number_list_jas
        elif counter==3:
            months_average=number_list_ond

        data = atNC.variables[plot_var][months_average, :, :];
        data[:,mask != 1] = np.nan;
        data= np.nanmean(data, axis=0);


        #find the maximum and minimum of data in the subplots so they are all
        #standardised and         
        data1 = atNC.variables[plot_var][number_list_jfm, :, :];
        data1[:,mask != 1] = np.nan;
        data1= np.nanmean(data1, axis=0);
        data2 = atNC.variables[plot_var][number_list_amj, :, :];
        data2[:,mask != 1] = np.nan;
        data2= np.nanmean(data2, axis=0);
        data3= atNC.variables[plot_var][number_list_jas, :, :];
        data3[:,mask != 1] = np.nan;
        data3= np.nanmean(data3, axis=0);        
        data4 = atNC.variables[plot_var][number_list_ond, :, :];
        data4[:,mask != 1] = np.nan;
        data4= np.nanmean(data4, axis=0);

        dataall_subplots=([data1,data2,data3,data4]);
        
        if plot_var == "pCO2":
            cmap_touse = plt.cm.get_cmap('inferno_r');
        else:
            cmap_touse = plt.cm.inferno;
            
        if plot_var == "SST":
            data=data-273.15;    
            
        # these cover the full range in time
        # maxvar = np.nanmax(varALL);
        # minvar = np.nanmin(varALL);
        
        # arbitary limits
        if region == "oceansoda_amazon_plume":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);#2100;
                minvar = np.nanmin(dataall_subplots);#1000; 
            elif plot_var == "AT":
                maxvar = np.nanmax(dataall_subplots);#2400;
                minvar = np.nanmin(dataall_subplots);#600; 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);#37;
                minvar = np.nanmin(dataall_subplots);#32; 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots);#30;
                minvar = np.nanmin(dataall_subplots);#22; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);#8.5;
                minvar = np.nanmin(dataall_subplots);#7.5; 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);#3000;
                minvar = np.nanmin(dataall_subplots);#200; 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);#0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);#0; 
            else:
                pass
                
        elif region == "oceansoda_congo":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "AT":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots)-273.15;
                minvar = np.nanmin(dataall_subplots)-273.15; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);#
                minvar = np.nanmin(dataall_subplots);# 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);
                minvar = np.nanmin(dataall_subplots); 
            else:
                pass
        elif region == "oceansoda_mediterranean":
            if plot_var == "DIC":
                maxvar = np.nanmax(dataall_subplots);#2300;
                minvar = np.nanmin(dataall_subplots);# 1950; 
            elif plot_var == "SSS":
                maxvar = np.nanmax(dataall_subplots);#37;
                minvar = np.nanmin(dataall_subplots);# 32; 
            elif plot_var == "SST":
                maxvar = np.nanmax(dataall_subplots);#25;
                minvar = np.nanmin(dataall_subplots);# 10; 
            elif plot_var == "pH_free": 
                maxvar = np.nanmax(dataall_subplots);#8.2;
                minvar = np.nanmin(dataall_subplots);# 8; 
            elif plot_var == "pCO2":
                maxvar = np.nanmax(dataall_subplots);#600;
                minvar = np.nanmin(dataall_subplots);# 200; 
            elif plot_var == "saturation_aragonite":
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);# 0; 
            elif plot_var == "saturation_calcite":                 
                maxvar = np.nanmax(dataall_subplots);#10;
                minvar = np.nanmin(dataall_subplots);# 0; 
            else:
                pass
        else:
            pass
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': 16};
        gl.ylabel_style = {'size': 16};
        f3_ax.text(-0.12, 0.55, 'Latitude $(^{\circ})$', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',fontsize=16,
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.13, 'Longitude $(^{\circ})$', va='bottom', ha='center',
            rotation='horizontal', rotation_mode='anchor',fontsize=16,
            transform=f3_ax.transAxes)
        
        if region == "oceansoda_amazon_plume":
            f3_ax.set_xlim(-74, -30);
            f3_ax.set_ylim(-4, 26);
        elif region == "oceansoda_congo":
            f3_ax.set_xlim(-3, 18);
            f3_ax.set_ylim(-11, 5);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_xlim(-6, 40);
            f3_ax.set_ylim(28, 48);
        else:
            pass
        
        f3_ax.coastlines();
        resol = '50m'  # use data at this scale
        land = cfeature.NaturalEarthFeature('physical', 'land', \
        scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
        f3_ax.add_feature(land, facecolor='beige')

        m = plt.cm.ScalarMappable(cmap=cmap_touse)
        m.set_array(plot_var)
        m.set_clim(minvar, maxvar)
        cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))  
        cb.ax.tick_params(labelsize=16)

        if plot_var == "DIC":
            cb.ax.set_title("DIC ($\mu mol \ kg^{-1}$)",fontsize=16);
        elif plot_var == "AT":
            cb.ax.set_title("TA ($\mu mol \ kg^{-1}$)",fontsize=16);
        elif plot_var == "SSS":
            cb.ax.set_title("Salinity (PSU)",fontsize=16);
        elif plot_var == "SST":
            cb.ax.set_title("SST (($^\circ$ C))",fontsize=16);
        elif plot_var == "pH_free":
            cb.ax.set_title("pH",fontsize=16);
        elif plot_var == "pCO2":
            cb.ax.set_title("pCO$_{2}$ (ppm)",fontsize=16);
        elif plot_var == "saturation_aragonite":
            cb.ax.set_title("$\Omega$ Aragonite",fontsize=16);
        elif plot_var == "saturation_calcite":
            cb.ax.set_title("$\Omega$ Calcite",fontsize=16);
        else:
            pass
        
        # #add a marker for the mouth of the Amazon
        # Mark some particular places with a small circle and a name label...
        # Define some test points with latitude and longitude coordinates.
                # #add a marker for the mouth of the Amazon
        if region == "oceansoda_amazon_plume":
            # Mark some particular places with a small circle and a name label...
            # Define some test points with latitude and longitude coordinates.
            city_data = [('Amazon River mouth', -1.455833,-48.503889 ), 
                          ('Orinoco River mouth', 8.616667, -62.250000),
                          ('Maroni River mouth', 5.745833, -53.968333)]
                            #('Essequibo River mouth', 7.033333, -58.450000)
            crs_latlon = ccrs.PlateCarree()
            # Place a single marker point and a text annotation at each place.
            for name, lat, lon in city_data:
            
                plt.plot(lon, lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                          markerfacecolor='black', markeredgecolor='white',transform=crs_latlon)
                # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
                # so for this one we transform the coordinates with a Cartopy call.
                at_x, at_y = f3_ax.projection.transform_point(lon, lat,src_crs=crs_latlon)
                if name =="Maroni River mouth":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-150, -40), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Orinoco River mouth":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-100, -30), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                else:
                # place_labels=plt.text(lon, lat,name,fontsize=15,weight='bold')
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-200, -10), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
        
        elif region == "oceansoda_congo":
            # Mark some particular places with a small circle and a name label...
            # Define some test points with latitude and longitude coordinates.
            city_data = [('Congo River', -6.060116, 12.494964 ), 
                          ('Niger Delta', 4.30996655147119, 6.22693906885239),
                          ('Ogooué River', -1.027200, 8.884800),
                          ('Sanaga River', 3.559338, 9.652175)]
            crs_latlon = ccrs.PlateCarree()
            # Place a single marker point and a text annotation at each place.
            for name, lat, lon in city_data:
            
                plt.plot(lon, lat, marker='o', markersize=7.0, markeredgewidth=2.5,
                          markerfacecolor='black', markeredgecolor='white',transform=crs_latlon)
                # NOTE: the "plt.annotate call" does not have a "transform=" keyword,
                # so for this one we transform the coordinates with a Cartopy call.
                at_x, at_y = f3_ax.projection.transform_point(lon, lat,src_crs=crs_latlon)
                if name =="Congo River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(-20, 50), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Niger Delta":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(130, -5), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))
                elif name =="Ogooué River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(50, 0), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5))   
                elif name =="Sanaga River":
                    plt.annotate(name, xy=(at_x, at_y), xytext=(50, -40), textcoords='offset points',
                        color='black', backgroundcolor='white', size='large',
                        arrowprops=dict(arrowstyle='->', color='white', linewidth=2.5)) 
                    
         
        # # add the slice as a line
        if region == "oceansoda_amazon_plume":
            ad_lat, ad_lon = 4, -52
            liv_lat, liv_lon = 23, -52
            x1, y1 = f3_ax.projection.transform_point(ad_lon, ad_lat,src_crs=crs_latlon)
            x2, y2 = f3_ax.projection.transform_point(liv_lon, liv_lat,src_crs=crs_latlon)
    
            f3_ax.plot([x1, x2], [y1, y2],
                      color='black', linewidth=1, marker='s', markersize=3, markerfacecolor='black',linestyle='--',transform=crs_latlon)
        elif region == "oceansoda_congo":
            ad_lat, ad_lon = -6, -1.5
            liv_lat, liv_lon = -6, 11.5
            x1, y1 = f3_ax.projection.transform_point(ad_lon, ad_lat,src_crs=crs_latlon)
            x2, y2 = f3_ax.projection.transform_point(liv_lon, liv_lat,src_crs=crs_latlon)
    
            f3_ax.plot([x1, x2], [y1, y2],
                      color='black', linewidth=1, marker='s', markersize=3, markerfacecolor='black',linestyle='--',transform=crs_latlon)
           
        
        if counter==0:
            f3_ax.set_title("JFM",fontsize=16);
        elif counter==1:
            f3_ax.set_title("AMJ",fontsize=16);
        elif counter==2:
            f3_ax.set_title("JAS",fontsize=16);
        elif counter==3:
            f3_ax.set_title("OND",fontsize=16);
            
        if counter==0:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(a)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==1:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(b)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==2:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(c)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
        elif counter==3:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(d)', transform=f3_ax.transAxes, fontsize=16, weight="bold");
            
        counter=counter+1;
             
    #check that folder exists to save images, otherwise create it
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    fig6.savefig(path.join("os_plots\{0}_seasonal".format(region), "{0}_monthly_averages.png".format(plot_var)));
    fig6.tight_layout()
    plt.close();
    del varALL;    
    
    
    #### Figure 1 RMSDe algo comparison
    
import pickle;
import pandas as pd
import numpy as np

        #### Load in Amazon TA pickle
with open("output/algo_metrics/Amazon_AT_finalscores_pickled.txt", "rb") as myFile:
    Amazon_at_dict = pickle.load(myFile)

#get a full list of the combinations
list_name=[];
for key in Amazon_at_dict:
    list_name.append(key)
    print(key)
    
#loop through combinations and get a full list of the unique algos 
#in each list and then append them
full_algos_list=pd.Series();

for loop in np.arange(0, len(Amazon_at_dict), 1):
    algos=Amazon_at_dict[list_name[loop]]['algorithm']
    full_algos_list=full_algos_list.append(algos)
    
unqiuealgos_at_amazon=full_algos_list.unique();


RMSD_matrix_AT_AMAZON = np.empty((len(unqiuealgos_at_amazon),len(Amazon_at_dict),))
RMSD_matrix_AT_AMAZON[:] = np.nan
for loop2 in np.arange(0, len(unqiuealgos_at_amazon), 1): # for each algorithm
    for loop3 in np.arange(0, len(Amazon_at_dict), 1): #loop through each combination 
        A=unqiuealgos_at_amazon[loop2]
        B=Amazon_at_dict[list_name[loop3]]['algorithm']
        B_unique_sorted, B_idx = np.unique(B, return_index=True)
        B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
        C=B_unique_sorted[B_in_A_bool]
        index=np.where(B_in_A_bool)[0]
        
        if not C: #if algorithm is not in that combination - skip
            print("Algo not run for that combination")
        else: #if algorithm is in that combination get the rmsd
            print("Algo run for this combination")
            RMSD_matrix_AT_AMAZON[loop2,loop3]=np.NAN;
            addmatrix=Amazon_at_dict[list_name[loop3]]['final_wrmsd'][index]
            if isinstance(addmatrix, np.float64):
                print("do nothing") 
                RMSD_matrix_AT_AMAZON[loop2,loop3]=np.NAN;
            else:
                RMSD_matrix_AT_AMAZON[loop2,loop3]=addmatrix;


unqiuealgos_at_amazon2=("(Astor et al. 2017)","(Astor et al. 2017)","(Astor et al. 2017)","(Astor et al. 2017)","(Astor et al. 2017)","(Brewer et al. 1995)","(Cai et al. 2010)","(Cooley et al. 2006)","(Goyet et al. 1998)","(Lee et al. 2006)","(Lefèvre et al. 2010)","(Millero et al. 1998)","(Sasse et al. 2013)","(Sasse et al. 2013)","(Ternon et al. 2000)","(Takahashi et al. 2014)")
                        


#### Load in Amazon DIC pickle 

with open("output/algo_metrics/Amazon_DIC_finalscores_pickled.txt", "rb") as myFile:
    Amazon_DIC_dict = pickle.load(myFile)

#get a full list of the combinations
list_name=[];
for key in Amazon_DIC_dict:
    list_name.append(key)
    print(key)
    
#loop through combinations and get a full list of the unique algos 
#in each list and then append them
full_algos_list=pd.Series();

for loop in np.arange(0, len(Amazon_DIC_dict), 1):
    algos=Amazon_DIC_dict[list_name[loop]]['algorithm']
    full_algos_list=full_algos_list.append(algos)
    
unqiuealgos_dic_amazon=full_algos_list.unique();


RMSD_matrix_DIC_AMAZON = np.empty((len(unqiuealgos_dic_amazon),len(Amazon_DIC_dict),))
RMSD_matrix_DIC_AMAZON[:] = np.nan
for loop2 in np.arange(0, len(unqiuealgos_dic_amazon), 1): # for each algorithm
    for loop3 in np.arange(0, len(Amazon_DIC_dict), 1): #loop through each combination 
        A=unqiuealgos_dic_amazon[loop2]
        B=Amazon_DIC_dict[list_name[loop3]]['algorithm']
        B_unique_sorted, B_idx = np.unique(B, return_index=True)
        B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
        C=B_unique_sorted[B_in_A_bool]
        index=np.where(B_in_A_bool)[0]
        
        if not C: #if algorithm is not in that combination - skip
            print("Algo not run for that combination")
        else: #if algorithm is in that combination get the rmsd
            print("Algo run for this combination")
            RMSD_matrix_DIC_AMAZON[loop2,loop3]=np.NAN;
            addmatrix=Amazon_DIC_dict[list_name[loop3]]['final_wrmsd'][index]
            if isinstance(addmatrix, np.float64):
                print("do nothing") 
                RMSD_matrix_DIC_AMAZON[loop2,loop3]=np.NAN;
            else:
                RMSD_matrix_DIC_AMAZON[loop2,loop3]=addmatrix;


        #### Load in Congo TA pickle 
with open("output/algo_metrics/Congo_AT_finalscores_pickled.txt", "rb") as myFile:
    Congo_at_dict = pickle.load(myFile)

unqiuealgos_dic_amazon2=("(Cooley et al. 2006)","(Lee et al. 2000)","(Lefèvre et al. 2010)","(Lefèvre et al. 2017)","(Ternon et al. 2000)","(Brewer et al. 1995)","(Sasse et al. 2013)","(Sasse et al. 2013)")


#get a full list of the combinations
list_name_Congo_AT=[];
for key in Congo_at_dict:
    list_name_Congo_AT.append(key)
    print(key)
    
#loop through combinations and get a full list of the unique algos 
#in each list and then append them
full_algos_list=pd.Series();

for loop in np.arange(0, len(Congo_at_dict), 1):
    algos=Congo_at_dict[list_name_Congo_AT[loop]]['algorithm']
    full_algos_list=full_algos_list.append(algos)
    
unqiuealgos_at_Congo=full_algos_list.unique();


RMSD_matrix_AT_Congo = np.empty((len(unqiuealgos_at_Congo),len(Congo_at_dict),))
RMSD_matrix_AT_Congo[:] = np.nan
for loop2 in np.arange(0, len(unqiuealgos_at_Congo), 1): # for each algorithm
    for loop3 in np.arange(0, len(Congo_at_dict), 1): #loop through each combination 
        A=unqiuealgos_at_Congo[loop2]
        B=Congo_at_dict[list_name_Congo_AT[loop3]]['algorithm']
        B_unique_sorted, B_idx = np.unique(B, return_index=True)
        B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
        C=B_unique_sorted[B_in_A_bool]
        index=np.where(B_in_A_bool)[0]
        
        if not C: #if algorithm is not in that combination - skip
            print("Algo not run for that combination")
            RMSD_matrix_AT_Congo[loop2,loop3]=np.NAN;
        else: #if algorithm is in that combination get the rmsd
            print("Algo run for this combination")
            RMSD_matrix_AT_Congo[loop2,loop3]=np.NAN;
            addmatrix=Congo_at_dict[list_name_Congo_AT[loop3]]['final_wrmsd'][index]
            if isinstance(addmatrix, np.float64):
                print("do nothing") 
                RMSD_matrix_AT_Congo[loop2,loop3]=np.NAN;
            else:
                RMSD_matrix_AT_Congo[loop2,loop3]=addmatrix;


        #### Load in Congo DIC pickle 

with open("output/algo_metrics/CONGO_DIC_finalscores_pickled.txt", "rb") as myFile:
    CONGO_DIC_dict = pickle.load(myFile)
    
unqiuealgos_at_Congo2=("(Goyet et al. 1998)","(Lee et al. 2006)","(Takahashi et al. 2014)");



#get a full list of the combinations
list_name_Congo_DIC=[];
for key in CONGO_DIC_dict:
    list_name_Congo_DIC.append(key)
    print(key)
    
#loop through combinations and get a full list of the unique algos 
#in each list and then append them
full_algos_list=pd.Series();

for loop in np.arange(0, len(CONGO_DIC_dict), 1):
    algos=CONGO_DIC_dict[list_name_Congo_DIC[loop]]['algorithm']
    full_algos_list=full_algos_list.append(algos)
    
unqiuealgos_dic_CONGO=full_algos_list.unique();


RMSD_matrix_DIC_CONGO = np.empty((len(unqiuealgos_dic_CONGO),len(CONGO_DIC_dict),))
RMSD_matrix_DIC_CONGO[:] = np.nan
for loop2 in np.arange(0, len(unqiuealgos_dic_CONGO), 1): # for each algorithm
    for loop3 in np.arange(0, len(CONGO_DIC_dict), 1): #loop through each combination 
        A=unqiuealgos_dic_CONGO[loop2]
        B=CONGO_DIC_dict[list_name_Congo_DIC[loop3]]['algorithm']
        B_unique_sorted, B_idx = np.unique(B, return_index=True)
        B_in_A_bool = np.in1d(B_unique_sorted, A, assume_unique=True)
        C=B_unique_sorted[B_in_A_bool]
        index=np.where(B_in_A_bool)[0]
        
        if not C: #if algorithm is not in that combination - skip
            print("Algo not run for that combination")
        else: #if algorithm is in that combination get the rmsd
            print("Algo run for this combination")
            RMSD_matrix_DIC_CONGO[loop2,loop3]=np.NAN;
            addmatrix=CONGO_DIC_dict[list_name_Congo_DIC[loop3]]['final_wrmsd'][index]
            if isinstance(addmatrix, np.float64):
                print("do nothing") 
                RMSD_matrix_DIC_CONGO[loop2,loop3]=np.NAN;
            else:
                RMSD_matrix_DIC_CONGO[loop2,loop3]=addmatrix;

unqiuealgos_dic_CONGO2=("(Bakker et al. 1999)","(Bakker et al. 1999)","(Bakker et al. 1999)","(Lee et al. 2000)","(Vangriesheim et al. 2009)","(Vangriesheim et al. 2009)","(Brewer et al. 1995)")   


#### now make the figure
figsize = (16,24);
fig1 = plt.figure(figsize=figsize);
gs = fig1.add_gridspec(46, 50)

f1_ax1 = fig1.add_subplot(gs[0:15, 5:35])
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,0]), 1),RMSD_matrix_AT_AMAZON[:,0],"r|",markersize=25, label=list_name[0]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,1]), 1),RMSD_matrix_AT_AMAZON[:,1],"b|", markersize=25,label=list_name[1]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,2]), 1),RMSD_matrix_AT_AMAZON[:,2],"k|", markersize=25,label=list_name[2]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,3]), 1),RMSD_matrix_AT_AMAZON[:,3],"rx",markersize=25, label=list_name[3]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,4]), 1),RMSD_matrix_AT_AMAZON[:,4],"bx",markersize=25, label=list_name[4]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,5]), 1),RMSD_matrix_AT_AMAZON[:,5],"kx", markersize=25,label=list_name[5]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,6]), 1),RMSD_matrix_AT_AMAZON[:,6],"r^", markersize=25, label=list_name[6]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,7]), 1),RMSD_matrix_AT_AMAZON[:,7],"b^", markersize=25,label=list_name[7]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,8]), 1),RMSD_matrix_AT_AMAZON[:,8],"k^", markersize=25,label=list_name[8]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,9]), 1),RMSD_matrix_AT_AMAZON[:,9],"r_", markersize=25,label=list_name[9]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,10]), 1),RMSD_matrix_AT_AMAZON[:,10],"b_",markersize=25, label=list_name[10]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,11]), 1),RMSD_matrix_AT_AMAZON[:,11],"k_",markersize=25, label=list_name[11]);
f1_ax1.set_xticks(np.arange(0, len(unqiuealgos_at_amazon), 1))
f1_ax1.set_xticklabels(list(unqiuealgos_at_amazon), rotation=90,fontsize=24)
f1_ax1.set_ylim(0, 100);
f1_ax1.set_ylabel("RMSDe",fontsize=24);
f1_ax1.set_title('Amazon TA',fontsize=24);
f1_ax1.legend(loc="upper left",ncol=1);
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);

f1_ax2 = fig1.add_subplot(gs[23:38, 5:35])
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,0]), 1),RMSD_matrix_DIC_AMAZON[:,0],"r|",markersize=25, label=list_name[0]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,1]), 1),RMSD_matrix_DIC_AMAZON[:,1],"b|", markersize=25,label=list_name[1]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,2]), 1),RMSD_matrix_DIC_AMAZON[:,2],"k|", markersize=25,label=list_name[2]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,3]), 1),RMSD_matrix_DIC_AMAZON[:,3],"rx",markersize=25, label=list_name[3]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,4]), 1),RMSD_matrix_DIC_AMAZON[:,4],"bx",markersize=25, label=list_name[4]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,5]), 1),RMSD_matrix_DIC_AMAZON[:,5],"kx", markersize=25,label=list_name[5]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,6]), 1),RMSD_matrix_DIC_AMAZON[:,6],"r^", markersize=25, label=list_name[6]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,7]), 1),RMSD_matrix_DIC_AMAZON[:,7],"b^", markersize=25,label=list_name[7]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,8]), 1),RMSD_matrix_DIC_AMAZON[:,8],"k^", markersize=25,label=list_name[8]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,9]), 1),RMSD_matrix_DIC_AMAZON[:,9],"r_", markersize=25,label=list_name[9]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,10]), 1),RMSD_matrix_DIC_AMAZON[:,10],"b_",markersize=25, label=list_name[10]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,11]), 1),RMSD_matrix_DIC_AMAZON[:,11],"k_",markersize=25, label=list_name[11]);
f1_ax2.set_xticks(np.arange(0, len(unqiuealgos_dic_amazon), 1))
f1_ax2.set_xticklabels(list(unqiuealgos_dic_amazon), rotation=90,fontsize=24)
f1_ax2.set_ylim(0, 100);
f1_ax2.set_ylabel("RMSDe",fontsize=24);
f1_ax2.set_title('Amazon DIC',fontsize=24);
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);

f1_ax3 = fig1.add_subplot(gs[0:15, 40:50])
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,0]), 1),RMSD_matrix_AT_Congo[:,0],"r|",markersize=25, label=list_name_Congo_AT[0]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,1]), 1),RMSD_matrix_AT_Congo[:,1],"b|", markersize=25,label=list_name_Congo_AT[1]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,2]), 1),RMSD_matrix_AT_Congo[:,2],"k|", markersize=25,label=list_name_Congo_AT[2]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,3]), 1),RMSD_matrix_AT_Congo[:,3],"rx",markersize=25, label=list_name_Congo_AT[3]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,4]), 1),RMSD_matrix_AT_Congo[:,4],"bx",markersize=25, label=list_name_Congo_AT[4]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,5]), 1),RMSD_matrix_AT_Congo[:,5],"kx", markersize=25,label=list_name_Congo_AT[5]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,6]), 1),RMSD_matrix_AT_Congo[:,6],"r_", markersize=25,label=list_name_Congo_AT[6]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,7]), 1),RMSD_matrix_AT_Congo[:,7],"b_", markersize=25,label=list_name_Congo_AT[7]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,8]), 1),RMSD_matrix_AT_Congo[:,8],"k_", markersize=25,label=list_name_Congo_AT[8]);
f1_ax3.set_xticks(np.arange(0, len(unqiuealgos_at_Congo), 1))
f1_ax3.set_xticklabels(list(unqiuealgos_at_Congo), rotation=90,fontsize=24)
f1_ax3.set_ylim(0, 100);
f1_ax3.set_title('Congo TA',fontsize=24);
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);

f1_ax4 = fig1.add_subplot(gs[23:38, 40:50])
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,0]), 1),RMSD_matrix_DIC_CONGO[:,0],"r|",markersize=25, label=list_name_Congo_AT[0]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,1]), 1),RMSD_matrix_DIC_CONGO[:,1],"b|", markersize=25,label=list_name_Congo_AT[1]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,2]), 1),RMSD_matrix_DIC_CONGO[:,2],"k|", markersize=25,label=list_name_Congo_AT[2]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,3]), 1),RMSD_matrix_DIC_CONGO[:,3],"rx",markersize=25, label=list_name_Congo_AT[3]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,4]), 1),RMSD_matrix_DIC_CONGO[:,4],"bx",markersize=25, label=list_name_Congo_AT[4]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,5]), 1),RMSD_matrix_DIC_CONGO[:,5],"kx", markersize=25,label=list_name_Congo_AT[5]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,6]), 1),RMSD_matrix_DIC_CONGO[:,6],"r_", markersize=25,label=list_name_Congo_AT[6]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,7]), 1),RMSD_matrix_DIC_CONGO[:,7],"b_", markersize=25,label=list_name_Congo_AT[7]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,8]), 1),RMSD_matrix_DIC_CONGO[:,8],"k_", markersize=25,label=list_name_Congo_AT[8]);
f1_ax4.set_xticks(np.arange(0, len(unqiuealgos_dic_CONGO), 1))
f1_ax4.set_xticklabels(list(unqiuealgos_dic_CONGO), rotation=90,fontsize=24)
f1_ax4.set_ylim(0, 100);
f1_ax4.set_title('Congo DIC',fontsize=24);
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);

fig1.savefig("os_plots\\Fig1_RMSDe.png");





list_name_short=("ESACCI_SSS-ESACCI_SST","ESACCI_SSS-CORA","ESACCI_SSS-OISST","CORA-ESACCI_SST","CORA-CORA","CORA-OISST","RSSSMAP-ESACCI_SST","RSSMAP-CORA","RSSSMAP-OISST","ISAS-ESACCI_SST","ISAS-CORA","ISAS-OISST")

color_teal = [18/255,150/255,155/255];
color_peach = [251/255,111/255,66/255];
color_firebreak = [178/255,34/255,34/255];

#### now make the figure
figsize = (24*1.2,26*1.2);
fig1 = plt.figure(figsize=figsize);
gs = fig1.add_gridspec(45, 50)

f1_ax1 = fig1.add_subplot(gs[0:16, 10:50])
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,0],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,0]), 1),"|", color =color_firebreak ,markersize=25,mew=10, label=list_name_short[0]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,1],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,1]), 1),"|", color =color_peach, markersize=25,mew=10,label=list_name_short[1]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,2],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,2]), 1),"|", color =color_teal, markersize=25,mew=10,label=list_name_short[2]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,3],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,3]), 1),"x", color =color_firebreak ,markersize=25,mew=6, label=list_name_short[3]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,4],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,4]), 1),"x", color =color_peach,markersize=25,mew=6, label=list_name_short[4]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,5],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,5]), 1),"x", color =color_teal, markersize=25,mew=6,label=list_name_short[5]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,6],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,6]), 1),"^", color =color_firebreak , markersize=25,mew=3, label=list_name_short[6]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,7],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,7]), 1),"^", color =color_peach, markersize=25,mew=3,label=list_name_short[7]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,8],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,8]), 1),"^", color =color_teal, markersize=25,mew=3,label=list_name_short[8]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,9],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,9]), 1),"_", color =color_firebreak , markersize=25,mew=10,label=list_name_short[9]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,10],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,10]), 1),"_", color =color_peach,markersize=25,mew=10, label=list_name_short[10]);
f1_ax1.plot(RMSD_matrix_AT_AMAZON[:,11],np.arange(0, len(RMSD_matrix_AT_AMAZON[:,11]), 1),"_", color =color_teal,markersize=25,mew=10, label=list_name_short[11]);
f1_ax1.set_yticks(np.arange(0, len(unqiuealgos_at_amazon), 1))
f1_ax1.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
f1_ax1.set_yticklabels(list(unqiuealgos_at_amazon2), fontsize=24)
f1_ax1.set_xlim(0, 100);
f1_ax1.set_xlabel("RMSDe ($\mu mol \ kg^{-1}$)",fontsize=24);
f1_ax1.set_title('Amazon TA',fontsize=24, weight="bold");
#f1_ax1.legend(loc="upper left",ncol=2);
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);
plt.grid(linestyle="-")
yyy=f1_ax1.get_ylim()
plt.ylim(yyy[0]-0.5, yyy[1]+0.5)
f1_ax1.text(-0.2, 1.06, '(a)', transform=f1_ax1.transAxes, fontsize=24, weight="bold");

f1_ax3 = fig1.add_subplot(gs[20:23, 10:50])
f1_ax3.plot(RMSD_matrix_AT_Congo[:,0],np.arange(0, len(RMSD_matrix_AT_Congo[:,0]), 1),"|", color =color_firebreak ,markersize=25,mew=10, label=list_name_Congo_AT[0]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,1],np.arange(0, len(RMSD_matrix_AT_Congo[:,1]), 1),"|", color =color_peach, markersize=25,mew=10,label=list_name_Congo_AT[1]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,2],np.arange(0, len(RMSD_matrix_AT_Congo[:,2]), 1),"|", color =color_teal, markersize=25,mew=10,label=list_name_Congo_AT[2]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,3],np.arange(0, len(RMSD_matrix_AT_Congo[:,3]), 1),"x", color =color_firebreak ,markersize=25,mew=6, label=list_name_Congo_AT[3]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,4],np.arange(0, len(RMSD_matrix_AT_Congo[:,4]), 1),"x", color =color_peach,markersize=25,mew=6, label=list_name_Congo_AT[4]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,5],np.arange(0, len(RMSD_matrix_AT_Congo[:,5]), 1),"x", color =color_teal, markersize=25,mew=10,label=list_name_Congo_AT[5]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,6],np.arange(0, len(RMSD_matrix_AT_Congo[:,6]), 1),"_", color =color_firebreak , markersize=25,mew=10,label=list_name_Congo_AT[6]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,7],np.arange(0, len(RMSD_matrix_AT_Congo[:,7]), 1),"_", color =color_peach, markersize=25,mew=10,label=list_name_Congo_AT[7]);
f1_ax3.plot(RMSD_matrix_AT_Congo[:,8],np.arange(0, len(RMSD_matrix_AT_Congo[:,8]), 1),"_", color =color_teal, markersize=25,mew=10,label=list_name_Congo_AT[8]);
f1_ax3.set_yticks(np.arange(0, len(unqiuealgos_at_Congo2), 1))
f1_ax3.set_yticklabels(list(unqiuealgos_at_Congo2),fontsize=24)
f1_ax3.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
f1_ax3.set_xlim(0, 100);
f1_ax3.set_title('Congo TA',fontsize=24, weight="bold");
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);
plt.grid(linestyle="-")
yyy=f1_ax3.get_ylim()
plt.ylim(yyy[0]-0.5, yyy[1]+0.5)
f1_ax3.set_xlabel("RMSDe ($\mu mol \ kg^{-1}$)",fontsize=24);
f1_ax3.text(-0.2, 1.06, '(b)', transform=f1_ax3.transAxes, fontsize=24, weight="bold");


f1_ax2 = fig1.add_subplot(gs[27:35, 10:50])
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,0],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,0]), 1),"|", color =color_firebreak ,markersize=25,mew=10, label=list_name[0]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,1],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,1]), 1),"|", color =color_peach, markersize=25,mew=10,label=list_name[1]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,2],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,2]), 1),"|", color =color_teal, markersize=25,mew=10,label=list_name[2]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,3],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,3]), 1),"x", color =color_firebreak ,markersize=25,mew=6, label=list_name[3]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,4],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,4]), 1),"x", color =color_peach,markersize=25,mew=6, label=list_name[4]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,5],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,5]), 1),"x", color =color_teal, markersize=25,mew=6,label=list_name[5]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,6],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,6]), 1),"^", color =color_firebreak , markersize=25,mew=3, label=list_name[6]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,7],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,7]), 1),"^", color =color_peach, markersize=25,mew=3,label=list_name[7]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,8],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,8]), 1),"^", color =color_teal, markersize=25,mew=3,label=list_name[8]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,9],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,9]), 1),"_", color =color_firebreak , markersize=25,mew=10,label=list_name[9]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,10],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,10]), 1),"_", color =color_peach,markersize=25,mew=10, label=list_name[10]);
f1_ax2.plot(RMSD_matrix_DIC_AMAZON[:,11],np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,11]), 1),"_", color =color_teal,markersize=25,mew=10, label=list_name[11]);
f1_ax2.set_yticks(np.arange(0, len(unqiuealgos_dic_amazon), 1))
f1_ax2.set_yticklabels(list(unqiuealgos_dic_amazon2),fontsize=24)
f1_ax2.set_xlim(0, 100);
f1_ax2.set_xlabel("RMSDe ($\mu mol \ kg^{-1}$)",fontsize=24);
f1_ax2.set_title('Amazon DIC',fontsize=24, weight="bold");
f1_ax2.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);
plt.grid(linestyle="-")
yyy=f1_ax2.get_ylim()
plt.ylim(yyy[0]-0.5, yyy[1]+0.5)
f1_ax2.text(-0.2, 1.06, '(c)', transform=f1_ax2.transAxes, fontsize=24, weight="bold");


f1_ax4 = fig1.add_subplot(gs[39:45, 10:50])
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,0],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,0]), 1),"|", color =color_firebreak ,markersize=25,mew=10, label=list_name_Congo_AT[0]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,1],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,1]), 1),"|", color =color_peach, markersize=25,mew=10,label=list_name_Congo_AT[1]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,2],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,2]), 1),"|", color =color_teal, markersize=25,mew=10,label=list_name_Congo_AT[2]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,3],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,3]), 1),"x", color =color_firebreak ,markersize=25,mew=6, label=list_name_Congo_AT[3]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,4],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,4]), 1),"x", color =color_peach,markersize=25,mew=6, label=list_name_Congo_AT[4]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,5],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,5]), 1),"x", color =color_teal, markersize=25,mew=6,label=list_name_Congo_AT[5]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,6],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,6]), 1),"_", color =color_firebreak , markersize=25,mew=10,label=list_name_Congo_AT[6]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,7],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,7]), 1),"_", color =color_peach, markersize=25,mew=10,label=list_name_Congo_AT[7]);
f1_ax4.plot(RMSD_matrix_DIC_CONGO[:,8],np.arange(0, len(RMSD_matrix_DIC_CONGO[:,8]), 1),"_", color =color_teal,markersize=25,mew=10,label=list_name_Congo_AT[8]);
f1_ax4.set_yticks(np.arange(0, len(unqiuealgos_dic_CONGO), 1))
f1_ax4.set_yticklabels(list(unqiuealgos_dic_CONGO2), fontsize=24)
f1_ax4.set_xlim(0, 100);
f1_ax4.set_title('Congo DIC',fontsize=24, weight="bold");
f1_ax4.set_xticks([0,10,20,30,40,50,60,70,80,90,100])
plt.xticks(fontsize=24);
plt.yticks(fontsize=24);
plt.grid(linestyle="--")
yyy=f1_ax4.get_ylim()
plt.ylim(yyy[0]-0.5, yyy[1]+0.5)
f1_ax4.set_xlabel("RMSDe ($\mu mol \ kg^{-1}$)",fontsize=24);
f1_ax4.text(-0.2, 1.06, '(d)', transform=f1_ax4.transAxes, fontsize=24, weight="bold");


#legend plotting

f1_ax1.legend(loc="center",bbox_to_anchor=(0.4, 1.25), fontsize=24,ncol=4);

fig1.savefig("os_plots\\Fig1_RMSDe.png");



import pickle 
#### Load in Amazon RADII pickle
with open("output/Amazon_radii.txt", "rb") as myFile:
    Radii_masks_dict = pickle.load(myFile)
    
atAmazonNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_amazon_plume", VAR="AT"));
lats = atAmazonNC.variables["lat"][:];
lons = atAmazonNC.variables["lon"][:];  
  
import matplotlib.pyplot as plt        

fig3, ax = plt.subplots()
f3_ax = plt.axes(projection=ccrs.PlateCarree())
gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
gl.xlabels_top = False;
gl.ylabels_right = False;
gl.xformatter = LONGITUDE_FORMATTER;
gl.yformatter = LATITUDE_FORMATTER;
gl.xlabel_style = {'size': 24};
gl.ylabel_style = {'size': 24};
f3_ax.text(-0.2, 0.55, 'Latitude $(^{\circ})$', va='bottom', ha='center',
    rotation='vertical', rotation_mode='anchor',fontsize=24,
    transform=f3_ax.transAxes)
f3_ax.text(0.5, -0.25, 'Longitude $(^{\circ})$', va='bottom', ha='center',
    rotation='horizontal', rotation_mode='anchor',fontsize=24,
    transform=f3_ax.transAxes)

   
f3_ax.set_xlim(-74, -30);
f3_ax.set_ylim(-4, 26);

f3_ax.coastlines();
resol = '50m'  # use data at this scale
land = cfeature.NaturalEarthFeature('physical', 'land', \
scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
f3_ax.add_feature(land, facecolor='beige')
for radii in Radii_masks_dict:   

    z=Radii_masks_dict[radii]
    z = np.ma.masked_array(z, z < 0.5)
    from numpy import random
    color_random=[random.rand()*255,random.rand()*255,random.rand()*255]
    
    discrete_colors = [(random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255),
    (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255), (random.rand()*255,random.rand()*255,random.rand()*255)]
    discrete_colors = [(r/255., g/255., b/255.) for r, g, b in discrete_colors]         

    my_colormap = colors.ListedColormap(discrete_colors)
    c = f3_ax.pcolor(lons,lats,z, vmin=0.5, cmap=my_colormap, vmax=1.5,)


fig3.tight_layout()
plt.show()

fig3.savefig("os_plots\\Fig10_Transect_radii.png");


#### single plots of SST and SAL for presentation
regions = ["oceansoda_amazon_plume"];#, "oceansoda_st_lawrence"];
video_vars = ["SSS","SST"];
for region in regions:
        #### Define the two netCDF files for each region
    dicNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="DIC"), 'r');
    atNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="AT"), 'r');
    time = get_datetimes(dicNC["time"][:]);
    
    
    #Maps of the variables
    for plot_var in video_vars:
        outputPath = "plots/";
        maskPath = "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data/osoda_region_masks_v2.nc";
        maskNC = Dataset(maskPath, 'r');
        
        dateToPlot3 = datetime(2015, 7, 1); #Which time slice to plot?
    
        datetoplotall=[dateToPlot3]
        
        mask = maskNC.variables[region][:];
        xmin = min(np.where(np.flipud(mask)==1)[1]);
        xmax = max(np.where(np.flipud(mask)==1)[1]);
        ymin = min(np.where(np.flipud(mask)==1)[0]);
        ymax = max(np.where(np.flipud(mask)==1)[0]);
        pad = 4; #size of xlim/ylim extents
        

        ticksize = 16;

        dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
        varALL = dicNC.variables["{0}".format(plot_var)][:];
        hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
        
        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        lats = dicNC.variables["lat"];
        lons = dicNC.variables["lon"];
        
        counter=0
        for subplot_loop_no in datetoplotall:
        
            timeIndex = get_time_index(datetoplotall[counter], firstIndexSecs = int(dicNC.variables["time"][0]));
            data = dicNC.variables[plot_var][timeIndex, :, :];
            data[mask != 1] = np.nan;
            date = dates[timeIndex]

            if plot_var == "pCO2":
                cmap_touse = plt.cm.get_cmap('inferno_r');
            else:
                cmap_touse = plt.cm.inferno;
                
            if plot_var == "SST":
                data=data-273.15;    

            # arbitary limits
            if region == "oceansoda_amazon_plume":
                if plot_var == "DIC":
                    maxvar = 2000;
                    minvar = 1400; 
                elif plot_var == "SSS":
                    maxvar = np.nanmax(data);
                    minvar = np.nanmin(data); 
                elif plot_var == "SST":
                    maxvar = np.nanmax(data);
                    minvar = np.nanmin(data); 

                else:
                    pass
            fig8, ax = plt.subplots(figsize=(14, 14))
            f8_ax = plt.axes(projection=ccrs.PlateCarree())           
            gl = f8_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
            contPlot1 = f8_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            gl.xlabels_top = False;
            gl.ylabels_right = False;
            gl.xformatter = LONGITUDE_FORMATTER;
            gl.yformatter = LATITUDE_FORMATTER;
            gl.xlabel_style = {'size': 24};
            gl.ylabel_style = {'size': 24};
            f8_ax.set_xlabel("Longitude($^\circ$)", fontsize=24);
            f8_ax.set_ylabel("Latitude($^\circ$)", fontsize=24);
            
            if region == "oceansoda_amazon_plume":
                f8_ax.set_xlim(-74, -30);
                f8_ax.set_ylim(-4, 26);
            else:
                pass
            
            f8_ax.coastlines();
            resol = '50m'  # use data at this scale
            land = cfeature.NaturalEarthFeature('physical', 'land', \
            scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land']);
            f8_ax.add_feature(land, facecolor='beige')
        
            m = plt.cm.ScalarMappable(cmap=cmap_touse)
            m.set_array(plot_var)
            m.set_clim(minvar, maxvar)
            cb = plt.colorbar(m); #, boundaries=np.linspace(0, 2, 6))   
            cb.ax.tick_params(labelsize=24)
            
            f8_ax.text(-0.05, 2.5, 'Latitude $(^{\circ})$', va='bottom', ha='center',
                rotation='vertical', rotation_mode='anchor',fontsize=24,
                transform=f3_ax.transAxes)
            f8_ax.text(1.5, 1, 'Longitude $(^{\circ})$', va='bottom', ha='center',
                rotation='horizontal', rotation_mode='anchor',fontsize=24,
                transform=f3_ax.transAxes)
            
            if plot_var == "DIC":
                cb.set_label("DIC ($\mu mol \ kg^{-1}$)", fontsize=24);
                #f8_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")), fontsize=24);  
            elif plot_var == "SSS":
                cb.set_label("Salinity (PSU)", fontsize=24);
                #f8_ax.set_title("Salinity for {0} {1}".format(date.year, format(date.month, "02d")), fontsize=24);
            elif plot_var == "SST":
                cb.set_label("SST ($^\circ$ C)", fontsize=24);
                #f8_ax.set_title("SST for {0} {1}".format(date.year, format(date.month, "02d")), fontsize=24);
            else:
                pass

            plt.show()
            fig8.savefig(path.join("os_plots\{0}_seasonal".format(region,plot_var), "{0}singlemonth.png".format(plot_var)));
                
        