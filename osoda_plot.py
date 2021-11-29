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
meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, years, annualMeanPlumeATAmazon, annualMeanNotPlumeATAmazon, annualMeanAllATAmazon = extract_var_means(atAmazonNC, "AT_pred");
dates = get_datetimes(atAmazonNC.variables["time"][:]);

# atMississippiNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mississippi", VAR="AT"));
# meanPlumeATMississippi, meanNotPlumeATMississippi, meanAllATMississippi, years, annualMeanPlumeATMississippi, annualMeanNotPlumeATMississippi, annualMeanAllATMississippi = extract_var_means(atMississippiNC, "AT_pred");

atCongoNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_congo", VAR="AT"));
meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo, years, annualMeanPlumeATCongo, annualMeanNotPlumeATCongo, annualMeanAllATCongo = extract_var_means(atCongoNC,"AT_pred");

# atStLawrenceNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_st_lawrence", VAR="AT"));
# meanPlumeATStLawrence, meanNotPlumeATStLawrence, meanAllATStLawrence, years, annualMeanPlumeATStLawrence, annualMeanNotPlumeATStLawrence, annualMeanAllATStLawrence = extract_var_means(atStLawrenceNC, "AT_pred");

atMediterraneanNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mediterranean", VAR="AT"));
meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean, years, annualMeanPlumeATMediterranean, annualMeanNotPlumeATMediterranean, annualMeanAllATMediterranean = extract_var_means(atMediterraneanNC, "AT_pred");

#min and max values for y scale:
#listAll = [meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, meanPlumeATMississippi, meanNotPlumeATMississippi, meanAllATMississippi, meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo, meanPlumeATStLawrence, meanNotPlumeATStLawrence, meanAllATStLawrence,meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean];
listAll = [meanPlumeATAmazon, meanNotPlumeATAmazon, meanAllATAmazon, meanPlumeATMediterranean, meanNotPlumeATMediterranean, meanAllATMediterranean, meanPlumeATCongo, meanNotPlumeATCongo, meanAllATCongo];

minY = np.nanmin(listAll)*0.95;
maxY = np.nanmax(listAll)*1.05;

#dic
dicAmazonNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_amazon_plume", VAR="DIC"));
meanPlumeDICAmazon, meanNotPlumeDICAmazon, meanAllDICAmazon, years, annualMeanPlumeDICAmazon, annualMeanNotPlumeDICAmazon, annualMeanAllDICAmazon= extract_var_means(dicAmazonNC, "DIC_pred");
dates = get_datetimes(dicAmazonNC.variables["time"][:]);

# dicMississippiNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mississippi", VAR="DIC"));
# meanPlumeDICMississippi, meanNotPlumeDICMississippi, meanAllDICMississippi, years, annualMeanPlumeDICMississippi, annualMeanNotPlumeDICMississippi, annualMeanAllDICMississippi = extract_var_means(dicMississippiNC, "DIC_pred");

dicCongoNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_congo", VAR="DIC"));
meanPlumeDICCongo, meanNotPlumeDICCongo, meanAllDICCongo,years, annualMeanPlumeDICCongo, annualMeanNotPlumeDICCongo, annualMeanAllDICCongo = extract_var_means(dicCongoNC, "DIC_pred");

# dicStLawrenceNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_st_lawrence", VAR="DIC"));
# meanPlumeDICStLawrence, meanNotPlumeDICStLawrence, meanAllDICStLawrence, years, annualMeanPlumeDICStLawrence, annualMeanNotPlumeDICStLawrence, annualMeanAllDICStLawrence = extract_var_means(dicStLawrenceNC, "DIC_pred");

dicMediterraneanNC = Dataset(inputTemplate.safe_substitute(REGION="oceansoda_mediterranean", VAR="DIC"));
meanPlumeDICMediterranean, meanNotPlumeDICMediterranean, meanAllDICMediterranean, years, annualMeanPlumeDICMediterranean, annualMeanNotPlumeDICMediterranean, annualMeanAllDICMediterranean = extract_var_means(dicMediterraneanNC, "DIC_pred");

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
video_vars = ["DIC_pred", "SSS","SST", "pH_free","pCO2", "saturation_aragonite", "saturation_calcite"];
create_animations="False";

for region in regions:
        #### Define the two netCDF files for each region
    dicNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="DIC"), 'r');
    atNC = Dataset(inputTemplate.safe_substitute(REGION=region, VAR="AT"), 'r');
    time = get_datetimes(dicNC["time"][:]);
    
        #### Calculate stats for all key variables
    #DIC
    meanPlumeDIC, meanNotPlumeDIC, meanAllDIC, _, _, _, _ = extract_var_means(dicNC, "DIC_pred");
    stdPlumeDIC, stdNotPlumeDIC, stdAllDIC, _, _, _, _ = extract_var_stds(dicNC, "DIC_pred");
    minPlumeDIC, minNotPlumeDIC, minAllDIC, _, _, _, _ = extract_var_mins(dicNC, "DIC_pred");
    maxPlumeDIC, maxNotPlumeDIC, maxAllDIC, _, _, _, _ = extract_var_maxs(dicNC, "DIC_pred");
    meanPlumeDICuncertainty, meanNotPlumeDICuncertainty, meanAllDICuncertainty, _, _, _, _ = extract_var_means(dicNC, "DIC_pred_combined_uncertainty");

    #AT
    meanPlumeAT, meanNotPlumeAT, meanAllAT, _, _, _, _ = extract_var_means(atNC, "AT_pred");
    stdPlumeAT, stdNotPlumeAT, stdAllAT, _, _, _, _ = extract_var_stds(atNC, "AT_pred");
    minPlumeAT, minNotPlumeAT, minAllAT, _, _, _, _ = extract_var_mins(atNC, "AT_pred");
    maxPlumeAT, maxNotPlumeAT, maxAllAT, _, _, _, _ = extract_var_maxs(atNC, "AT_pred");
    meanPlumeATuncertainty, meanNotPlumeATuncertainty, meanAllATuncertainty, _, _, _, _ = extract_var_means(atNC, "AT_pred_combined_uncertainty");
    
    #pH
    meanPlumepH_dic, meanNotPlumepH_dic, meanAllpH_dic, _, _, _, _ = extract_var_means(dicNC, "pH");
    stdPlumepH_dic, stdNotPlumepH_dic, stdAllpH_dic, _, _, _, _ = extract_var_stds(dicNC, "pH");
    minPlumepH_dic, minNotPlumepH_dic, minAllpH_dic, _, _, _, _ = extract_var_mins(dicNC, "pH");
    maxPlumepH_dic, maxNotPlumepH_dic, maxAllpH_dic, _, _, _, _ = extract_var_maxs(dicNC, "pH");
    meanPlumepHuncertainty_dic, meanNotPlumepHuncertainty_dic, meanAllpHuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "u_pH");
    
    meanPlumepH_at, meanNotPlumepH_at, meanAllpH_at, _, _, _, _ = extract_var_means(atNC, "pH");
    stdPlumepH_at, stdNotPlumepH_at, stdAllpH_at, _, _, _, _ = extract_var_stds(atNC, "pH");
    minPlumepH_at, minNotPlumepH_at, minAllpH_at, _, _, _, _ = extract_var_mins(atNC, "pH");
    maxPlumepH_at, maxNotPlumepH_at, maxAllpH_at, _, _, _, _ = extract_var_maxs(atNC, "pH");
    meanPlumepHuncertainty_at, meanNotPlumepHuncertainty_at, meanAllpHuncertainty_at, _, _, _, _ = extract_var_means(atNC, "u_pH");

    #hydrogen_free
    meanPlumehydrogen_free_dic, meanNotPlumehydrogen_free_dic, meanAllhydrogen_free_dic, _, _, _, _ = extract_var_means(dicNC, "hydrogen_free");
    stdPlumehydrogen_free_dic, stdNotPlumehydrogen_free_dic, stdAllhydrogen_free_dic, _, _, _, _ = extract_var_stds(dicNC, "hydrogen_free");
    minPlumehydrogen_free_dic, minNotPlumehydrogen_free_dic, minAllhydrogen_free_dic, _, _, _, _ = extract_var_mins(dicNC, "hydrogen_free");
    maxPlumehydrogen_free_dic, maxNotPlumehydrogen_free_dic, maxAllhydrogen_free_dic, _, _, _, _ = extract_var_maxs(dicNC, "hydrogen_free");
    meanPlumehydrogen_freeuncertainty_dic, meanNotPlumehydrogen_freeuncertainty_dic, meanAllhydrogen_freeuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "u_hydrogen_free");
    
    meanPlumehydrogen_free_at, meanNotPlumehydrogen_free_at, meanAllhydrogen_free_at, _, _, _, _ = extract_var_means(atNC, "hydrogen_free");
    stdPlumehydrogen_free_at, stdNotPlumehydrogen_free_at, stdAllhydrogen_free_at, _, _, _, _ = extract_var_stds(atNC, "hydrogen_free");
    minPlumehydrogen_free_at, minNotPlumehydrogen_free_at, minAllhydrogen_free_at, _, _, _, _ = extract_var_mins(atNC, "hydrogen_free");
    maxPlumehydrogen_free_at, maxNotPlumehydrogen_free_at, maxAllhydrogen_free_at, _, _, _, _ = extract_var_maxs(atNC, "hydrogen_free");
    meanPlumehydrogen_freeuncertainty_at, meanNotPlumehydrogen_freeuncertainty_at, meanAllhydrogen_freeuncertainty_at, _, _, _, _ = extract_var_means(atNC, "u_hydrogen_free");

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
    meanPlumepco2uncertainty_dic, meanNotPlumepco2uncertainty_dic, meanAllpco2uncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "u_pCO2");

    meanPlumepco2_at, meanNotPlumepco2_at, meanAllpco2_at, _, _, _, _ = extract_var_means(atNC, "pCO2");
    stdPlumepco2_at, stdNotPlumepco2_at, stdAllpco2_at, _, _, _, _ = extract_var_stds(atNC, "pCO2");
    minPlumepco2_at, minNotPlumepco2_at, minAllpco2_at, _, _, _, _ = extract_var_mins(atNC, "pCO2");
    maxPlumepco2_at, maxNotPlumepco2_at, maxAllpco2_at, _, _, _, _ = extract_var_maxs(atNC, "pCO2");
    meanPlumepco2uncertainty_at, meanNotPlumepco2uncertainty_at, meanAllpco2uncertainty_at, _, _, _, _ = extract_var_means(atNC, "u_pCO2");
    
    #Aragonite saturation state
    meanPlumeOmegaAragonite_dic, meanNotPlumeOmegaAragonite_dic, meanAllOmegaAragonite_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_aragonite");
    stdPlumeOmegaAragonite_dic, stdNotPlumeOmegaAragonite_dic, stdAllOmegaAragonite_dic, _, _, _, _ = extract_var_stds(dicNC, "saturation_aragonite");
    minPlumeOmegaAragonite_dic, minNotPlumeOmegaAragonite_dic, minAllOmegaAragonite_dic, _, _, _, _ = extract_var_mins(dicNC, "saturation_aragonite");
    maxPlumeOmegaAragonite_dic, maxNotPlumeOmegaAragonite_dic, maxAllOmegaAragonite_dic, _, _, _, _ = extract_var_maxs(dicNC, "saturation_aragonite");
    meanPlumeOmegaAragoniteuncertainty_dic, meanNotPlumeOmegaAragoniteuncertainty_dic, meanAllOmegaAragoniteuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "u_saturation_aragonite");

    meanPlumeOmegaAragonite_at, meanNotPlumeOmegaAragonite_at, meanAllOmegaAragonite_at, _, _, _, _ = extract_var_means(atNC, "saturation_aragonite");
    stdPlumeOmegaAragonite_at, stdNotPlumeOmegaAragonite_at, stdAllOmegaAragonite_at, _, _, _, _ = extract_var_stds(atNC, "saturation_aragonite");
    minPlumeOmegaAragonite_at, minNotPlumeOmegaAragonite_at, minAllOmegaAragonite_at, _, _, _, _ = extract_var_mins(atNC, "saturation_aragonite");
    maxPlumeOmegaAragonite_at, maxNotPlumeOmegaAragonite_at, maxAllOmegaAragonite_at, _, _, _, _ = extract_var_maxs(atNC, "saturation_aragonite");
    meanPlumeOmegaAragoniteuncertainty_at, meanNotPlumeOmegaAragoniteuncertainty_at, meanAllOmegaAragoniteuncertainty_at, _, _, _, _ = extract_var_means(atNC, "u_saturation_aragonite");

    #Calcite saturation state
    meanPlumeOmegaCalcite_dic, meanNotPlumeOmegaCalcite_dic, meanAllOmegaCalcite_dic, _, _, _, _ = extract_var_means(dicNC, "saturation_calcite");
    stdPlumeOmegaCalcite_dic, stdNotPlumeOmegaCalcite_dic, stdAllOmegaCalcite_dic, _, _, _, _ = extract_var_stds(dicNC, "saturation_calcite");
    minPlumeOmegaCalcite_dic, minNotPlumeOmegaCalcite_dic, minAllOmegaCalcite_dic, _, _, _, _ = extract_var_mins(dicNC, "saturation_calcite");
    maxPlumeOmegaCalcite_dic, maxNotPlumeOmegaCalcite_dic, maxAllOmegaCalcite_dic, _, _, _, _ = extract_var_maxs(dicNC, "saturation_calcite");
    meanPlumeOmegaCalciteuncertainty_dic, meanNotPlumeOmegaCalciteuncertainty_dic, meanAllOmegaCalciteuncertainty_dic, _, _, _, _ = extract_var_means(dicNC, "u_saturation_calcite");

    meanPlumeOmegaCalcite_at, meanNotPlumeOmegaCalcite_at, meanAllOmegaCalcite_at, _, _, _, _ = extract_var_means(atNC, "saturation_calcite");
    stdPlumeOmegaCalcite_at, stdNotPlumeOmegaCalcite_at, stdAllOmegaCalcite_at, _, _, _, _ = extract_var_stds(atNC, "saturation_calcite");
    minPlumeOmegaCalcite_at, minNotPlumeOmegaCalcite_at, minAllOmegaCalcite_at, _, _, _, _ = extract_var_mins(atNC, "saturation_calcite");
    maxPlumeOmegaCalcite_at, maxNotPlumeOmegaCalcite_at, maxAllOmegaCalcite_at, _, _, _, _ = extract_var_maxs(atNC, "saturation_calcite");
    meanPlumeOmegaCalciteuncertainty_at, meanNotPlumeOmegaCalciteuncertainty_at, meanAllOmegaCalciteuncertainty_at, _, _, _, _ = extract_var_means(atNC, "u_saturation_calcite");

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
    with open('uncert_summary_{0}.csv'.format(region), 'w') as f:
        for key in summary_dict.keys():
            f.write("%s,%s\n"%(key,summary_dict[key]))



        #### Formatted table of summary stats
    field_names = ['Timeperiod', 'mean', 'std','min','max','uncertainty']
    dp_to_round=2;
    TA = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllAT),dp_to_round), 'std': np.round(np.std(stdAllAT),dp_to_round), 'min': np.round(np.min(minAllAT),dp_to_round), 'max':np.round( np.max(maxAllAT),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllAT[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllAT[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllAT[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllAT[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllAT[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllAT[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllAT[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllATuncertainty[number_list_ond]),dp_to_round)},
    ]
    
    DIC = [
    {'Timeperiod':'Annual' , 'mean': np.round(np.mean(meanAllDIC),dp_to_round), 'std': np.round(np.std(stdAllDIC),dp_to_round), 'min': np.round(np.min(minAllDIC),dp_to_round), 'max': np.round(np.max(maxAllDIC),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty),dp_to_round)},
    {'Timeperiod':'JFM', 'mean': np.round(np.mean(meanAllDIC[number_list_jfm]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_jfm]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_jfm]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_jfm]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_jfm]),dp_to_round)},
    {'Timeperiod':'AMJ' , 'mean': np.round(np.mean(meanAllDIC[number_list_amj]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_amj]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_amj]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_amj]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_amj]),dp_to_round)},
    {'Timeperiod':'JAS', 'mean': np.round(np.mean(meanAllDIC[number_list_jas]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_jas]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_jas]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_jas]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_jas]),dp_to_round)},
    {'Timeperiod':'OND' , 'mean': np.round(np.mean(meanAllDIC[number_list_ond]),dp_to_round), 'std': np.round(np.std(stdAllDIC[number_list_ond]),dp_to_round), 'min': np.round(np.min(minAllDIC[number_list_ond]),dp_to_round), 'max': np.round(np.max(maxAllDIC[number_list_ond]),dp_to_round), 'uncertainty': np.round(np.mean(meanAllDICuncertainty[number_list_ond]),dp_to_round)},
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
    with open('Gridded_data_summary_stats_{0}.csv'.format(region), 'w',newline='') as csvfile:
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
                if video_var == "DIC_pred":
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
                if video_var == "DIC_pred":
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
                if video_var == "DIC_pred":
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
                
                if video_var == "DIC_pred":
                    cb.set_label("DIC $(\mu mol^{-1})$");
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
    
    #### Figure 4 - Timeseries plots
    
    #plotting parameters
    fontsize = 12;
    figsizex = 9.0;
    figsizey = 2.0*5;
    figsize = (figsizex, figsizey);
    
    plt.figure(figsize=figsize);
    plt.title("carbonate_timeseries_{0}".format(region))

    ax11 = plt.subplot(6,1,1);
    ax11.plot(time, meanPlumeDIC, 'r--');#, label="plume");
    ax11.plot(time, meanNotPlumeDIC, 'r:');#, label="non-plume");
    ax11.plot(time, meanAllDIC, 'r');#, label="whole region");
    ax11.set_ylabel("DIC\n$(\mu mol \ kg^{-1})$", fontsize=fontsize);

    ax21 = plt.subplot(6,1,2);
    ax21.plot(time, meanPlumeAT, 'b--', label="plume");
    ax21.plot(time, meanNotPlumeAT, 'b:', label="non-plume");
    ax21.plot(time, meanAllAT, 'b', label="whole region");
    ax21.set_ylabel("TA\n$(\mu mol \ kg^{-1})$", fontsize=fontsize);
    
    ax31 = plt.subplot(6,1,3);
    ax31.plot(time, meanPlumepH_dic, 'r--', label="plume");
    ax31.plot(time, meanNotPlumepH_dic, 'r:', label="non-plume");
    ax31.plot(time, meanAllpH_dic, 'r', label="whole region");
    
    ax31.plot(time,meanPlumepH_at , 'b--', label="plume");
    ax31.plot(time,meanNotPlumepH_at, 'b:', label="non-plume");
    ax31.plot(time,meanAllpH_at , 'b', label="whole region");
    ax31.set_ylabel("pH", fontsize=fontsize);
    
    ax41 = plt.subplot(6,1,4);
    ax41.plot(time, meanPlumepco2_dic, 'r--', label="plume");
    ax41.plot(time, meanNotPlumepco2_dic, 'r:', label="non-plume");
    ax41.plot(time,meanAllpco2_dic , 'r', label="whole region");
    
    ax41.plot(time, meanPlumepco2_at , 'b--', label="plume");
    ax41.plot(time, meanNotPlumepco2_at, 'b:', label="non-plume");
    ax41.plot(time, meanAllpco2_at, 'b', label="whole region");
    ax41.set_ylabel("pCO$_{2}$ \n (ppm)", fontsize=fontsize);
    
    ax51 = plt.subplot(6,1,5);
    ax51.plot(time,meanPlumeOmegaAragonite_dic , 'r--', label="plume");
    ax51.plot(time, meanNotPlumeOmegaAragonite_dic, 'r:', label="non-plume");
    ax51.plot(time, meanAllOmegaAragonite_dic, 'r', label="whole region");
    
    ax51.plot(time,  meanPlumeOmegaAragonite_at, 'b--', label="plume");
    ax51.plot(time, meanNotPlumeOmegaAragonite_at, 'b:', label="non-plume");
    ax51.plot(time, meanAllOmegaAragonite_at, 'b', label="whole region");
    ax51.set_ylabel("Aragonite \n saturation \n state", fontsize=fontsize);
    
    ax61 = plt.subplot(6,1,6);
    ax61.plot(time,meanPlumeOmegaCalcite_dic , 'r--', label="plume");
    ax61.plot(time, meanNotPlumeOmegaCalcite_dic, 'r:', label="non-plume");
    ax61.plot(time, meanAllOmegaCalcite_dic, 'r', label="whole region");
    
    ax61.plot(time,  meanPlumeOmegaCalcite_at, 'b--', label="plume");
    ax61.plot(time, meanNotPlumeOmegaCalcite_at, 'b:', label="non-plume");
    ax61.plot(time, meanAllOmegaCalcite_at, 'b', label="whole region");
    ax61.set_ylabel("Calcite\n saturation \n state", fontsize=fontsize);
    ax61.set_xlabel("Time (years)", fontsize=fontsize);

    #legend plotting
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k--', label="Plume");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k:', label="Non-plume");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'k', label="Whole region");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'rs', label="Using best SST & SSS inputs from DIC algorithm");
    ax11.plot([np.nan, np.nan], [np.nan, np.nan], 'bs', label="Using best SST & SSS inputs from TA algorithm");
    ax11.legend(loc="center",bbox_to_anchor=(0.5, 1.5), fontsize=fontsize,ncol=2);
    
    # plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

    plt.savefig("os_plots\\fig4_carbonate_timeseries_{0}.png".format(region));
    plt.savefig("os_plots\\fig4_carbonate_timeseries_{0}.pdf".format(region));    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
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
                if plot_var == "DIC_pred":
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
                if plot_var == "DIC_pred":
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
                if plot_var == "DIC_pred":
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
        
            if plot_var == "DIC_pred":
                cb.set_label("DIC $(\mu mol^{-1})$");
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
        hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);

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
                if plot_var == "DIC_pred":
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
                if plot_var == "DIC_pred":
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
                if plot_var == "DIC_pred":
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
        
            if plot_var == "DIC_pred":
                cb.set_label("DIC $(\mu mol^{-1})$");
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
    
    
    hovmoll_vars = ["DIC_pred", "AT_pred", "pH_free","pCO2", "saturation_aragonite", "saturation_calcite"];
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
        
        if plot_var == "AT_pred":
            dates = [get_datetime(int(timeSecs)) for timeSecs in atNC.variables["time"][:]];
            varALL = atNC.variables["{0}".format(plot_var)][:];
            hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);
        else:
            dates = [get_datetime(int(timeSecs)) for timeSecs in dicNC.variables["time"][:]];
            varALL = dicNC.variables["{0}".format(plot_var)][:];
            hasData = np.array([np.any(varALL[t,:,:].mask==False) for t in range(0, varALL.shape[0])]);


        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
        
        lats = dicNC.variables["lat"];
        lons = dicNC.variables["lon"];

        

        if region == "oceansoda_congo": #latitude slice
            if plot_var == "AT_pred":
                data = atNC.variables[plot_var][:, 90+lat_to_slice,: ]; 
                data=data.T
                data.filled(np.nan) 
            else:
                data = dicNC.variables[plot_var][:,90+lat_to_slice,: ];
                data=data.T
                data.filled(np.nan) 
        else: #longitude slice
            if plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            if plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
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
            if plot_var == "DIC_pred":
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
        
        
        f3_ax = fig42.add_subplot(gs[counter_hov:counter_hov+1,0:1])
        
        if region == "oceansoda_congo":
            contPlot1 = f3_ax.pcolor(dates,lons,  data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            f3_ax.set_xlabel("Time", fontsize=labelsize);
            f3_ax.set_ylabel("Longitude($^\circ$E)", fontsize=labelsize);
        else:
            contPlot1 = f3_ax.pcolor(dates,lats,  data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
            f3_ax.set_xlabel("Time", fontsize=labelsize);
            f3_ax.set_ylabel("Latitude($^\circ$N)", fontsize=labelsize);
        
        if region == "oceansoda_mediterranean":
            f3_ax.set_title("Longitude slice ({0}$^\circ$E)".format(long_to_slice*-1), fontsize=labelsize);
        elif region == "oceansoda_congo":
            f3_ax.set_title("Latitude slice ({0}$^\circ$S)".format(lat_to_slice), fontsize=labelsize);
        else:
            f3_ax.set_title("Longitude slice ({0}$^\circ$W)".format(long_to_slice), fontsize=labelsize);

        
        x=np.where(hasData)[0];
        firsttimeind=x[0];
        lastimeind=x[-1];

        if region == "oceansoda_amazon_plume":
            f3_ax.set_ylim(0, 24);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        elif region == "oceansoda_congo":
            f3_ax.set_ylim(-3, 12);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        elif region == "oceansoda_mediterranean":
            f3_ax.set_ylim(28, 48);
            f3_ax.set_xlim(dates[firsttimeind], dates[lastimeind]);
        else:
            pass
        cb = plt.colorbar(contPlot1)
    
        if plot_var == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
        elif plot_var == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
 
        counter_hov=counter_hov+1;
        
        #check that folder exists to save images, otherwise create it
        from pathlib import Path
        Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
            
    fig42.savefig(path.join("os_plots\{0}_seasonal".format(region), "_hov_multivar_all_vars.png"));
    plt.close();
    del varALL;        


    #### DIC and TA in 2015
    plot_var="DIC_pred"
    plot_var2="AT_pred"
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
            if plot_var == "DIC_pred":
                maxvar = 2000;
                minvar = 1400; 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
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
            if plot_var == "DIC_pred":
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

    
        if plot_var == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
            f3_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
        elif plot_var == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
            if plot_var2 == "DIC_pred":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT_pred":
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
            if plot_var2 == "DIC_pred":
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
            if plot_var2 == "DIC_pred":
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
        if plot_var2 == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
            f3_ax.set_title("DIC for {0} {1}".format(date.year, format(date.month, "02d")));  
        elif plot_var2 == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
    plot_var="DIC_pred"
    plot_var2="AT_pred"
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
            if plot_var == "DIC_pred":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
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

        if plot_var == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
        elif plot_var == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
            if plot_var2 == "DIC_pred":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var2 == "AT_pred":
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
            if plot_var2 == "DIC_pred":
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
            if plot_var2 == "DIC_pred":
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
        if plot_var2 == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
        elif plot_var2 == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
    plot_var="DIC_pred"
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
    
    fig6 = plt.figure(figsize=(15,10))
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
            if plot_var == "DIC_pred":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
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
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'Latitude $^({\circ})$', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'Longitude $(^{\circ})$', va='bottom', ha='center',
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

        if plot_var == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
        elif plot_var == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
            
        if counter==0:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(a)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==1:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(b)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==2:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(c)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==3:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(d)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
            
        counter=counter+1;
            
    #check that folder exists to save images, otherwise create it
    from pathlib import Path
    Path("os_plots\{0}_seasonal".format(region)).mkdir(parents=True, exist_ok=True)
    
    fig6.savefig(path.join("os_plots\{0}_seasonal".format(region), "{0}_monthly_averages.png".format(plot_var)));
    fig6.tight_layout()
    plt.close();
    del varALL;     
    
    #### Figure 2b - TA mmm
    plot_var="AT_pred"
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
    
    fig6 = plt.figure(figsize=(15,10))
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
            if plot_var == "DIC_pred":
                maxvar = 2100;
                minvar = 1000; 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
                maxvar = np.nanmax(varALL);
                minvar = np.nanmin(varALL); 
            elif plot_var == "AT_pred":
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
            if plot_var == "DIC_pred":
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
        
        f3_ax = fig6.add_subplot(gs[counter:counter+1],projection=ccrs.PlateCarree())
        gl = f3_ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.25, linestyle='--');
        contPlot1 = f3_ax.pcolormesh(lons, lats, data, vmin=minvar, vmax=maxvar, cmap=cmap_touse);
        gl.xlabels_top = False;
        gl.ylabels_right = False;
        gl.xformatter = LONGITUDE_FORMATTER;
        gl.yformatter = LATITUDE_FORMATTER;
        gl.xlabel_style = {'size': ticksize};
        gl.ylabel_style = {'size': ticksize};
        f3_ax.text(-0.14, 0.55, 'Latitude $(^{\circ})$', va='bottom', ha='center',
            rotation='vertical', rotation_mode='anchor',
            transform=f3_ax.transAxes)
        f3_ax.text(0.5, -0.2, 'Longitude $(^{\circ})$', va='bottom', ha='center',
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

        if plot_var == "DIC_pred":
            cb.set_label("DIC $(\mu mol^{-1})$");
        elif plot_var == "AT_pred":
            cb.set_label("TA $(\mu mol^{-1})$");
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
            
        if counter==0:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(a)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==1:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(b)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==2:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(c)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
        elif counter==3:
            f3_ax_txt = f3_ax.text(-0.03, 1.06, '(d)', transform=f3_ax.transAxes, fontsize='x-large', weight="bold");
            
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
with open("Amazon_AT_fs_pickled.txt", "rb") as myFile:
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

        #### Load in Amazon DIC pickle 

with open("Amazon_DIC_fs_pickled.txt", "rb") as myFile:
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
with open("Congo_AT_fs_pickled.txt", "rb") as myFile:
    Congo_at_dict = pickle.load(myFile)

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

with open("CONGO_DIC_fs_pickled.txt", "rb") as myFile:
    CONGO_DIC_dict = pickle.load(myFile)

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
    
    
        #### now make the figure
figsize = (16,24);
fig1 = plt.figure(figsize=figsize);
gs = fig1.add_gridspec(46, 50)

f1_ax1 = fig1.add_subplot(gs[0:15, 5:35])
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,0]), 1),RMSD_matrix_AT_AMAZON[:,0],"r|",markersize=10, label=list_name[0]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,1]), 1),RMSD_matrix_AT_AMAZON[:,1],"b|", markersize=10,label=list_name[1]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,2]), 1),RMSD_matrix_AT_AMAZON[:,2],"y|", markersize=10,label=list_name[2]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,3]), 1),RMSD_matrix_AT_AMAZON[:,3],"rx",markersize=10, label=list_name[3]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,4]), 1),RMSD_matrix_AT_AMAZON[:,4],"bx",markersize=10, label=list_name[4]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,5]), 1),RMSD_matrix_AT_AMAZON[:,5],"yx", markersize=10,label=list_name[5]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,6]), 1),RMSD_matrix_AT_AMAZON[:,6],"r^", markersize=10, label=list_name[6]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,7]), 1),RMSD_matrix_AT_AMAZON[:,7],"b^", markersize=10,label=list_name[7]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,8]), 1),RMSD_matrix_AT_AMAZON[:,8],"y^", markersize=10,label=list_name[8]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,9]), 1),RMSD_matrix_AT_AMAZON[:,9],"r_", markersize=10,label=list_name[9]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,10]), 1),RMSD_matrix_AT_AMAZON[:,10],"b_",markersize=10, label=list_name[10]);
f1_ax1.plot(np.arange(0, len(RMSD_matrix_AT_AMAZON[:,11]), 1),RMSD_matrix_AT_AMAZON[:,11],"y_",markersize=10, label=list_name[11]);
f1_ax1.set_xticks(np.arange(0, len(unqiuealgos_at_amazon), 1))
f1_ax1.set_xticklabels(list(unqiuealgos_at_amazon), rotation=90,fontsize=12)
f1_ax1.set_ylim(0, 100);
f1_ax1.set_ylabel("RMSDe",fontsize=12);
f1_ax1.set_title('Amazon TA',fontsize=12);
f1_ax1.legend(loc="upper left",ncol=1);

f1_ax2 = fig1.add_subplot(gs[23:38, 5:35])
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,0]), 1),RMSD_matrix_DIC_AMAZON[:,0],"r|",markersize=10, label=list_name[0]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,1]), 1),RMSD_matrix_DIC_AMAZON[:,1],"b|", markersize=10,label=list_name[1]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,2]), 1),RMSD_matrix_DIC_AMAZON[:,2],"y|", markersize=10,label=list_name[2]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,3]), 1),RMSD_matrix_DIC_AMAZON[:,3],"rx",markersize=10, label=list_name[3]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,4]), 1),RMSD_matrix_DIC_AMAZON[:,4],"bx",markersize=10, label=list_name[4]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,5]), 1),RMSD_matrix_DIC_AMAZON[:,5],"yx", markersize=10,label=list_name[5]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,6]), 1),RMSD_matrix_DIC_AMAZON[:,6],"r^", markersize=10, label=list_name[6]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,7]), 1),RMSD_matrix_DIC_AMAZON[:,7],"b^", markersize=10,label=list_name[7]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,8]), 1),RMSD_matrix_DIC_AMAZON[:,8],"y^", markersize=10,label=list_name[8]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,9]), 1),RMSD_matrix_DIC_AMAZON[:,9],"r_", markersize=10,label=list_name[9]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,10]), 1),RMSD_matrix_DIC_AMAZON[:,10],"b_",markersize=10, label=list_name[10]);
f1_ax2.plot(np.arange(0, len(RMSD_matrix_DIC_AMAZON[:,11]), 1),RMSD_matrix_DIC_AMAZON[:,11],"y_",markersize=10, label=list_name[11]);
f1_ax2.set_xticks(np.arange(0, len(unqiuealgos_dic_amazon), 1))
f1_ax2.set_xticklabels(list(unqiuealgos_dic_amazon), rotation=90,fontsize=12)
f1_ax2.set_ylim(0, 100);
f1_ax2.set_ylabel("RMSDe",fontsize=12);
f1_ax2.set_title('Amazon DIC',fontsize=12);


f1_ax3 = fig1.add_subplot(gs[0:15, 40:50])
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,0]), 1),RMSD_matrix_AT_Congo[:,0],"r|",markersize=10, label=list_name_Congo_AT[0]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,1]), 1),RMSD_matrix_AT_Congo[:,1],"b|", markersize=10,label=list_name_Congo_AT[1]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,2]), 1),RMSD_matrix_AT_Congo[:,2],"y|", markersize=10,label=list_name_Congo_AT[2]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,3]), 1),RMSD_matrix_AT_Congo[:,3],"rx",markersize=10, label=list_name_Congo_AT[3]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,4]), 1),RMSD_matrix_AT_Congo[:,4],"bx",markersize=10, label=list_name_Congo_AT[4]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,5]), 1),RMSD_matrix_AT_Congo[:,5],"yx", markersize=10,label=list_name_Congo_AT[5]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,6]), 1),RMSD_matrix_AT_Congo[:,6],"r_", markersize=10,label=list_name_Congo_AT[6]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,7]), 1),RMSD_matrix_AT_Congo[:,7],"b_", markersize=10,label=list_name_Congo_AT[7]);
f1_ax3.plot(np.arange(0, len(RMSD_matrix_AT_Congo[:,8]), 1),RMSD_matrix_AT_Congo[:,8],"y_", markersize=10,label=list_name_Congo_AT[8]);
f1_ax3.set_xticks(np.arange(0, len(unqiuealgos_at_Congo), 1))
f1_ax3.set_xticklabels(list(unqiuealgos_at_Congo), rotation=90,fontsize=12)
f1_ax3.set_ylim(0, 100);
f1_ax3.set_title('Congo TA',fontsize=12);


f1_ax4 = fig1.add_subplot(gs[23:38, 40:50])
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,0]), 1),RMSD_matrix_DIC_CONGO[:,0],"r|",markersize=10, label=list_name_Congo_AT[0]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,1]), 1),RMSD_matrix_DIC_CONGO[:,1],"b|", markersize=10,label=list_name_Congo_AT[1]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,2]), 1),RMSD_matrix_DIC_CONGO[:,2],"y|", markersize=10,label=list_name_Congo_AT[2]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,3]), 1),RMSD_matrix_DIC_CONGO[:,3],"rx",markersize=10, label=list_name_Congo_AT[3]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,4]), 1),RMSD_matrix_DIC_CONGO[:,4],"bx",markersize=10, label=list_name_Congo_AT[4]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,5]), 1),RMSD_matrix_DIC_CONGO[:,5],"yx", markersize=10,label=list_name_Congo_AT[5]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,6]), 1),RMSD_matrix_DIC_CONGO[:,6],"r_", markersize=10,label=list_name_Congo_AT[6]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,7]), 1),RMSD_matrix_DIC_CONGO[:,7],"b_", markersize=10,label=list_name_Congo_AT[7]);
f1_ax4.plot(np.arange(0, len(RMSD_matrix_DIC_CONGO[:,8]), 1),RMSD_matrix_DIC_CONGO[:,8],"y_", markersize=10,label=list_name_Congo_AT[8]);
f1_ax4.set_xticks(np.arange(0, len(unqiuealgos_dic_CONGO), 1))
f1_ax4.set_xticklabels(list(unqiuealgos_dic_CONGO), rotation=90,fontsize=12)
f1_ax4.set_ylim(0, 100);
f1_ax4.set_title('Congo DIC',fontsize=12);


plt.savefig("os_plots\\Fig1_RMSDe.png");










