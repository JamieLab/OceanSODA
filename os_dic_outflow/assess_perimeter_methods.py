#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 12:17:44 2020

@author: verwirrt
"""

import pandas as pd;
import numpy as np;
from string import Template;
import matplotlib.pyplot as plt;
import calendar;
import datetime;
from .preprocess_discharge_data import load_discharge_data;

#returns river discharge in kg month-1
def get_discharge(df, year, month):
    row = df.loc[(df["date"] == datetime.datetime(year, month, 1))].iloc[0];
    discharge = row["monthly_discharge"]; #in m^3 month-1
    dischargeKg = discharge*1000.0; #in kg month-1: 1000 kg in 1 m^3 fresh water
    return dischargeKg;

#returns mean discharge for each month in kg month-1
def get_monthly_mean_discharge(df):
    meanDischarges = [];
    for month in range(1, 13):
        currentMonthRows = df.loc[df["date"].dt.month == month];
        meanDischarge = np.nanmean(currentMonthRows["monthly_discharge"]) * 1000.0; #means * 1000 to convert from m^3 month-1 to kg month-1
        meanDischarges.append(meanDischarge);
    return np.array(meanDischarges);

#returns DIC river outflow (and SD) in tg C month-1
def get_outflow(df, year, month):
    row = df.loc[(df["date"] == datetime.datetime(year, month, 1))].iloc[0];
    return row["dic_outflow"], row["dic_outflow_sd"];

def g_to_Tg(val):
    return val / (10.0**12);

def get_radii(dfKeys):
    radii = [];
    for key in dfKeys:
        if "radius" in key:
            radii.append(int(key.split("radius")[1]));
    radii = np.unique(radii);
    return radii;
            
def replace_0s_with_nans_in_dataframe(df):
    cols = [col for col in df.keys() if "dic_outflow" in col];
    df[cols] = df[cols].replace({0:np.nan});
    return df;

useLong = True;
DATASET_STR = "_min_range" if useLong else "";

region = "oceansoda_amazon_plume";

dicTimeseriesPathTemplate = Template("../output/${REGION}${DATASET_STR}/monthly_timeseries_${REGION}.csv");
dicAnnualPathTemplate = Template("../output/${REGION}${DATASET_STR}/interyear_timeseries_${REGION}.csv");

dicTimeseriesDF = pd.read_csv(dicTimeseriesPathTemplate.safe_substitute(REGION=region, DATASET_STR=DATASET_STR), parse_dates=["date"]);
dicTimeseriesDF = replace_0s_with_nans_in_dataframe(dicTimeseriesDF); #Remove 0s and replace with NaNs, because this 0s are where the plume misses the perimeter
dicAnnualDF = pd.read_csv(dicAnnualPathTemplate.safe_substitute(REGION=region, DATASET_STR=DATASET_STR));

radii = get_radii(dicAnnualDF.keys());
monthDates = [d.to_pydatetime() for d in pd.date_range("1960-01-01", "2019-11-01", freq="MS")];
meanMonthlyDischarge = load_discharge_data(monthDates, region);
meanMonthlyDischarge = meanMonthlyDischarge[(meanMonthlyDischarge["date"] >= datetime.datetime(1991, 1, 1)) & (meanMonthlyDischarge["date"] <= datetime.datetime(2019,12,31))];
meanAnnualDischarge = get_monthly_mean_discharge(meanMonthlyDischarge);



#####Plot annual estimates for all radii
fig = plt.figure(figsize=(10,6));
plt.fill_between(dicAnnualDF["month_name"], dicAnnualDF["dic_outflow_mean"]-dicAnnualDF["dic_outflow_mean_sd"], dicAnnualDF["dic_outflow_mean"]+dicAnnualDF["dic_outflow_mean_sd"], alpha=0.4, color='b', linewidth=0);
plt.plot(dicAnnualDF["dic_outflow_mean"], 'b', label="Mean DIC outflow", linewidth=2);
plt.plot(dicAnnualDF["dic_outflow_median"], 'b:', label="Median DIC outflow", linewidth=2);
plt.ylabel("DIC outflow (Tg C month-1)");
#Add individual radius estimates
for radius in radii:
    if radius in [min(radii), max(radii)]:
        plt.plot(dicAnnualDF["dic_outflow_radius"+str(radius)], 'k', label="radius=%d"%radius, alpha=1.0-(radius/(max(radii)+2)));
    else: #no label
        plt.plot(dicAnnualDF["dic_outflow_radius"+str(radius)], 'k', alpha=1.0-(radius/max(radii)));
plt.legend(loc=0);



#####Plot monthly estimates for all radii
hasMonthData = np.isfinite(dicTimeseriesDF["dic_outflow_mean"])
plt.subplots(2, 1, figsize=(10,6));
plt.subplot(2,1,1);
plt.fill_between(np.arange(0, len(dicTimeseriesDF["date"])), dicTimeseriesDF["dic_outflow_mean"]-dicTimeseriesDF["dic_outflow_mean_sd"], dicTimeseriesDF["dic_outflow_mean"]+dicTimeseriesDF["dic_outflow_mean_sd"], alpha=0.4, color='b', linewidth=0);
plt.plot(dicTimeseriesDF["dic_outflow_mean"], 'b', label="Mean DIC outflow", linewidth=2);
plt.plot(dicTimeseriesDF["dic_outflow_median"], 'b:', label="Median DIC outflow", linewidth=2);
plt.ylabel("DIC outflow (Tg C month-1)");
for radius in radii:
    plt.plot(dicTimeseriesDF["dic_outflow_radius"+str(radius)], 'k', alpha=0.1); #, label="radius=%d"%radius
plt.legend(loc=0);
plt.subplot(2,1,2);
plt.plot(meanMonthlyDischarge["date"].values[hasMonthData[:-1]], meanMonthlyDischarge["monthly_discharge"].values[hasMonthData[:-1]], 'k');
plt.ylabel("discharge (m3 month-1)");



######Calculate annual mean DIC outflows
dicAnnualEstimate = np.nansum(dicAnnualDF["dic_outflow_mean"]);
dicAnnualEstimateSD = np.sqrt(np.nansum(dicAnnualDF["dic_outflow_mean_sd"]**2));
print("Annual mean DIC outflow:", format(dicAnnualEstimate, ".4"), "(+/-"+format(dicAnnualEstimateSD, ".2")+") TgC yr-1");


def recalculate_mean_from_specific_radii(dicDF, radiiList):
    colsSelected = ["dic_outflow_radius"+str(radius) for radius in radiiList];
    
    subset = dicDF[colsSelected];
    meanData = np.nanmean(subset.values, axis=1);
    return meanData;


######Compare with different number of radii
numRadiiList = range(1, len(radii)+1); #[3, 5, 10, 15, 20, 25];
plt.subplots(2,1, figsize=(10,7));
plt.subplot(2,1,1);
plt.ylabel("DIC outflow (Tg C month-1)");
annualMeansDiffRadii = {};
for numRadii in numRadiiList:
    alpha = 1.0 - (numRadii/(max(numRadiiList)+2));
    radiiList = radii[0:numRadii];
    meanData = recalculate_mean_from_specific_radii(dicAnnualDF, radiiList);
    annualMeansDiffRadii[numRadii] = meanData.copy();
    if numRadii in [min(numRadiiList), max(numRadiiList)]:
        plt.plot(dicAnnualDF["month_name"], meanData, 'k', label="Num radii: "+str(numRadii), alpha=alpha);
    else: #Plot without label
        plt.plot(dicAnnualDF["month_name"], meanData, 'k', alpha=alpha);
plt.legend(loc=0);

plt.subplot(2,1,2);
plt.ylabel("DIC outflow (Tg C month-1)");
for numRadii in numRadiiList:
    alpha = 1.0 - (numRadii/(max(numRadiiList)+2));
    radiiList = radii[0:numRadii];
    meanData = recalculate_mean_from_specific_radii(dicTimeseriesDF, radiiList);
    plt.plot(dicTimeseriesDF["date"], meanData, 'k', label="Num radii: "+str(numRadii), alpha=alpha);
#plt.legend(loc=0);
plt.tight_layout();


######Flip the plot arround, plot single months at different radii to see if we can see convergence
plt.figure(figsize=(10,5));
for month in range(1, 13):
    monthStr = calendar.month_abbr[month];
    dataForMonth = [annualMeansDiffRadii[numRadii][month-1] for numRadii in annualMeansDiffRadii.keys()];
    colour = plt.cm.jet(month/12);
    plt.plot(numRadiiList, dataForMonth, 'g', label=monthStr, color=colour);#, alpha=(month+2)/14);

plt.ylabel("Month mean DIC outflow (TgC month-1)");
plt.xlabel("Number of radii used to calculate mean");
plt.legend(loc=0);
plt.tight_layout();








