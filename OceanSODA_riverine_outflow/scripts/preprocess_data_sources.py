#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 15:13:11 2020

@author: verwirrt
"""

import pandas as pd;
import numpy as np;
from os import path;
from datetime import timedelta;


def load_discharge_data(carbonateDates, region):
    #resample river discharge data to monthly totals spanning the same data range as the carbonate data
    #Assumes riverOutput is total daily discharge in m^3
    #dischargeColName is the column in the provided riverDischarge dataframe
    #scaleFactor allows conversion e.g. from discharge per second to discharge per day (must be set correctly to ensure each original discharge amount of correct when scaled to the full original time period)
    #returns a dictionary of dataframes for each ocean soda river region
    #   Each data frame contains the "date", "monthly_discharge" (m^3) and "monthly_discharge_sd" (standard deviation, m^3)
    def resample_to_carbonate_dates(riverDischarge, dischargeColName, carbonateDates, dateColName="date"):
        output = pd.DataFrame();
        output["date"] = carbonateDates;
        output.index = carbonateDates;
        
        #Calculate monthly discharge by resampling for month, and then trimming data range to match the carbonate output
        monthlyOutflow = riverDischarge.resample("M", label="left", loffset="1D", on=dateColName).sum();
        #monthlyOutflow[dischargeColName] = monthlyOutflow[dischargeColName] * scaleFactor; #convert from per second to per day (daily values already summed)
        inDateRange = (monthlyOutflow.index >= carbonateDates[0]) & (monthlyOutflow.index <= carbonateDates[-1]); #remove dates that have no carbonate data for
        output["monthly_discharge"] = monthlyOutflow[dischargeColName][inDateRange];
        
        
        #Calculate monthly discharge uncertainty
        monthlyOutflowSD = riverDischarge.resample("M", label="left", loffset="1D", on=dateColName).std();
        #monthlyOutflowSD[dischargeColName] = monthlyOutflowSD[dischargeColName] * scaleFactor; #convert from per second to per day (daily values already summed)
        output["monthly_discharge_sd"] = monthlyOutflowSD[dischargeColName][inDateRange];
        
        return output;
    
    
    if region == "oceansoda_amazon_plume":
        ######################
        # process amazon data
        ######################
        #Read data file
        amazonDischargePath = "../data/obidos/17050001_debits.csv";
        riverDischarge = pd.read_csv(amazonDischargePath, parse_dates=["date"], dayfirst=True);
        
        #Convert to daily discharge in m^3
        riverDischarge["discharge"] = riverDischarge["discharge"] * (60*60*24); #Scale from m^3 per second discharge to m^3 per day
        
        #Resample to calculate monthly data and trim off any unrequired dates
        riverDischarge = resample_to_carbonate_dates(riverDischarge, "discharge", carbonateDates); 
        
        #finished processing, store in dictionary
        return riverDischarge;
    
    
    
    if region == "oceansoda_congo":
        #####################
        # process congo data
        #####################
        #Read data file
        congoDischargePath = "../data/congo_brazzaville/50800000_debits.csv";
        riverDischarge = pd.read_csv(congoDischargePath, parse_dates=["date"], dayfirst=True);
        
        #Convert to daily discharge in m^3
        riverDischarge["valeur"] = riverDischarge["valeur"] * (60*60*24); #Scale from m^3 per second discharge to m^3 per day
        
        #Resample to calculate monthly data and trim off any unrequired dates
        riverDischarge = resample_to_carbonate_dates(riverDischarge, "valeur", carbonateDates); 
        
        #finished processing, store in dictionary
        return riverDischarge;
    
    
    
    if region == "oceansoda_mississippi":
        ###########################
        # process mississippi data
        ###########################
        #Read data files
        mississippiDischargePath = "../data/mississippi_tarbert_landing/tarbert_landing_1980_2020.csv";
        riverDischarge = pd.read_csv(mississippiDischargePath, parse_dates=["Date / Time"], dayfirst=True);
        
        #remove missing values
        riverDischarge["Flow (CFS)"][riverDischarge["Flow (CFS)"]=="M"] = np.nan;
        
        #Convert to daily discharge in m^3
        riverDischarge["discharge"] = riverDischarge["Flow (CFS)"].astype(float) * (60*60*24) / 35.3147; #Scale from cubic feet per second discharge to m^3 per day (35.3147 cubic feet in a cubic metre)
        
        #Resample to calculate monthly data and trim off any unrequired dates
        riverDischarge = resample_to_carbonate_dates(riverDischarge, "discharge", carbonateDates, dateColName="Date / Time"); 
        
        
        #finished processing, store in dictionary
        return riverDischarge;
    
    
    
    if region == "oceansoda_st_lawrence":
        ###########################
        # process st lawrence data
        ###########################
        #Read data files
        stLawrenceDischargePathHistoric = "../data/st_lawrence/Historic_Daily__Apr-30-2020_10_34_18AM.csv";
        riverDischargeHistoric = pd.read_csv(stLawrenceDischargePathHistoric, parse_dates=["Date"], dayfirst=True);
        stLawrenceDischargePathRecent = "../data/st_lawrence/Recent_5min_02OA016_QR_Apr-30-2020_10_58_35AM_cols_renamed.csv";
        riverDischargeRecent = pd.read_csv(stLawrenceDischargePathRecent, parse_dates=["Date"], dayfirst=True);
        
        ###combine historic and recent
        #riverDischargeRecent["Date2"] = riverDischargeRecent["Date"].copy();
        riverDischargeRecentDailyMean = riverDischargeRecent.resample("D", label="left", loffset="1D", on="Date").mean(); #calculate daily means from recent data
        riverDischargeRecentDailyMean.insert(2, "Date", riverDischargeRecentDailyMean.index-timedelta(days=1));
        riverDischargeRecentDailyMean.index = range(len(riverDischargeHistoric), len(riverDischargeHistoric)+len(riverDischargeRecentDailyMean));
        #combined dataframes
        riverDischarge = pd.concat([riverDischargeHistoric, riverDischargeRecentDailyMean]);
        
        #Convert to daily discharge in m^3
        riverDischarge["Value"] = riverDischarge["Value"] * (60*60*24); #Scale from m^3 per second discharge to m^3 per day
        
        #Resample to calculate monthly data and trim off any unrequired dates
        riverDischarge = resample_to_carbonate_dates(riverDischarge, "Value", carbonateDates, dateColName="Date"); 
        
        #finished processing, store in dictionary
        return riverDischarge;






