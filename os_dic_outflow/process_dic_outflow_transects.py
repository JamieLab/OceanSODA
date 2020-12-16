#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 09:28:25 2020

@author: tom holding
"""

import pandas as pd;
import numpy as np;
import matplotlib.pyplot as plt;
import datetime;
from os import path, makedirs;
from scipy.stats import pearsonr;

from .preprocess_discharge_data import load_discharge_data;


def g_to_Tg(val):
    return val / (10.0**12);

#returns river discharge in kg month-1 for a specific month in a specific year
def get_discharge_at_month_year(df, year, month):
    row = df.loc[(df["date"] == datetime.datetime(year, month, 1))].iloc[0];
    discharge = row["monthly_discharge"]; #in m^3 month-1
    dischargeKg = discharge*1000.0; #in kg month-1: 1000 kg in 1 m^3 fresh water
    return dischargeKg;

#returns the number of radii included in a radiiData dataset
#Assumes radii start from 1 in increment by 1 to n. Returns n.
def get_num_radii(radiiData):
    numRadii = 0;
    for key in radiiData.keys():
        if "dic_outflow_radius" in key:
            numRadii+=1;
    return numRadii;

#returns the mean DIC outflow using a set of radii (the set is assumed to include 1:numRadii inclusive)
def calculate_mean_of_radii_data_dic_outflow(radiiData, numRadii):
    radii = range(1, numRadii+1);
    colNames = ["dic_outflow_radius"+str(radius) for radius in radii];
    colNamesSD = ["dic_outflow_sd_radius"+str(radius) for radius in radii];
    
    radiiSubset = radiiData[colNames];
    radiiSubset.values[radiiSubset == 0] = np.nan; #remove any 0s (as these are where the circle perimeter does not intersect the plume at all)
    radiiSubsetSD = radiiData[colNamesSD];
    
    radiiSetMean = np.nanmean(radiiSubset.values, axis=1);
    radiiSetSD = np.sqrt(np.nansum(radiiSubsetSD.values**2, axis=1))/np.sum(np.isnan(radiiSubsetSD.values)==False, axis=1);
    
    return radiiSetMean, radiiSetSD;


#Choose the maximum number of radii to use such that the calculated mean DIC outflow peaks at same time as discharge
#   dischargeData: mean discharge for each month
#   radiiDataSeasonal: seasonal monthly DIC outflow for each radii (e.g. one value for each month of the year / radii combination)
def determine_num_radii_with_peak_alignment(dischargeData, radiiDataSeasonal):
    meanSeasonalDischarge = dischargeData.groupby([dischargeData.date.dt.month])["monthly_discharge"].mean().values;
    maxDischargeIndex = np.argmax(meanSeasonalDischarge);
    
    selectedNumRadii = 0;
    distanceBetweenPeaks = 12;
    
    maxNumRadii = get_num_radii(radiiDataSeasonal);
    for numRadii in range(1, maxNumRadii+1):
        #Calculate mean DIC outflow using numRadii radii
        radiiSetMeanSeasonalDischarge, radiiSetSDSeasonalDischarge = calculate_mean_of_radii_data_dic_outflow(radiiDataSeasonal, numRadii);
        
        #Find maximum discharge month
        maxOutflowIndex = np.argmax(radiiSetMeanSeasonalDischarge);
        
        #If the current numRadii reduce or maintain the distance to peak, use it
        if np.abs(maxOutflowIndex - maxDischargeIndex) <= distanceBetweenPeaks:
            selectedNumRadii = numRadii;
            distanceBetweenPeaks = np.abs(maxOutflowIndex - maxDischargeIndex);
        
#        #Compare with maxDischargeIndex and if the same, update selectedNumRadii
#        if maxOutflowIndex == maxDischargeIndex:
#            selectedNumRadii = numRadii;
    
    return selectedNumRadii, distanceBetweenPeaks;


def determine_num_radii_with_correlation(dischargeData, radiiDataSeasonal):
    meanSeasonalDischarge = dischargeData.groupby([dischargeData.date.dt.month])["monthly_discharge"].mean().values;
    
    selectedNumRadii = 0;
    bestCorrelationCoefficient = 0;
    maxNumRadii = get_num_radii(radiiDataSeasonal);
    for numRadii in range(1, maxNumRadii+1):
        #Calculate mean DIC outflow using numRadii radii
        radiiSetMeanSeasonalDischarge, radiiSetSDSeasonalDischarge = calculate_mean_of_radii_data_dic_outflow(radiiDataSeasonal, numRadii);
        
        #calculate correlation coefficient
        newCorrelationCoefficient, p = pearsonr(meanSeasonalDischarge, radiiSetSDSeasonalDischarge);
        #print(newCorrelationCoefficient, bestCorrelationCoefficient);
        
        #Compare to current best correlation coefficient and update as necessary
        if newCorrelationCoefficient > bestCorrelationCoefficient:
            bestCorrelationCoefficient = newCorrelationCoefficient;
            selectedNumRadii = numRadii;
    
    return selectedNumRadii, bestCorrelationCoefficient;


#Plot to compare alignment of temporal patterns in discharge and calculated DIC outflow using the peak alignment and the correlation methods
# for choosing number of radii
def peak_align_diagnostic_plot(region, dischargeDataMonths, radiiData, numRadii_amp, numRadii_corr, offset, correlationCoefficient):
    #Calculate seasonal means
    meanDischarge = dischargeDataMonths.groupby([dischargeDataMonths.date.dt.month])["monthly_discharge"].mean().values;
    dicOutflowMean_amp, dicOutflowSD_amp = calculate_mean_of_radii_data_dic_outflow(radiiData, numRadii_amp);
    dicOutflowMean_corr, dicOutflowSD_corr = calculate_mean_of_radii_data_dic_outflow(radiiData, numRadii_corr);
    #normalise
    meanDischarge /= np.nanmax(meanDischarge);
    dicOutflowMean_amp /= np.nanmax(dicOutflowMean_amp);
    dicOutflowMean_corr /= np.nanmax(dicOutflowMean_corr)
    #plot
    plt.figure(figsize=(10,10));
    plt.subplot(2,1,1);
    plt.title(region+" peak alignment diagnostic plot");
    plt.plot(meanDischarge, 'k', label="discharge (normalised)");
    plt.plot(dicOutflowMean_amp, 'r', label="CT outflow (normalised, peak amplitude)\nUsing {0} radii, peak offset = {1}".format(numRadii_amp, offset));
    plt.legend(loc=0);
    plt.subplot(2,1,2);
    plt.plot(meanDischarge, 'k', label="discharge (normalised)");
    plt.plot(dicOutflowMean_corr, 'r', label="CT outflow (normalised, correlation)\nUsing {0} radii, correlation coefficient = {1}".format(numRadii_corr, correlationCoefficient));
    plt.legend(loc=0);
    plt.tight_layout();



#inputDataRootTemplate: Template (REGION) defining the directory containing the monthly and interyear time series for each radii
def process_dic_outflow_transects(inputDataRootTemplate, regions, outputDirectory):
    fontSize = 10;
    legendSize = 10;
    figsizex = 6.5;
    figsizey = 3.0;
    figureSize = (figsizex, figsizey);
    
    if path.exists(outputDirectory) == False:
        makedirs(outputDirectory);
    plotOutputPath = outputDirectory;

    summaryData = {};
    
    #region = "oceansoda_amazon_plume";
    for region in regions:
        print("Processing DIC outflow for "+region+"...");
        
        #Load previously calculated DIC outflow data for each transect
        radiiDataMonthsPath = path.join(inputDataRootTemplate.safe_substitute(REGION=region), "monthly_timeseries_"+region+".csv");
        radiiDataMonths = pd.read_csv(radiiDataMonthsPath, parse_dates=["date"]);
        radiiDataSeasonalPath = path.join(inputDataRootTemplate.safe_substitute(REGION=region), "interyear_timeseries_"+region+".csv");
        radiiDataSeasonal = pd.read_csv(radiiDataSeasonalPath);
        
        monthDates = [d.to_pydatetime() for d in pd.date_range("1960-01-01", "2019-11-01", freq="MS")];
        dischargeMonths = load_discharge_data(monthDates, region); #discharge in m^3 month-1
        dischargeSeasonal = dischargeMonths.groupby([dischargeMonths.date.dt.month])["monthly_discharge"].mean().values; #discharge in m^3 month-1
        
        
        #Run using the discharge and DIC outflow peak amplitude matching to choose number of radii to use
        numRadii_amp, distanceBetweenPeaks = determine_num_radii_with_peak_alignment(dischargeMonths, radiiDataSeasonal);
        print("\tnumRadii (peak  applitude):", numRadii_amp, "distanceBeterrnPeaks:", distanceBetweenPeaks);
        radiiSetMean_amp, radiiSetSD_amp = calculate_mean_of_radii_data_dic_outflow(radiiDataMonths, numRadii_amp);
        radiiSetSeasonalMean_amp, radiiSetSeasonalSD_amp = calculate_mean_of_radii_data_dic_outflow(radiiDataSeasonal, numRadii_amp);
        summaryData[region+"_num_radii_amp"] = numRadii_amp;
        summaryData[region+"_distance_between_peaks"] = distanceBetweenPeaks;
        
        #Run using the discharge and DIC outflow maximum correlation to choose number of radii to use
        numRadii_corr, correlationCoefficient = determine_num_radii_with_correlation(dischargeMonths, radiiDataSeasonal);
        print("\tnumRadii (correlation):", numRadii_corr, "correlationCoefficient:", correlationCoefficient);
        radiiSetMean_corr, radiiSetSD_corr = calculate_mean_of_radii_data_dic_outflow(radiiDataMonths, numRadii_corr);
        radiiSetSeasonalMean_corr, radiiSetSeasonalSD_corr = calculate_mean_of_radii_data_dic_outflow(radiiDataSeasonal, numRadii_corr);
        summaryData[region+"_num_radii_corr"] = numRadii_corr;
        summaryData[region+"_correlation_coefficient"] = correlationCoefficient;
        
        #Diagnostic plot (useful to check that the method worked)
        if (distanceBetweenPeaks != 0):
            print("\n*** WARNING: When selecting num radii using peak alignment the distance between discharge and DIC outflow peak could not be reduced to 0. Minimum distances was {0}.".format(distanceBetweenPeaks));
        peak_align_diagnostic_plot(region, dischargeMonths, radiiDataSeasonal, numRadii_amp, numRadii_corr, distanceBetweenPeaks, correlationCoefficient);
        plt.savefig(path.join(plotOutputPath, region+"_peak_diagnostic.png"));
        
        
        ##################Annual values
        summaryData[region+"_annual_DIC_outflow_amp"] = np.sum(radiiSetSeasonalMean_amp); #Tg C yr-1
        summaryData[region+"_annual_DIC_outflow_uncert_amp"] = np.sqrt(np.sum(radiiSetSeasonalSD_amp**2)); #Tg C yr-1
        summaryData[region+"_annual_DIC_outflow_corr"] = np.sum(radiiSetSeasonalMean_corr); #Tg C yr-1
        summaryData[region+"_annual_DIC_outflow_uncert_corr"] = np.sqrt(np.sum(radiiSetSeasonalSD_corr**2)); #Tg C yr-1
        
        ##################Interannual variation
        radiiDataYears = radiiDataMonths.groupby(radiiDataMonths.date.dt.year).sum();
        radiiDataYearsMean_amp, radiiDataYearsSD_amp = calculate_mean_of_radii_data_dic_outflow(radiiDataYears, numRadii_amp);
        summaryData[region+"_annual_DIC_outflow_SD_amp"] = np.nanstd(radiiDataYearsMean_amp);
        summaryData[region+"_annual_DIC_outflow_CoV_amp"] = np.nanstd(radiiDataYearsMean_amp) / np.nanmean(radiiDataYearsMean_amp);
        radiiDataYearsMean_corr, radiiDataYearsSD_corr = calculate_mean_of_radii_data_dic_outflow(radiiDataYears, numRadii_corr);
        summaryData[region+"_annual_DIC_outflow_SD_corr"] = np.nanstd(radiiDataYearsMean_corr);
        summaryData[region+"_annual_DIC_outflow_CoV_corr"] = np.nanstd(radiiDataYearsMean_corr) / np.nanmean(radiiDataYearsMean_corr);
        
        
        ##################Plot seasonal values
        plt.figure(figsize=figureSize);
        plt.fill_between(radiiDataSeasonal["month_name"], radiiSetSeasonalMean_amp-radiiSetSeasonalSD_amp, radiiSetSeasonalMean_amp+radiiSetSeasonalSD_amp, color='r', alpha=0.4, linewidth=0);
        plt.plot(radiiSetSeasonalMean_amp, 'r', linewidth=2, label="peak alignment method, $n_r={0}$".format(numRadii_amp));
        plt.fill_between(radiiDataSeasonal["month_name"], radiiSetSeasonalMean_corr-radiiSetSeasonalSD_corr, radiiSetSeasonalMean_corr+radiiSetSeasonalSD_corr, color='b', alpha=0.4, linewidth=0);
        plt.plot(radiiSetSeasonalMean_corr, 'b', linewidth=2, label="correlation method, $n_r={0}$".format( numRadii_corr));
        plt.ylabel("CT outflow (Tg C month$^{-1}$)", fontsize=fontSize);
        
        
        #Add comparative estimates to the plots based on region
        if region == "oceansoda_amazon_plume":
            ###Add Richey et al 1990 and 1991 eestimates
            richeyDICConcentration1 = 500.0; #umol kg-1, two values to span the range quoted in Richey et al 1991
            richeyDICConcentration2 = 600.0; #umol kg-1, two values to span the range quoted in Richey et al 1991
            richeyDICConcentration1grams = richeyDICConcentration1 * 12.0107 / 1000000; #g C kg-1
            richeyDICConcentration2grams = richeyDICConcentration2 * 12.0107 / 1000000; #g C kg-1
            richeyDICConcentration1grams = richeyDICConcentration1grams * 1000; #g C m-3 #1000 kg in 1m^3
            richeyDICConcentration2grams = richeyDICConcentration2grams * 1000; #g C m-3 #1000 kg in 1m^3
            #dates = pd.date_range('1982-01-01','1984-12-01', freq='MS').to_pydatetime();
            #discharge = np.array([get_discharge_at_month_year(dischargeMonths, d.year, d.month) for d in dates]).reshape(3, 12).T;
            dicOutflow1 = g_to_Tg(dischargeSeasonal*richeyDICConcentration1grams); #Tg C month-1
            dicOutflow2 = g_to_Tg(dischargeSeasonal*richeyDICConcentration2grams); #Tg C month-1
            plt.plot(range(0,12), dicOutflow2, label="Richey et al (2000) (upper)");
            plt.plot(range(0,12), dicOutflow1, label="Richey et al (2000) (lower)");
            
            ###Add druffel et al 2005 results (single cruise, riverine end-member extrapolation from salinity relationship)
            druffelDICConcentration = 363.0; #umol kg-1, see section 4.1
            druffelDICConcentrationGrams = druffelDICConcentration * 12.0107 / 1000000; #g C kg-1
            druffelDICConcentrationGrams = druffelDICConcentrationGrams  * 1000; #g C m-3 #1000 kg in 1m^3
            #Druffel 1991 November
            #discharge = get_discharge(dfAllDischarge, 1991, 11); #kg month-1
            dicOutflowNov = g_to_Tg(dischargeSeasonal[11-1]*druffelDICConcentrationGrams); #Tg C month-1
            #Druffel 1991 December
            #discharge = get_discharge(dfAllDischarge, 1991, 12); #kg month-1
            dicOutflowDec = g_to_Tg(dischargeSeasonal[12-1]*druffelDICConcentrationGrams); #Tg C month-1
            #plot
            plt.scatter([11-1, 12-1], [dicOutflowNov, dicOutflowDec], label="Druffel et al 2005");
    
        
        plt.legend(loc=0, fontsize=legendSize);
        plt.tick_params(labelsize=fontSize)
        plt.tight_layout();
        plt.savefig(path.join(plotOutputPath, region+"_seasonal.pdf"));
        plt.savefig(path.join(plotOutputPath, region+"_seasonal.png"));
        #plt.close();
        
        
        ##################Plot monthly values
        date = [str(d.year)+"-"+str(format(d.month, "02d")) for d in radiiDataMonths["date"]]; #x vals to plot
        #date = [datetime.datetime(d.year, d.month, 1) for d in radiiDataMonths["date"]];
        plt.figure(figsize=figureSize);
        plt.fill_between(date, radiiSetMean_amp-radiiSetSD_amp, radiiSetMean_amp+radiiSetSD_amp, color='r', alpha=0.4, linewidth=0);
        plt.plot(date, radiiSetMean_amp, 'r', linewidth=2, label="peak alignment method, $n_r={0}$".format(numRadii_amp));
        plt.fill_between(date, radiiSetMean_corr-radiiSetSD_corr, radiiSetMean_corr+radiiSetSD_corr, color='b', alpha=0.4, linewidth=0);
        plt.plot(date, radiiSetMean_corr, 'b', linewidth=2, label="correlation method, $n_r={0}$".format(numRadii_corr));
        
        
        wHasData = np.where(np.isfinite(radiiSetMean_amp));
        plt.xticks(np.arange(wHasData[0][0], wHasData[0][-1], 12*2));
        plt.ylabel("CT outflow (Tg C month$^{-1}$)", fontsize=fontSize);
        plt.tick_params(labelsize=fontSize);
        # plt.ylim(0, 20);
        
        # ## Add second axis with yearly values
        # ax = plt.twinx();
        # years = [datetime.datetime(y, 7, 1) for y in radiiDataYears.index]
        # ax.plot(years, radiiDataYearsMean_amp, 'r-', linewidth=2, label="annual mean (peak alignment)"); #do as scatter?
        # ax.fill_between(years, radiiDataYearsMean_amp+radiiDataYearsSD_amp, radiiDataYearsMean_amp-radiiDataYearsSD_amp, color='r', alpha=0.2, linewidth=0);#, hatch="x");
        # #error bar for radiiDataYearsSD_amp
        # ax.plot(years, radiiDataYearsMean_corr, 'b-', linewidth=2, label="annual mean (maximum correlation)"); #do as scatter?
        # ax.fill_between(years, radiiDataYearsMean_corr+radiiDataYearsSD_corr, radiiDataYearsMean_corr-radiiDataYearsSD_corr, color='b', alpha=0.2, linewidth=0);#, hatch="x");
        # ax.set_ylabel("CT outflow (Tg C year$^{-1}$)");
        # ax.set_ylim(0, 60);
        
        plt.legend(loc=0, fontsize=legendSize);
        plt.tight_layout();
        plt.savefig(path.join(plotOutputPath, region+"_monthly.pdf"));
        plt.savefig(path.join(plotOutputPath, region+"_monthly.png"));
        #plt.close();
    
        
        ################Plot yearly time series
        plt.figure(figsize=figureSize);
        years = [datetime.datetime(y, 7, 1) for y in radiiDataYears.index]
        plt.plot(years, radiiDataYearsMean_amp, 'r', label="annual mean (peak alignment)"); #do as scatter?
        plt.fill_between(years, radiiDataYearsMean_amp+radiiDataYearsSD_amp, radiiDataYearsMean_amp-radiiDataYearsSD_amp, color='r', alpha=0.4, linewidth=0);
        #error bar for radiiDataYearsSD_amp
        plt.plot(years, radiiDataYearsMean_corr, 'b', label="annual mean (maximum correlation)"); #do as scatter?
        plt.fill_between(years, radiiDataYearsMean_corr+radiiDataYearsSD_corr, radiiDataYearsMean_corr-radiiDataYearsSD_corr, color='b', alpha=0.4, linewidth=0);
        
        #wHasData = np.where(np.isfinite(radiiSetMean_amp));
        #plt.xticks(np.arange(wHasData[0][0], wHasData[0][-1], 12*2));
        plt.ylabel("CT outflow (Tg C yr$^{-1}$)", fontsize=fontSize);
        plt.legend(loc=0, fontsize=legendSize);
        plt.tick_params(labelsize=fontSize);
        plt.tight_layout();
        plt.savefig(path.join(plotOutputPath, region+"_yearly.pdf"));
        plt.savefig(path.join(plotOutputPath, region+"_yearly.png"));
        
    #Output some summary stats at the end
    for region in regions:
        print(region);
        print("\tMean annual DIC outflow (amp): ", summaryData[region+"_annual_DIC_outflow_amp"], "+/-", summaryData[region+"_annual_DIC_outflow_uncert_amp"]);
        print("\tnumRadii (amp):", summaryData[region+"_num_radii_amp"], "distance between peaks:", summaryData[region+"_distance_between_peaks"]);
        print("\tYearly mean SD (amp): ", summaryData[region+"_annual_DIC_outflow_SD_amp"]);
        print("\tYearly mean CoV (amp): ", summaryData[region+"_annual_DIC_outflow_CoV_amp"]); #coefficient of variation
        print("");
        print("\tMean annual DIC outflow (corr): ", summaryData[region+"_annual_DIC_outflow_corr"], "+/-", summaryData[region+"_annual_DIC_outflow_uncert_corr"]);
        print("\tnumRadii (corr):", summaryData[region+"_num_radii_corr"], "correlation coefficient:", summaryData[region+"_correlation_coefficient"]);
        print("\tYearly mean SD (corr): ", summaryData[region+"_annual_DIC_outflow_SD_corr"]);
        print("\tYearly mean CoV (corr): ", summaryData[region+"_annual_DIC_outflow_CoV_corr"]); #coefficient of variation
