#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#####
#
# This version uses a band of 'virtual stations' to estimate DIC outflow.
#

from netCDF4 import Dataset;
import matplotlib.pyplot as plt;
import numpy as np;
import pandas as pd;
from os import path, makedirs;
from datetime import datetime, timedelta;
from skimage.draw import circle_perimeter;

import os_dic_outflow.analysis_tools as analysis_tools;
import os_dic_outflow.preprocess_discharge_data as preprocess_discharge_data;


###Create a circle perimeter that acts as our bondary for counting DIC that 'leaves' the river mouth region
def create_perimeter_mask(shape, perimeterCentre, perimeterRadius):
    #Get coordinates of a circular perimeter
    circleMask = np.zeros(shape, dtype=int);
    cx, cy = circle_perimeter(perimeterCentre[0], perimeterCentre[1], perimeterRadius);
    circleMask[cx, cy] = 1;
    return circleMask;

#assumes negative lat is top of image (i.e. upside down globe)
def latlon_to_grid_indices(lat, lon, latRes, lonRes):
    ilat = int(90.0//latRes)+int(lat//latRes);
    ilon = int(180.0//lonRes)+int(lon//lonRes);
    return ilat, ilon;


#simple plot used to help debug / verify correct running
def debug_plot(savePath, sss, perimeterMask, plumeMask, dateStr, region, closeAfter=True):
    if region == "oceansoda_amazon_plume":
        xlims = (104, 154);
        ylims = (60, 100);
    elif region == "oceansoda_congo":
        xlims = (176, 196);
        ylims = (85, 100);
    elif region == "oceansoda_mississippi":
        xlims = (80, 99);
        ylims = (57, 66);
    else:
        xlims = (0, 360);
        ylims = (0, 180);
    
    plt.subplots(2, 1, figsize=(10,20));
    plt.subplot(2,1,1);
    plt.title(dateStr+"\nplume highlighted");
    plt.imshow(np.flipud(sss-(perimeterMask*20)));
    plt.colorbar();
    plt.xlim(xlims[0], xlims[1]);
    plt.ylim(ylims[0], ylims[1]);
    plt.clim(25, 38);
    
    plt.subplot(2,1,2);
    plt.title("Plume mask+perimeter")
    tmpMatrix = plumeMask.astype(float);
    tmpMatrix[np.where(np.isfinite(sss)==False)] = np.nan;
    tmpMatrix+=(perimeterMask*2);
    plt.imshow(np.flipud(tmpMatrix));
    #plt.colorbar();
    plt.xlim(xlims[0], xlims[1]);
    plt.ylim(ylims[0], ylims[1]);
 
    plt.savefig(savePath);
    if closeAfter:
        plt.close();



# Calculates DIC outflow from a river using circulate virtual transects for a list of radii.
#carbonateParametersTemplate: Template (REGION, LATRES, LONRES, OUTPUTVAR) containing gridded carbonate parameter time series data
#outputDirectoryTemplate: Template (region) for root output path (analysis output and plots written here)
#regionMaskPath: path to OceanSODA region mask definition netCDF file
#gridAreasPath: path to pre-computed lon/lat grid areas
#perimeterRadii: list / iterator of radii to use to create perimeters / virtual transects
#numSamples: How many SSS and DIC samples to use (for uncertainty analysis)
def calculate_dic_outflow_from_circle_transects(carbonateParametersTemplate, outputDirectoryTemplate, regions, regionMaskPath, gridAreasPath, perimeterRadii=range(1, 25), numSamples=100, verbose=True):
    regionMaskNC = Dataset(regionMaskPath, 'r');
    gridAreas = np.genfromtxt(gridAreasPath, delimiter=","); #manually read as newer pyproj isn't compatible with qt5 dependencies
    
    #Convenient definitions
    plumeSalinityThreshold = 35.0; #Below this salinity is considered plume
    lonRes = latRes = 1;
    makePlots = True; #Create and save plots?
    closePlots = True; #Close plots automatically to avoid cluttering screen
    #Serves as the centre point for circle perimeters
    riverMouthCoords = {"oceansoda_amazon_plume": latlon_to_grid_indices(0.0, -50.0, latRes, lonRes),
                        "oceansoda_congo": latlon_to_grid_indices(-6.0, 12.0, latRes, lonRes),
                        "oceansoda_mississippi": latlon_to_grid_indices(29.0, -89.0, latRes, lonRes),
                        "oceansoda_st_lawrence": latlon_to_grid_indices(49.0, -61.0, latRes, lonRes),
                        };
    
    zeroPlumes=[];
    skipped = [];
    for region in regions: #["oceansoda_amazon_plume"]:#regions:
        #Setup region specific output directories
        baseOutputPath = outputDirectoryTemplate.safe_substitute(REGION=region);
        if path.exists(baseOutputPath) == False:
            makedirs(baseOutputPath);
        basePlotPath = baseOutputPath;
        
        #Read input data sets
        #Carbonate parameters
        try:
            carbonatePath = carbonateParametersTemplate.safe_substitute(REGION=region, LATRES=latRes, LONRES=lonRes, OUTPUTVAR="DIC");
            carbonateParameters = Dataset(carbonatePath);
        except FileNotFoundError:
            print("Skipping", region, "because no gridded carbonate parameter file found for it.");
            skipped.append((region, carbonatePath));
            continue;
        
        #mask data
        regionMask = regionMaskNC.variables[region][:];

    #### create_perimeter_radii_masks
        #Virtual perimeter for measuring plume DIC throughput:
        #Create a circle perimeter that acts as our bondary for counting DIC that 'leaves' the river mouth region
        perimeterMasks = {perimeterRadius: create_perimeter_mask((int(180/latRes), int(360/lonRes)), riverMouthCoords[region], perimeterRadius) for perimeterRadius in perimeterRadii};
        
        #Store these dictionaries at pickle txt files for plotting
        #They contain the final scores data as tables
        import pickle
        with open("output/Amazon_radii.txt", "wb") as myFile:
            pickle.dump(perimeterMasks, myFile)
        
    #### Time data
        numTimePoints = len(carbonateParameters.variables["time"]);
        carbonateParamsTimeRange = (datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][0])), datetime(1980, 1, 1) + timedelta(seconds=float(carbonateParameters.variables["time"][-1])));
        carbonateDates = analysis_tools.create_date_month_array(carbonateParamsTimeRange[0], carbonateParamsTimeRange[1]);
        
    #### River discharge data
        monthlyDischarge = preprocess_discharge_data.load_discharge_data(carbonateDates, region);
        
        
    #### Create netCDF file to store results - this will be updated with each iteration
        outputPathNetCDF = path.join(baseOutputPath, "dic_outflow_"+region+".nc");
        ncout = analysis_tools.create_netCDF_file(outputPathNetCDF, carbonateParameters, numSamples, regionMask, gridAreas);
        
        
        #Create a data frame to contain all the values used in the calculation
    #### River discharge. carbonate data is monthly averages, so calculate monthly means of river discharge.
        monthlyDF = pd.DataFrame();
        monthlyDF["date"] = carbonateDates;
        monthlyDF.index = carbonateDates;
        monthlyDF["plume_surface_area"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_surface_area_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_mean_thickness"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_mean_thickness_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_volume"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_volume_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_total_dic"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_total_dic_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_total_riverine_dic"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["plume_total_riverine_dic_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["river_discharge"] = monthlyDischarge["monthly_discharge"];
        monthlyDF["river_discharge_sd"] = monthlyDischarge["monthly_discharge_sd"];
        #Add individual values for each radius
        for perimeterRadius in perimeterRadii:
            monthlyDF["dic_outflow_radius"+str(perimeterRadius)] = np.full((len(carbonateDates)), np.nan, dtype=float);
            monthlyDF["dic_outflow_sd_radius"+str(perimeterRadius)] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow_mean"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow_mean_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow_median"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        monthlyDF["dic_outflow_most_points"] = np.full((len(carbonateDates)), 0, dtype=int);
        monthlyDF["dic_outflow_most_points_estimate"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        #monthlyDF["dic_outflow_median_sd"] = np.full((len(carbonateDates)), np.nan, dtype=float);
        
        
    #### loop through months - start
            # Split calculations by time point because memory usage is too large to store many samples for all time points at once
        mostPointsRadii = [];
        for t in range(0, numTimePoints):
            currentDatetime = carbonateDates[t];
            if verbose:
                print(region, t+1, "of", numTimePoints, "("+str(currentDatetime)+")");
            
            #### extract data for this time point #SSS in PSU, DIC in umol kg-1
            dic, DIC_Combined_uncertainty_dueto_RMSD_and_Bias, sss, SSS_uncertainty = analysis_tools.extract_data(carbonateParameters, regionMask, t=t);
            
            #Don't perform computations if everything is nan
            if np.all(np.isfinite(dic)==False):
                continue;
            
            
            #Generate DIC and SSS samplestime point
            sssSamples = np.random.normal(sss, SSS_uncertainty, size=(numSamples,)+sss.shape);
            dicSamples = np.random.normal(dic, DIC_Combined_uncertainty_dueto_RMSD_and_Bias, size=(numSamples,)+dic.shape);
            

            #### Calculate plume masks (1==inside plume)
            plumeMask = analysis_tools.calculate_plume_mask(sss, plumeSalinityThreshold=plumeSalinityThreshold);
            plumeMaskSamples = analysis_tools.calculate_plume_mask(sssSamples, plumeSalinityThreshold=plumeSalinityThreshold);
                
            if np.sum(plumeMask) == 0:
                zeroPlumes.append(carbonateDates[t]);
            
            #### Calculate plume surface area (m^2)
            surfaceArea, griddedSurfaceArea = analysis_tools.calculate_plume_surface_area(plumeMask, gridAreas); #surfaceArea in m^2
            surfaceAreaSamples, griddedSurfaceAreaSamples = analysis_tools.calculate_plume_surface_area(plumeMaskSamples, gridAreas, calculateForSamples=True);
            monthlyDF.loc[currentDatetime, "plume_surface_area"] = surfaceArea;
            monthlyDF.loc[currentDatetime, "plume_surface_area_sd"] = np.std(surfaceAreaSamples);
            
            
            #### Calculate plume thickness(m)
            meanPlumeThickness, griddedPlumeThickness = analysis_tools.plume_thickness_coles2013(sss, plumeMask, useModelled=False, sampleUncertainty=False); #thickness in metres (m)
            meanPlumeThicknessSamples, griddedPlumeThicknessSamples = analysis_tools.plume_thickness_coles2013(sssSamples, plumeMaskSamples, useModelled=False, sampleUncertainty=True);
            monthlyDF.loc[currentDatetime, "plume_mean_thickness"] = meanPlumeThickness;
            monthlyDF.loc[currentDatetime, "plume_mean_thickness_sd"] = np.std(meanPlumeThicknessSamples);
            
            
            #### Calculate plume volume (m^3)
            plumeVolume, griddedPlumeVolume = analysis_tools.calculate_plume_volume(griddedSurfaceArea, griddedPlumeThickness, calculateForSamples=False); #volumes in m^3
            plumeVolumeSamples, griddedPlumeVolumeSamples = analysis_tools.calculate_plume_volume(griddedSurfaceAreaSamples, griddedPlumeThicknessSamples, calculateForSamples=True);
            monthlyDF.loc[currentDatetime, "plume_volume"] = plumeVolume;
            monthlyDF.loc[currentDatetime, "plume_volume_sd"] = np.std(plumeVolumeSamples);
            
            
            #### Calculate gridded mean plume DIC, in the process also calculate meanplume salinity (umol kg-1)
            #This uses relationship between mean plume salinity and sss in Hu 2004 (PS = 0.881*SSS+4.352, see section 3.4 of https://doi.org/10.1016/j.dsr2.2004.04.001)
            #and assumes DIC is perfectly conservative with salinity
            #I.e. the proportion mean plume salinity to surface salinity is assumed to equal mean plume DIC / surface DIC, thus rearranging gives:
            #   meanPlumeDIC = surfaceDIC * (meanPlumeSalinity / surfaceSalinity)
            #mean plume DIC given in umol kg-1, mean plume SSS given in PSU
            griddedMeanDICConc, griddedMeanSSS = analysis_tools.calculate_mean_dic_sss(sss, dic, plumeMask, interceptUncertaintyRatio=0.0, slopeUncertaintyRatio=0.0, calculateForSamples=False); #DIC in umol kg-1, salinity in PSU
            griddedMeanDICConcSamples, griddedMeanSSSSamples = analysis_tools.calculate_mean_dic_sss(sssSamples, dicSamples, plumeMaskSamples, interceptUncertaintyRatio=0.1, slopeUncertaintyRatio=0.1, calculateForSamples=True);
            
            
            #### calculate total plume DIC (mols C)
            totalDIC, griddedTotalDIC = analysis_tools.calculate_total_plume_dic(griddedMeanDICConc, griddedPlumeVolume, calculateForSamples=False); #DIC in mol
            totalDICSamples, griddedTotalDICSamples = analysis_tools.calculate_total_plume_dic(griddedMeanDICConcSamples, griddedPlumeVolumeSamples, calculateForSamples=True);
            monthlyDF.loc[currentDatetime, "plume_total_dic"] = totalDIC;
            monthlyDF.loc[currentDatetime, "plume_total_dic_sd"] = np.std(totalDICSamples);
            
            
            #### calculate river originating DIC estimate within the plume (mols C)
            #Assumes a linear interpolation between mean salinity of 0 (100% of DIC originates from river) to salinity 35 (0% DIC from river)
            totalRiverineDIC, griddedTotalRiverineDIC = analysis_tools.interpolate_riverine_dic(griddedTotalDIC, griddedMeanSSS, calculateForSamples=False); #mols
            totalRiverineDICSamples, griddedTotalRiverineDICSamples = analysis_tools.interpolate_riverine_dic(griddedTotalDICSamples, griddedMeanSSSSamples, calculateForSamples=True); #mols
            monthlyDF.loc[currentDatetime, "plume_total_riverine_dic"] = totalRiverineDIC;
            monthlyDF.loc[currentDatetime, "plume_total_riverine_dic_sd"] = np.std(totalRiverineDICSamples);
            
            
            #### loop over each perimeter radius and create a composite estimate
            dicOutflowEstimates = []; #one for each perimeter radius
            dicOutflowSDEstimates = []; #one for each perimeter radius
            mostPoints = 0;
            mostPointsRadius = None;
            for perimeterRadius in perimeterRadii:
                perimeterMask = perimeterMasks[perimeterRadius];
                #calculate the amount of riverine DIC on the perimeter
                plumeOnPerimeterMask = np.zeros(sss.shape, dtype=int);
                plumeOnPerimeterMask[np.where(plumeMask & perimeterMask)] = 1; #grid cells which are in the plume and on perimeter = 1
                
                totalPerimeterDIC = np.nansum(griddedTotalDIC[plumeOnPerimeterMask==1]);
                #again with samples
                plumeOnPerimeterMaskSamples = np.zeros(sssSamples.shape, dtype=int);
                wPlumeOnPerimeterMaskLocations = np.where(plumeMaskSamples & np.broadcast_to(perimeterMask, plumeMaskSamples.shape));
                plumeOnPerimeterMaskSamples[wPlumeOnPerimeterMaskLocations] = 1; #grid cells which are in the plume and on perimeter = 1
                sampleGroups = wPlumeOnPerimeterMaskLocations[0]; #this is the indices of dimension 0 (i.e. which SSS sample it corresponds to)
                totalPerimeterDICSamples = np.bincount(sampleGroups, weights=griddedTotalDICSamples[plumeOnPerimeterMaskSamples==1]); #sum within SSS sample groups. Equivalent to nansum(a, axis=0) but we couldn't use that because we indexed with np.where which resulted in a 1d array
                
                
                #### get perimeter DIC concentrations
                perimeterDICConc = np.full(sss.shape, np.nan, dtype=float);
                perimeterDICConc[plumeOnPerimeterMask==1] = griddedMeanDICConc[plumeOnPerimeterMask==1]; #umol kg-1
                #get just the riverine concentration
                #Note that summing the riverinePerimeterDICConcs works, because the linear interpolation takes away the ocean component of the concentration, leaving only riverine concentration that has been diluted by the volume of the plume in the grid cell.
                #Summing thus integrates over volume (or, assuming the perimeter is of equal age, back through time to the point at which it was discharged from the river mouth)
                #Resulting DIC concentration units: umol kg-1
                sumTotalRiverinePerimeterDICConc, riverinePerimeterDICConc = analysis_tools.interpolate_riverine_dic(perimeterDICConc, griddedMeanSSS, calculateForSamples=False);
                #repeat the process for samples
                perimeterDICConcSamples = np.full(sssSamples.shape, np.nan, dtype=float);
                perimeterDICConcSamples[plumeOnPerimeterMaskSamples==1] = griddedMeanDICConcSamples[plumeOnPerimeterMaskSamples==1]; #umol kg-1
                sumTotalRiverinePerimeterDICConcSamples, riverinePerimeterDICConcSamples = analysis_tools.interpolate_riverine_dic(perimeterDICConcSamples, griddedMeanSSSSamples, calculateForSamples=True);
                
                
                #### Calculate the DIC otuflow from the river
                monthlyDischargeKg = monthlyDischarge.loc[currentDatetime]["monthly_discharge"] * 1000; #Convert discharge from m^3 month-1 to kg month-1
                dicOutflow = sumTotalRiverinePerimeterDICConc * monthlyDischargeKg; #concentration: umol month-1
                dicOutflowMol = dicOutflow/1000000; #umol month-1 to mol month-1
                dicOutflowTGrams = analysis_tools.mol_to_TgC(dicOutflowMol);
                #repeat with samples
                monthlyDischargeKgSamples = np.random.normal(monthlyDischarge.loc[currentDatetime]["monthly_discharge"], monthlyDischarge.loc[currentDatetime]["monthly_discharge_sd"], size=numSamples) * 1000; #Multiply by 1000 to convert discharge from m^3 month-1 to kg month-1
                dicOutflowSamples = sumTotalRiverinePerimeterDICConcSamples * monthlyDischargeKgSamples; #concentration: umol month-1
                dicOutflowMolSamples = dicOutflowSamples/1000000; #umol month-1 to mol month-1
                dicOutflowTGramsSamples = analysis_tools.mol_to_TgC(dicOutflowMolSamples);
                
                #### Store the estimates for this perimeter radius
                monthlyDF.loc[currentDatetime, "dic_outflow_radius"+str(perimeterRadius)] = dicOutflowTGrams;
                monthlyDF.loc[currentDatetime, "dic_outflow_sd_radius"+str(perimeterRadius)] = np.std(dicOutflowTGramsSamples);
                dicOutflowEstimates.append(dicOutflowTGrams);
                dicOutflowSDEstimates.append(np.std(dicOutflowTGramsSamples));
                
                numPoints = np.nansum(plumeOnPerimeterMask);
                if numPoints > mostPoints:
                    mostPoints = numPoints;
                    mostPointsRadius = perimeterRadius;
                    mostPointsEstimate = dicOutflowTGrams;
    

            
            #### Calculate a composite permeter radius
            #remove any 0s
            dicOutflowEstimates = np.array([val for val in dicOutflowEstimates if (val != 0) and (np.isfinite(val))]);
            dicOutflowSDEstimates = np.array([val for val in dicOutflowSDEstimates if (val != 0) and (np.isfinite(val))]);
            #Means
            monthlyDF.loc[currentDatetime, "dic_outflow_mean"] = np.nanmean(dicOutflowEstimates);
            dicOutflowMeanSD = np.sqrt(np.nansum(dicOutflowSDEstimates**2))/np.sum(np.isfinite(dicOutflowSDEstimates));
            monthlyDF.loc[currentDatetime, "dic_outflow_mean_sd"] = dicOutflowMeanSD;
            #Medians
            monthlyDF.loc[currentDatetime, "dic_outflow_median"] = np.nanmedian(dicOutflowEstimates);
            #Number of points found/used
            monthlyDF.loc[currentDatetime, "dic_outflow_most_points"] = mostPoints;
            monthlyDF.loc[currentDatetime, "dic_outflow_most_points_estimate"] = mostPointsEstimate;
            mostPointsRadii.append(mostPointsRadius);
            
            
            #### Update output netCDF file for this time point
            griddedPlumeVolumeSD = np.nanstd(griddedPlumeVolumeSamples, axis=0);
            griddedTotalDICSD = np.nanstd(griddedTotalDICSamples, axis=0);
            analysis_tools.update_gridded_time_point_netCDF(ncout, t, griddedPlumeVolume, griddedPlumeVolumeSD, griddedTotalDIC, griddedTotalDICSD, plumeMask);
    
        
            
        #### Calculate inter-annual values for each month        
        interyearDF = analysis_tools.calculate_inter_year_monthly_means(monthlyDF);
        
        #### calculate mean annual values
        annualDF = analysis_tools.calculate_annual_values(interyearDF);
        
        
        #### Write time series to netCDF file
        analysis_tools.write_timeseries_to_netCDF(ncout, monthlyDF["plume_volume"].values, monthlyDF["plume_volume_sd"].values,
                                   monthlyDF["plume_total_dic"].values, monthlyDF["plume_total_dic_sd"].values,
                                   monthlyDF["plume_total_riverine_dic"].values, monthlyDF["plume_total_riverine_dic_sd"].values,
                                   monthlyDF["dic_outflow"].values, monthlyDF["dic_outflow_sd"].values,
                                   monthlyDF["river_discharge"].values, monthlyDF["river_discharge_sd"].values,
                                   interyearDF["dic_outflow"].values, interyearDF["dic_outflow_sd"].values);
        print("Written netCDF file output to: " + baseOutputPath);
        
        ncout.close();
        
        #Write csv files
        if path.exists(baseOutputPath) == False:
            makedirs(baseOutputPath);
        monthlyPath = path.join(baseOutputPath, "monthly_timeseries_"+region+".csv");
        monthlyDF.to_csv(monthlyPath, sep=",");
        interyearPath = path.join(baseOutputPath, "interyear_timeseries_"+region+".csv");
        interyearDF.to_csv(interyearPath, sep=",");
        annualPath = path.join(baseOutputPath, "annual_means_"+region+".csv");
        annualDF.insert(0, "region", region);
        annualDF.to_csv(annualPath, sep=",");
        print("Written csv file output to:" + baseOutputPath);
    
        
        #### plotting
        if makePlots == True:
            #Visualise plume volume area over time
            outputPath = path.join(basePlotPath, "monthly_plume_volume.pdf");
            analysis_tools.plot_with_uncertainty(monthlyDF, "plume_volume", "date", r"plume volume ($m^3$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise total plume DIC content over time
            outputPath = path.join(basePlotPath, "monthly_plume_dic.pdf");
            analysis_tools.plot_with_uncertainty(monthlyDF, "plume_total_dic", "date", r"DIC content of plume ($mols C$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise total plume DIC content over time
            outputPath = path.join(basePlotPath, "monthly_plume_riverine_dic.pdf");
            analysis_tools.plot_with_uncertainty(monthlyDF, "plume_total_riverine_dic", "date", r"riverine DIC content of plume ($mols C$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise DIC discharge over time
            outputPath = path.join(basePlotPath, "monthly_dic_discharge.pdf");
            analysis_tools.plot_with_uncertainty(monthlyDF, "dic_outflow", "date", r"riverine DIC discharge ($TgC$)", x=monthlyDF["date"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise mean (inter-year) monthly plume volume
            outputPath = path.join(basePlotPath, "interyear_plume_volume.pdf");
            analysis_tools.plot_with_uncertainty(interyearDF, "plume_volume", "", r"mean plume volume ($m^3$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise mean (inter-year) monthly total plume DIC
            outputPath = path.join(basePlotPath, "interyear_plume_dic.pdf");
            analysis_tools.plot_with_uncertainty(interyearDF, "plume_total_dic", "", r"mean DIC content of plume ($mols C$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
            
            #Visualise mean (inter-year) monthly DIC discharge
            outputPath = path.join(basePlotPath, "interyear_dic_discharge.pdf");
            analysis_tools.plot_with_uncertainty(interyearDF, "dic_outflow", "", r"mean riverine DIC discharge ($TgC$)", x=interyearDF["month_name"], outputPath=outputPath, closeAfter=closePlots);
            
    
    if len(skipped) != 0:
        print("Some regions were skipped due to missing gridded carbonate parameter files:");
        for tup in skipped:
            print("\t"+tup[0]+":", tup[1]);



# plt.figure();
# plt.plot(monthlyDF["dic_outflow_most_points"]); plt.title("dic_outflow_most_points");
# plt.figure();
# plt.plot(monthlyDF["dic_outflow_most_points_estimate"]); plt.title("dic_outflow_most_points_estimate");


