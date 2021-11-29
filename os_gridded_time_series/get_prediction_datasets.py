#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 00:23:35 2020

Downloads and processes all prediction data sets for calculating
gridded time series

@author: Tom Holding
"""

import argparse;
from string import Template;
from os import path, makedirs;
import numpy as np;
from netCDF4 import Dataset;
import shutil;
import urllib.request as request;
from contextlib import closing;
import ssl;
#import urllib;
import calendar;
from datetime import datetime, timedelta;
import tarfile;

#For downloads requiring authentication use requests directly.
import requests as rrequests;
#Disable insecure ssl warnings
from requests.packages.urllib3.exceptions import InsecureRequestWarning
rrequests.packages.urllib3.disable_warnings(InsecureRequestWarning)

from ftplib import FTP; #for FTP servers requiring authentication

def download_all_woa_nutrients(destinationRoot):
    def _do_woa_download(urlTemplate, destinationTemplate, name):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            print("Downloading WOA", name, "for month:", monthStr+"...");
        
            #generate source url and destination path
            url = urlTemplate.safe_substitute(MM=monthStr);
            destination = destinationTemplate.safe_substitute(MM=monthStr);
            
            #Make directory if it doesn't exist
            if path.exists(path.dirname(destination)) == False:
                makedirs(path.dirname(destination));
            
            #Do not overwrite any data
            if path.exists(destination) == True:
                print("Skipping", name, "month:", monthStr, "to prevent overwriting existing file.");
                continue;
    
            #Try the download
            try:
                with closing(request.urlopen(url)) as req:
                    with open(destination, 'wb') as dest:
                        shutil.copyfileobj(req, dest);

            except Exception as e:
                print("Error downloading", name, monthStr);
                print(e);
    
    _do_woa_download(Template("ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/WOA18/DATA/nitrate/netcdf/all/1.00/woa18_all_n${MM}_01.nc"),
                     Template(path.join(destinationRoot, "WOA_nitrate/woa18_all_n${MM}_01.nc")),
                     "nitrate");
    
    _do_woa_download(Template("ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/WOA18/DATA/phosphate/netcdf/all/1.00/woa18_all_p${MM}_01.nc"),
                     Template(path.join(destinationRoot, "WOA_phosphate/woa18_all_p${MM}_01.nc")),
                     "phosphate");
    
    _do_woa_download(Template("ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/WOA18/DATA/silicate/netcdf/all/1.00/woa18_all_i${MM}_01.nc"),
                     Template(path.join(destinationRoot, "WOA_silicate/woa18_all_i${MM}_01.nc")),
                     "silicate");
    
    _do_woa_download(Template("ftp://ftp.nodc.noaa.gov/pub/data.nodc/woa/WOA18/DATA/oxygen/netcdf/all/1.00/woa18_all_o${MM}_01.nc"),
                     Template(path.join(destinationRoot, "WOA_dissolved_oxygen/woa18_all_o${MM}_01.nc")),
                     "dissolved oxygen");


def process_all_woa_nutrients(dataRoot, downloadedRoot, oceanMaskPath):
    def get_uncertainty_value(varName, oceanName):
        #Data from  and https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol4.pdf
        uncertainties = {"dissolved_oxygen": {"atlantic": 9.8, "pacific": 11.5, "indian": 10.6}, #table 6 of https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol3.pdf
                         "nitrate": {"atlantic": 1.6, "pacific": 1.9, "indian": 1.7}, #table 6 of https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol4.pdf
                         "phosphate": {"atlantic": 0.11, "pacific": 0.13, "indian": 0.12}, #table 6 of https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol4.pdf
                         "silicate": {"atlantic": 3.1, "pacific": 3.5, "indian": 13.0}}; #table 6 of https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol4.pdf
        
        return uncertainties[varName][oceanName]


    downloadedFileTemplate = Template(path.join(downloadedRoot, "WOA_${VARNAME}/woa18_all_${VARCODE}${MM}_01.nc"));
    processedFileTemplate = Template(path.join(dataRoot, "WOA_${VARNAME}/woa18_all_${VARCODE}${MM}_processed.nc"));
    
    varNameCodeTuples = [("dissolved_oxygen", "o"),
                         ("nitrate", "n"),
                         ("phosphate", "p"),
                         ("silicate", "i")];

    #read ocean mask file
    oceanMask = np.flipud(Dataset(oceanMaskPath, 'r').variables["sea-mask"][0,:,:]);

    #ocean mask values
    ATLANTIC=30;
    PACIFIC=70;
    INDIAN=50;
    
    
    for varName, varCode in varNameCodeTuples:
        print("Processing WOA", varName);
        
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            downloadedFilePath = downloadedFileTemplate.safe_substitute(VARNAME=varName, VARCODE=varCode, MM=monthStr);
            
            #check path exists, if not, make it
            processedFilePath = processedFileTemplate.safe_substitute(VARNAME=varName, VARCODE=varCode, MM=monthStr);
            if path.exists(path.dirname(processedFilePath)) == False:
                makedirs(path.dirname(processedFilePath));
            
            #check file doesn't already exist, if so, report error
            if path.exists(processedFilePath) == True:
                print("File already exists for WOA {0} {1}. Skipping to avoid overwriting existing file.".format(varName, monthStr));
                continue;
    
            #copy downloaded file and open it to append uncertainty data to
            shutil.copyfile(downloadedFilePath, processedFilePath)
            nc = Dataset(processedFilePath, 'a');
            
            #Calculate gridded uncertainty
            griddedUncertainty = np.full((180, 360), np.nan);
            griddedUncertainty[oceanMask==ATLANTIC] = get_uncertainty_value(varName, "atlantic");
            griddedUncertainty[oceanMask==PACIFIC] = get_uncertainty_value(varName, "pacific");
            griddedUncertainty[oceanMask==INDIAN] = get_uncertainty_value(varName, "indian");
            
            #Add the uncertainty data to the NC file
            var = nc.createVariable(varCode+"_uncertainty", float, ("time", "lat", "lon"));
            var.long_name = "Basin scale surface (0-500 metre) uncertainty as reported by Table 6 in https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol3.pdf and https://data.nodc.noaa.gov/woa/WOA18/DOC/woa18_vol4.pdf"
            var.units = "micromoles_per_kilogram"
            var.grid_mapping = "crs";
            var.coordinates = "time lat lon";
            var[:] = griddedUncertainty;
            
            #close the file
            nc.close();


def download_rss_smap_sss(destinationRoot, startYear, endYear):
    urlTemplate = Template("https://podaac-tools.jpl.nasa.gov/drive/files/allData/smap/L3/RSS/V4/monthly/SCI/${YYYY}/RSS_smap_SSS_L3_monthly_${YYYY}_${MM}_FNL_v04.0.nc");
    destinationTemplate = Template(path.join(destinationRoot, "RSS_SMAP_SSS/RSS_smap_SSS_L3_monthly_${YYYY}_${MM}_FNL_v04.0.nc"));
    
    #Check start year makes sense
    if startYear < 2015:
        print("Starting RSS SMAP SSS download at 2015 because there is no data before this year.");
        startYear = 2015;
    
    for year in range(startYear, endYear+1):
        for imonth in range(0, 12):
            if (year == 2015) & (imonth < 3):
                continue; #no data for Jan-Mar in 2015.
            
            monthStr = format(imonth+1, "02d");
            print("Downloading", year, monthStr+"...");
            
            url = urlTemplate.safe_substitute(YYYY=year, MM=monthStr);
            destination = destinationTemplate.safe_substitute(YYYY=year, MM=monthStr);
            
            #Create directory if it doesn't exist
            if path.exists(path.dirname(destination)) == False:
                makedirs(path.dirname(destination));
            
            #If file already exists, don't overwrite it
            if path.exists(destination) == True:
                print("Skipping download of RSS SMAP SSS for", year, monthStr, "to avoid overwriting existing file.");
                continue;
            
            #Try the download
            try:
                request = rrequests.get(url, auth=("tholding", "OAUlJNofE8plHl9ql4f"), verify=False, stream=True);
                request.raw.decode_content = True;
                with open(destination, 'wb') as f:
                    shutil.copyfileobj(request.raw, f);
            except Exception as e:
                print("Error downloading RSS SMAP SSS:", year, monthStr);
                print(e);

def process_rss_smap_sss(dataRoot, downloadedRoot, startYear, stopYear):
    #Rebin 2D data to a specified shape.
    #Converts to a 4D array, then takes the mean in first and last dimension to give newly binned array.
    def rebin(m, shape):
        sh = (shape[0], m.shape[0]//shape[0], shape[1], m.shape[1]//shape[1]);
        return m.reshape(sh).mean(-1).mean(1);
    
    def rebin_sum(m, shape):
        sh = (shape[0], m.shape[0]//shape[0], shape[1], m.shape[1]//shape[1]);
        return m.reshape(sh).sum(-1).sum(1);
    
    #No data before 2015
    if startYear < 2015:
        startYear = 2015;
    
    
    downloadedTemplate = Template(path.join(downloadedRoot, "RSS_SMAP_SSS/RSS_smap_SSS_L3_monthly_${YYYY}_${MM}_FNL_v04.0.nc"));
    processedTemplate = Template(path.join(dataRoot, "RSS_SMAP_SSS/RSS_smap_SSS_L3_monthly_${YYYY}_${MM}_FNL_v04.0_processed.nc"));
    outputRes = 1.0;
    
    
    skipped = [];
    for year in range(startYear, stopYear+1):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            try:
                inputNC = Dataset(downloadedTemplate.safe_substitute(YYYY=year, MM=monthStr), 'r');
            except FileNotFoundError:
                skipped.append((year, monthStr));
                continue;

            #make directory if it doesn't exist
            outputPath = processedTemplate.safe_substitute(YYYY=year, MM=monthStr);
            if path.exists(path.dirname(outputPath)) == False:
                makedirs(path.dirname(outputPath));
            if path.exists(outputPath) == True:
                skipped.append((year, monthStr));
                continue;
            
            print("Processing RSS SMAP SSS data for ", year, monthStr);
                
            sss = inputNC["sss_smap"][:];
            sssErr = inputNC["sss_smap_uncertainty"][:];
            
            #resample
            sssRebinned = rebin(sss, (int(180/outputRes), int(360/outputRes)));
            sssCounts = rebin_sum(sss.mask==False, (int(180/outputRes), int(360/outputRes)));
            sssErrs = rebin_sum(sssErr**2, (int(180/outputRes), int(360/outputRes))); #sum of squares
            sssErrs = np.sqrt(sssErrs) / sssCounts; #square root and divide by n
            
            #roll by 180 degrees longitude
            sssRebinned = np.roll(sssRebinned, int(180/outputRes), axis=1);
            sssCounts = np.roll(sssCounts, int(180/outputRes), axis=1);
            sssErrs = np.roll(sssErrs, int(180/outputRes), axis=1);
            
            ###Write output file
            ncout = Dataset(outputPath, 'w');
            
            ncout.createDimension("lat", int(180/outputRes));
            ncout.createDimension("lon", int(360/outputRes));
            
            #dimension variables
            var = ncout.createVariable("lat", float, ("lat",));
            var.units = "lat (degrees North)";
            var[:] = np.arange(-90, 90, outputRes)+(0.5*outputRes);
            var = ncout.createVariable("lon", float, ("lon",));
            var.units = "lon (degrees East)";
            var[:] = np.arange(-180, 180, outputRes)+(0.5*outputRes);
            
            #data variables
            var = ncout.createVariable("sss", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Mean monthly sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = sssRebinned;
            
            var = ncout.createVariable("sss_err", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Propagated uncertainty in the resampled monthly mean sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = sssErrs;
            
            var = ncout.createVariable("sss_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original RSS SMAP salinity data that were used to calculate resampled grid cell.";
            var[:] = sssCounts;
            
            ncout.close();
    
    print("Number RSS SMAP SSS data files skipped (not processed):", len(skipped));
    for val in skipped:
        print("\t", val);


def download_oisst_sst(downloadRoot, startYear, endYear):
    #OLD URL: https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/avhrr-only/201701/avhrr-only/201504/avhrr-only-v2.20150401.nc
    #OLD2 URL: https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/202005/oisst-avhrr-v02r01.20200501.nc
    #NEW URL: https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/${YYYY}${MM}/oisst-avhrr-v02r01.${YYYY}${MM}${DD}.nc
    #wget -r -A.nc https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/avhrr-only/201601/

    sourceTemplate = Template("https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/${YYYY}${MM}/oisst-avhrr-v02r01.${YYYY}${MM}${DD}.nc");
    #sourceTemplate = Template("https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/access/avhrr-only/${YYYY}${MM}/avhrr-only-v2.${YYYY}${MM}${DD}.nc");
    destinationTemplate = Template(path.join(downloadRoot, "OISST_SST/avhrr-only-v2.${YYYY}${MM}${DD}.nc"));

    notDownloadedError = [];
    notDownloadedSkipped = [];
    downloadedFiles = [];
    
    #sslContext = ssl._create_unverified_context(); #Ignore SSL verification
    for year in range(startYear, endYear+1):
        for month in range(1, 13):            
            monthStr = format(month, "02d");
            for day in range(1, calendar.monthrange(year, month)[1]+1):
                dayStr = format(day, "02d");
                url = sourceTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr);
                
                destination = destinationTemplate.safe_substitute(YYYY=year, MM=monthStr, DD=dayStr);
                
                #Create destination directory if it doesn't already exist
                if path.exists(path.dirname(destination)) == False:
                    makedirs(path.dirname(destination));
                    
                if path.exists(destination) == False:
                    print("Downloading:", path.basename(url));
                    try:
                        request = rrequests.get(url, verify=False, stream=True);
                        request.raw.decode_content = True;
                        with open(destination, 'wb') as f:
                            shutil.copyfileobj(request.raw, f);
                    except Exception as e:
                        print(type(e), e.args);
                        notDownloadedError.append((url, e));
                else: #Skip file, because it already seems to be there...
                    notDownloadedSkipped.append(url);
    
    print("Completed with:\n\t", len(downloadedFiles), "file(s) downloaded\n\t", len(notDownloadedError), "file(s) skipped due to error downloading.\n\t", len(notDownloadedSkipped), "file(s) skipped to prevent overwriting existing files.");


def process_oisst_sst(dataRoot, downloadedRoot,referenceFilename, startYear, endYear):
    def get_year_month_pairs(startYear, startMonth, endYear, endMonth):
        pairs = [];
        for year in range(startYear, endYear+1):
            for month in range(1, 13):
                if year == startYear and month < startMonth:
                    continue;
                if year == endYear and month > endMonth:
                    continue;
                
                pairs.append( (year, month) );
        return pairs;
    
    #Rebin 2D data to a specified shape.
    #Converts to a 4D array, then takes the mean in first and last dimension to give newly binned array.
    def rebin(m, shape):
        sh = (shape[0], m.shape[0]//shape[0], shape[1], m.shape[1]//shape[1]);
        return m.reshape(sh).mean(-1).mean(1);
    
    def rebin_sum(m, shape):
        sh = (shape[0], m.shape[0]//shape[0], shape[1], m.shape[1]//shape[1]);
        return m.reshape(sh).sum(-1).sum(1);
    
    def resample_reynolds(stopYear, stopMonth, startYear=1981, startMonth=9, lonResolution=1.0, latResolution=1.0, sourceTemplate=Template("downloaded_files/avhrr-only-v2.${YYYY}${MM}${DD}.nc"), destinationRootDirectory=""):
        if ((360.0 / lonResolution) % 1.0 != 0) or ((180.0 / latResolution) % 1.0 != 0):
            raise ValueError("Longitude or latitude resolution must be an exact multiple of longitude or latitude (respectively).");
        if lonResolution < 0.25 or latResolution < 0.25:
            raise ValueError("Longitude or latitude resolution cannot be smaller than 0.25 (downloaded data resolution).");
        
        #List to support multiple alternate file name templates (they changed the format half way through)
        if type(sourceTemplate) is not list:
            sourceTemplate = [sourceTemplate];
        
        #referenceFilename = path.join(path.dirname(reynoldspathe), "REFERENCE_FILE_FOR_METADATA-REYNOLDS.nc"); #Path to netCDF reference file.
        referenceNc = Dataset(referenceFilename, 'r');
        
        yearMonthsToRun = get_year_month_pairs(startYear, startMonth, stopYear, stopMonth); #inclusive
        
        filesSkipped = [];
        filesResampled = [];
        for (year, month) in yearMonthsToRun:
            currentOutputDir = path.join(downloadedRoot, "reynolds_avhrr_only_monthly_resampled_"+str(lonResolution)+"x"+str(latResolution), str(year));
            
            if path.exists(currentOutputDir) == False:
                makedirs(currentOutputDir);
            
            monthStr = format(month, "02d");
            
            curOutputpath = path.join(currentOutputDir, str(year)+monthStr+"01_OCF-SST-GLO-1M-100-REYNOLDS_"+str(lonResolution)+"x"+str(latResolution)+".nc");
            if path.exists(curOutputpath) == True: #Don't overwrite existing files.
                #print("Skipping year/month", year, monthStr)
                filesSkipped.append(curOutputpath);
                continue;
            
            print("Processing year/month:", year, month)
            sstVals = np.ma.empty((calendar.monthrange(year, month)[1], int(180/latResolution), int(360/lonResolution)));
            sstCountVals = np.ma.empty((calendar.monthrange(year, month)[1], int(180/latResolution), int(360/lonResolution)));
            sstErrVals = np.ma.empty((calendar.monthrange(year, month)[1], int(180/latResolution), int(360/lonResolution)));
            iceVals = np.ma.empty((calendar.monthrange(year, month)[1], int(180/latResolution), int(360/lonResolution)));
    
            for day in range(1, calendar.monthrange(year, month)[1]+1):
                dayStr = format(day, "02d");
                #try each supplied sourceTemplate until we find one that matches a file
                dailyInputNc = None;
                for filePathTemplate in sourceTemplate:
                    if path.exists(filePathTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr)):
                        dailyInputNc = Dataset(filePathTemplate.safe_substitute(YYYY=str(year), MM=monthStr, DD=dayStr), 'r');
                        break;
                if dailyInputNc is None: #No matching template was found
                    print(" *** WARNING: No sourceTemplate matched for ", year, monthStr, dayStr);
                
                sst = dailyInputNc.variables["sst"][0,0,:,:];
                sstRebinned = rebin(sst, (int(180/latResolution), int(360/lonResolution)));
                sstVals[day-1,:,:] = sstRebinned;
    
                sstCounts = sst.mask==False;
                sstCounts = rebin_sum(sstCounts, (int(180/latResolution), int(360/lonResolution)));
                hasData = np.where(sstCounts > 0);
                sstCountVals[day-1,:,:] = sstCounts;
    
                sstErrSquared = dailyInputNc.variables["err"][0,0,:,:]**2;
                sstErrRebinned = np.sqrt(rebin_sum(sstErrSquared, (int(180/latResolution), int(360/lonResolution))));
                sstErrRebinned[hasData] /= sstCounts[hasData];
                sstErrVals[day-1,:,:] = sstErrRebinned;
                #allErrs[hasData] += errs[hasData]**2; #sum errors as variance
                
                ice = dailyInputNc.variables["ice"][0,0,:,:];
                iceRebinned = rebin(ice, (int(180/latResolution), int(360/lonResolution)));
                iceVals[day-1,:,:] = iceRebinned;
        
        
            #Calculate data for monthly output file
            #SST
            means = sstVals.mean(0);
            counts = sstCountVals.sum(0);
            errs = np.sqrt(np.sum(sstErrVals**2, axis=0) / sstVals.count(0));
            
            means = np.roll(means, -int(180/latResolution), axis=1);
            counts = np.roll(counts, -int(180/latResolution), axis=1);
            errs = np.roll(errs, -int(180/latResolution), axis=1);
            
            #ICE
            meansIce = iceVals.mean(0);
            meansIce = np.roll(meansIce, -int(180/latResolution), axis=1);
        
            
            #Write to netCDF
            newnc = Dataset(curOutputpath, 'w');
            newnc.createDimension(u'time', 1); #create dimensions
            newnc.createDimension(u'lat', int(180/latResolution));
            newnc.createDimension(u'lon', int(360/lonResolution));
        
            refTime = datetime(1981, 1, 1, 0, 0, 0);
            curTime = datetime(year, month, 1, 0, 0, 0);
            td = curTime-refTime;
            timeVar = newnc.createVariable("time", "float64", (u"time",));
            timeVar.setncatts({k: referenceNc.variables["time"].getncattr(k) for k in referenceNc.variables["time"].ncattrs()}); # Copy variable attributes
            timeVar[:] = [td.total_seconds()];
        
            lon = np.arange(-180.0+(0.5*lonResolution), 180, lonResolution); #-179.5 to 179.5 in steps of 1
            lonVar = newnc.createVariable("lon", "float64", (u"lon",));
            lonVar.setncatts({k: referenceNc.variables["lon"].getncattr(k) for k in referenceNc.variables["lon"].ncattrs()}); # Copy variable attributes
            lonVar[:] = lon;
            
            lat = np.arange(-90.0+(0.5*latResolution), 90, latResolution); #-89.5 to 89.5 in steps of 1
            latVar = newnc.createVariable("lat", "float64", (u"lat",));
            latVar.setncatts({k: referenceNc.variables["lat"].getncattr(k) for k in referenceNc.variables["lat"].ncattrs()}); # Copy variable attributes
            latVar[:] = lat;
            
            #SST
            meanVar = newnc.createVariable("sst_mean", "float32", (u"time", u"lat", u"lon"));
            meanVar.setncatts({k: referenceNc.variables["sst_mean"].getncattr(k) for k in referenceNc.variables["sst_mean"].ncattrs()}); # Copy variable attributes
            meanVar[:] = np.ma.expand_dims(means, axis=0);
            
            countVar = newnc.createVariable("sst_count", "uint32", (u"time", u"lat", u"lon"));
            countVar.setncatts({k: referenceNc.variables["sst_count"].getncattr(k) for k in referenceNc.variables["sst_count"].ncattrs()}); # Copy variable attributes
            countVar[:] = np.ma.expand_dims(counts, axis=0);
            
            stddevVar = newnc.createVariable("sst_stddev", "float32", (u"time", u"lat", u"lon"));
            stddevVar.setncatts({k: referenceNc.variables["sst_stddev"].getncattr(k) for k in referenceNc.variables["sst_stddev"].ncattrs()}); # Copy variable attributes
            stddevVar[:] = np.ma.expand_dims(errs, axis=0);
            
            #ICE
            meanIce = newnc.createVariable("ice", "float32", (u"time", u"lat", u"lon"));
            meanIce.valid_min = 0.0;
            meanIce.valid_max = 100.0;
            meanIce.units = "percentage";
            meanIce.long_name = "Sea ice concentration";
            meanIce.fill_value = -999.0;
            meanIce[:] = np.ma.expand_dims(meansIce, axis=0);
            
            newnc.close();
            filesResampled.append(curOutputpath);
        
        print("Completed with:\n\t", len(filesResampled), "file(s) resampled\n\t", len(filesSkipped), "file(s) skipped to prevent overwriting existing files.");
    
    
    #Set start year and month to prevent starting before data exists
    if startYear < 1981: #No data before 1981
        startYear = 1981;
    if startYear == 1981:
        startMonth = 9; #No data before September 1981
    else:
        startMonth = 1;
    #Run the wrapped function
    sourceTemplate=Template(path.join(downloadedRoot, "OISST_SST", "avhrr-only-v2.${YYYY}${MM}${DD}.nc"));
    destinationRootDirectory=path.join(dataRoot, "OISST_SST");
    resample_reynolds(stopYear=endYear, stopMonth=12, startYear=startYear, startMonth=startMonth, sourceTemplate=sourceTemplate, destinationRootDirectory=destinationRootDirectory);


def download_isas_sss_sst(dataRoot):
    url = "https://www.seanoe.org/data/00412/52367/data/53158.tar.gz"; #One single 25GB file
    destination = path.join(dataRoot, "downloaded", "ISAS_SSS_SST", "53158.tar.gz");
    
    if path.exists(destination):
        print("Skipping ISAS SSS SST download to avoid overwriting existing file.");
        return;
        
    if path.exists(path.dirname(destination)) == False:
        makedirs(path.dirname(destination));
    
    print("Downloading ISAS SSS SST (large file):", path.basename(url));
    try:
        request = rrequests.get(url, verify=False, stream=True);
        request.raw.decode_content = True;
        with open(destination, 'wb') as f:
            shutil.copyfileobj(request.raw, f);
    except Exception as e:
        print("Error downloading ISAS SSS SST");
        print(type(e), e.args);

def uncompress_isas_sss_sst(downloadedRoot):
    tarFilePath = path.join(downloadedRoot, "ISAS_SSS_SST", "53158.tar.gz");
    destinationPath = path.join(downloadedRoot, "ISAS_SSS_SST", "uncompressed");
    
    #Check if an uncompressed file already exists, if it does assume it's already been uncompressed
    if path.exists(path.join(destinationPath, "2002", "ISAS15_DM_20020115_fld_PSAL.nc")) == True:
        print("Skipping uncompression of ISAS SSS SST to avoid overwriting existing files...");
        return;
    
    print("Uncompressing ISAS SSS SST (large file)...");
    try:
        tar = tarfile.open(tarFilePath);
        tar.extractall(path=destinationPath);
        tar.close();
    except FileNotFoundError as e:
        print("ISAS SSS SST: Could not untar downloaded file:\n", e);
        raise e;


def process_isas_sss_sst(dataRoot, downloadedRoot, startYear, endYear):
    def process_slice(valData, errData, outputRes=1.0):
        newGrid = np.full((180, 360), np.nan, dtype=float);
        newGridCount = np.zeros((180, 360), dtype=float);
        newGridErr = np.full((180, 360), np.nan, dtype=float);
        for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
            for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                if iCoordMeshes[ilat, ilon] is not None:
                    newGrid[ilat, ilon] = np.mean(valData[iCoordMeshes[ilat, ilon]]);
                    newGridCount[ilat, ilon] = len(iCoordMeshes[ilat, ilon][0]) * len(iCoordMeshes[ilat, ilon][0][0]);
                    newGridErr[ilat, ilon] = np.sum(errData[iCoordMeshes[ilat, ilon]]**2);
        newGridErr = newGridErr / newGridCount;
        
        return newGrid, newGridCount, newGridErr;
    
    downloadedFileTemplate = Template(path.join(downloadedRoot, "ISAS_SSS_SST", "uncompressed", "${YYYY}/ISAS15_DM_${YYYY}${MM}15_fld_${VARNAME}.nc"));
    varNameTemperature = "TEMP";
    varNameSalinity = "PSAL";
    processedFileTemplate = Template(path.join(dataRoot, "ISAS_SSS_SST", "${YYYY}/ISAS15_DM_${YYYY}_${MM}_processed.nc"));
    outputRes = 1.0;
    
    #Constrain years to temporal range of the dataset
    if startYear < 2002:
        startYear = 2002;
    if endYear > 2015:
        endYear = 2015;
    years = range(startYear, endYear+1);
    
    iCoordMeshes = None; #Initially this is None but will be calculated exactly once.
    
    for year in years:
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            
            #Create directory if it doens't already exist
            processedFilePath = processedFileTemplate.safe_substitute(YYYY=year, MM=monthStr);
            if path.exists(path.dirname(processedFilePath)) == False:
                makedirs(path.dirname(processedFilePath));
            
            #Check to see if file processed output file already exists, don't overwrite
            if path.exists(processedFilePath) == True:
                print("Skipping ISAS SSS SST at", year, monthStr, "to avoid overwriting existing files.");
                continue;
            
            print("Processing", year, monthStr+"...");
    
            temperaturePath = downloadedFileTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME=varNameTemperature);
            tempNC = Dataset(temperaturePath, 'r');
            tempData = tempNC.variables["TEMP"][0,0,:,:];
            tempErrData = tempNC.variables["TEMP_ERR"][0,0,:,:];
            ISASLats = tempNC.variables["latitude"][:];
            ISASLons = tempNC.variables["longitude"][:];
            
            salinityPath = downloadedFileTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME=varNameSalinity);
            salNC = Dataset(salinityPath, 'r');
            salData = salNC.variables["PSAL"][0,0,:,:];
            salErrData = salNC.variables["PSAL_ERR"][0,0,:,:];
            
            #Calculate binning information
            if iCoordMeshes is None: #Only do this once because it's computationally expensive but the same for all months
                print("Calculating grid cell mapping...");
                iCoordMeshes = np.full((180, 360), None, dtype=object);
                for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
                    print("Grid cell mapping for latitude", lat);
                    for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                        wlat = np.where((ISASLats >= lat) & (ISASLats < (lat+outputRes)));
                        wlon = np.where((ISASLons >= lon) & (ISASLons< (lon+outputRes)));
                        
                        if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                            iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
            
            #Process a slice
            newTemp, newTempCount, newTempErr = process_slice(tempData, tempErrData);
            newSal, newSalCount, newSalErr = process_slice(salData, salErrData);
            
            #Output to new netCDF file
            ncout = Dataset(processedFilePath, 'w');
            
            ncout.createDimension("lat", int(180/outputRes));
            ncout.createDimension("lon", int(360/outputRes));
            
            #dimension variables
            var = ncout.createVariable("lat", float, ("lat",));
            var.units = "lat (degrees North)";
            var[:] = np.arange(-90, 90, outputRes)+(0.5*outputRes);
            var = ncout.createVariable("lon", float, ("lon",));
            var.units = "lon (degrees East)";
            var[:] = np.arange(-180, 180, outputRes)+(0.5*outputRes);
            
            #data variables
            var = ncout.createVariable("TEMP", float, ("lat", "lon"));
            var.units = "Degrees Celcius";
            var.long_name = "Mean monthly sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newTemp;
            
            var = ncout.createVariable("TEMP_err", float, ("lat", "lon"));
            var.units = "Degrees Celcius";
            var.long_name = "Uncertainty in the monthly mean sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newTempErr;
            
            var = ncout.createVariable("TEMP_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original ISAS temperature data that were used to calculate resampled grid cell.";
            var[:] = newTempCount;
            
            var = ncout.createVariable("PSAL", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Mean monthly sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = newSal;
            
            var = ncout.createVariable("PSAL_err", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Uncertainty in the monthly mean sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = newSalErr;
            
            var = ncout.createVariable("PSAL_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original ISAS salinity data that were used to calculate resampled grid cell.";
            var[:] = newSalCount;
            
            ncout.close();
            

def download_cci_ostia_sst(downloadedRoot, startYear, endYear):
    import cdsapi;
    c = cdsapi.Client();
    
    outputTemplate = Template(path.join(downloadedRoot, "ESACCI_SST_OSTIA", "CCI_OSTIA_SST_${YYYY}_${MM}_${DD}.tar.gz"));
    
    if startYear < 1981:
        startYear = 1981;
    if startYear == 1981:
        startiMonth = 8;
    else:
        startiMonth = 0;
    if endYear > 2016:
        endYear = 2016;
        
    endiMonth = 11;
    
    
    for year in range(startYear, endYear+1)[::-1]:
        for imonth in range(0, 12):
            #Ignore missing months from start and end of range
            if (year == startYear) and (imonth < startiMonth):
                continue;
            if (year == endYear) and (imonth > endiMonth):
                continue;
            
            monthStr = format(imonth+1, "02d");
            
            for day in range(1, calendar.monthrange(year, imonth+1)[1]+1):
                dayStr = format(day, "02d");
                
                outputPath = outputTemplate.safe_substitute(YYYY=year, MM=monthStr, DD=dayStr);
                if path.exists(path.dirname(outputPath)) == False:
                    makedirs(path.dirname(outputPath));
                
                if path.exists(outputPath):
                    print("Skipping ESACCI SST OSTIA for:", year, monthStr, dayStr, "to avoid overwritting existing files");
                    continue;
                else:
                    print("Downloading RSACCI SST OSTIA for:", year, monthStr, dayStr, "...");
                
                c.retrieve(
                    "satellite-sea-surface-temperature-ensemble-product",
                    {
                        "variable": "all",
                        "format": "tgz",
                        "day": dayStr,
                        "month": monthStr,
                        "year": str(year),
                    },
                    outputPath);

def uncompress_cci_ostia_sst(downloadedRoot, startYear, endYear):
    tarPathTemplate = Template(path.join(downloadedRoot, "ESACCI_SST_OSTIA", "CCI_OSTIA_SST_${YYYY}_${MM}_${DD}.tar.gz"));
    untarPathTemplate = Template(path.join(downloadedRoot, "CCI_OSTIA_SST", "netCDF", "${YYYY}"));
    
    if startYear < 1981:
        startYear = 1981;
    if startYear == 1981:
        startiMonth = 8;
    else:
        startiMonth = 0;
    if endYear > 2016:
        endYear = 2016;
        
    endiMonth = 11;
    
    missingFiles = [];
    for year in range(startYear, endYear+1):
        for imonth in range(0, 12):
            #Ignore missing months from start and end of range
            if (year == startYear) and (imonth < startiMonth):
                continue;
            if (year == endYear) and (imonth > endiMonth):
                continue;
            
            monthStr = format(imonth+1, "02d");
            
            for day in range(1, calendar.monthrange(year, imonth+1)[1]+1):
                dayStr = format(day, "02d");
                tarFilePath = tarPathTemplate.safe_substitute(YYYY=year, MM=monthStr, DD=dayStr);
                untarFilePath = untarPathTemplate.safe_substitute(YYYY=year, MM=monthStr, DD=dayStr);
                
                if path.exists(path.dirname(untarFilePath)) == False:
                    makedirs(path.dirname(untarFilePath));
                if path.exists(untarFilePath) == False:
                    print("Skipping CCI OSTIA SST uncompression to prevent overwriting existing files for:", year, monthStr, dayStr);
                    continue;
                
                print("Uncompressing CCI OSTIA SST", year, monthStr, dayStr+"...");
                #Try the uncompression
                try:
                    tar = tarfile.open(tarFilePath);
                    tar.extractall(path=untarFilePath);
                    tar.close();
                except FileNotFoundError:
                    missingFiles.append((year, monthStr, dayStr));
                    continue;
    
    print("Finished with", len(missingFiles), "missing files.");
    for v in missingFiles:
        print("\t", v);


def process_cci_ostia_sst(dataRoot, downloadedRoot, startYear, endYear):
    def process_slice(valData, errData, countData, outputRes=1.0):
        newGrid = np.full((180, 360), np.nan, dtype=float);
        newGridCount = np.zeros((180, 360), dtype=float);
        newGridErr = np.full((180, 360), np.nan, dtype=float);
        for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
            for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                if iCoordMeshes[ilat, ilon] is not None:
                    newGrid[ilat, ilon] = np.mean(valData[iCoordMeshes[ilat, ilon]]);
                    newGridCount[ilat, ilon] = np.sum(countData[iCoordMeshes[ilat, ilon]]);
                    #newGridCount[ilat, ilon] = len(iCoordMeshes[ilat, ilon][0]) * len(iCoordMeshes[ilat, ilon][0][0]);
                    newGridErr[ilat, ilon] = np.sqrt(np.sum(errData[iCoordMeshes[ilat, ilon]]**2));
        
        newGridErr[newGridCount!=0] = newGridErr[newGridCount!=0] / newGridCount[newGridCount!=0];
        
        return newGrid, newGridCount, newGridErr;
    
    
    downloadedFileTemplate = Template(path.join(downloadedRoot, "CCI_OSTIA_SST", "netCDF", "${YYYY}/${YYYY}${MM}${DD}120000-ESACCI-L4_GHRSST-SST-GMPE-GLOB_CDR2.0-v02.0-fv01.0.nc"));
    processedFileTemplate = Template(path.join(dataRoot, "ESACCI_SST_OSTIA", "processed", "${YYYY}/ESACCI-SST-OSTIA-LT-v02.0-fv01.1_${YYYY}_${MM}_processed.nc"));
    outputRes = 1.0;
    
    #Bound temporal range to data set range
    if startYear < 1981:
        startYear = 1981;
    if endYear > 2016:
        endYear = 2016;
    overwrite = False;
    
    iCoordMeshes = None; #Initially this is None but will be calculated exactly once. It's a long calculation so doesn't want to be repeated.
    skipped = [];
    for year in range(startYear, endYear+1):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            print("Processing", year, monthStr+"...");
            daysInMonth = calendar.monthrange(year, imonth+1)[1];
            
            processedFilePath = processedFileTemplate.safe_substitute(YYYY=year, MM=monthStr);
            if (path.exists(processedFilePath) == False) and (overwrite == True):
                print("Skipping", year, monthStr, "to prevent overwriting existing file.");
                continue;
            
            #Store data for each day of the month
            monthVals = np.empty((daysInMonth, 720, 1440), dtype=float);
            monthValsErr = np.empty((daysInMonth, 720, 1440), dtype=float);
    
            #Read input data
            for iday in range(0, daysInMonth):
                dayStr = format(iday+1, "02d");
                cciPath = downloadedFileTemplate.safe_substitute(YYYY=year, MM=monthStr, DD=dayStr);
                if path.exists(cciPath) == False:
                    skipped.append((year, monthStr));
                    break;
                cciNC = Dataset(cciPath, 'r');
                monthVals[iday,:,:] = cciNC.variables["analysed_sst"][0,:,:];
                monthValsErr[iday,:,:] = cciNC.variables["standard_deviation"][0,:,:];
            
            #If any of the days were skipped, then move to the next month
            if (year, monthStr) in skipped:
                print("Skipping CCI OSTIA SST file for:", year, monthStr);
                continue;
            
            #Calculate means for the month from daily values
            monthVals[monthVals == -32768.0] = np.nan;
            monthVals = np.nanmean(monthVals, axis=0);
            monthValsErr[monthValsErr == -32768.0] = np.nan;
            monthValsCounts = np.sum(np.isfinite(monthValsErr), axis=0);
            monthValsErr = np.sqrt(np.nansum(monthValsErr**2, axis=0));
            monthValsErr[monthValsCounts!=0] = monthValsErr[monthValsCounts!=0] / monthValsCounts[monthValsCounts!=0];
            
            
            #Calculate binning information
            if iCoordMeshes is None: #Only do this once because it's computationally expensive but the same for all months
                CCILats = cciNC.variables["lat"][:];
                CCILons = cciNC.variables["lon"][:];
                print("Calculating grid cell mapping...");
                iCoordMeshes = np.full((180, 360), None, dtype=object);
                for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
                    print("Grid cell mapping for latitude", lat);
                    for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                        wlat = np.where((CCILats >= lat) & (CCILats < (lat+outputRes)));
                        wlon = np.where((CCILons >= lon) & (CCILons< (lon+outputRes)));
                        
                        if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                            iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
            
            
            #Process a slice
            newVals, newCountCount, newValsErr = process_slice(monthVals, monthValsErr, monthValsCounts);
            
            #Output to new netCDF file
            if path.exists(path.dirname(processedFilePath)) == False:
                makedirs(path.dirname(processedFilePath));
            ncout = Dataset(processedFilePath, 'w');
            
            ncout.createDimension("lat", int(180/outputRes));
            ncout.createDimension("lon", int(360/outputRes));
            
            #dimension variables
            var = ncout.createVariable("lat", float, ("lat",));
            var.units = "lat (degrees North)";
            var[:] = np.arange(-90, 90, outputRes)+(0.5*outputRes);
            var = ncout.createVariable("lon", float, ("lon",));
            var.units = "lon (degrees East)";
            var[:] = np.arange(-180, 180, outputRes)+(0.5*outputRes);
            
            #data variables
            var = ncout.createVariable("sst", float, ("lat", "lon"));
            var.units = "Degrees Kelvin";
            var.long_name = "Mean monthly sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newVals;
            
            var = ncout.createVariable("sst_err", float, ("lat", "lon"));
            var.units = "Degrees Kelvin";
            var.long_name = "Uncertainty in the monthly mean sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newValsErr;
            
            var = ncout.createVariable("sst_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original CCI SST-OSTIA data that were used to calculate resampled grid cell.";
            var[:] = newCountCount;
            
            ncout.close();


def download_esacci_sss_smos(downloadRoot, startYear, endYear):
    urlTemplate = Template("ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sea_surface_salinity/data/v01.8/30days/${YYYY}/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_Monthly_CENTRED_15Day_25km-${YYYY}${MM}15-fv1.8.nc");
    destinationTemplate = Template(path.join(downloadRoot, "ESACCI_SSS_SMOS", "${YYYY}/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_Monthly_CENTRED_15Day_25km-${YYYY}${MM}15-fv1.8.nc"));

    if startYear < 2010:
        startYear = 2010;
    
    for year in range(startYear, endYear+1):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            print("Downloading ESACCI SSS SMOS for:", year, monthStr+"...");
        
            #generate source url and destination path
            url = urlTemplate.safe_substitute(YYYY=year, MM=monthStr);
            destination = destinationTemplate.safe_substitute(YYYY=year, MM=monthStr);
            
            #Make directory if it doesn't exist
            if path.exists(path.dirname(destination)) == False:
                makedirs(path.dirname(destination));
            
            #Do not overwrite any data
            if path.exists(destination) == True:
                print("Skipping ESACCI SSS SMOS for:", year, monthStr, "to prevent overwriting existing file.");
                continue;
    
            #Try the download
            try:
                with closing(request.urlopen(url)) as req:
                    with open(destination, 'wb') as dest:
                        shutil.copyfileobj(req, dest);
    
            except Exception as e:
                print("Error downloading ESACCI SSS SMOS for:", year, monthStr);
                print(e);


def process_esacci_sss_smos(dataRoot, downloadedRoot, startYear, endYear):
    def process_slice(valData, errData, outputRes=1.0):
        newGrid = np.full((180, 360), np.nan, dtype=float);
        newGridCount = np.zeros((180, 360), dtype=float);
        newGridErr = np.full((180, 360), np.nan, dtype=float);
        for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
            for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                if iCoordMeshes[ilat, ilon] is not None:
                    newGrid[ilat, ilon] = np.mean(valData[iCoordMeshes[ilat, ilon]]);
                    newGridCount[ilat, ilon] = len(iCoordMeshes[ilat, ilon][0]) * len(iCoordMeshes[ilat, ilon][0][0]);
                    newGridErr[ilat, ilon] = np.sqrt(np.sum(errData[iCoordMeshes[ilat, ilon]]**2));
        
        newGridErr[newGridCount!=0] = newGridErr[newGridCount!=0] / newGridCount[newGridCount!=0];
        
        return newGrid, newGridCount, newGridErr;
    
    downloadedFileTemplate = Template(path.join(downloadedRoot, "ESACCI_SSS_SMOS", "${YYYY}/ESACCI-SEASURFACESALINITY-L4-SSS-MERGED_OI_Monthly_CENTRED_15Day_25km-${YYYY}${MM}15-fv1.8.nc"));
    processedFileTemplate = Template(path.join(dataRoot, "ESACCI_SSS_SMOS", "processed", "${YYYY}/ESACCI-SEASURFACESALINITY-L4-CENTRED15Day_${YYYY}_${MM}_processed.nc"));
    outputRes = 1.0;
    
    if startYear < 2010:
        startYear = 2010;
    
    iCoordMeshes = None; #Initially this is None but will be calculated exactly once.
    for year in range(startYear, endYear+1):
        for imonth in range(0, 12):
            monthStr = format(imonth+1, "02d");
            
            processedFilePath = processedFileTemplate.safe_substitute(YYYY=year, MM=monthStr);
            if path.exists(processedFilePath) == True:
                print("Skipping ESACCI SSS SMOS for:", year, monthStr, "(to avoid overwriting existing files).");
                continue;
            if path.exists(path.dirname(processedFilePath)) == False:
                makedirs(path.dirname(processedFilePath));
            
            
            print("Processing ESACCI SSS SMOS for:", year, monthStr+"...");
    
            #Read input data
            cciPath = downloadedFileTemplate.safe_substitute(YYYY=year, MM=monthStr);
            if path.exists(cciPath) == False:
                print("Skipping ESACCI SSS SMOS for:", year, monthStr, "(no raw downloaded file found).");
                continue;
            cciNC = Dataset(cciPath, 'r');
            salData = cciNC.variables["sss"][0,:,:];
            salErrData = cciNC.variables["sss_random_error"][0,:,:];
            CCILats = cciNC.variables["lat"][:];
            CCILons = cciNC.variables["lon"][:];
            
            #Calculate binning information
            if iCoordMeshes is None: #Only do this once because it's computationally expensive but the same for all months
                print("Calculating grid cell mapping...");
                iCoordMeshes = np.full((180, 360), None, dtype=object);
                for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
                    #print("Grid cell mapping for latitude", lat);
                    for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                        wlat = np.where((CCILats >= lat) & (CCILats < (lat+outputRes)));
                        wlon = np.where((CCILons >= lon) & (CCILons< (lon+outputRes)));
                        
                        if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                            iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
            
            #Process a slice
            newSal, newSalCount, newSalErr = process_slice(salData, salErrData);
            
            #Output to new netCDF file
            ncout = Dataset(processedFilePath, 'w');
            
            ncout.createDimension("lat", int(180/outputRes));
            ncout.createDimension("lon", int(360/outputRes));
            
            #dimension variables
            var = ncout.createVariable("lat", float, ("lat",));
            var.units = "lat (degrees North)";
            var[:] = np.arange(-90, 90, outputRes)+(0.5*outputRes);
            var = ncout.createVariable("lon", float, ("lon",));
            var.units = "lon (degrees East)";
            var[:] = np.arange(-180, 180, outputRes)+(0.5*outputRes);
            
            #data variables
            var = ncout.createVariable("sss", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Mean monthly sea surface salinity resamples to a 1x1 degree spatial resolution";
            var[:] = newSal;
            
            var = ncout.createVariable("sss_err", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Uncertainty in the monthly mean sea surface salinity resamples to a 1x1 degree spatial resolution";
            var[:] = newSalErr;
            
            var = ncout.createVariable("sss_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original CCI salinity data that were used to calculate resampled grid cell.";
            var[:] = newSalCount;
            
            ncout.close();


def download_cora_sss_sst(downloadRoot, startYear, endYear, ftpUser, ftpPass):
    urlTemplate = Template("Core/INSITU_GLO_TS_OA_REP_OBSERVATIONS_013_002_b/CORIOLIS-GLOBAL-CORA-OBS_FULL_TIME_SERIE/field/${YYYY}/OA_CORA5.2_${YYYY}${MM}15_fld_${VARNAME}.nc");
    destinationTemplate = Template(path.join(downloadRoot, "CORA_SSS_SST", "${YYYY}/OA_CORA5.2_${YYYY}${MM}15_fld_${VARNAME}.nc"))
    #ftp://tholding@my.cmems-du.eu/Core/INSITU_GLO_TS_OA_REP_OBSERVATIONS_013_002_b/CORIOLIS-GLOBAL-CORA-OBS_FULL_TIME_SERIE/field/1990/OA_CORA5.2_19900115_fld_PSAL.nc
    varNames = ["PSAL", "TEMP"];
    
    if startYear < 1990:
        startYear = 1990;
    
    #Open ftp connection and log in
    ftp = FTP("my.cmems-du.eu");
    ftp.login(ftpUser, ftpPass);
    
    for year in range(startYear, endYear+1):
        for month in range(1, 13):
            monthStr = format(month, "02d");
            for varName in varNames:
                url = urlTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME=varName, UNAME=ftpUser, PASS=ftpPass);
                destination = destinationTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME=varName);
                
                if path.exists(destination) == True:
                    print("Skipping download of CORA SSS SST for: {0} {1} {2} to avoid overwriting existing file.".format(varName, year, monthStr));
                    continue;
                if path.exists(path.dirname(destination)) == False:
                    makedirs(path.dirname(destination));
                
                #Try the download
                try:
                    print("Downloading CORA SSS SST {0} {1} {2} (large file)".format(varName, year, monthStr));
                    ftp.retrbinary("RETR "+url, open(destination, 'wb').write);
                except Exception as e:
                    print("Error downloading CORA SSS SST for: {0} {1} {2}".format(varName, year, monthStr));
                    print(e);
                    
    ftp.close();
    
    
def process_cora_sss_sst(dataRoot, downloadRoot, startYear, endYear):
    def convert_date(inputTime, baseDate):
        currentDate = baseDate + timedelta(days=int(inputTime));
        return currentDate.year, currentDate.month-1, format(currentDate.month, "02d");
    
    def process_slice(valData, errData, iCoordMeshes, outputRes=1.0):
        newGrid = np.full((180, 360), np.nan, dtype=float);
        newGridCount = np.zeros((180, 360), dtype=float);
        newGridErr = np.full((180, 360), np.nan, dtype=float);
        for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
            for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                if iCoordMeshes[ilat, ilon] is not None:
                    newGrid[ilat, ilon] = np.mean(valData[iCoordMeshes[ilat, ilon]]);
                    newGridCount[ilat, ilon] = len(iCoordMeshes[ilat, ilon][0]) * len(iCoordMeshes[ilat, ilon][0][0]);
                    newGridErr[ilat, ilon] = np.sum(errData[iCoordMeshes[ilat, ilon]]**2);
        newGridErr = newGridErr / newGridCount;
        
        return newGrid, newGridCount, newGridErr;
    
    downloadedTemplate = Template(path.join(downloadRoot, "CORA_SSS_SST", "${YYYY}/OA_CORA5.2_${YYYY}${MM}15_fld_${VARNAME}.nc"));
    processedTemplate = Template(path.join(dataRoot, "CORA_SSS_SST", "${YYYY}/OA_CORA5.2_${YYYY}_${MM}_processed.nc"));
    outputRes = 1.0;
    
    if startYear < 1990:
        startYear = 1990;
    
    iCoordMeshes = None; #Initially this is None but will be calculated exactly once.
    
    for year in range(startYear, endYear+1):
        for month in range(1, 13):
            monthStr = format(month, "02d");
            processedPath = processedTemplate.safe_substitute(YYYY=year, MM=monthStr);
            
            #check for existing files / directory structure
            if path.exists(processedPath) == True:
                print("Skipping CORA SSS SST to avoid overwriting existing files for:", year, monthStr);
                continue;
            
            #open both PSAL and TEMP files#
            try:
                psalPath = downloadedTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME="PSAL");
                psalNC = Dataset(psalPath, 'r');
                tempPath = downloadedTemplate.safe_substitute(YYYY=year, MM=monthStr, VARNAME="TEMP");
                tempNC = Dataset(tempPath, 'r');
            except FileNotFoundError as e:
                print("Error with CORA SSS SST for", year, month, "(couldn't find downloaded file:");
                print(e);
                continue;

            print("Processing CORA SSS SST for:", year, monthStr);            
            #Extract and process data
            #inputTime = psalNC.variables["time"][:];
            inputLats = psalNC.variables["latitude"][:];
            inputLons = psalNC.variables["longitude"][:];
            salData = psalNC.variables["PSAL"][0, 0, :, :];
            salErrData = psalNC.variables["PSAL_ERR"][0, 0, :, :];
            tempData = tempNC.variables["TEMP"][0, 0, :, :];
            tempErrData = tempNC.variables["TEMP_ERR"][0, 0, :, :];
            
            #Calculate binning information
            if iCoordMeshes is None: #Only do this once because it's computationally expensive but the same for all months
                print("Calculating grid cell mapping...");
                iCoordMeshes = np.full((180, 360), None, dtype=object);
                for ilat, lat in enumerate(np.arange(-90, 90, outputRes)):
                    for ilon, lon in enumerate(np.arange(-180, 180, outputRes)):
                        wlat = np.where((inputLats >= lat) & (inputLats < (lat+outputRes)));
                        wlon = np.where((inputLons >= lon) & (inputLons< (lon+outputRes)));
                        
                        if (len(wlat[0]) > 0) & (len(wlon[0]) > 0):
                            iCoordMeshes[ilat, ilon] = np.meshgrid(wlat[0], wlon[0]);
            
            #bin mean and error data for salinity and temperature
            newTemp, newTempCount, newTempErr = process_slice(tempData, tempErrData, iCoordMeshes);
            newSal, newSalCount, newSalErr = process_slice(salData, salErrData, iCoordMeshes);
            
            #Create output directory structure, if needed
            if path.exists(path.dirname(processedPath)) == False:
                makedirs(path.dirname(processedPath));
                
            #Create and write to output netCDF file
            ncout = Dataset(processedPath, 'w');
            
            ncout.createDimension("lat", int(180/outputRes));
            ncout.createDimension("lon", int(360/outputRes));
            
            #dimension variables
            var = ncout.createVariable("lat", float, ("lat",));
            var.units = "lat (degrees North)";
            var[:] = np.arange(-90, 90, outputRes)+(0.5*outputRes);
            var = ncout.createVariable("lon", float, ("lon",));
            var.units = "lon (degrees East)";
            var[:] = np.arange(-180, 180, outputRes)+(0.5*outputRes);
            
            #data variables
            var = ncout.createVariable("TEMP", float, ("lat", "lon"));
            var.units = "Degrees Celcius";
            var.long_name = "Mean monthly sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newTemp;
            
            var = ncout.createVariable("TEMP_err", float, ("lat", "lon"));
            var.units = "Degrees Celcius";
            var.long_name = "Uncertainty in the monthly mean sea surface temperature resampled to a 1x1 degree spatial resolution";
            var[:] = newTempErr;
            
            var = ncout.createVariable("TEMP_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original OA_CORA5.2 temperature data that were used to calculate resampled grid cell.";
            var[:] = newTempCount;
            
            var = ncout.createVariable("PSAL", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Mean monthly sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = newSal;
            
            var = ncout.createVariable("PSAL_err", float, ("lat", "lon"));
            var.units = "PSU";
            var.long_name = "Uncertainty in the monthly mean sea surface salinity resampled to a 1x1 degree spatial resolution";
            var[:] = newSalErr;
            
            var = ncout.createVariable("PSAL_count", int, ("lat", "lon"));
            var.units = "count (integer)";
            var.long_name = "Number of samples from the original OA_CORA5.2 salinity data that were used to calculate resampled grid cell.";
            var[:] = newSalCount;
            
            ncout.close();
            psalNC.close();
            tempNC.close();
            

if __name__ == "__main__":

    startYear = 2018; #clArgs.start_year; #year to try to start downloading data for
    endYear = 2020; #clArgs.end_year; #last year to try to download data for

    dataRoot='K:\downloaded\data\predictiondatasets'
    downloadedRoot='K:\downloaded\data\predictiondatasets'
    oceanMaskPath= "C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data\\World_Seas-IHO-mask.nc"
    destinationRoot='K:\downloaded\data\predictiondatasets'
    downloadRoot='K:\downloaded\data\predictiondatasets'
    referenceFilename='C:\\Users\\rps207\\Documents\\Python\\2021-Oceansoda_pre_github\\OceanSODA-master\\aux_data\\reference_data\\REFERENCE_FILE_FOR_METADATA-REYNOLDS.nc'
    destinationRootDirectory='K:\downloaded\data\predictiondatasets'
    
    # download_all_woa_nutrients(destinationRoot)
    # process_all_woa_nutrients(dataRoot, downloadedRoot, oceanMaskPath)
    
    # download_rss_smap_sss(dataRoot, 2015, 2021)
    # process_rss_smap_sss(dataRoot, downloadedRoot, 2015, 2021)
    
    download_oisst_sst(downloadRoot, 1981, 2020)
    process_oisst_sst(dataRoot, downloadedRoot,referenceFilename,1981, 2020)
    
    # #note that the files are no longer hosted there for the 2015 dataset
    # #so i was forced to get them directly from Tom. 
    # # download_isas_sss_sst(dataRoot)
    # # uncompress_isas_sss_sst(downloadedRoot)
    # # process_isas_sss_sst(dataRoot, downloadedRoot, startYear=2015, endYear=2021);
    
    # download_cci_ostia_sst(downloadedRoot, 1981, 2016)
    # uncompress_cci_ostia_sst(downloadedRoot, 1981, 2016)
    # process_cci_ostia_sst(dataRoot, downloadedRoot, 1981, 2016)
    
    # download_esacci_sss_smos(downloadedRoot, startYear=2010, endYear=2021);
    # process_esacci_sss_smos(dataRoot, downloadedRoot, startYear=2010, endYear=2021);
    
    # ftpUser='rsims';
    # ftpPass='CMEMS_Richard_2021';
    # download_cora_sss_sst(downloadRoot, 1990, 2021, ftpUser, ftpPass)
    # process_cora_sss_sst(dataRoot, downloadRoot, startYear=1990, endYear=2021)
