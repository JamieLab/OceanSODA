#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:23:46 2020

@author: verwirrt
"""

import numpy as np;
from netCDF4 import Dataset;

nc = Dataset("algorithms/LeeEtAl2006.nc", 'r');
L06names = nc.variables['regionName'][:];
L06coefs = nc.variables['coefficients'][:];
L06tempRange = nc.variables['temperatureRange'][:];
L06salRange = nc.variables['salinityRange'][:];
L06rms2 = nc.variables['rmsds'][:] ** 2;
L06region = nc.variables['regionNo'][:];
L06backup = nc.variables['backupRegionNo'][:];
L06backup2 = nc.variables['secondBackupRegionNo'][:];


ilat = None; jlon = None; #These are the grid-rounded lon/lat values for each row where complete data exists
lat = lon = None; #These are the lat/lon of each complete row
data = matchupData = None; #Original matchup dataset
yearPos = None; #index of the year column in data/matchupData
atPos = None; #Index of AT column
sssPos = None;
sstPos = None;
LeeEtAlStartYear = 1999;
verbose=None;

def Lee06_at(w0, atVariance): # never refit
    '''calculate error in Lee et al 2006 At algorithm'''
    regions = L06region[ilat, jlon] #Creates a list of region numbers #(one for each row)
    backup = L06backup[ilat, jlon] #Creates a list of region numbers #(one for each row)
    backup2 = L06backup2[ilat, jlon] #Creates a list of region numbers #(one for each row)
    
    #store list of indices of non-zero rows
    if 'LeeEtAlStartYear' in globals():
        validFit = np.nonzero(np.logical_and(
           data[w0, yearPos] >= LeeEtAlStartYear,
           L06coefs[regions, 0] > -990.))[0]
    else:
        validFit = np.nonzero(L06coefs[regions, 0] > -990.)[0]
    w = w0[validFit] #select row numbers only
    nw = len(w) #number of data points
    if nw == 0:
        print('No Lee06 data, max year =', data[w0, yearPos].max())
        return []
    
    #subset all the data with just valid rows
    lats = lat[validFit]
    lons = lon[validFit]
    regions = regions[validFit]
    backup = backup[validFit]
    backup2 = backup2[validFit]
    talk = data[w, atPos]
    ssss = data[w, sssPos]
    ssts = data[w, sstPos]
    
    #sanity check
    if ssts.min() < -99.:
        raise ValueError('Lee06 low SST',ssts.min())
        
    #For each row, index = indices of the current row
    for index in range(nw): # check for region changes
        
        ##### Update region based on various SST and SSS conditions
        #if SST > threshold, set region to backup region
        if ssts[index] > L06tempRange[regions[index], 1]:
            regions[index] = backup[index]
        #if SST < another threshold, run some other logic to find the correct region (equation) and update region
        elif ssts[index] < L06tempRange[regions[index], 0]:
            if regions[index] in [1,2]:
                regions[index] = backup[index]
                if ssts[index] < L06tempRange[regions[index], 0]:
                    regions[index] = backup2[index]
        if regions[index] == 2 and ssts[index] > L06tempRange[1, 0] and ssss[index] > L06salRange[2, 1]:
            regions[index] = 1
            
    #all rows that are within their region specific temperature and salinity ranges
    inRange = np.logical_and.reduce((ssts >= L06tempRange[regions, 0],
       ssts <= L06tempRange[regions, 1], ssss >= L06salRange[regions, 0],
       ssss <= L06salRange[regions, 1]))
    
    #subset validFit, to give only rows which are complete and within their region's SSS and SST ranges
    validFit = validFit[inRange]
    w = w[inRange]
    regions = regions[inRange]
    talk = talk[inRange]
    ssts = ssts[inRange]
    ssss = ssss[inRange]
    lons = lons[inRange]
    #subset data to be both valid and in range (already subsetted based on inRange values previously...)
    
    #For each row, get the coefficients used for the relevant equation
    coefs = L06coefs[regions, :]
    
    #For each row get it's respective RMSD
    fitRmsd2s = L06rms2[regions]
    
    #Calculate the model total alkalinity
    modelTalk = (coefs[:, 0] + coefs[:, 1] * (ssss - 35) +
       coefs[:, 2] * (ssss - 35) ** 2 +
       coefs[:, 3] * (ssts - coefs[:, 4]) +
       coefs[:, 5] * (ssts - coefs[:, 4]) ** 2 +
       coefs[:, 6] * (ssts - coefs[:, 4]) * lons)
    
    # NB should we subtract in situ variance from rms2?
    #Calculate the squared errors using the row by row list of (from original fit) RMSDs
    if hasattr(atVariance, '__len__'):
        squaredErrors = fitRmsd2s + atVariance[validFit] # add SST and SSS errors if known
    else:
        squaredErrors = fitRmsd2s + atVariance # add SST and SSS errors if known
    
    return printFit('AT', 'Lee06', [sssName, sstName], modelTalk, talk,
       yearMonth.min(), yearMonth.max(), squaredErrors, w, nTalk,
       [coefs, ssss, ssts, talk, modelTalk])