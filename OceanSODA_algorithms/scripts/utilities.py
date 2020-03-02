#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 10:49:55 2020

Misc. utilities used for loading and processing data

@author: tom holding
"""

from netCDF4 import Dataset;
import numpy as np;
import pandas as pd;
import rpy2.robjects as ro;


#Reads yearly matchup database netCDF files and concatinates them into a single dataframe.
#   years: and iterable of integer years
#   settings: the global settings dictionary
def load_matchup_to_dataframe(settings, commonNames=None):
    #Concatinate each year
    dfList = [];
    
    if commonNames == None:
        commonNames = settings["columnMap"].keys()
    
    for year in settings["years"]:
        matchupNC = Dataset(settings["matchupDatasetTemplate"].safe_substitute(YYYY=year), 'r');
        
        #Create a pandas dataframe from the netCDF file
        df = pd.DataFrame();
        for commonName in commonNames:
            try:
                df[commonName] = matchupNC[settings["columnMap"][commonName]][:];
            except IndexError:
                print("Missing data: ", year, commonName, settings["columnMap"][commonName]);
                #print("Year", year, "has been skipped.");
        dfList.append(df);
    matchupData = pd.concat(dfList, ignore_index=True);
    
    #Convert date from time in seconds cince 1980-01-01 to a pd.datetime object
    matchupData["date"] = convert_time_to_date(matchupData["date"]);
    
    return matchupData;

#Converts time in seconds to data
#   dt: time delta in seconds from the base date (pandas series)
#   baseDate: base date from which the change in time (dt) is from
def convert_time_to_date(dt, baseDate=pd.datetime(1980, 1, 1)):
    dates = baseDate + pd.to_timedelta(dt, 'S');
    return dates;


#Given a dataframe, returns a subset containing only complete rows.
def subset_complete_rows(df):
    selectedRows = np.all(df.notna(), axis=1)
    subset = df.loc[selectedRows];
    return subset;


#subset from a gridded netCDF mask. Makes no assumptions as to shape of mask
def subset_from_mask(data, maskNC, maskName, maskValue=1):
    #Read the mask and lon/lat dimension values
    mask = maskNC.variables[maskName][:];
    lats = maskNC.variables["lat"][:];
    lons = maskNC.variables["lon"][:];
    
    #find index of closest mask lats value to each row
    a, b = np.meshgrid(lats, data["lat"]);
    latIndices = np.abs(a-b).argmin(axis=1);
    a, b = np.meshgrid(lons, data["lon"]); #and again for longitudes
    lonIndices = np.abs(a-b).argmin(axis=1);
    
    subset = data.iloc[mask[(latIndices, lonIndices)] == maskValue];
    return subset;


##subset from a square gridded netCDF mask
#def old_subset_from_mask(data, maskNC, maskName):
#    #Read the mask and lon/lat dimension values
#    mask = maskNC.variables[maskName][:];
#    lats = maskNC.variables["lat"][:];
#    lons = maskNC.variables["lon"][:];
#    
#    #find the min and max longitude based on the mask
#    latRange = np.any(mask, axis=1);
#    minLat = np.floor(lats[latRange].min());
#    maxLat = np.floor(lats[latRange].max());
#    
#    #And again for longitude
#    lonRange = np.any(mask, axis=0);
#    minLon = np.floor(lons[lonRange].min());
#    maxLon = np.floor(lons[lonRange].max());
#    
#    #subset data and return
#    subset = data[(data["lat"]>=minLat) & (data["lat"]<=maxLat) &
#                  (data["lon"]>=minLon) & (data["lon"]<=maxLon)];
#    return subset;

#Subset a dataframe by including any points which fall inside any boxes formed by a list of lon and lat ranges
def subset_from_inclusive_coord_list(includedLons, includedLats, data):
    if len(includedLons) != len(includedLats):
        raise ValueError("There must be an equal number of lon and lat ranges");
    
    #If empty, default to global (no regional subsetting needed)
    if len(includedLons) == 0:
        return data;
    
    toKeep = np.full((len(data,)), False, dtype=bool);
    for i in range(len(includedLons)): #for each pair of lon/lat ranges, keep any datapoints that fall into this range.
        lonMin = includedLons[i][0];
        lonMax = includedLons[i][1]
        latMin = includedLats[i][0];
        latMax = includedLats[i][1];
        toKeep = toKeep | ((data["lon"] >= lonMin) & (data["lon"] <= lonMax) & (data["lat"] >= latMin) & (data["lat"] <= latMax));
    
    subset = data.loc[toKeep];
    return subset;




#Adapted from: https://github.com/iga202/Pathfinders_SeaCarb
def run_sea_carb(tflag, atData, dicData, sssData, sstData, verbose = False):
    '''Creates an R function to arrange data and run the carb function within the SeaCarb package.
    It assumes data are in a python dictionary, with the optional inputs defined by variable name
    At present this only takes the pHscale switch
    Ian Ashton, 22/07/2015
    Added k1k2, kf, ks and b switches
    Peter Land 18/11/15'''
    
    ro.r('library(seacarb)');
    #Perturb the data here (or make a sister function that perturbs the data to be called if necessary.
    ro.r('''
    f<-function(flags,val1,val2,S,T,Patm,P,Pt,Sit,k1k2,kf,ks,b,pHs){
        flags = as.numeric(matrix(data = flags, nrow = length(flags), ncol = 1))
        val1 = as.numeric(matrix(data = val1, nrow = length(val1), ncol = 1))
        val2 = as.numeric(matrix(data = val2, nrow = length(val2), ncol = 1))
        S = as.numeric(matrix(data = S, nrow = length(S), ncol = 1))
        T = as.numeric(matrix(data = T, nrow = length(T), ncol = 1))
        out = carb(flags,val1,val2,S,T,Patm,P,Pt,Sit,k1k2=k1k2,kf=kf,ks=ks,b=b,pHscale=pHs)
        return(out)
     }
     ''');
    
    rCarb = ro.r['f']#define rCarb as above function
    nd = len(indata[var1name])
    out = {}
    fails = []
    if 1: #try:
        res = rCarb(tflag, indata[var1name], indata[var2name], indata['SSS'], indata['SST'], 0, 0, "x", "x", "d", "u74", "T") # Call R function
        for l,v in res.items(): # For each parameter in the output, convert into Python dictionary, out
            out[l] = []
            for j in range(len(v)):
                if not type(v[j]) in [float, int]:
                    raise ValueError(res.items(),l,v,j,v[j])
                out[l].append(v[j])
    if 0: #except:
        out = {}
        nfails = 0
        for index in np.xrange(nd):
            try:
                res = rCarb(tflag, indata[var1name][index],
                    indata[var2name][index], indata['S'][index],
                    indata['T'][index], indata['Patm'], indata['P'],
                    indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],
                    indata['ks'], indata['b'], indata['pHscale']) # Call R function
                for l,v in res.items(): # For each parameter in the output, convert into Python dictionary, out
                    if l in out:
                        out[l].append(v[0])
                    else:
                        out[l] = [v[0]]
                    if not type(v[0]) in [float, int]:
                        raise ValueError(res.items(), l, v)
            except:
                #traceback.print_exc()
                nfails+=1
                fails.append(index)
                for key in out:
                    out[key].append(-99999)
       
                print(nfails, 'fails out of', nd)
                # print('GOOD ->',indata[var1name][index-10],indata[var2name][index-10],indata['S'][index-10],indata['T'][index-10], indata['Patm'], indata['P'],indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],indata['ks'], indata['b'], indata['pHscale'])
                # print('BAD ->',indata[var1name][index],indata[var2name][index],indata['S'][index],indata['T'][index], indata['Patm'], indata['P'],indata['Pt'], indata['Sit'], indata['k1k2'], indata['kf'],indata['ks'], indata['b'], indata['pHscale'])
            if verbose:# or nfails == nd:
                print(nfails, 'fails out of', nd);
    for l in out:
        out[l] = np.array(out[l])
    
    return(out,fails)

