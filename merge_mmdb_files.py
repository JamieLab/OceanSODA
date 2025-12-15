#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Created on 8 October 2025

# Merging bonus MMDB files with the orginial files
# @author: Daniel Ford

import numpy as np;
from netCDF4 import Dataset
import os
import geopandas as gpd
from shapely.geometry import Point
from math import cos, sin, asin, sqrt, radians

def calc_distance(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees):
    from: https://stackoverflow.com/questions/4913349/haversine-formula-in-python-bearing-and-distance-between-two-gps-points/4913653#4913653
    """
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) # convert decimal degrees to radians
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2  #haversine formula
    c = 2 * asin(sqrt(a))
    km = 6371 * c
    return km

def calc_distance_to_coastline(longitude,latitude,coastline):
    target_coordinate=Point(longitude,latitude )
    return coastline.distance(target_coordinate).values[0]

def distance_degrees_to_kilometers(distance,coord=[0,0]):
    coord_plus=[c+distance for c in coord]
    coord_minus=[c-distance for c in coord]
    return (calc_distance(*coord,*coord_plus)+calc_distance(*coord,*coord_minus))*0.5

def calc_distance_to_coastline_km(longitude,latitude,coastline):
    target_coordinate=Point(longitude,latitude )

    return distance_degrees_to_kilometers(coastline.distance(target_coordinate).values[0],[longitude,latitude])

def merge_mmdb(mmdb,mmdb_extra,start_yr = 1980,end_yr=2024,extra_name = ''):
    for i in range(start_yr,end_yr+1):
        #Check if mmdb file exists
        mmdb_file = mmdb.replace('%Y',str(i))
        mm = os.path.exists(mmdb_file)
        #Check if extra file exists
        mmdb_file_extra = mmdb_extra.replace('%Y',str(i))
        mm_e = os.path.exists(mmdb_file_extra)
        #If both exist merge the extra netcdf variables into the mmdb with the extra name.
        if mm & mm_e:
            print('Both files exist - copying variables for year: '+str(i))
            c = Dataset(mmdb_file,'a')
            d = Dataset(mmdb_file_extra,'r')
            keys = list(d.variables.keys())
            keys.remove('lat'); keys.remove('lon'); keys.remove('region_id'); keys.remove('time');
            for j in range(len(keys)):
                print('Copying: ' + keys[j])
                if extra_name+keys[j] not in list(c.variables.keys()):
                    print('Variable doesnt already exist...')
                    var = c.createVariable(extra_name+keys[j],'double',('region_id'),fill_value=d.variables[keys[j]]._FillValue)
                data = d.variables[keys[j]][:]
                c.variables[extra_name+keys[j]][:] = data
                diction = d[keys[j]].__dict__
                diction.pop('_FillValue')
                c.variables[extra_name+keys[j]].setncatts(diction)
            c.close()
            d.close()
        else:
            print('Both files dont exist - skipping year: ' +str(i))
                #var._FillValue = np.nan

def add_bathymetry(mmdb,bathymetry_file,start_yr=1980,end_yr=2024):
    c = Dataset(gebco_bathymetry_file,'r')
    lon = np.array(c['lon'])
    lat = np.array(c['lat'])
    depth = np.array(c['elevation'])
    c.close()

    for i in range(start_yr,end_yr+1):
        mmdb_file = mmdb.replace('%Y',str(i))
        mm = os.path.exists(mmdb_file)
        if mm:
            print('Adding bathymetry information for year: '+str(i))
            c = Dataset(mmdb_file,'a')
            lon_m = np.array(c['lon'])
            lat_m = np.array(c['lat'])
            coast = np.zeros((len(lon_m))); coast[:] = np.nan

            for j in range(len(coast)):
                g = np.where(np.abs(lon_m[j] - lon) == np.min(np.abs(lon_m[j] - lon)))[0]
                f = np.where(np.abs(lat_m[j] - lat) == np.min(np.abs(lat_m[j] - lat)))[0]
                coast[j] = depth[f[0],g[0]]

            if 'elevation' not in list(c.variables.keys()):
                print('Variable doesnt already exist...')
                var = c.createVariable('elevation','double',('region_id'),fill_value=np.nan)

            c.variables['elevation'][:] = coast
            c.variables['elevation'].long_name = 'Elevation from GEBCO'
            c.variables['elevation'].Units = 'm'
            c.variables['elevation'].direction = 'Negative is below sea level'
            c.variables['elevation'].source = gebco_bathymetry_file
            c.close()

def add_coastline(mmdb,coastline_file,start_yr=1980,end_yr=2024):
    coastline_data = gpd.read_file(natural_earth_coastline_file)
    coastline = gpd.GeoSeries(coastline_data.geometry.unary_union)

    for i in range(start_yr,end_yr+1):
        mmdb_file = mmdb.replace('%Y',str(i))
        mm = os.path.exists(mmdb_file)
        if mm:
            print('Adding coastline information for year: '+str(i))
            c = Dataset(mmdb_file,'a')
            lon_m = np.array(c['lon'])
            lat_m = np.array(c['lat'])
            coast = np.zeros((len(lon_m))); coast[:] = np.nan

            for j in range(len(coast)):
                coast[j] = calc_distance_to_coastline_km(lon_m[j],lat_m[j],coastline)


            if 'distance_to_coast' not in list(c.variables.keys()):
                print('Variable doesnt already exist...')
                var = c.createVariable('distance_to_coast','double',('region_id'),fill_value=np.nan)

            c.variables['distance_to_coast'][:] = coast
            c.variables['distance_to_coast'].long_name = 'Distance to coastline with respect to land mask'
            c.variables['distance_to_coast'].Units = 'km'
            c.variables['distance_to_coast'].source = coastline_file
            c.close()

mmdb_file_struct = 'F:/OceanSODA/data/matchup_datasets/oceansoda/OCEANSODA-MMDB-%Y-fv01.nc'
natural_earth_coastline_file = './data/coastline/ne_10m_coastline/ne_10m_coastline.shp'
gebco_bathymetry_file = 'F:/Data/Bathymetry/GEBCO_2023.nc'

# extra_file_struct = 'F:/OceanSODA/data/matchup_datasets/oceansoda/bonus/ethz-monthly-1deg/OCEANSODA-MMDB-%Y-fv01.nc'
# extra = 'ETHZ_OHOA_month_'
#
# merge_mmdb(mmdb_file_struct,extra_file_struct,extra_name = extra)
#
#
# extra_file_struct = 'F:/OceanSODA/data/matchup_datasets/oceansoda/bonus/pml-hrcoastal-monthly-1deg/OCEANSODA-MMDB-%Y-fv01.nc'
# extra = 'PML_OHOA_month_'
#
# merge_mmdb(mmdb_file_struct,extra_file_struct,extra_name = extra)
#
#
# extra_file_struct = 'F:/OceanSODA/data/matchup_datasets/oceansoda/bonus/pml-hrcoastal/OCEANSODA-MMDB-%Y-fv01.nc'
# extra = 'PML_OHOA_'
#
# merge_mmdb(mmdb_file_struct,extra_file_struct,extra_name = extra)

# add_bathymetry(mmdb_file_struct,gebco_bathymetry_file)
add_coastline(mmdb_file_struct,natural_earth_coastline_file)
