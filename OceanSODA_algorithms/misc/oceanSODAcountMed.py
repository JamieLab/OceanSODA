#! /usr/bin/env python
'''utility to read in OceanSODA data and count Mediterranean measurements'''

from __future__ import print_function
import numpy as np, math, csv, sys, os, shutil, argparse, operator, cPickle
from datetime import datetime, timedelta
from calendar import timegm  # returns seconds since 1/1/1970 00:00
from scipy.interpolate import interp1d
from ftplib import FTP
from netCDF4 import Dataset
from glob import glob
from point import Point, gcd, angles, project
from region import regionGlobals, Region
from fCO2topCO2 import ftopCO2 #fCO2topCO2
from xCO2topCO2 import xCO2topCO2
from rhoT import rhoT
from matplotlib import pyplot as plt

wMed = np.zeros((180, 360))
eMed = np.zeros((180, 360))
with Dataset('IHO_World_Seas.nc') as nc:
    for sea in 'Alboran Sea,Balearic Sea,Mediterranean Sea - Western Basin,'\
            'Tyrrhenian Sea'.split(','):  # Alboran includes a snippet of E Atlantic
        data = nc.variables[sea][0]
        wMed += data
    for sea in 'Aegean Sea,Adriatic Sea,'\
            'Mediterranean Sea - Eastern Basin'.split(','):
        data = nc.variables[sea][0]
        eMed += data

wMed[eMed > wMed] = 0.  # make mutually exclusive
eMed[wMed > eMed] = 0.
kg2bar = 9.81 / 101325.
firstYear = 0
lastYear = np.inf
verbYear = 0
defaultHalfwidth = timedelta(5.).total_seconds()
defaultRadius = 50.
carbonateParams = 'PCO2W,DPCO2,AT,DIC,PH'.split(',')
nMedW = {param:0L for param in carbonateParams}
nMedE = {param:0L for param in carbonateParams}

def addMed(param, line, pos):
    '''check whether lat, lon = param[pos] is in the Med, if so add to param counter'''
    lat, lon = line[pos[0]], line[pos[1]]
    if len(lat) == 0: return
    if lat[-1] in'NS':
        if lat[-1] == 'S': lat = '-' + lat
        lat = lat[:-1]
    if lon[-1] in'EW':
        if lon[-1] == 'W': lon = '-' + lon
        lon = lon[:-1]
    ilat, jlon = int(90. - float(lat)), int(float(lon) + 180.) % 360
    try:
        if wMed[ilat, jlon]:
            nMedW[param] += 1
        elif eMed[ilat, jlon]:
            nMedE[param] += 1
    except Exception as e: raise ValueError(e, lat, lon, param, pos, line)

nParams = len(carbonateParams)
SOCATfile = 'insitu/SOCAT/SOCATv2019_reanalysed_subskin.tsv'
SOCAToutFile = 'insitu/SOCAT/SOCATv2019_reanalysed_subskinFiltered.csv'
GLODAPv2file = 'insitu/GLODAPv2/GLODAPv2.2019_Merged_Master_File.csv'
GLODAPv2outFile = 'insitu/GLODAPv2/GLODAPv2.2019filtered.csv'
LDEOfile = 'insitu/LDEO/LDEO_Database_V2018.csv'
LDEOoutFile = 'insitu/LDEO/LDEOfiltered.csv'
AMTfiles = sorted(glob('insitu/AMT/AMT*_carbonate_*.csv'))
AMTpCO2files = sorted(glob('insitu/AMT/pCO2_AMT*.csv'))
AMToutFile = 'insitu/AMT/AMTdata.csv'
AMTpCO2outFile = 'insitu/AMT/AMTpCO2data.csv'
ARGOoutFile = 'insitu/BIO-ARGO/bioArgoData.csv'
yearlyDir = 'insitu/yearlyData/'
regionFile = 'insitu/OceanSODAregions.csv'
outHeader = 'DATETIME,LAT,LON,DIST,DEPTH,SAMPLEDEPTH,SST,SSS,PCO2W,PCO2WR,'\
    'PCO2W_QC,PCO2W_ERROR,PCO2WRT,DPCO2,DPCO2R,DPCO2_ERROR,AT,AT_QC,AT_ERROR,'\
    'DIC,DIC_QC,DIC_ERROR,PH,PH_QC,PH_ERROR'.split(',')
carbPos = {param:outHeader.index(param) for param in carbonateParams}
regionsHeader = 'REGION,DATETIME,LAT,LON,N,INNER DATETIME,LAT,LON,HALFWIDTH(S),RADIUS(KM)'.split(',')
headerPos = {name:i for i, name in enumerate(outHeader)}
nOut = len(outHeader)
inSituErrors = {'AT':2.5, 'DIC':2.5, 'PH':.02, 'PCO2W':2., 'DPCO2':np.sqrt(8.)} # nominal 'state-of-art' errors
inSituErrors = {'AT':6., 'DIC':4., 'PH':.02, 'PCO2W':2., 'DPCO2':np.sqrt(8.)} # GLODAPv2 AT and DIC uncertainties
inSituErrorRatios = {'AT':.005, 'DIC':.005, 'PCO2W':.03} # nominal 'state-of-art' errors
inSituSquaredErrors = {}
for item in inSituErrors:
   inSituSquaredErrors[item] = inSituErrors[item] ** 2
PCO2errors = {'A':2., 'B':2., 'C':5., 'D':5., '':5.}
# LDEO quoted average uncertainty is 2.5, but this is MUCH less rigorous than SOCAT!
# Precautionary value to avoid preferring LDEO over SOCAT when weighting
LDEOerror = 5.
useErrorRatios = True
warnDiff = {'LAT':np.inf, 'LON':np.inf, 'DIST':1., 'DEPTH':1., 'SAMPLEDEPTH':1.,
    'SST':.5, 'SSS':.5, 'PCO2W':2., 'DPCO2':2., 'AT': 1., 'AT_QC':.01,
    'DIC':1., 'DIC_QC':.01, 'PH':.01, 'PH_QC':.01, 'PH_ERROR':np.inf}
failDiff = {'LAT':np.inf, 'LON':np.inf, 'DIST':10., 'DEPTH':10., 'SAMPLEDEPTH':5.,
    'SST':1., 'SSS':1., 'PCO2W':5., 'DPCO2':5., 'AT': 5., 'AT_QC':.01,
    'DIC':5., 'DIC_QC':.01, 'PH':.01, 'PH_QC':.01, 'PH_ERROR':np.inf}

parser = argparse.ArgumentParser(
    description = 'download/read OceanSODA data and create matchups database')
parser.add_argument('-R','--read', action = 'store_true',
   help="read all original data, download latest ARGO")
parser.add_argument('-S','--SOCAT', action = 'store_true',
   help="read SOCAT original data")
parser.add_argument('-L','--LDEO', action = 'store_true',
   help="read LDEO original data (NB will also read SOCAT)")
parser.add_argument('-G','--GLODAP', action = 'store_true',
   help="read GLODAPv2018 original data")
parser.add_argument('-M','--AMT', action = 'store_true',
   help="read ARGO original data")
parser.add_argument('-A','--ARGO', action = 'store_true',
   help="read ARGO original data")
parser.add_argument('-U','--updateARGO', action = 'store_true',
   help="download latest ARGO original data")
parser.add_argument('-d','--divide', action = 'store_true',
   help="divide original data into years")
parser.add_argument('-u','--uncertainty', action = 'store_true',
   help="calculate uncertainties")
parser.add_argument('-s','--sort', action = 'store_true',
   help="sort original data by datetime, lat and lon")
parser.add_argument('-c','--collate', action = 'store_true',
   help="collate sorted original data, combining coincident lines")
parser.add_argument('-r','--reconstruct', action = 'store_true',
   help="reconstruct missing carbonate data in each line")
parser.add_argument('-f','--find', action = 'store_true',
   help="find regions")
parser.add_argument('-n','--noclobber', action = 'store_true',
   help="update existing regions")
parser.add_argument('-t','--test', action = 'store_true',
   help="run in test mode with no output")
parser.add_argument('-v', '--verb', action = 'store_true',
    help="verbose output")
args = parser.parse_args()
verb = args.verb
reading = args.read
splitting = args.divide
sorting = args.sort
uncertainty = args.uncertainty
collating = args.collate
reconstructing = args.reconstruct
finding = args.find
noclobber = args.noclobber
test = args.test
if not (reading or splitting or sorting or uncertainty or collating or finding
        or args.GLODAP or args.LDEO or args.SOCAT or args.AMT or args.ARGO or
        args.updateARGO):  # nothing = everything
    reading = True
    splitting = True
    sorting = True
    uncertainty = True
    collating = True
    finding = True
GLODAP = reading or args.GLODAP
LDEO = reading or args.LDEO
SOCAT = reading or args.SOCAT or LDEO
AMT = reading or args.AMT
updateARGO = args.updateARGO
ARGO = reading or args.ARGO or updateARGO
# might make these command line params at some point
radiusKm, halfwidthSeconds = regionGlobals()
relSecond = 1. / halfwidthSeconds
relKm = 1. / radiusKm

def fillDate(myDate):
    '''zero-pads dates such as 5/5/1964 to 05/05/1964 to keep strptime happy'''
    d, m, y = myDate.split('/')
    if int(d) <= 9:
        d = '0' + d
    if int(m) <= 9:
        m = '0' + m
    return '/'.join([d, m, y])

def resid(x):
    x = float(x)
    return timedelta(x - int(x))

start = datetime.now()
print('started', start)
if SOCAT:
    # read in SOCAT data
    print('SOCAT')
    with open(SOCATfile, 'U') as f:
        reader = csv.reader(f, delimiter = '	')
        second = False
        for header in reader:
            if len(header) == 0: continue
            if header[0] == 'Expocode':
                if second: break
                second = True
        depthPos = header.index('sample_depth [m]')
        pCO2pos = header.index('pCO2_reanalysed [uatm]')
        QCpos = header.index('QC_Flag')
        QCpos2 = header.index('fCO2rec_flag')
        if LDEO:
            SOCATpos = []  # sort SOCAT by date/lat/lon for comparison with LDEO
        nBadSOCAT = 0L
        for line in reader:
            for pos, item in enumerate(line):
                if item.lower() == 'nan': line[pos] = ''
            if line[depthPos] != '' and float(line[depthPos]) > 10.:
                if verb: print(line[depthPos])
                continue
            try:
                f = float(line[pCO2pos])
                if f < 0.:
                    nBadSOCAT += 1L
                    continue
            except ValueError:
                nBadSOCAT += 1L
                continue
            if LDEO:  # include QC E+ in LDEO matchups to eliminate bad data from LDEO
                y, m, d, h, M = [int(item) for item in line[4:9]]
                S, lon, lat = [float(item) for item in line[9:12]]
                try:
                    if S == 60.:  # WTF?
                        d = datetime(y, m, d, h, M) + timedelta(seconds = 60)
                    else:
                       d = datetime(y, m, d, h, M, int(S), int((S - int(S)) *  1e6))
                except Exception as e:
                    raise ValueError(e,y,m,d,h,M,S)
                SOCATpos += [[timegm(d.timetuple()), lat, lon]]
            if line[QCpos] not in 'ABCD':
                if verb: print(line[QCpos])
                continue
            if line[QCpos2] != '2': raise ValueError(line[QCpos], line[QCpos2], line)
            addMed('PCO2W', line, [11, 10])
    if LDEO:
        nSOCAT = len(SOCATpos)
        mSOCAT = nSOCAT - 1
        SOCATpos.sort(key = operator.itemgetter(*[0, 1, 2]))
        SOCATpos = np.transpose(SOCATpos)
    print(nBadSOCAT, 'bad SOCAT data')
    print('W Med', nMedW)
    print('E Med', nMedE)

if GLODAP:
    # read in GLODAPv2 data
    print('GLODAPv2', datetime.now() - start)
    with open(GLODAPv2file) as f:
        reader = csv.reader(f)
        header = reader.next()
        tco2Pos = header.index('tco2')
        talkPos = header.index('talk')
        phPos = header.index('phts25p0')
        depthPos = header.index('depth')
        nBadGLODAP = 0L
        for line in reader:
            for pos, item in enumerate(line):
                if item == '-9999': line[pos] = ''
            if line[depthPos] == '':
                nBadGLODAP += 1
                continue
            if float(line[depthPos]) > 10.:
                continue
            bad = True
            for param, pos in zip(['AT', 'DIC', 'PH'], [tco2Pos, talkPos, phPos]):
                try:
                    f = float(line[pos])
                    if f >= 0:
                        bad = False
                        addMed(param, line, [8, 9])
                except ValueError:
                    continue
            if bad:
                nBadGLODAP += 1L
                continue
    print(nBadGLODAP, 'bad GLODAP data')
    print('W Med', nMedW)
    print('E Med', nMedE)

if LDEO:
    # read in LDEO data
    print('LDEO', datetime.now() - start)
    with open(LDEOfile) as f:
        reader = csv.reader(f)
        header = reader.next()
        header = reader.next()
        pCO2pos = header.index('PCO2_SST')
        f = interp1d(SOCATpos[0], xrange(nSOCAT))
        dtMax = 1.  # find out how far SOCAT and LDEO can differ and still be the same point
        distMax = .5
        nLines = 0
        nWonkyDates = 0
        nBadLDEO = 0L
        for line in reader:
            if verb:
                if nLines % 100000 == 0:
                    print('line', nLines, datetime.now() - start)
                nLines += 1
            try:
                x = float(line[pCO2pos])
                if x < 0.:
                    if nBadLDEO == 0: print('bad pCO2', line[pCO2pos])
                    nBadLDEO += 1L
                    continue
            except ValueError:
                if nBadLDEO == 0: print('bad pCO2', line[pCO2pos])
                nBadLDEO += 1L
                continue
            lat, lon = [float(x) for x in line[2:4]]
            dateStrs = line[4].split('/')
            year = int(dateStrs[2])
            jd0 = float(line[5]) - 1.
            try:
                oldD = datetime(year, int(dateStrs[0]), int(dateStrs[1]))
                goodDate = True
            except ValueError:
                goodDate = False
            d = datetime(year, 1, 1) + timedelta(jd0)
            dD = (d - oldD).days
            if goodDate and abs(dD) > 1:
                goodDate = False
                if dD > 185:
                    year += 1
                    d = datetime(year, 1, 1) + timedelta(jd0)
                if dD < -185:
                    year -= 1
                    d = datetime(year, 1, 1) + timedelta(jd0)
            newDate = '/'.join([str(d.month), str(d.day), str(d.year)])
            if not goodDate:
                if verb and line[4][0] != '0' and dateStrs[1] != '0':
                    print(line[4], newDate)
                line[4] = newDate  # overwrite bad date
                nWonkyDates += 1
            myDate = timegm((d).timetuple())
            if myDate < SOCATpos[0, 0] or myDate > SOCATpos[0, -1]:
                addMed('PCO2W', line, [2, 3])  # outside SOCAT date range
                continue
            printing = newDate == '10/23/1957'
            pos0 = int(round(f(myDate)))
            dt = abs(SOCATpos[0, pos0] - myDate)
            if printing: print(line, pos0, SOCATpos[:, pos0], myDate, dt)
            if dt > 30.:  # time not found in SOCAT
                addMed('PCO2W', line, [2, 3])
                continue
            pos = pos0  # now find nearest SOCAT data at the same time
            minDist = 1.  # minimum distance found in km
            down = True
            while True:
                try:
                    dist = gcd([lat, lon], SOCATpos[1:3, pos], latLon = True)
                except Exception as e:
                    raise ValueError(e, pos, pos0, SOCATpos.shape, nSOCAT,
                        SOCATpos[0,[0,-1]], nLines, line)
                if printing: print(pos, [lat, lon], SOCATpos[1:3, pos], down, dist, minDist, minPos, minDt)
                if dist < minDist:
                    minDist = dist
                    minPos = pos
                    minDt = dt
                    if minDt > 30.: raise ValueError(minDist,minPos,minDt)
                if down:
                    if pos not in [0, pos0]:
                        dt = abs(SOCATpos[0, pos - 1] - myDate)
                    if pos == 0 or dt > 30.:
                        down = False
                        pos = pos0
                        if pos >= mSOCAT: break
                        dt = abs(SOCATpos[0, pos + 1] - myDate)
                        continue
                    pos -= 1
                else:
                    if pos >= mSOCAT: break
                    dt = abs(SOCATpos[0, pos + 1] - myDate)
                    if dt > 30.: break
                    pos += 1
            if minDist == 1.:  # position not found in SOCAT
                addMed('PCO2W', line, [2, 3])
                continue
            if minDt > dtMax or minDist > distMax:
                if verb: print(myDate, lat, lon, SOCATpos[:, minPos], minDt, minDist)
                dtMax = max(dtMax, minDt)
                distMax = max(distMax, minDist)
    print(nWonkyDates, 'wonky dates')
    print(nBadLDEO, 'bad LDEO data')
    print('W Med', nMedW)
    print('E Med', nMedE)

if AMT:
    # read in AMT data
    print('AMT', datetime.now() - start)
    nBadAMT = 0L
    for AMTfile in AMTfiles:
        with open(AMTfile) as f:
            reader = csv.reader(f)
            while True:
                header = reader.next()
                if header[0] == 'CRUISE': break
            tco2Pos = header.index('DIC_umol_kg')
            talkPos = header.index('TA_umol_kg')
            phPos = header.index('pH in situ')  # total scale
            depthPos = header.index('Depth')
            blankLines = False
            for line in reader:
                if line[0] == '':
                    blankLines = True
                    continue
                if blankLines: raise ValueError()
                for pos, item in enumerate(line):
                    if item in ['L', 'ND']: line[pos] = ''
                if line[depthPos] == '':
                    nBadAMT += 1
                    if verb: print('bad AMT depth', AMTfile, line[depthPos])
                    continue
                if float(line[depthPos]) > 10.: continue
                bad = True
                for pos in [tco2Pos, talkPos, phPos]:
                    try:
                        f = float(line[pos])
                        if f >= 0: bad = False
                    except ValueError:
                        continue
                if bad:
                    nBadAMT += 1L
                    if verb: print('bad AMT', AMTfile, line)
                    continue
                for param, pos in zip(['AT', 'DIC', 'PH'], [tco2Pos, talkPos, phPos]):
                    addMed(param, line, [5, 6])
        if verb: print(nBadAMT, AMTfile)
    print(nBadAMT, 'bad AMT data')
    nBadAMT = 0L
    for AMTfile in AMTpCO2files:
        with open(AMTfile) as f:
            reader = csv.reader(f)
            while True:
                header = reader.next()
                if header[0] == 'JD_GMT': break
            pCO2pos = header.index('pCO2_sw[uatm]')
            nBadAMT = 0L
            for line in reader:
                try:
                    f = float(line[pCO2pos])
                    if f < 0.:
                        nBadAMT += 1L
                        continue
                except ValueError:
                    nBadAMT += 1L
                    continue
                addMed('PCO2W', line, [3, 4])
    print(nBadAMT, 'bad AMT pCO2 data')
    print('W Med', nMedW)
    print('E Med', nMedE)

if ARGO:
    # read in BIO-ARGO data
    print('ARGO', datetime.now() - start)
    with open(ARGOoutFile) as f:
        reader = csv.reader(f)#, delimiter = ';')
        header = reader.next()
        for line in reader:
            addMed('PH', line, [2, 3])

    print('BIO-ARGO')
    print('W Med', nMedW)
    print('E Med', nMedE)

print('finished filtering', datetime.now() - start)

