# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 12:42:46 2015

@author: thecakeisalie

Maintained by Daniel Hueholt
Last updated: 5/23/2019
"""
import glob
import os
import numpy as np
import re
from datetime import datetime, timedelta
#import pymap3d
import netCDF4 as nc
from scipy.ndimage.filters import minimum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from math import radians, cos, sin, asin, atan2, sqrt
import pandas as pd
import nexradaws

def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)

def get_AWSlist(start_datetime, end_datetime, radar):
    """
    Function to get a list of AWS files in a given range for a given radar

    :param start_datetime: datetime of start time
    :param end_datetime: datetime of end time
    :param radar: string of radar
    :return: list of AWS file objects
    """

    # establish connection
    conn = nexradaws.NexradAwsInterface()

    # get list of available files
    availscans = conn.get_avail_scans_in_range(start_datetime, end_datetime, radar)

    # remove files that aren't readable
    availscans = [file for file in availscans if 'MDM' not in file.create_filepath('', True)[1]]

    # convert to list of filepaths to read in with pyart
    files = [scan.create_filepath(basepath='s3://noaa-nexrad-level2/', keep_aws_structure=True)[1] for scan in availscans]
    files = [path.replace("\\", '/') for path in files]

    return files


def get_azimuth(radar, sweepnum):
    """
    DESCIPTION: For RHI's, extracts the azimuth angle from the radar object.
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    sweepnum = Natural number, usually 0 or 1 (for up sweep and down sweep)
    
    OUTPUTS:
    azi = The azimuth angle for the secifies sweep number.
    """
    if sweepnum==0:
        azi = radar.azimuth['data'][0]
    else:
        azi = radar.azimuth['data'][np.size(radar.azimuth['data'])-1]
    return azi

def get_filelist(inpath, wildcard, savefile):
    """
    DESCRIPTION: Generates a list of files that contain the specified wildcard 
        file their file name that in the specified path.
        
    INPUTS:
    inpath = A string that specifies the path to the desired files. 
    wildcard = A string that sets the phrase which is common to all desired
        files. The files in the inpath directory will be sorted using this 
        variable.
    savefile = A boolean value where True will save off a text files containing
        the names of the desired files separated by a newline, "\n".
        
    OUPUTS:
    filelist = A list of desired file names.
    """
    # Determine all files in inpath that begin with wildcard
    os.chdir(inpath)
    ident = '*' + wildcard
    identifier = ident + '*'
    filenamelist = glob.glob(identifier)
    
    # Sort the list
    filelist = sorted(filenamelist)
    
    if savefile == True:
        # Save the list
        name = "filelist_%s.txt" % wildcard
        txtfile = open(name, "w")
        for item in filelist:
            txtfile.write("%s\n" % item)
        txtfile.close()
        
    return filelist

def get_filelist_range(inpath, starttime, endtime, dataType):
    
    if dataType in ('stitched', 'subset', 'interp', 'level2'):
        dateFmt = '%Y%m%d_%H%M%S'
        regStr = '\d{8}_\d{6}' # YYYYMMDD_hhmmss
    elif dataType == 'imgs':
        dateFmt = '%Y%m%d%H%M%S'
        regStr = '\d{14}'
    elif dataType == 'era5':
        dateFmt = '%Y%m%d'
        regStr = '\d{8}' # YYYYMMDD
    else:
        print("Invalid data type")
        
    # Convert start and end times to datetime
    start_datetime = datetime.strptime(starttime, '%Y%m%d%H%M%S')
    end_datetime   = datetime.strptime(endtime,   '%Y%m%d%H%M%S')
    
    # Find string of 4 capital letters to identify radar (just to make sure we're getting the right files)
    #radar = re.search('[A-Z]{4}', inpath).group()
    
    # List files
    allFiles = glob.glob(inpath + '*' + start_datetime.strftime('%Y') + '*')
    
    if len(allFiles) == 0:
        filelist = []
    else:
        allFiles = sorted(allFiles) # sort files
        regex = re.compile(regStr) # extract date from files
        fileTimes = [regex.search(file).group() for file in allFiles]
        fileDatetimes = [datetime.strptime(time, dateFmt) for time in fileTimes] # convert to datetime
        
        # Find indices
        start_idx = nearest_ind(fileDatetimes, start_datetime)
        end_idx = nearest_ind(fileDatetimes, end_datetime)
        
        # if the first file is more than an hour from the start time
        if start_idx==end_idx: # ((fileDatetimes[start_idx]-start_datetime).total_seconds())/3600 > 1 or 
            filelist = []
        else:
            # Get subset and check if any are repeated
            subsetTimes = fileDatetimes[start_idx:end_idx+1]
            if (len(subsetTimes)-len(set(subsetTimes))) > 50:
                repeatFlag = True
            else:
                repeatFlag = False
            
            subsetList = allFiles[start_idx:end_idx+1]
            exts = ['.gz', '.Z']#, '_9', '_10', '_7', '_8']
            
            if all(file.endswith(tuple(exts)) for file in subsetList):
                filelist = [x for x in subsetList if not x.endswith('.Z')]
            else:
                if repeatFlag:
                    filelist = [x for x in subsetList if x.endswith('.gz')]
                else:
                    filelist = subsetList
            
    return filelist

def getDatetimes(filelist, dataType):

    if dataType in ('stitched', 'subset', 'interp', 'level2'):
        dateFmt = '%Y%m%d_%H%M%S'
        regStr = '\d{8}_\d{6}' # YYYYMMDD_hhmmss
    elif dataType == 'imgs':
        dateFmt = '%Y%m%d%H%M%S'
        regStr = '\d{14}'
    elif dataType == 'era5':
        dateFmt = '%Y%m%d'
        regStr = '\d{8}' # YYYYMMDD
    else:
        print("Invalid data type")
    
    # List files
    allFiles = sorted(filelist) # sort files
    regex = re.compile(regStr) # extract date from files
    fileTimes = [regex.search(file).group() for file in allFiles]
    fileDatetimes = [datetime.strptime(time, dateFmt) for time in fileTimes] # convert to datetime
    
    return fileDatetimes

        
def azi_calculator(azi_lines,max_length):
    """
    DESCRIPTION: Uses trigonometry to calculates coordinates to overlay RHI azimuths
        on a PPI plot.
        
    INPUTS:
        azi_lines: array containing azimuths to be plotted (in degrees, e.g. [134,224])
        max_length: hypoteneuse. For a radar PPI image, this will be the x_lim/y_lim
        
    OUTPUTS:
        x_c: array containing x-coordinates
        y_c: array containing y-coordinates
    
    See Master_plotter for an example of how this is used in practice.
    
    """
    x_c = np.empty(np.size(azi_lines))
    y_c = np.empty(np.size(azi_lines))
    
    for count in range(0,np.size(azi_lines)):
        azi = azi_lines[count]
        math_deg = (450 - azi)%360
        azi_rad = np.deg2rad(math_deg)
        
        x_c[count] = np.cos(azi_rad)*max_length
        y_c[count] = np.sin(azi_rad)*max_length
#    for count in range(0,np.size(azi_lines)):
#        azi = azi_lines[count]
#        if azi > 0 and azi <= 90:
#                azi_rad = np.deg2rad(azi)
#                x_c[count] = np.cos(azi_rad)*max_length
#                y_c[count] = np.sin(azi_rad)*max_length
#        elif azi > 270 and azi <= 360:
#                azi_rad = np.deg2rad(azi-90)
#                x_c[count] = np.cos(azi_rad)*max_length
#                y_c[count] = -np.sin(azi_rad)*max_length
#        elif azi > 180 and azi <= 270:
#                azi_rad = np.deg2rad(azi-180)
#                x_c[count] = -np.cos(azi_rad)*max_length
#                y_c[count] = -np.sin(azi_rad)*max_length
#        elif azi > 90 and azi <= 180:
#                azi_rad = np.deg2rad(azi-270)
#                x_c[count] = -np.cos(azi_rad)*max_length
#                y_c[count] = np.sin(azi_rad)*max_length
            
    return x_c,y_c

def grid_selector(case):
    
    if case == 'NE':
        
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(700/dx),(500/dx)+1,1); ypts = (dy*1000)*np.arange(-(540/dy),(660/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 
    
    elif case in ('NC_IMPACTS', 'IMPACTS_NC'):
    
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(300/dx),(400/dx)+1,1); ypts = (dy*1000)*np.arange(-(300/dy),(400/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case in ('Midwest_IMPACTS', 'IMPACTS_Midwest', 'IMPACTS_Ohio'):
        
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(700/dx),(500/dx)+1,1); ypts = (dy*1000)*np.arange(-(500/dy),(500/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case == 'NC':
        
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(550/dx),(350/dx)+1,1); ypts = (dy*1000)*np.arange(-(500/dy),(400/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case == 'NC_Alt':
        
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(700/dx),(300/dx)+1,1); ypts = (dy*1000)*np.arange(-(450/dy),(450/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case == 'Waves_NE':
        
        dx = dy = 0.5
        xpts = (dx*1000)*np.arange(-(700/dx),(500/dx)+1,1); ypts = (dy*1000)*np.arange(-(540/dy),(660/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case == 'individual':
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(200/dx),(200/dx)+1,1); ypts = (dy*1000)*np.arange(-(200/dy),(200/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 

    elif case == 'ind_highres':
        dx = dy = 0.5
        xpts = (dx*1000)*np.arange(-(200/dx),(200/dx)+1,1); ypts = (dy*1000)*np.arange(-(200/dy),(200/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 
    
    elif case == 'FINESST':
        
        dx = dy = 2
        xpts = (dx*1000)*np.arange(-(500/dx),(700/dx)+1,1); ypts = (dy*1000)*np.arange(-(400/dy),(850/dy)+1,1); zpts = np.zeros([1,len(xpts)]); 
    else:
        
        raise ValueError('Case not recognized')
        
    return xpts, ypts, zpts

def datetime_range(start,end,delta):
    
    saveTimes = []
    
    saveTimes.append(start)
    
    current = start
    
    while current < end:
        
        #yield current
        
        current+=delta
        
        saveTimes.append(current)

    return saveTimes

def radarCase(case):
    
    if case in ('NE', 'Waves_NE'):
        radars = ['KOKX', 'KBOX', 'KDIX', 'KDOX', 'KENX', 'KGYX',
                  'KTYX', 'KBGM', 'KBUF', 'KLWX', 'KCCX', 'KCXX']
        
    elif case in ('NC_IMPACTS', 'IMPACTS_NC', 'NC', 'NC_Alt'):
        radars = ['KRAX', 'KMHX', 'KAKQ', 'KLTX', 'KGSP', 'KFCX']
        
    elif case in ('Midwest_IMPACTS', 'IMPACTS_Midwest'):
        radars = ['KIND', 'KDVN', 'KILX', 'KEAX', 'KLSX',
                  'KLOT', 'KVWX', 'KILN', 'KIWX']

    elif case in ('IMPACTS_Ohio'):
        radars = ['KILN', 'KIND', 'KVWX', 'KLVX', 'KJKL', 
                  'KCLE', 'KILX', 'KRLX', 'KIWX', 'KLOT']
        
    else:
        raise ValueError('Case not recognized')
        
    return radars

def getERA5pres(era5Files, latLims, lonLims):
    
    timeArr = []
    latArr  = []
    lonArr  = []
    presArr = []
    
    # Loop through era5 files:
    for iera, eraFile in enumerate(era5Files):
        
        # read in file
        eraNC = nc.Dataset(eraFile)
        # extract variables
        eraTimes = eraNC.variables['time']
        eraLon   = eraNC.variables['longitude'][:]
        eraLat   = eraNC.variables['latitude'][:]
        # convert to datetime
        eraDatetime = list(nc.num2date(eraTimes[:], eraTimes.units, eraTimes.calendar))
        timeArr.extend(eraDatetime)
        # get limits for analysis
        minLat = nearest_ind(eraLat, latLims[0]); maxLat = nearest_ind(eraLat, latLims[1])
        minLon = nearest_ind(eraLon, lonLims[0]); maxLon = nearest_ind(eraLon, lonLims[1])
        lonSub = eraLon[minLon:maxLon+1]; latSub = eraLat[maxLat:minLat+1]
        # create empty arrays to store data
        mPresLat = np.empty_like(eraDatetime)
        mPresLon = np.empty_like(eraDatetime)
        mPres    = np.empty_like(eraDatetime)
        
        # Loop through each time
        for itime, time in enumerate(eraDatetime):
        
            # subset of Pressure array
            eraMSLPsub = eraNC.variables['msl'][itime,maxLat:minLat+1,minLon:maxLon+1]
            # where is lowest pressure
            mPresInd = np.where(eraMSLPsub == np.min(eraMSLPsub))
            # if there are multiple points, choose closest to previous
            if len(mPresInd[0]) == 1:
                mPresLat[itime] = mPresInd[0][0] 
                mPresLon[itime] = mPresInd[1][0]
            elif len(mPresInd[0]) != 1:
                prevLat = mPresLat[itime-1]; prevLon = mPresLon[itime-1]
                nearestInd = nearest_ind(mPresInd[0], prevLat)
                mPresLat[itime] = mPresInd[0][nearestInd]
                mPresLon[itime] = mPresInd[1][nearestInd]
            # store minimum pressure
            mPres[itime] = np.min(eraMSLPsub)
            
            latArr.append(latSub[mPresLat[itime]])
            lonArr.append(lonSub[mPresLon[itime]])
            presArr.append(mPres[itime])
           
    # empty list for storing track2 info
    datetimeArr = [datetime.strptime(itime.strftime(), itime.format) for itime in timeArr]
    fullLatList = []
    fullLonList = []
    fullTimeList = []
    
    # loop through all points
    for i in np.arange(len(latArr)-1):
        # get list of times every 5 minutes
        times = datetime_range(datetimeArr[i], datetimeArr[i+1], timedelta(minutes=5))
        # compute distance between 2 points every 5 minutes
        latList, lonList = pymap3d.vincenty.track2(latArr[i], lonArr[i], latArr[i+1], lonArr[i+1], pymap3d.ellipsoid.Ellipsoid(model='wgs84'), npts=len(times))
        # save
        fullTimeList.extend(times[:-1])
        fullLatList.extend(latList[:-1])
        fullLonList.extend(lonList[:-1])
    
    # Append last time to arrays
    fullTimeList.append(datetimeArr[-1]); fullLatList.append(latArr[-1]); fullLonList.append(lonArr[-1])
    # convert to numpy arrays
    fullTimeList = np.array(fullTimeList); fullLatList = np.array(fullLatList); fullLonList = np.array(fullLonList)
    # convert longitudes to all same type
    fullLonList = np.array([x if x<0 else x-360 for x in fullLonList])

    return fullTimeList, fullLatList, fullLonList

def detectLocalMinima(image):
    
    neighborhood = generate_binary_structure(2,2)
    
    local_min = minimum_filter(image, 15)==image
    
    background = (image==0)
    
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)
    
    detected_peaks = local_min ^ eroded_background
    
    return detected_peaks

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula
    dlon = lon2-lon1
    dlat = lat2-lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2*asin(sqrt(a))
    # radius of earth 6371 km
    km = 6371 * c
    
    return km

def bearingCalc(lon1, lat1, lon2, lat2):
    
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    
    dL = lon2 - lon1
    
    x = cos(lat2) * sin(dL)
    y = (cos(lat1) * sin(lat2)) - (sin(lat1)*cos(lat2)*cos(dL))
    
    bearing = atan2(x, y)
    
    return np.rad2deg(bearing)

# def trackPoints(track):
#
#     # set non trackable points to NaN
#     track.loc[track['confidence']==0, ['lat', 'lon', 'pres']] = np.nan
#
#         # empty list for storing track2 info
#     datetimeArr = list(pd.to_datetime(track['Time']))
#     fullLatList = []
#     fullLonList = []
#     fullTimeList = []
#
#     # loop through all points
#     for i in np.arange(len(track)-1):
#         # get list of times every 5 minutes
#         times = datetime_range(datetimeArr[i], datetimeArr[i+1], timedelta(minutes=5))
#         # compute distance between 2 points every 5 minutes
#         latList, lonList = pymap3d.vincenty.track2(track['lat'].iloc[i], track['lon'].iloc[i],
#                                                    track['lat'].iloc[i+1], track['lon'].iloc[i+1],
#                                                    pymap3d.ellipsoid.Ellipsoid(model='wgs84'), npts=len(times))
#         # save
#         fullTimeList.extend(times[:-1])
#         fullLatList.extend(latList[:-1])
#         fullLonList.extend(lonList[:-1])
#
#     # Append last time to arrays
#     fullTimeList.append(datetimeArr[-1]); fullLatList.append(track['lat'].iloc[-1]); fullLonList.append(track['lon'].iloc[-1])
#     # convert to numpy arrays
#     fullTimeList = np.array(fullTimeList); fullLatList = np.array(fullLatList); fullLonList = np.array(fullLonList)
#     # convert longitudes to all same type
#     fullLonList = np.array([x if x<0 else x-360 for x in fullLonList])
#
#     dict = {'time':fullTimeList, 'lat': fullLatList, 'lon':fullLonList}
#
#     fullTrack = pd.DataFrame(dict)
#
#     return fullTrack
    
def getWindow(array, lat, lon, center, buffer):
    
    xlim, ylim = array.shape
    
    xMin = center[0] - buffer
    xMax = center[0] + buffer
    
    yMin = center[1] - buffer
    yMax = center[1] + buffer
    
    if xMin < 0: xMin = 0
    if yMin < 0: yMin = 0
    if xMax > (xlim-1): xMax = xlim-1
    if yMax > (ylim-1): yMax = ylim-1
    
    arraySub = array[xMin:xMax, yMin:yMax]
    
    latSub = lat[xMin:xMax]
    lonSub = lon[yMin:yMax]
    
    return arraySub, latSub, lonSub
    
def getStormDates(year='all'):
    
    # .csv file with storm dates
    date_csv = "/home/disk/zathras/ltomkins/python/mckenzie_list_1996_2021.csv"

    # read in file and get start times, add 1 day to get end times
    dates = pd.read_csv(date_csv, header=None, names=['start'], parse_dates=['start'], date_parser=pd.to_datetime)
    
    if year != 'all':
        dates = dates[dates['start'].dt.year == year]
        
    dates['end'] = dates.start + pd.Timedelta(days=1)
    
    # get list of start and end times for era5 fetcher
    storm_times = []
    
    for index, storm in dates.iterrows():
        storm_times.append([storm['start'].to_pydatetime(), storm['end'].to_pydatetime()])
        
    return storm_times
