# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 13:53:15 2020

Contains various functions for specific manipulation of PyART Radar objects.

Version Date: 13 May 2020
Author: Laura Tomkins
North Carolina State University
Undergraduate Research Assistant at Environment Analytics

NEXRAD RHI version
"""



import pyart
import numpy as np
from skimage import morphology

#%%
def nearest_ind(items, pivot):
    time_diff = np.abs([date - pivot for date in items])
    return time_diff.argmin(0)

def getDataSweep(radar, field, tilt):
    
    """
    DESCRIPTION: Returns the lowest tilt of a specified field which has data. 
    New NEXRAD scan strategies have some data for some fields in some titls but not all.
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    field = Names of field in the radar object
    tilt = Tilt to use to filter radar object.
    
    OUTPUTS:
    itilt = Index of the lowest tilt for the field with data.
        
    """
   
    elevAngles = radar.fixed_angle['data']
    lowTilts = np.where(np.abs(elevAngles-tilt) == np.min(np.abs(elevAngles-tilt)))[0]
    
    dataTilts = []

    for itilt in np.arange(len(lowTilts)):
        
        tiltSlice = radar.get_slice(lowTilts[itilt])

        data = radar.fields[field]['data'][tiltSlice] 
        
        if ~data.mask.all():     # if NOT all data is masked, append, if all data is masked continue        
            
            dataTilts.append(itilt)
        
    return dataTilts

def returnLowest(radar, field):
    
    """
    DESCRIPTION: Returns the lowest tilt of a specified field which has data. 
    New NEXRAD scan strategies have some data for some fields in some titls but not all.
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    field = Names of field in the radar object
    tilt = Tilt to use to filter radar object.
    
    OUTPUTS:
    itilt = Index of the lowest tilt for the field with data.
        
    """
   
    elevAngles = radar.fixed_angle['data']

    for itilt in np.arange(0, len(elevAngles)):
        
        tiltSlice = radar.get_slice(itilt)

        data = radar.fields[field]['data'][tiltSlice] 
        
        if ~data.mask.all():     # if NOT all data is masked, append, if all data is masked continue        
            
            break
        
    return itilt

def returnAll(radar, field):
    
    """
    DESCRIPTION: Returns the all tilts of a specified field which has data. 
    New NEXRAD scan strategies have some data for some fields in some titls but not all.
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    field = Names of field in the radar object
    
    OUTPUTS:
    itilt = Index of the lowest tilt for the field with data.
        
    """
   
    elevAngles = radar.fixed_angle['data']
    
    goodTilts = []
    goodAngs = []

    for itilt in np.arange(0, len(elevAngles)):
        
        tiltSlice = radar.get_slice(itilt)

        data = radar.fields[field]['data'][tiltSlice] 
        
        if ~data.mask.all():     # if NOT all data is masked, append, if all data is masked continue        
            
            goodTilts.append(itilt)
            goodAngs.append(elevAngles[itilt])
        
    return goodTilts, goodAngs

def removeDuplicate(angles, tilts):
    
    dupTilts = []
    for idx, ang in enumerate(angles):
        if len(np.where(angles == ang)[0]) > 1:
            dupTilts.append(idx)
            
    goodTilts = tilts.copy()
    
    if len(dupTilts) > 0:
        del goodTilts[dupTilts[1]]
    
    return goodTilts

def extractSweeps(radar, tilts):
    
    """
    DESCRIPTION: Returns the lowest tilt of a specified field which has data. 
    New NEXRAD scan strategies have some data for some fields in some titls but not all.
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    field = Names of field in the radar object
    tilt = Tilt to use to filter radar object.
    
    OUTPUTS:
    itilt = Index of the lowest tilt for the field with data.
        
    """
    tiltsExtract = []
    elevAngles = radar.fixed_angle['data']
    
    for tilt in tilts:
        
        ind = nearest_ind(elevAngles, tilt)
        tiltsExtract.append(ind)
        
    radar = radar.extract_sweeps(tiltsExtract)
    
    return radar

def trimRange(radar, limit):#), lvl2_flag=False):
    
    """
    DESCRIPTION: Trims range of radar fields to a specified limit
    
    INPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions.
    limit = Range limit in kilometers
    
    OUTPUTS:
    radar = A python object structure that contains radar information. Created
        by PyART in one of the pyart.io.read functions. Fields will stop at the range limit.
        
    """

    radar_range = radar.range['data']
    range_tile = np.tile(radar_range.transpose(), (radar.nrays, 1))
    
    range_mask = np.ma.masked_greater(range_tile, limit)
    
    radar.range['data'] = np.ma.masked_greater(radar_range, limit)
    
    fields = list(radar.fields.keys())
    
    for ifield in fields:
        
        radar.fields[ifield]['data'] = np.ma.masked_where(np.ma.getmask(range_mask), radar.fields[ifield]['data'])
        
       
    return radar

def trimUR(radar, lvl2_flag=True):
    
    radar_range = radar.range['data']
    range_tile = np.tile(radar_range.transpose(), (radar.nrays, 1))

    if lvl2_flag:        
        vel_field = 'velocity'
    else:         
        vel_field = 'dealised_velocity'
    
    vel_tilt = getDataSweep(radar, vel_field, 0.5)
    
    for itilt in vel_tilt:
        
        vel_slice = radar.get_slice(itilt)
        vel_rmax = np.min(radar.instrument_parameters['unambiguous_range']['data'][vel_slice])
        range_mask = np.ma.masked_greater(range_tile, vel_rmax)
        
        radar.fields[vel_field]['data'][vel_slice] = np.ma.masked_where(np.ma.getmask(range_mask[vel_slice]), radar.fields[vel_field]['data'][vel_slice])

    return radar

def removeSmallObjects(radar, fields):
    
    for field in fields:
        
        for sweep in np.arange(radar.nsweeps):
            
            temp_var = radar.fields[field]['data'][radar.get_slice(sweep)]
            
            var_int = ~temp_var.mask
            var_filter = morphology.remove_small_objects(var_int, 128, connectivity=8)
            new_mask = ~var_filter
            
            radar.fields[field]['data'][radar.get_slice(sweep)].mask = new_mask
            
    return radar


def calcWaves(radarA, radarB):
    

    # set arrays for adding data to
    waves_full = radarB.fields["dealiased_velocity"]['data'].copy()
    diff_full = radarB.fields["dealiased_velocity"]['data'].copy()
    
    for itilt in np.arange(radarB.nsweeps):
        
        print('Calculating waves for sweep ' + str(itilt))
        
#    if index == 1: # for first file
#        
#        wavesA_full = radarA.fields["dealiased_velocity"]['data'].copy()
#        diffA_full = radarA.fields["dealiased_velocity"]['data'].copy()
#    
#    for itilt in both_tilts: # loop through tilts
        
        # get azimuths and range for each file
        aziA = radarA.get_azimuth(itilt); aziB = radarB.get_azimuth(itilt)
        distA = radarA.range['data']; distB = radarB.range['data']
        
        # if files have same number of azimuths
        if np.shape(aziA) == np.shape(aziB):
            
            # if files have same range
            if np.shape(distA) == np.shape(distB):
                
                # get first slice (tilt with data)
                sliceA = radarA.get_slice(itilt)
                sliceB = radarB.get_slice(itilt)
                
                # get dealiased velocity and shift so zero azimuth is first
                dvelA = radarA.fields["dealiased_velocity"]['data'][sliceA]
                dvelA_shift = np.roll(dvelA, -np.where(aziA == np.amin(aziA))[0], axis=0) # negative to roll forward
                #aziA_shift = np.roll(aziA, -np.where(aziA == np.amin(aziA))[0])
                
                dvelB = radarB.fields["dealiased_velocity"]['data'][sliceB]                
                dvelB_shift = np.roll(dvelB, -np.where(aziB == np.amin(aziB))[0], axis=0)
                #aziB_shift = np.roll(aziB, -np.where(aziB == np.amin(aziB))[0])
                
                # subtract tilts
                differenceB = abs(dvelB_shift) - abs(dvelA_shift)
                
                # saving raw difference
                diff_shift = np.roll(differenceB, np.where(aziB == np.amin(aziB))[0], axis=0) # positive to roll backward
                diff_mask = np.ma.masked_where(np.ma.getmask(dvelB), diff_shift)
                
                diff_full[sliceB] = diff_mask
                
                # make waves a binary field
                waves = differenceB
                waves[waves >  -1.0] = 0
                waves[waves <= -1.0] = 1
                                    
                # before we save the waves field we have to rotate it back to original file
                waves_shift = np.roll(waves, np.where(aziB == np.amin(aziB))[0], axis=0) # positive to roll backward
                waves_mask = np.ma.masked_where(np.ma.getmask(dvelB), waves_shift)
                
                waves_full[sliceB] = waves_mask   
                
                # Handling waves in first file
#                if index == 1:
#                    
#                    diffA_shift = np.roll(differenceB, np.where(aziA == np.amin(aziA))[0], axis=0) # positive to roll backward
#                    diffA_mask = np.ma.masked_where(np.ma.getmask(dvelA), diffA_shift)
#                    
#                    diffA_full[sliceA] = diffA_mask
#
#                    wavesA_shift = np.roll(waves, np.where(aziA == np.amin(aziA))[0], axis=0) # positive to roll backward
#                    wavesA_mask = np.ma.masked_where(np.ma.getmask(dvelA), wavesA_shift)
#                    
#                    wavesA_full[sliceA] = wavesA_mask                     
                
            else:
                
                print("Files do not have same range")
                
        else:
            
            print("Files do not have same number of azimuths")
            
    # if files have a mismatch number of tilts, add some masked data to the field
#    if len(velB_tilts) > len(both_tilts):
#        
#       unable = list(set(velB_tilts) - set(both_tilts)) 
#       
#       for jtilt in unable:
#           
#           jslice = radarB.get_slice(jtilt)
#           waves_full[jslice] = np.ma.masked_all(np.shape(radarB.fields["dealiased_velocity"]['data'][jslice]))
            
    # add field to radar object
    radarB.add_field_like('dealiased_velocity','waves', waves_full, True)
    radarB.add_field_like('dealiased_velocity','diff', diff_full, True)
    
    return radarB

# def RHIwaves(radarA, radarB):
#
#     # Find what tilts to use -- need to check that both are the same? In this case, yes.
#     velA_tilts, velA_angs = returnAll(radarA, 'velocity')
#     refA_tilts, refA_angs = returnAll(radarA, 'reflectivity')
#     bothA_tilts = list(set(velA_tilts).intersection(refA_tilts))
#     bothA_angs = list(set(velA_angs).intersection(refA_angs))
#
#     velB_tilts, velB_angs = returnAll(radarB, 'velocity')
#     refB_tilts, refB_angs = returnAll(radarB, 'reflectivity')
#     bothB_tilts = list(set(velB_tilts).intersection(refB_tilts))
#     bothB_angs = list(set(velB_angs).intersection(refB_angs))
#
#     both_tilts = list(set(bothA_tilts).intersection(bothB_tilts))
#     both_angs = list(set(bothA_angs).intersection(bothB_angs))
#
#     # Extract tilts
#     radarA = radarA.extract_sweeps(both_tilts)
#     radarB = radarB.extract_sweeps(both_tilts)
#
#     # Trim range
#     radarA = trimUR(radarA)#, 127000)
#     radarB = trimUR(radarB)#, 127000)
#
#     # filter velocity by reflectivity before dealiasing
#     radarA.fields['velocity']['data'] = np.ma.masked_where(radarA.fields['reflectivity']['data'] < 0, radarA.fields['velocity']['data'])
#     radarB.fields['velocity']['data'] = np.ma.masked_where(radarB.fields['reflectivity']['data'] < 0, radarB.fields['velocity']['data'])
#
#     # Remove small objects
#     radarA = removeSmallObjects(radarA, ['velocity'])
#     radarB = removeSmallObjects(radarB, ['velocity'])
#
#     # Dealias velocity
#     radarA = qc_fun.dealias(radarA, 'velocity', "dealiased_velocity")
#     radarB = qc_fun.dealias(radarB, 'velocity', "dealiased_velocity")
#
#     # Calculate waves
#     radarB = calcWaves(radarA, radarB)
#
#     # Filter waves
#     waves_full = radarB.fields['waves']['data']
#     waves_filtered = np.empty_like(waves_full)
#     starts = radarB.sweep_start_ray_index['data']
#     ends = radarB.sweep_end_ray_index['data']
#
#     for start, end in zip(starts, ends):
#
#         waves_subset = waves_full[start:end]
#         roll_factor = (end-start+1)/2
#         waves_subset_rotated = np.roll(waves_subset, int(roll_factor), axis=0)
#
#         waves_subset_filtered = morphology.remove_small_objects(waves_subset.astype(bool), 32, connectivity=4)
#         waves_subset_rotated_filtered = morphology.remove_small_objects(waves_subset_rotated.astype(bool), 32, connectivity=4)
#
#         waves_subset_rotated_filtered_rotated = np.roll(waves_subset_rotated_filtered, -int(roll_factor), axis=0)
#
#         waves_combined = np.logical_and(waves_subset_filtered, waves_subset_rotated_filtered_rotated)
#
#         waves_filtered[start:end] = waves_combined.astype(np.float32)
#
#     radarB.fields['waves']['data'] = waves_filtered
#
#     return radarB
