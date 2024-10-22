#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 11:13:17 2024

@author: hannah

this script calcualtes the elevation angle between surface and HST for each HST visit 

THIS SCRIPT IS NOT AUTOMATED and requires editing between visits
"""

import numpy as np
import spiceypy as spice
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
#list(colormaps)

# leap seconds kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")
# this one is for iau - planetary constants
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has coordinate systems in it
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")
# kernal for planet locations
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/de436s.bsp")
# hst kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/hst.bsp")

# mac
root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/'

radius_e = 71942 + 240
radius_p = 66854 + 240

# dataframes
visit_times = pd.read_csv(root_folder+'dataframes/visit_times.csv')

# grab all visit time data
visit_04 = visit_times['Visit_04'].to_numpy()
visit_05 = visit_times['Visit_05'].to_numpy()
visit_08 = visit_times['Visit_08'].to_numpy()
visit_09 = visit_times['Visit_09'].to_numpy()
visit_10 = visit_times['Visit_10'].to_numpy()
visit_11 = visit_times['Visit_11'].to_numpy()
visit_12 = visit_times['Visit_12'].to_numpy()
visit_13 = visit_times['Visit_13'].to_numpy()
visit_15 = visit_times['Visit_15'].to_numpy()
visit_16 = visit_times['Visit_16'].to_numpy()
visit_17 = visit_times['Visit_17'].to_numpy()
visit_18 = visit_times['Visit_18'].to_numpy()
visit_19 = visit_times['Visit_19'].to_numpy()
visit_20 = visit_times['Visit_20'].to_numpy()
visit_21 = visit_times['Visit_21'].to_numpy()


def elevAngles(radius_e, radius_p, visit):
    import numpy as np
    # need to furnsh some kernal? meta kernal probably?

    '''
    VISIT
    '''
    ettime = []
    # convert visit into ephermous time to be able to compare with iono time
    for p in range(len(visit)):
        times = spice.str2et(visit[p])
        ettime.append(times)
        
    
    radius_e_p = np.array([radius_e, radius_e, radius_p])

    
    # longitudes
    import numpy as np
    lons = np.linspace(0, 360, num=1441)
    lons = lons[:-1]
    lons = lons/180*np.pi
    
    # latitudes
    lats = np.linspace(-90, 90, num=721)
    lats = lats[:-1]
    lats = lats/180*np.pi


    # store angles
    elev_angles = np.full((len(lons),len(lats),len(ettime)), np.nan)
    #breakpoint()
    
    # vector for HST -> Jupiter
    vec = np.array([spice.spkez(spice.bodn2c('HST'), et, 'IAU_JUPITER','NONE',spice.bodn2c('JUPITER BARYCENTER'))[0][:3] for et in ettime]) #[0][:3]
    #positions, lighttimes = spice.spkpos('HST',ettime,'IAU_JUPITER','NONE','JUPITER BARYCENTER')
    
    # calculate viewing angle from HST for each long/lat
    for lon_angle in range(len(lons)):
        if not lon_angle%10:
            print('Progress:',int(lon_angle/len(lons)*100),'%')
        for lat_angle in range(len(lats)):
            lon = lons[lon_angle]
            lat = lats[lat_angle]
            sfnm = spice.surfnm(radius_e_p[0],radius_e_p[1],radius_e_p[2], spice.latrec(radius_e,lon,lat))
            for et_angle in range(len(ettime)):
                elev_angles[lon_angle, lat_angle, et_angle] = spice.vsep(vec[et_angle,:], sfnm)
    
    return (np.pi/2 - elev_angles),ettime,vec #-(elev_angles - np.pi/2)




visit = visit_05
angles, times, vecs = elevAngles(radius_e, radius_p, visit)
angles_deg = angles*180/np.pi
np.savez(root_folder+'elev_angles/elevation_angles_visit_05_test.npz',angles)
