#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 17:05:52 2023

@author: hannah

function to interpolate the 60-second magnetometer data down to 30 second averaged

works by expanding array of epermarous times and filling out
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import spiceypy as spice
import scipy.interpolate as interp
import glob
import datetime as dt
import matplotlib.dates as mdates
import math

# leap seconds kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")
# this one is for iau - planetary constants
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has cooridinate systems in it
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")

# load in big 1 minute averaged dataframe
mag_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/mag_df.csv')

mag_df_min = mag_df['minute'].to_numpy()

# acc_time = []
# for i in range(len(mag_df_min)):
#     if 
#     diff = mag_df_min[i+1]-mag_df_min[i]

df = pd.DataFrame(mag_df_min)
df = df.assign(Min=mag_df_min)
df['Accumulation_Time']= df['Min'].diff()
    

# convert to empherous time in order to interpolate
utc = mag_df['UTC'].to_numpy()
et_time = spice.str2et(utc)

# create extra spaces for new data
array_len = (len(et_time)) * 2
# now expand time array for interpolation
new_time_array = np.linspace(min(et_time), max(et_time), num=array_len)

# first things first - check accumulation time
fdoy = mag_df['fdoy'].to_numpy()

# function to check accumulation time
def accumulation_time(fdoy):
    accumulation_time = []
    for i in range(len(fdoy)-1):
        diff = fdoy[i+1] - fdoy[i]
        diff = diff * 24 *60
        if diff > 1:
            #print(diff)
            frac, intNum = math.modf(diff)
            diff = 60 + (frac * 60)
            #print(diff)
        else:
            diff = diff*60
        accumulation_time.append(diff)
        
    times = []
    # need to make this into for loop as can only convert one thing into date time at a time
    for p in range(len(mag_df['UTC'])):
        timetest = dt.datetime.strptime(mag_df['UTC'][p],'%Y-%jT%H:%M:%S.%f')
        times.append(timetest)
        
    # fig = plt.figure(figsize=(12,4))
    # ax = plt.subplot(1,1,1)  
    # ax.plot(times[:-1], accumulation_time,'.',markersize=0.1)
    # ax.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    # ax.xaxis.set_major_locator(plt.MaxNLocator(15))
    # ax.set_xlabel('DOY')
    # ax.set_ylabel('Accumulation Time (s)')
    # ax.set_ylim(0,120)
    
# # call funtction
# accumulation_time(fdoy)

# for interpolation need Br, Bt, Bn
Br = mag_df['Br'].to_numpy()
Bt = mag_df['Bt'].to_numpy()
Bn = mag_df['Bn'].to_numpy()

times = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(mag_df['UTC'])):
    timetest = dt.datetime.strptime(mag_df['UTC'][p],'%Y-%jT%H:%M:%S.%f')
    times.append(timetest)

new_fdoy = np.linspace(min(fdoy),max(fdoy),num=array_len)
interpolated_Br = np.interp(new_time_array, et_time, Br)
interpolated_Bt = np.interp(new_time_array, et_time, Bt)
interpolated_Bn = np.interp(new_time_array, et_time, Bn)
#interpolated_fdoy = np.interp(new_time_array, et_time, fdoy)

norm_time = spice.et2utc(new_time_array, "ISOD", 5)

# # quick check to see interpolated and expanded fdoy are equal
# fig = plt.figure(figsize=(12,4))
# ax = plt.subplot(1,1,1)  
# ax.plot(new_time_array, new_fdoy, '.',markersize=10.0)
# ax.plot(new_time_array,interpolated_fdoy,'.',markersize=1)
# #ax.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
# #ax.xaxis.set_major_locator(plt.MaxNLocator(15))
# ax.set_xlabel('DOY')
# ax.set_ylabel('fdoy (s)')

# now need to put these all into a dataframe....
mag_30sec_df = pd.DataFrame()
mag_30sec_df = mag_30sec_df.assign(Juno_Time=norm_time)
mag_30sec_df = mag_30sec_df.assign(Fractional_DOY=new_fdoy)
mag_30sec_df = mag_30sec_df.assign(ET_Time=new_time_array)
mag_30sec_df = mag_30sec_df.assign(Br=interpolated_Br)
mag_30sec_df = mag_30sec_df.assign(Bt=interpolated_Bt)
mag_30sec_df = mag_30sec_df.assign(Bn=interpolated_Bn)

# can save out file
mag_30sec_df.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/mag_30sec_df.csv',index=False)
