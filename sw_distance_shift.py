#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 15:20:33 2023

@author: hannah
"""

# need to turn this into a module or copy it over?

import pandas as pd
import numpy as np
import joy_model_python as jmp # import joy_model.py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime as dt

sw_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/solar_wind_data.csv',delimiter=',')
juno_loc = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/juno_position_df.csv',delimiter=',')

ram_pressure_array = sw_df['RAM_PRESSURE_PROTONS_NPA'].to_numpy() # pressure 
pressure_uncertainty_array = sw_df['RAM_PRESSURE_PROTONS_NPA_UNCERTAINTY'].to_numpy() # pressure uncertainty

sw_speed = sw_df[['V_KMPS']].to_numpy()

X_juno = juno_loc['X'].to_numpy()
Y_juno = juno_loc['Y'].to_numpy()
#time_juno = juno_loc['UTC'][:-9].to_numpy()


times = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(sw_df['UTC'])):
    timetest = dt.datetime.strptime(sw_df['UTC'][p],'%Y-%jT%H:%M:%S.%f')
    times.append(timetest)

pressure = ram_pressure_array
nose_locs_mp, nose_locs_bs = jmp.multi_nose(pressure)

nose_locs_bs_km = []
for j in range(len(nose_locs_bs)):
    nose_loc = nose_locs_bs[j] * 71492
    nose_locs_bs_km.append(nose_loc)
    
X_juno_RJ = []
Y_juno_RJ = []
pos_juno_RJ = []
for j in range(len(X_juno)):
    toX_RJ = X_juno[j]/71492
    X_juno_RJ.append(toX_RJ)
    toY_RJ = Y_juno[j]/71492
    pos_Y = abs(toY_RJ)
    Y_juno_RJ.append(toY_RJ)
    pos_juno_RJ.append(pos_Y)

sw_angle = np.deg2rad(11)

#delta_x = []
distance_x = []
for i in range(len(pressure[0:54785])):
    delta = np.abs(Y_juno[i])*np.tan(sw_angle)
    #delta_x.append(delta)
    dx = nose_locs_bs_km[i] - X_juno[i] + delta
    distance_x.append(dx)
    
# deltas = []
# distance_x = []
# for i in range(len(pressure[0:54785])):
#     delta = np.abs(Y_juno_RJ[i])*np.tan(sw_angle)
#     deltas.append(delta)
#     dx = nose_locs_bs[i] - X_juno_RJ[i] + delta
#     distance_x.append(dx)
    
time_x = []
# this needs to be in seconds - km / kms-1 will = s
for j in range(len(pressure[0:54785])):
    shift = distance_x[j]/sw_speed[j]
    time_x.append(shift)
    # need to use distance in km NOT RJ
    
fig = plt.figure(figsize=(10,6))
ax = plt.subplot(1,1,1)
ax.plot(times[0:54785],distance_x)
ax.plot(times[0:54785],X_juno)
ax.plot(times[0:54785],Y_juno)
ax.plot(times[0:54785],nose_locs_bs_km[0:54785])
#ax.plot(times[0:54785],deltas)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
ax.xaxis.set_major_locator(plt.MaxNLocator(15))
ax.set_xlabel('DOY')
ax.set_ylabel('Distance Shift (RJ)')
# ax.set_ylabel('Distance Shift (km) x10$^7$')

