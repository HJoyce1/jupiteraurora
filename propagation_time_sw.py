#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 17:33:22 2024

@author: hannah
"""

# load in relevant modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates

# sw_end = 54785 # final data point for pre-DOY 155 data (before accumulation time issues)
                    # will probably need to factor this inot module

def propagation_time(velocity,X_juno,Y_juno,juno_time,mp_loc,bs_loc,key): # juno_time in utc isod

    RJ = 71492 # jupiter radius in km
    #breakpoint()
    times = []
    if key == 'compilier':
        # #need to make this into for loop as can only convert one thing into date time at a time
        for p in range(len(juno_time)):
            timetest = dt.datetime.strptime(juno_time[p],'%Y-%m-%d %H:%M:%S.%f') # for complier
            times.append(timetest)
        
    else:
        for p in range(len(juno_time)):
            timetest = dt.datetime.strptime(juno_time[p],'%Y-%jT%H:%M:%S.%f') # for splitter
            times.append(timetest)
        
    times = np.array(times)
    #breakpoint()
    
    # print(velocity[0])
    # print(X_juno[0]/RJ)
    # print(Y_juno[0]/RJ)
    # print(mp_loc[0]/RJ)
    # print(bs_loc[0]/RJ)
        
    ''' time shift to bow shock '''
    sw_angle = np.deg2rad(11)
    
    delta_x = []
    distance_shift = []
    for i in range(len(velocity)):
        delta = np.abs(Y_juno[i])*np.tan(sw_angle)
        delta_x.append(delta)
        dx = bs_loc[i] - X_juno[i] + delta
        distance_shift.append(dx)
    
    time_shift = []
    # this needs to be in seconds - km / kms-1 will = s
    for j in range(len(velocity)):
        shift = distance_shift[j]/velocity[j]
        time_shift.append(shift)
        # need to use distance in km NOT RJ
    #distance_shift1 = np.array(distance_shift)
   # print(min(distance_shift1/RJ))   
    #print('time shift',min(time_shift))
    #print(min(delta_x))
    #print(min(distance_shift))
    
    
    ''' magnetosheath '''
    # need to calculate difference between mp and bs location for distance
    diff_dist = []
    for k in range(len(mp_loc)):
        bs_mp = bs_loc[k] - mp_loc[k]
        diff_dist.append(bs_mp)
        
    # now need a linear decrease for each point from 0.26 of sw velocity to 30kms-1
    # calculate deceleration via suvat and get time there?
    # v**2 = u**2+2as (distance = s, u = initial speed, v = final speed)
    # rearrange: a = v**2 - u**2 / 2s
    #acceleration = []
    time_mp_bs = []
    end_speed = 30
    for l in range(len(velocity)):
        start_speed = velocity[l]*0.26
        speed_diff = start_speed - end_speed #end_speed - start_speed
        #time_to_mp = -(diff_dist[l]/speed_diff) * np.log(end_speed/start_speed)
        time_to_mp = (diff_dist[l]/speed_diff) * np.log(start_speed/end_speed)
        time_mp_bs.append(time_to_mp)
    
    #     speed_diff = end_speed**2 - start_speed**2
    #     acceleration = (end_speed**2 - start_speed**2)/(2*diff_dist[l])

    #     if acceleration >= 0:
    #         print(f"Warning: unexpected acceleration at index {l}, a = {acceleration}, u = {start_speed}, v = {end_speed}")

    #     # time w/ accleration = v - u / a
    #     time_to_mp = (end_speed - start_speed)/acceleration
    #     time_mp_bs.append(time_to_mp)
    #print('min magnetosheath',np.min(time_mp_bs))
        
    
    ''' magnetopause to ionosphere '''
    # now use the magic formula from Cowley and Bunce 2003b
    # 26*(magnetopause_location/20)**2.7
    to_iono_time = []
    for m in range(len(mp_loc)):
        mp_to_RJ = mp_loc[m]/RJ
        time_iono = 26*((mp_to_RJ/20)**2.7)
        to_iono_time.append(time_iono)
    #print('to_iono',min(to_iono_time))
    
        
    ''' overall time '''
    overall_time = []
    time_before_shift = []
    for n in range(len(time_shift)):
        # time calc is bs to mp time + iono time - time shift
        time_sum = time_mp_bs[n] + to_iono_time[n] - time_shift[n]
        # for checking time w/out tims shift
        time_pos = time_mp_bs[n] + to_iono_time[n]
        overall_time.append(time_sum)
        time_before_shift.append(time_pos)
    
    overall_time_array = np.array(overall_time)
    
    fig = plt.figure(figsize=(15,30))
    ax1 = plt.subplot(8,1,1)
    ax1.plot(times, overall_time_array/3600,'.')
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
    ax1.yaxis.set_major_locator(plt.MaxNLocator(20))
    ax1.xaxis.set_major_locator(plt.MaxNLocator(15))
    ax1.set_xlabel('Time (DOY)')
    ax1.set_ylabel('Time Taken to Reach Ionosphere (s)')
    
    ''' Ionospheric Time '''
    iono_time = []
    for p in range(len(times)):
        time_calc = times[p] + dt.timedelta(seconds=overall_time[p])
        iono_time.append(time_calc)
        
    return times, iono_time,overall_time,time_shift,time_mp_bs,to_iono_time, delta_x, distance_shift

'''
# ------- Dataframe  ------

# making a dataframe - want to make a big data frame with everything in it
travel_times = pd.DataFrame()


# assign data to dataframe
travel_times = travel_times.assign(Juno_SW_Detection_Time=times[0:sw_end])
travel_times = travel_times.assign(Time_Impacts_Ionosphere=iono_time)
travel_times = travel_times.assign(Total_Travel_Time=overall_time)
travel_times = travel_times.assign(Time_Shift_To_Bow_Shock=time_shift)
travel_times = travel_times.assign(Magnetosheath_Travel_Time=time_mp_bs)
travel_times = travel_times.assign(Magnetopause_To_Ionosphere_Time=to_iono_time)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
# travel_times.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# export dataframe to read into other files
travel_times.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/travel_times_sw_df.csv',index=False)
