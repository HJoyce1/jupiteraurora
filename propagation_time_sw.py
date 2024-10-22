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

def propagation_time(pressure,velocity,X_juno,Y_juno,juno_time,mp_loc,bs_loc): # juno_time in utc isod

    RJ = 71492 # jupiter radius in km
    
    times = []
    # need to make this into for loop as can only convert one thing into date time at a time
    for p in range(len(juno_time)):
        timetest = dt.datetime.strptime(juno_time[p],'%Y-%jT%H:%M:%S.%f')
        times.append(timetest)
        
    times_array = np.array(times)
        
    ''' time shift to bow shock '''
    sw_angle = np.deg2rad(11)
    
    #delta_x = []
    distance_shift = []
    for i in range(len(pressure)):
        delta = np.abs(Y_juno[i])*np.tan(sw_angle)
        #delta_x.append(delta)
        dx = bs_loc[i] - X_juno[i] + delta
        distance_shift.append(dx)
        
    time_shift = []
    # this needs to be in seconds - km / kms-1 will = s
    for j in range(len(pressure)):
        shift = distance_shift[j]/velocity[j]
        time_shift.append(shift)
        # need to use distance in km NOT RJ
    distance_shift1 = np.array(distance_shift)
    print(min(distance_shift1/RJ))   
    print('time shift',min(time_shift))
    
    
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
    acceleration = []
    time_mp_bs = []
    end_speed = 30
    for l in range(len(velocity)):
        start_speed = velocity[l]*0.26
        speed_diff = end_speed**2 - start_speed**2
        speed_diff_over_diff_dist = speed_diff/(2*diff_dist[l])
        acceleration.append(speed_diff_over_diff_dist)
    
        # time w/ accleration = v - u / a
        time_to_mp = end_speed - start_speed/speed_diff_over_diff_dist
        time_mp_bs.append(time_to_mp)
    print('min magnetosheath',min(time_mp_bs))
        
    
    ''' magnetopause to ionosphere '''
    # now use the magic formula from Cowley and Bunce 2003b
    # 26*(magnetopause_location/20)**2.7
    to_iono_time = []
    for m in range(len(mp_loc)):
        mp_to_RJ = mp_loc[m]/RJ
        time_iono = 26*((mp_to_RJ/20)**2.7)
        to_iono_time.append(time_iono)
    print('to_iono',min(to_iono_time))
    
        
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
        
    # negative_overall = (np.where(time_before_shift < min(time_shift)))
    # negative_overall = list(negative_overall)
    # negative_overall = np.array(negative_overall).T
    
    # dates_negative = []
    # speeds_negative = []
    # bss_negative = []
    # juno_loc_Y = []
    # juno_loc_X = []
    # for v in negative_overall:
    #     #print(v)
    #     dates = times_array[v]
    #     dates_negative.append(dates)
        
    #     #breakpoint()
    #     speeds = np.array(velocity)
    #     speeds = velocity[v]
    #     speeds_negative.append(speeds)
        
    #     bss = np.array(bs_loc)
    #     bss = bs_loc[v]
    #     bss_negative.append(bss)
        
    #     x_juno = X_juno[v]/RJ
    #     juno_loc_X.append(x_juno)      
        
    #     y_juno = Y_juno[v]/RJ
    #     juno_loc_Y.append(y_juno)
        
    # dates_negative = np.array(dates_negative)#.T
    # bss_negative = np.array(bss_negative)#.T
    # speeds_negative = np.array(speeds_negative)#.T
    # juno_negative_X = np.array(juno_loc_X)#.T
    # juno_negative_Y = np.array(juno_loc_Y)#.T
    
    # print('min bs for negative time', min(bss_negative))
    # print('max bs for negative time', max(bss_negative))
    # print('min speed for negative time',min(speeds_negative))
    # print('max speed for negative time',max(speeds_negative))
    
    # print('max X for negative time',max(juno_negative_X))
    # print('max Y for negative time',max(juno_negative_Y))
    # print('min X for negative time',min(juno_negative_X))
    # print('min Y for negative time',min(juno_negative_Y))
        
        
    # #print(max(overall_time))
    # #print(min(overall_time))
    
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
        
    return times, iono_time,overall_time,time_shift,time_mp_bs,to_iono_time
