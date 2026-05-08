#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  5 15:31:29 2025

@author: hannah

visualisiation test of coordinate transform of magnetic field (IMF) data into Jupiter frame
this is done within the all_sw_data_compiler but can be visualised here

"""
    
import spiceypy as spice
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates
import pandas as pd

plotting = 'yes'

# leap seconds kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")
# this one is for iau - planetary constants
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has cooridinate systems in it
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/temp.tf")
# juno location
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160522_160729_160909.bsp")
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160312_160522_160614.bsp")

# load in big 1 minute averaged dataframe
mag_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/mag_df.csv')

new_mag_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/mag_30sec_df.csv')

mag_df_min = mag_df['minute'].to_numpy()

df = pd.DataFrame(mag_df_min)
df = df.assign(Min=mag_df_min)

utc = mag_df['UTC'].to_numpy()
et_time = spice.str2et(utc)

utc_datetimes = [datetime.strptime(spice.et2utc(et, 'ISOC', 0), '%Y-%m-%dT%H:%M:%S') for et in et_time]

# Extract Day of Year from each datetime
doy = [dt.timetuple().tm_yday for dt in utc_datetimes]

Br = mag_df['Br'].to_numpy()
Bt = mag_df['Bt'].to_numpy()
Bn = mag_df['Bn'].to_numpy()

# Build RTN magnetic field vectors
B_vectors_rtn = np.stack([Br, Bt, Bn], axis=1)  # shape (N, 3)

# Get transformation matrices from RTN to IAU_JUPITER
transform_matrices = np.array([spice.pxform('JUNO_SUN_EQU_RTN', 'JUNO_JM', t) for t in et_time])  # shape (N, 3, 3)

# Transform each B vector using the rotation matrix at each time
B_vectors = np.einsum('nij,nj->ni', transform_matrices, B_vectors_rtn)  # shape (N, 3)

# Extract components in JM frame
Bx = B_vectors[:, 0] #* -1
By = B_vectors[:, 1] #* -1
Bz = B_vectors[:, 2] 

clock = []
for j in range(len(Bn)):
# use arctan2 electric boogaloo
# 0 is Bz+, +/-180 is Bz-, +90 is By+, -90 is By-
    theta_c = np.arctan2(By[j],Bz[j])
    theta_c_deg = np.degrees(theta_c)
    clock.append(theta_c_deg)

if plotting == 'yes':    
    # Plotting
    plt.figure(figsize=(20, 4))
    plt.plot(utc_datetimes, Bx)
    plt.plot(utc_datetimes, By)
    plt.plot(utc_datetimes, Bz)
    # Format the x-axis
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gcf().autofmt_xdate()
    plt.ylabel('Magnetic Field Strength (nT)')
    plt.xlabel('UTC Date')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    plt.figure(figsize=(50, 4))
    plt.scatter(utc_datetimes, clock, s=1)
    import matplotlib.ticker as ticker
    
    # Custom formatter for DOY labels
    def doy_formatter(x, pos):
        try:
            date = mdates.num2date(x)  # convert from matplotlib float to datetime
            return f"{date.timetuple().tm_yday}"
        except:
            return ''
    
    # Set major ticks (e.g., one per day)
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(doy_formatter))
    plt.gcf().autofmt_xdate()
    plt.ylabel('Clock Angle (Degrees)')
    plt.xlabel('DOY')
    plt.tight_layout()
    
    saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/new_clock_2025_jm.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    plt.show()


if plotting == 'yes':
    Bz_pos = []
    Bz_neg = []
    By_pos = []
    By_neg = []
    
    for angle in range(len(clock)):
        if clock[angle] > -45 and clock[angle] < 45:
            Bz_pos.append(clock[angle])
        
        elif clock[angle] > 45 and clock[angle] < 135:
            By_pos.append(clock[angle])
    
        elif clock[angle] < -45 and clock[angle] > -135:
            By_neg.append(clock[angle])
            
        else:
            Bz_neg.append(clock[angle])
            
    
            
    clock_dist = np.array([len(Bz_pos),len(By_pos),len(Bz_neg),len(By_neg)])
    clock_labels = ['+ Bz','+ By  ','- Bz  ',' - By']
    
    plt.pie(clock_dist, labels = clock_labels)
    plt.show() 

# create extra spaces for new data
array_len = (len(et_time)) * 2
fdoy = mag_df['fdoy'].to_numpy()
# now expand time array for interpolation
new_time_array = np.linspace(min(et_time), max(et_time), num=array_len)
new_fdoy = np.linspace(min(fdoy),max(fdoy),num=array_len)
interpolated_Bx = np.interp(new_time_array, et_time, Bx)
interpolated_By = np.interp(new_time_array, et_time, By)
interpolated_Bz = np.interp(new_time_array, et_time, Bz)

norm_time = spice.et2utc(new_time_array, "ISOD", 5)

norm_times = [datetime.strptime(spice.et2utc(et, 'ISOC', 0), '%Y-%m-%dT%H:%M:%S') for et in new_time_array]

clock2 = []
for j in range(len(interpolated_Bz)):
# use arctan2 electric boogaloo
# 0 is Bz+, +/-180 is Bz-, +90 is By+, -90 is By-
    theta_c2 = np.arctan2(interpolated_By[j],interpolated_Bz[j])
    theta_c_deg2 = np.degrees(theta_c2)
    clock2.append(theta_c_deg2)


if plotting == 'yes':
    plt.figure(figsize=(50, 4))
    plt.scatter(norm_times, clock2, s=1)
    import matplotlib.ticker as ticker
    
    # Custom formatter for DOY labels
    def doy_formatter(x, pos):
        try:
            date = mdates.num2date(x)  # convert from matplotlib float to datetime
            return f"{date.timetuple().tm_yday}"
        except:
            return ''
    
    # Set major ticks (e.g., one per day)
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    plt.gca().xaxis.set_major_formatter(ticker.FuncFormatter(doy_formatter))
    # Format the x-axis
    #plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=1))
    #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gcf().autofmt_xdate()
    plt.ylabel('Clock Angle (Degrees)')
    plt.xlabel('DOY')
    plt.tight_layout()
    
    saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/new_clock_2025_jm.jpg')
    plt.savefig(saveloc,bbox_inches='tight',dpi=400)
    plt.show()
