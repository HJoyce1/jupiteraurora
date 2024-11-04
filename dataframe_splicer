"""
Created on Fri Oct 25 16:28:34 2024

@author: hannah
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import spiceypy as spice

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'
    
# leap seconds kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")

# constant parameters
RJ = 71492 # jupiter radius in km
sw_end = 54785 # final data point for relevant solar wind section

# load in dataframes
sw_df = pd.read_csv(root_folder+'solar_wind_data.csv')
juno_loc = pd.read_csv(root_folder+'juno_position_df.csv',delimiter=',')
mag_df = pd.read_csv(root_folder+'mag_30sec_df.csv')

mag_df = mag_df.rename(columns={"Juno_Time":"UTC"})

UTC_sw = sw_df['UTC'].to_numpy()

UTC_mag = mag_df['UTC'].to_numpy()

times_sw = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(UTC_sw)):
    timetest = dt.datetime.strptime(UTC_sw[p],'%Y-%jT%H:%M:%S.%f')
    times_sw.append(timetest)
    
times_mag = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(UTC_mag)):
    timetest = dt.datetime.strptime(UTC_mag[p],'%Y-%jT%H:%M:%S.%f')
    times_mag.append(timetest)
    
mag_df = mag_df.assign(Time=times_mag)
sw_df = sw_df.assign(Time=times_sw)

# --------

new_df = pd.merge_asof(sw_df, mag_df, on="Time", tolerance=pd.Timedelta("30s"))

new_df.to_csv(root_folder+'merged_sw_mag.csv',index=False)
