"""
Created on Wed Dec 13 16:23:34 2023

@author: hannah

this script collects all the different solar wind data and combines it into one large dataframe 
for convinenece - there are two different versions saved out, organised in terms of juno
detection time or time impacts ionosphere
"""

# load in relevant modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import spiceypy as spice
import mag_data_calculations as mdc 
import joy_model_python as jmp
import propagation_time_sw as ptsw
import reconnection_voltage as rv


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


# grab arrays needed from mag dataframe
Br = mag_df['Br'].to_numpy()
Bn = mag_df['Bn'].to_numpy()
Bt = mag_df['Bt'].to_numpy()
juno_time_mag = mag_df['Fractional_DOY'].to_numpy() 

# trim mag data to only use from 2 mins 30 secs in to match sw data
Br = Br[5:sw_end+5]
Bn = Bn[5:sw_end+5]
Bt = Bt[5:sw_end+5]
juno_time_mag = juno_time_mag[5:]#sw_end+5]

# grab  arrays from sw dataframe
# pressure
pressure_array = sw_df['RAM_PRESSURE_PROTONS_NPA'].to_numpy() # pressure 
pressure_uncertainty_array = sw_df['RAM_PRESSURE_PROTONS_NPA_UNCERTAINTY'].to_numpy() # pressure uncertainty - need to use this at some point
# speed
sw_speed_array = sw_df['V_KMPS'].to_numpy()
sw_uncertainty_array = sw_df['V_KMPS_UNCERTAINTY'].to_numpy() # uncertainty in speed, should use at some point / factor into time diff?
# probably don't want this ^ but want another uncertainty array based on travel time

# trim sw data to decrease run time
pressure_array = pressure_array[0:sw_end]
sw_speed_array = sw_speed_array[0:sw_end]

nose_locs_mp, nose_locs_bs = jmp.multi_nose(pressure_array)

nose_locs_mp_km = []
nose_locs_bs_km = []
# loop through magnetopause and bow shock locations to convert to km
for i in range(len(nose_locs_mp)):
    # magnetopause
    mp =  nose_locs_mp[i]*RJ
    nose_locs_mp_km.append(mp)
    # bowshock
    bs =  nose_locs_bs[i]*RJ
    nose_locs_bs_km.append(bs)

# juno position arrays
X_juno = juno_loc['X'].to_numpy()
Y_juno = juno_loc['Y'].to_numpy()
X_juno_RJ = juno_loc['XRJ'].to_numpy()
Y_juno_RJ = juno_loc['YRJ'].to_numpy()

fig = plt.figure(figsize=(20,5))
ax1 = plt.subplot(1,1,1)
ax1.scatter(np.linspace(0,sw_end,len(nose_locs_mp_km)),nose_locs_mp_km,s=0.01)
ax1.scatter(np.linspace(0,sw_end,len(nose_locs_mp_km)),nose_locs_bs_km,s=0.01)
#ax1.scatter(nose_locs_mp_km, nose_locs_bs_km, s=0.001)

# save plot
saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/bs_vs_mp_full.jpg')
plt.savefig(saveloc,bbox_inches='tight',dpi=400)


# --------- travel times -------

time_utc = sw_df['UTC']

juno_time_sw,iono_time,tot_travel_time,time_shift,magnetosheath_t_time,ionosphere_t_time = ptsw.propagation_time(pressure_array,sw_speed_array,X_juno,Y_juno,time_utc[0:sw_end],nose_locs_mp_km[0:sw_end],nose_locs_bs_km[0:sw_end])

# get et times for ionotimes
#iono_time_string = np.array_str(iono_time)
iono_time_string = [str(item) for item in iono_time]
iono_time_str = np.array(iono_time_string)

ettimes = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(iono_time_str)):
    times = spice.str2et(iono_time_str[p])
    ettimes.append(times)
    
# ----- mag data calculations ------

# some calculations
clock_angle, clock_angle_range = mdc.clock_angle_calculator(Bn, Bt,'no','no')
B_perp = mdc.B_perp_calculator(Bn, Bt, 'no')


# ------- reconnection voltages ------


LL, HL_pos, HL_neg = rv.reconnection_voltages(clock_angle, B_perp, sw_speed_array, nose_locs_mp_km[0:sw_end],juno_time_sw)


# ------- save out to dataframe ------

# put in new dataframe
juno_data_df = pd.DataFrame()

# add data to dataframe 
juno_data_df = juno_data_df.assign(Juno_Detection_Time=juno_time_sw)
juno_data_df = juno_data_df.assign(Juno_Location_X=X_juno) # juno x
juno_data_df = juno_data_df.assign(Juno_Location_X_RJ=X_juno_RJ)
juno_data_df = juno_data_df.assign(Juno_Location_Y=Y_juno) # juno x
juno_data_df = juno_data_df.assign(Juno_Location_Y_RJ=Y_juno_RJ)

juno_data_df = juno_data_df.assign(Time_Impacts_Ionosphere=iono_time)
juno_data_df = juno_data_df.assign(ET_Iono_Time=ettimes)
juno_data_df = juno_data_df.assign(Total_Travel_Time=tot_travel_time)
juno_data_df = juno_data_df.assign(Time_Shift_To_Bow_Shock=time_shift)
juno_data_df = juno_data_df.assign(Magnetosheath_Travel_Time=magnetosheath_t_time)
juno_data_df = juno_data_df.assign(Magnetopause_To_Ionosphere_Time=ionosphere_t_time)

juno_data_df = juno_data_df.assign(SW_Pressure=pressure_array[0:sw_end]) # pressure
juno_data_df = juno_data_df.assign(SW_Velocity=sw_speed_array[0:sw_end]) # velocity
juno_data_df = juno_data_df.assign(Br=Br)# br
juno_data_df = juno_data_df.assign(Bt=Bt) # bt
juno_data_df = juno_data_df.assign(Bn=Bn) # bn
juno_data_df = juno_data_df.assign(Clock_Angle=clock_angle) # clock angle
juno_data_df = juno_data_df.assign(Clock_Angle_Error=clock_angle_range)
juno_data_df = juno_data_df.assign(B_Perp=B_perp) # b perp
juno_data_df = juno_data_df.assign(Bow_Shock_Stand_Off=nose_locs_bs_km[0:sw_end])#[0:sw_end]) # bs stand off
juno_data_df = juno_data_df.assign(Bow_Shock_Stand_Off_RJ=nose_locs_bs[0:sw_end])#[0:sw_end]) 
juno_data_df = juno_data_df.assign(Magnetopause_Stand_Off=nose_locs_mp_km[0:sw_end])#[0:sw_end]) # mp stand off
juno_data_df = juno_data_df.assign(Magnetopause_Stand_Off_RJ=nose_locs_mp[0:sw_end])#[0:sw_end]) 
juno_data_df = juno_data_df.assign(Low_Latitude_Reconnection_Voltage=LL)
juno_data_df = juno_data_df.assign(High_Latitude_Reconnection_Voltage_BY_POS=HL_pos)
juno_data_df = juno_data_df.assign(High_Latitude_Reconnection_Voltage_BY_NEG=HL_neg)


# export dataframe to read into other files
juno_data_df.to_csv(root_folder+'juno_data_big_df_realtime_updated.csv',index=False)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
juno_data_df.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# export SORTED dataframe to read into other files
juno_data_df.to_csv(root_folder+'juno_data_big_df_ionotime_updated.csv',index=False)
