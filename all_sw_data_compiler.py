"""
Created on Fri Oct 25 16:57:27 2024

@author: hannah
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
# this one is for iau - planetary constants
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has cooridinate systems in it
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_jj_jm.tf") # custom coordinate frames in jupiter frame, thanks to Rob Wilson
# juno location
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160522_160729_160909.bsp")
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160312_160522_160614.bsp")

# constant parameters
RJ = 71492 # jupiter radius in km
sw_end = 54785 # final data point for relevant solar wind section

# load in dataframes
sw_df = pd.read_csv(root_folder+'merged_sw_mag_v3.csv')
juno_loc = pd.read_csv(root_folder+'juno_position_df.csv',delimiter=',')
#travel_times = pd.read_csv(root_folder+'travel_times_sw_df.csv')

# grab  arrays from sw dataframe
# pressure
pressure_array = sw_df['RAM_PRESSURE_PROTONS_NPA'].to_numpy() # pressure 
#pressure_uncertainty_array = sw_df['RAM_PRESSURE_PROTONS_NPA_UNCERTAINTY'].to_numpy() # pressure uncertainty - need to use this at some point
# speed
sw_speed_array = sw_df['V_KMPS'].to_numpy()
#sw_uncertainty_array = sw_df['V_KMPS_UNCERTAINTY'].to_numpy() # uncertainty in speed, should use at some point / factor into time diff?
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
    
mp_loc_km = np.array(nose_locs_mp_km)

# juno position arrays
X_juno = juno_loc['X'].to_numpy()
Y_juno = juno_loc['Y'].to_numpy()
X_juno_RJ = juno_loc['XRJ'].to_numpy()
Y_juno_RJ = juno_loc['YRJ'].to_numpy()


# --------- travel times -------

time_utc = sw_df['Time']

juno_time_sw,iono_time,tot_travel_time,time_shift,magnetosheath_t_time,ionosphere_t_time, delta, dx = ptsw.propagation_time(sw_speed_array,X_juno,Y_juno,time_utc[0:sw_end],nose_locs_mp_km[0:sw_end],nose_locs_bs_km[0:sw_end],'compilier')

# get et times for ionotimes
#iono_time_string = np.array_str(iono_time)
iono_time_string = [str(item) for item in iono_time]
iono_time_str = np.array(iono_time_string)

ettimes = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(iono_time_str)):
    times = spice.str2et(iono_time_str[p])
    ettimes.append(times)
    
# ---- coordinate transform ------

Br = sw_df['Br'].to_numpy()
Bt = sw_df['Bt'].to_numpy()
Bn = sw_df['Bn'].to_numpy()

Br = Br[0:sw_end]
Bt = Bt[0:sw_end]
Bn = Bn[0:sw_end]

# Build RTN magnetic field vectors
B_vectors_rtn = np.stack([Br, Bt, Bn], axis=1)  # shape (N, 3)

# Get transformation matrices from RTN to IAU_JUPITER
transform_matrices = np.array([spice.pxform('JUNO_SUN_EQU_RTN', 'JUNO_JM', t) for t in ettimes])  # shape (N, 3, 3)

# Transform each B vector using the rotation matrix at each time
B_vectors = np.einsum('nij,nj->ni', transform_matrices, B_vectors_rtn)  # shape (N, 3)

# Extract components in JM frame
Bx = B_vectors[:, 0] #* -1
By = B_vectors[:, 1] #* -1
Bz = B_vectors[:, 2] 

# clock = []
# for j in range(len(Bn)):
# # use arctan2 electric boogaloo
# # 0 is Bz+, +/-180 is Bz-, +90 is By+, -90 is By-
#     theta_c = np.arctan2(By[j],Bz[j])
#     theta_c_deg = np.degrees(theta_c)
#     clock.append(theta_c_deg)
    
# ----- mag data calculations ------

# some calculations
clock_angle = mdc.clock_angle_calculator(Bz, By,'no','no')
B_perp = mdc.B_perp_calculator(Bz, By, 'no')


# ------- reconnection voltages ------


LL_rec, HL_rec_pos, HL_rec_neg = rv.reconnection_voltages(clock_angle, B_perp, sw_speed_array,nose_locs_mp_km[0:sw_end],juno_time_sw)

# needs Br, Bn, Bt, velocity, pressure, magnetopause location in km, b_perpendicular, clock angle & location (ie low, neg_By, pos_By - dtermine orienation of planet field vs clock angle)
rec_gersh = rv.reconnection_gershman(Bx,By,Bz,sw_speed_array,pressure_array,mp_loc_km[0:sw_end],B_perp,clock_angle)

KH_gersh = rv.kelvin_helmholtz(sw_speed_array,pressure_array,mp_loc_km[0:sw_end])
KH_dawn, KH_dusk = rv.kelvin_helmholtz_dd(sw_speed_array,pressure_array,mp_loc_km[0:sw_end])
    
times = juno_time_sw

# ------- optional plot checker ------
# fig = plt.figure(figsize=(30,30))
# ax = plt.subplot(1,1,1)
# ax.plot(times,(KH/1e9), '.', markersize=0.2) #times[:end_first]
# #ax.set_yscale('log')
# ax.set_xlabel('Time (DOY)')
# ax.set_ylabel('Kelvin Helmholtz Power (GW)')
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%j'))
# ax.xaxis.set_major_locator(plt.MaxNLocator(15))

# saveloc = ('/Users/hannah/OneDrive - Lancaster University/aurora/test_plot.jpg')
# plt.savefig(saveloc,bbox_inches='tight',dpi=400)

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
juno_data_df = juno_data_df.assign(delta=delta)
juno_data_df = juno_data_df.assign(dx=dx)

juno_data_df = juno_data_df.assign(SW_Pressure=pressure_array[0:sw_end]) # pressure
juno_data_df = juno_data_df.assign(SW_Velocity=sw_speed_array[0:sw_end]) # velocity
juno_data_df = juno_data_df.assign(Br=Br[0:sw_end])# br
juno_data_df = juno_data_df.assign(Bt=Bt[0:sw_end]) # bt
juno_data_df = juno_data_df.assign(Bn=Bn[0:sw_end]) # bn
juno_data_df = juno_data_df.assign(Bx=Bx)# br
juno_data_df = juno_data_df.assign(By=By) # bt
juno_data_df = juno_data_df.assign(Bz=Bz) # bn
juno_data_df = juno_data_df.assign(Clock_Angle=clock_angle) # clock angle
juno_data_df = juno_data_df.assign(B_Perp=B_perp) # b perp
juno_data_df = juno_data_df.assign(Bow_Shock_Stand_Off=nose_locs_bs_km[0:sw_end])#[0:sw_end]) # bs stand off
juno_data_df = juno_data_df.assign(Bow_Shock_Stand_Off_RJ=nose_locs_bs[0:sw_end])#[0:sw_end]) 
juno_data_df = juno_data_df.assign(Magnetopause_Stand_Off=nose_locs_mp_km[0:sw_end])#[0:sw_end]) # mp stand off
juno_data_df = juno_data_df.assign(Magnetopause_Stand_Off_RJ=nose_locs_mp[0:sw_end])#[0:sw_end]) 
juno_data_df = juno_data_df.assign(Low_Latitude_Reconnection_Voltage=LL_rec)
juno_data_df = juno_data_df.assign(High_Latitude_Reconnection_Voltage_BY_POS=HL_rec_pos)
juno_data_df = juno_data_df.assign(High_Latitude_Reconnection_Voltage_BY_NEG=HL_rec_neg)
juno_data_df = juno_data_df.assign(Gershman_Reconnection_Power=rec_gersh)
juno_data_df = juno_data_df.assign(Kelvin_Helmholtz_Power=KH_gersh)
juno_data_df = juno_data_df.assign(Kelvin_Helmholtz_Power_Dawn=KH_dawn)
juno_data_df = juno_data_df.assign(Kelvin_Helmholtz_Power_Dusk=KH_dusk)


# # export dataframe to read into other files
juno_data_df.to_csv(root_folder+'juno_data_big_df_realtime_test.csv',index=False)

# # rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
juno_data_df.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
juno_data_df.to_csv(root_folder+'juno_data_big_df_ionotime_test.csv',index=False)
