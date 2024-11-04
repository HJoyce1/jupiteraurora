"""
Created on Fri Oct 18 14:15:48 2024

@author: hannah
"""


# load in relevant modules
import pandas as pd
import spiceypy as spice
import numpy as np
import datetime as dt
from datetime import datetime
import mag_data_calculations as mdc
import joy_model_python as jmp
import propagation_time_sw as ptsw
import reconnection_voltage as rv

# reconnection voltages will not work until mag data and sw plasma data is equal in length

# mac
root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'
# windows
# root_folder = ''

# leap seconds kernal - need this for ephemerous time
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")


RJ = 71492

sw_data = pd.read_csv(root_folder+'merged_sw_mag.csv')
# mag_df = pd.read_csv(root_folder+'mag_30sec_df.csv')
juno_loc_155 = pd.read_csv(root_folder+'juno_position_155_df.csv')
juno_loc_156_plus = pd.read_csv(root_folder+'juno_position_156_plus_df.csv')
juno_loc_175 = pd.read_csv(root_folder+'juno_position_175_df.csv')

X_juno_155 = juno_loc_155['X'].to_numpy()
Y_juno_155 = juno_loc_155['Y'].to_numpy()
X_juno_RJ_155 = juno_loc_155['XRJ'].to_numpy()
Y_juno_RJ_155 = juno_loc_155['YRJ'].to_numpy()

X_juno_156_plus = juno_loc_156_plus['X'].to_numpy()
Y_juno_156_plus = juno_loc_156_plus['Y'].to_numpy()
X_juno_RJ_156_plus = juno_loc_156_plus['XRJ'].to_numpy()
Y_juno_RJ_156_plus = juno_loc_156_plus['YRJ'].to_numpy()

X_juno_175 = juno_loc_175['X'].to_numpy()
Y_juno_175 = juno_loc_175['Y'].to_numpy()
X_juno_RJ_175 = juno_loc_175['XRJ'].to_numpy()
Y_juno_RJ_175 = juno_loc_175['YRJ'].to_numpy()


times = sw_data['Time'].to_numpy()
#mag_times = _df['ET_Time'].to_numpy()

pressure = sw_data['RAM_PRESSURE_PROTONS_NPA'].to_numpy()
velocity = sw_data['V_KMPS'].to_numpy()

Br = sw_data['Br'].to_numpy()
Bn = sw_data['Bn'].to_numpy()
Bt = sw_data['Bt'].to_numpy()


junotime = []
for p in range(len(times)):
    timetest = dt.datetime.strptime(times[p],'%Y-%m-%d %H:%M:%S.%f')
    junotime.append(timetest)

junotime = np.array(junotime)

# not sure why tf I have to do this but get this error if convert to et directly from junotime:
#    AttributeError: 'datetime.datetime' object has no attribute 'encode'
juno_time = pd.DataFrame()
juno_time = juno_time.assign(Juno_Time_SW=junotime)
juno_time.to_csv(root_folder+'juno_times_sw_df.csv',index=False)

all_juno_times = pd.read_csv(root_folder+'juno_times_sw_df.csv')
juno_times = all_juno_times['Juno_Time_SW'].to_numpy()

et_time = spice.str2et(juno_times)



def between_idxs(dataset,first,last):
    # empty arrays to fill
    idx = []
    between_times = []
    # loop through dataset
    for count, times in enumerate(dataset):
        if(times >= first and times <= last):
            # record index number
            idx.append(count)
            # record value
            between_times.append(times)
    return between_times, idx


def index_compare(index, variable):
    index_values = []
    for i  in index:
        index_values.append(variable[i])
    return index_values


def time_diff(d1,d2):
    time_diff = []
    for i in range(len(d1)):
        d1 = datetime.strptime(d1, '%Y-%jT%H:%M:%S.%f')
        d2 = datetime.strptime(d2, '%Y-%jT%H:%M:%S.%f')
        diff = abs((d2 - d1))
        time_diff.append(diff)
    return time_diff



# ---------- visit 25 --------
doy155_start = "2016-06-03 01:00:00"
doy155_end = "2016-06-04 00:00:00"

first_155 = spice.str2et(doy155_start)
last_155 = spice.str2et(doy155_end)

time_155, idx_155 = between_idxs(et_time,first_155,last_155) # indexes for sw df only


time_155_utc = spice.et2utc(time_155,'ISOD',6)

# grab data
pressure_155 = index_compare(idx_155, pressure)
velocity_155 = index_compare(idx_155, velocity)

Br_155 = index_compare(idx_155,Br)
Bt_155 = index_compare(idx_155,Bt)
Bn_155 = index_compare(idx_155,Bn)


mp_loc_155, bs_loc_155 = jmp.multi_nose(pressure_155)
clock_angle_155, clock_angle_range = mdc.clock_angle_calculator(Bn_155, Bt_155,'no','no')
B_perp_155 = mdc.B_perp_calculator(Bn_155, Bt_155, 'no')

mp_loc_155_km = []
bs_loc_155_km = []
# loop through magnetopause and bow shock locations to convert to km
for i in range(len(mp_loc_155)):
    # magnetopause
    mp = mp_loc_155[i]*RJ
    mp_loc_155_km.append(mp)
    # bowshock
    bs = bs_loc_155[i]*RJ
    bs_loc_155_km.append(bs)

times_155,iono_time_155,overall_time_155,time_shift_155,time_mp_bs_155,to_iono_time_155 = ptsw.propagation_time(pressure_155,velocity_155,X_juno_155,Y_juno_155,time_155_utc,mp_loc_155_km,bs_loc_155_km)
LL_155, HL_pos_155, HL_neg_155 = rv.reconnection_voltages(clock_angle_155,B_perp_155,velocity_155,mp_loc_155_km,times_155)

doy_155_sw = pd.DataFrame()
doy_155_sw = doy_155_sw.assign(UTC=times_155)
doy_155_sw = doy_155_sw.assign(Pressure=pressure_155)
doy_155_sw = doy_155_sw.assign(Velocity=velocity_155)
doy_155_sw = doy_155_sw.assign(Juno_X_Distance=X_juno_155)
doy_155_sw = doy_155_sw.assign(Juno_X_Distance_RJ=X_juno_RJ_155)
doy_155_sw = doy_155_sw.assign(Juno_Y_Distance=Y_juno_RJ_155)
doy_155_sw = doy_155_sw.assign(Total_Travel_Time=overall_time_155)
doy_155_sw = doy_155_sw.assign(Time_Shift_To_Bow_Shock=time_shift_155)
doy_155_sw = doy_155_sw.assign(Magnetosheath_Travel_Time=time_mp_bs_155)
doy_155_sw = doy_155_sw.assign(Magnetopause_To_Ionosphere=to_iono_time_155)
doy_155_sw = doy_155_sw.assign(Time_Impacts_Ionosphere=iono_time_155)
doy_155_sw = doy_155_sw.assign(Bow_Shock_Location=bs_loc_155)
doy_155_sw = doy_155_sw.assign(Bow_Shock_Location_km=bs_loc_155_km)
doy_155_sw = doy_155_sw.assign(Magnetopause_Location=mp_loc_155)
doy_155_sw = doy_155_sw.assign(Magnetopause_Location_km=mp_loc_155_km)
doy_155_sw = doy_155_sw.assign(Br=Br_155)
doy_155_sw = doy_155_sw.assign(Bt=Bt_155)
doy_155_sw = doy_155_sw.assign(Bn=Bn_155)
doy_155_sw = doy_155_sw.assign(Clock_Angle=clock_angle_155)
doy_155_sw = doy_155_sw.assign(B_perp=B_perp_155)
doy_155_sw = doy_155_sw.assign(Low_Latitude_Reconnection=LL_155)
doy_155_sw = doy_155_sw.assign(High_Latitude_Reconnection_By_Pos=HL_pos_155)
doy_155_sw = doy_155_sw.assign(High_Latitude_Reconnection_By_Neg=HL_neg_155)

# # export dataframe to read into other files
doy_155_sw.to_csv(root_folder+'juno_data_doy_155.csv',index=False)

# # # rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
doy_155_sw.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # # export SORTED dataframe to read into other files
doy_155_sw.to_csv(root_folder+'juno_data_doy_155_ionotime.csv',index=False)


# -------- visits 26,27,28 -------

doy156_start = "2016-06-04 00:00:00"
doy158_end = "2016-06-07 00:00:00"

first_156 = spice.str2et(doy156_start)
last_158 = spice.str2et(doy158_end)

time_156_plus, idx_156_plus = between_idxs(et_time,first_156,last_158) # indexes for sw df only

time_156_plus_utc = spice.et2utc(time_156_plus,'ISOD',6)

# grab data
pressure_156_plus = index_compare(idx_156_plus, pressure)
velocity_156_plus = index_compare(idx_156_plus, velocity)

Br_156_plus = index_compare(idx_156_plus,Br)
Bt_156_plus = index_compare(idx_156_plus,Bt)
Bn_156_plus = index_compare(idx_156_plus,Bn)


mp_loc_156_plus, bs_loc_156_plus = jmp.multi_nose(pressure_156_plus)
clock_angle_156_plus, clock_angle_range = mdc.clock_angle_calculator(Bn_156_plus, Bt_156_plus,'no','no')
B_perp_156_plus = mdc.B_perp_calculator(Bn_156_plus, Bt_156_plus, 'no')

mp_loc_156_plus_km = []
bs_loc_156_plus_km = []
# loop through magnetopause and bow shock locations to convert to km
for i in range(len(mp_loc_156_plus)):
    # magnetopause
    mp = mp_loc_156_plus[i]*RJ
    mp_loc_156_plus_km.append(mp)
    # bowshock
    bs = bs_loc_156_plus[i]*RJ
    bs_loc_156_plus_km.append(bs)

times_156_plus,iono_time_156_plus,overall_time_156_plus,time_shift_156_plus,time_mp_bs_156_plus,to_iono_time_156_plus = ptsw.propagation_time(pressure_156_plus,velocity_156_plus,X_juno_156_plus,Y_juno_156_plus,time_156_plus_utc,mp_loc_156_plus_km,bs_loc_156_plus_km)
LL_156_plus, HL_pos_156_plus, HL_neg_156_plus = rv.reconnection_voltages(clock_angle_156_plus,B_perp_156_plus,velocity_156_plus,mp_loc_156_plus_km,times_156_plus)

doy_156_plus_sw = pd.DataFrame()
doy_156_plus_sw = doy_156_plus_sw.assign(UTC=times_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Pressure=pressure_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Velocity=velocity_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Juno_X_Distance=X_juno_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Juno_X_Distance_RJ=X_juno_RJ_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Juno_Y_Distance=Y_juno_RJ_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Total_Travel_Time=overall_time_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Time_Shift_To_Bow_Shock=time_shift_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Magnetosheath_Travel_Time=time_mp_bs_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Magnetopause_To_Ionosphere=to_iono_time_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Time_Impacts_Ionosphere=iono_time_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Bow_Shock_Location=bs_loc_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Bow_Shock_Location_km=bs_loc_156_plus_km)
doy_156_plus_sw = doy_156_plus_sw.assign(Magnetopause_Location=mp_loc_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Magnetopause_Location_km=mp_loc_156_plus_km)
doy_156_plus_sw = doy_156_plus_sw.assign(Br=Br_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Bt=Bt_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Bn=Bn_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Clock_Angle=clock_angle_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(B_perp=B_perp_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(Low_Latitude_Reconnection=LL_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(High_Latitude_Reconnection_By_Pos=HL_pos_156_plus)
doy_156_plus_sw = doy_156_plus_sw.assign(High_Latitude_Reconnection_By_Neg=HL_neg_156_plus)

# # export dataframe to read into other files
doy_156_plus_sw.to_csv(root_folder+'juno_data_doy_156_plus.csv',index=False)

# # rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
doy_156_plus_sw.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
doy_156_plus_sw.to_csv(root_folder+'juno_data_doy_156_plus_ionotime.csv',index=False)


# -------- visits 34 & 35 -----------

doy175_start = "2016-06-23 00:00:00"
doy176_end = "2016-06-24 08:15:29"

first_175 = spice.str2et(doy175_start)
last_176 = spice.str2et(doy176_end)

time_175, idx_175 = between_idxs(et_time,first_175,last_176) # indexes for sw df only

time_175_utc = spice.et2utc(time_175,'ISOD',6)

# grab data
pressure_175 = index_compare(idx_175, pressure)
velocity_175 = index_compare(idx_175, velocity)

Br_175 = index_compare(idx_175,Br)
Bt_175 = index_compare(idx_175,Bt)
Bn_175 = index_compare(idx_175,Bn)

mp_loc_175, bs_loc_175 = jmp.multi_nose(pressure_175)
clock_angle_175, clock_angle_range = mdc.clock_angle_calculator(Bn_175, Bt_175,'no','no')
B_perp_175 = mdc.B_perp_calculator(Bn_175, Bt_175, 'no')

mp_loc_175_km = []
bs_loc_175_km = []
# loop through magnetopause and bow shock locations to convert to km
for i in range(len(mp_loc_175)):
    # magnetopause
    mp = mp_loc_175[i]*RJ
    mp_loc_175_km.append(mp)
    # bowshock
    bs = bs_loc_175[i]*RJ
    bs_loc_175_km.append(bs)

times_175,iono_time_175,overall_time_175,time_shift_175,time_mp_bs_175,to_iono_time_175 = ptsw.propagation_time(pressure_175,velocity_175,X_juno_175,Y_juno_175,time_175_utc,mp_loc_175_km,bs_loc_175_km)
LL_175, HL_pos_175, HL_neg_175 = rv.reconnection_voltages(clock_angle_175,B_perp_175,velocity_175,mp_loc_175_km,times_175)

doy_175_sw = pd.DataFrame()
doy_175_sw = doy_175_sw.assign(UTC=times_175)
doy_175_sw = doy_175_sw.assign(Pressure=pressure_175)
doy_175_sw = doy_175_sw.assign(Velocity=velocity_175)
doy_175_sw = doy_175_sw.assign(Juno_X_Distance=X_juno_175)
doy_175_sw = doy_175_sw.assign(Juno_X_Distance_RJ=X_juno_RJ_175)
doy_175_sw = doy_175_sw.assign(Juno_Y_Distance=Y_juno_RJ_175)
doy_175_sw = doy_175_sw.assign(Total_Travel_Time=overall_time_175)
doy_175_sw = doy_175_sw.assign(Time_Shift_To_Bow_Shock=time_shift_175)
doy_175_sw = doy_175_sw.assign(Magnetosheath_Travel_Time=time_mp_bs_175)
doy_175_sw = doy_175_sw.assign(Magnetopause_To_Ionosphere=to_iono_time_175)
doy_175_sw = doy_175_sw.assign(Time_Impacts_Ionosphere=iono_time_175)
doy_175_sw = doy_175_sw.assign(Bow_Shock_Location=bs_loc_175)
doy_175_sw = doy_175_sw.assign(Bow_Shock_Location_km=bs_loc_175_km)
doy_175_sw = doy_175_sw.assign(Magnetopause_Location=mp_loc_175)
doy_175_sw = doy_175_sw.assign(Magnetopause_Location_km=mp_loc_175_km)
doy_175_sw = doy_175_sw.assign(Br=Br_175)
doy_175_sw = doy_175_sw.assign(Bt=Bt_175)
doy_175_sw = doy_175_sw.assign(Bn=Bn_175)
doy_175_sw = doy_175_sw.assign(Clock_Angle=clock_angle_175)
doy_175_sw = doy_175_sw.assign(B_perp=B_perp_175)
doy_175_sw = doy_175_sw.assign(Low_Latitude_Reconnection=LL_175)
doy_175_sw = doy_175_sw.assign(High_Latitude_Reconnection_By_Pos=HL_pos_175)
doy_175_sw = doy_175_sw.assign(High_Latitude_Reconnection_By_Neg=HL_neg_175)

# # export dataframe to read into other files
doy_175_sw.to_csv(root_folder+'juno_data_doy_175.csv',index=False)

# # rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
doy_175_sw.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
doy_175_sw.to_csv(root_folder+'juno_data_doy_175_ionotime.csv',index=False)
