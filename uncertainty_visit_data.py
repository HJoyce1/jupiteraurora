"""
Created on Tue Sep  3 14:53:02 2024

@author: hannah

this one deals with new ionosphere times if 10% error is added to total travel time from
juno to jupiter
"""

import pandas as pd
import spiceypy as spice
import numpy as np


error='50'

# leap seconds kernal - need this for ephemerous time
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

# dataframes
visit_times = pd.read_csv(f'{root_folder}visit_times.csv')

df = pd.read_csv(f'{root_folder}juno_data_big_df_ionotime_new_errors_{error}.csv')


df_155 = pd.read_csv(f'{root_folder}juno_data_doy_155_ionotime_errors_{error}.csv')
df_156_plus = pd.read_csv(f'{root_folder}juno_data_doy_156_plus_ionotime_errors_{error}.csv')
df_175 = pd.read_csv(f'{root_folder}juno_data_doy_175_ionotime_errors_{error}.csv')


# grab all visit time data
visit_01 = visit_times['Visit_01'].to_numpy()
visit_02 = visit_times['Visit_02'].to_numpy()
visit_03 = visit_times['Visit_03'].to_numpy()
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
visit_23 = visit_times['Visit_23'].to_numpy()
visit_24 = visit_times['Visit_24'].to_numpy()
visit_25 = visit_times['Visit_25'].to_numpy()
visit_26 = visit_times['Visit_26'].to_numpy()
visit_27 = visit_times['Visit_27'].to_numpy()
visit_28 = visit_times['Visit_28'].to_numpy()
visit_34 = visit_times['Visit_34'].to_numpy()
visit_35 = visit_times['Visit_35'].to_numpy()


# ------------ 136 to 154 ----------------


pressure_array = df['SW_Pressure'].to_numpy()
velocity_array = df['SW_Velocity'].to_numpy()
clock_angle = df['Clock_Angle'].to_numpy()
clock_angle_err = df['Clock_Angle_Error'].to_numpy()
b_perp = df['B_Perp'].to_numpy()
Br = df['Br'].to_numpy()
Bt = df['Bt'].to_numpy()
Bn = df['Bn'].to_numpy()
juno_time = df['Juno_Detection_Time'].to_numpy()
bs_loc = df['Bow_Shock_Stand_Off_RJ'].to_numpy()
mp_loc = df['Magnetopause_Stand_Off_RJ'].to_numpy()
LL_rec = df['Low_Latitude_Reconnection_Voltage'].to_numpy()
HL_rec_neg = df['High_Latitude_Reconnection_Voltage_BY_NEG'].to_numpy()
HL_rec_pos = df['High_Latitude_Reconnection_Voltage_BY_POS'].to_numpy()


# ------------

plus_error = df['Ionosphere_Times_Plus_Error'].to_numpy()
minus_error = df['Ionosphere_Times_Minus_Error'].to_numpy()

plus_error_et = []
for p in range(len(plus_error)):
    times = spice.str2et(plus_error[p])
    plus_error_et.append(times)
    
minus_error_et = []
for p in range(len(minus_error)):
    times = spice.str2et(minus_error[p])
    minus_error_et.append(times)
    
plus_error_et = np.array(plus_error_et)
minus_error_et = np.array(minus_error_et)

# --------------- 155 ----------------

pressure_array_155 = df_155['Pressure'].to_numpy()
velocity_array_155 = df_155['Velocity'].to_numpy()
clock_angle_155 = df_155['Clock_Angle'].to_numpy()
b_perp_155 = df_155['B_perp'].to_numpy()
Br_155 = df_155['Br'].to_numpy()
Bt_155 = df_155['Bt'].to_numpy()
Bn_155 = df_155['Bn'].to_numpy()
juno_time_155 = df_155['UTC'].to_numpy()
bs_loc_155 = df_155['Bow_Shock_Location'].to_numpy()
mp_loc_155 = df_155['Magnetopause_Location'].to_numpy()
LL_rec_155 = df_155['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_155 = df_155['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_155 = df_155['High_Latitude_Reconnection_By_Pos'].to_numpy()


# ------------

plus_error_155 = df_155['Ionosphere_Times_Plus_Error'].to_numpy()
minus_error_155 = df_155['Ionosphere_Times_Minus_Error'].to_numpy()

plus_error_et_155 = []
for p in range(len(plus_error_155)):
    times = spice.str2et(plus_error_155[p])
    plus_error_et_155.append(times)
    
minus_error_et_155 = []
for p in range(len(minus_error_155)):
    times = spice.str2et(minus_error_155[p])
    minus_error_et_155.append(times)
    
plus_error_et_155 = np.array(plus_error_et_155)
minus_error_et_155 = np.array(minus_error_et_155)


# --------------- 156 ----------------

pressure_array_156_plus = df_156_plus['Pressure'].to_numpy()
velocity_array_156_plus = df_156_plus['Velocity'].to_numpy()
clock_angle_156_plus = df_156_plus['Clock_Angle'].to_numpy()
b_perp_156_plus = df_156_plus['B_perp'].to_numpy()
Br_156_plus = df_156_plus['Br'].to_numpy()
Bt_156_plus = df_156_plus['Bt'].to_numpy()
Bn_156_plus = df_156_plus['Bn'].to_numpy()
juno_time_156_plus = df_156_plus['UTC'].to_numpy()
bs_loc_156_plus = df_156_plus['Bow_Shock_Location'].to_numpy()
mp_loc_156_plus = df_156_plus['Magnetopause_Location'].to_numpy()
LL_rec_156_plus = df_156_plus['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_156_plus = df_156_plus['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_156_plus = df_156_plus['High_Latitude_Reconnection_By_Pos'].to_numpy()


# ------------

plus_error_156_plus = df_156_plus['Ionosphere_Times_Plus_Error'].to_numpy()
minus_error_156_plus = df_156_plus['Ionosphere_Times_Minus_Error'].to_numpy()

plus_error_et_156_plus = []
for p in range(len(plus_error_156_plus)):
    times = spice.str2et(plus_error_156_plus[p])
    plus_error_et_156_plus.append(times)
    
minus_error_et_156_plus = []
for p in range(len(minus_error_156_plus)):
    times = spice.str2et(minus_error_156_plus[p])
    minus_error_et_156_plus.append(times)
    
plus_error_et_156_plus = np.array(plus_error_et_156_plus)
minus_error_et_156_plus = np.array(minus_error_et_156_plus)


# --------------- 175 ----------------

pressure_array_175 = df_175['Pressure'].to_numpy()
velocity_array_175 = df_175['Velocity'].to_numpy()
clock_angle_175 = df_175['Clock_Angle'].to_numpy()
b_perp_175 = df_175['B_perp'].to_numpy()
Br_175 = df_175['Br'].to_numpy()
Bt_175 = df_175['Bt'].to_numpy()
Bn_175 = df_175['Bn'].to_numpy()
juno_time_175 = df_175['UTC'].to_numpy()
bs_loc_175 = df_175['Bow_Shock_Location'].to_numpy()
mp_loc_175 = df_175['Magnetopause_Location'].to_numpy()
LL_rec_175 = df_175['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_175 = df_175['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_175 = df_175['High_Latitude_Reconnection_By_Pos'].to_numpy()


# ------------

plus_error_175 = df_175['Ionosphere_Times_Plus_Error'].to_numpy()
minus_error_175 = df_175['Ionosphere_Times_Minus_Error'].to_numpy()

plus_error_et_175 = []
for p in range(len(plus_error_175)):
    times = spice.str2et(plus_error_175[p])
    plus_error_et_175.append(times)
    
minus_error_et_175 = []
for p in range(len(minus_error_175)):
    times = spice.str2et(minus_error_175[p])
    minus_error_et_175.append(times)
    
plus_error_et_175 = np.array(plus_error_et_175)
minus_error_et_175 = np.array(minus_error_et_175)


# -----------

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

# -----------

'''
VISIT 01
'''
ettime_01 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_01)):
    times = spice.str2et(visit_01[p])
    ettime_01.append(times)

# grab first and last timesteps of the visit
first_01 = ettime_01[0]
last_01 = ettime_01[-1]

# run function to determine times and index values of those times for the visit
times_plus_01, idx_plus_01 = between_idxs(plus_error_et,first_01,last_01)
times_minus_01, idx_minus_01 = between_idxs(minus_error_et,first_01,last_01)


# ---------  + % error ------------

iono_time_plus_01 = plus_error[idx_plus_01]
pressure_plus_01 = pressure_array[idx_plus_01]
clock_plus_01 = clock_angle[idx_plus_01]
clock_err_plus_01 = clock_angle_err[idx_plus_01]
LL_plus_01 = LL_rec[idx_plus_01]
HL_plus_01_pos = HL_rec_pos[idx_plus_01] 
HL_plus_01_neg = HL_rec_neg[idx_plus_01]
b_perp_plus_01 = b_perp[idx_plus_01]



# --------- - % error ------------

iono_time_minus_01 = minus_error[idx_minus_01]
pressure_minus_01 = pressure_array[idx_minus_01]
clock_minus_01 = clock_angle[idx_minus_01]
clock_err_minus_01 = clock_angle_err[idx_minus_01]
LL_minus_01 = LL_rec[idx_minus_01]
HL_minus_01_pos = HL_rec_pos[idx_minus_01]
HL_minus_01_neg = HL_rec_neg[idx_minus_01]
b_perp_minus_01 = b_perp[idx_minus_01]


# -------- dataframe -----------

# save data out for use in other modules
df_01_p = pd.DataFrame()
df_01_m = pd.DataFrame()

df_01_p = df_01_p.assign(PLUS_Big_DF_Index=idx_plus_01)
df_01_p = df_01_p.assign(PLUS_Ionosphere_Time=iono_time_plus_01)
df_01_p = df_01_p.assign(PLUS_Pressure=pressure_plus_01)
df_01_p = df_01_p.assign(PLUS_Clock_Angle=clock_plus_01)
df_01_p = df_01_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_01)
df_01_p = df_01_p.assign(PLUS_LL_rec_V=LL_plus_01)
df_01_p = df_01_p.assign(PLUS_HL_rec_V_pos=HL_plus_01_pos)
df_01_p = df_01_p.assign(PLUS_HL_rec_V_neg=HL_plus_01_neg)
df_01_p = df_01_p.assign(PLUS_B_Perp=b_perp_plus_01)
df_01_p.to_csv(f'{root_folder}visit_01_times_plus_error_{error}_new.csv',index=False)


df_01_m = df_01_m.assign(MINUS_Big_DF_Index=idx_minus_01)
df_01_m = df_01_m.assign(MINUS_Ionosphere_Time=iono_time_minus_01)
df_01_m = df_01_m.assign(MINUS_Pressure=pressure_minus_01)
df_01_m = df_01_m.assign(MINUS_Clock_Angle=clock_minus_01)
df_01_m = df_01_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_01)
df_01_m = df_01_m.assign(MINUS_LL_rec_V=LL_minus_01)
df_01_m = df_01_m.assign(MINUS_HL_rec_V_pos=HL_minus_01_pos)
df_01_m = df_01_m.assign(MINUS_HL_rec_V_neg=HL_minus_01_neg)
df_01_m = df_01_m.assign(MINUS_B_Perp=b_perp_minus_01)
df_01_m.to_csv(f'{root_folder}visit_01_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 02
'''
ettime_02 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_02)):
    times = spice.str2et(visit_02[p])
    ettime_02.append(times)

# grab first and last timesteps of the visit
first_02 = ettime_02[0]
last_02 = ettime_02[-1]

# run function to determine times and index values of those times for the visit
times_plus_02, idx_plus_02 = between_idxs(plus_error_et,first_02,last_02)
times_minus_02, idx_minus_02 = between_idxs(minus_error_et,first_02,last_02)


# ---------  + % error ------------

iono_time_plus_02 = plus_error[idx_plus_02]
pressure_plus_02 = pressure_array[idx_plus_02]
clock_plus_02 = clock_angle[idx_plus_02]
clock_err_plus_02 = clock_angle_err[idx_plus_02]
LL_plus_02 = LL_rec[idx_plus_02]
HL_plus_02_pos = HL_rec_pos[idx_plus_02] 
HL_plus_02_neg = HL_rec_neg[idx_plus_02]
b_perp_plus_02 = b_perp[idx_plus_02]



# --------- - % error ------------

iono_time_minus_02 = minus_error[idx_minus_02]
pressure_minus_02 = pressure_array[idx_minus_02]
clock_minus_02 = clock_angle[idx_minus_02]
clock_err_minus_02 = clock_angle_err[idx_minus_02]
LL_minus_02 = LL_rec[idx_minus_02]
HL_minus_02_pos = HL_rec_pos[idx_minus_02]
HL_minus_02_neg = HL_rec_neg[idx_minus_02]
b_perp_minus_02 = b_perp[idx_minus_02]


# -------- dataframe -----------

# save data out for use in other modules
df_02_p = pd.DataFrame()
df_02_m = pd.DataFrame()

df_02_p = df_02_p.assign(PLUS_Big_DF_Index=idx_plus_02)
df_02_p = df_02_p.assign(PLUS_Ionosphere_Time=iono_time_plus_02)
df_02_p = df_02_p.assign(PLUS_Pressure=pressure_plus_02)
df_02_p = df_02_p.assign(PLUS_Clock_Angle=clock_plus_02)
df_02_p = df_02_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_02)
df_02_p = df_02_p.assign(PLUS_LL_rec_V=LL_plus_02)
df_02_p = df_02_p.assign(PLUS_HL_rec_V_pos=HL_plus_02_pos)
df_02_p = df_02_p.assign(PLUS_HL_rec_V_neg=HL_plus_02_neg)
df_02_p = df_02_p.assign(PLUS_B_Perp=b_perp_plus_02)
df_02_p.to_csv(f'{root_folder}visit_02_times_plus_error_{error}_new.csv',index=False)


df_02_m = df_02_m.assign(MINUS_Big_DF_Index=idx_minus_02)
df_02_m = df_02_m.assign(MINUS_Ionosphere_Time=iono_time_minus_02)
df_02_m = df_02_m.assign(MINUS_Pressure=pressure_minus_02)
df_02_m = df_02_m.assign(MINUS_Clock_Angle=clock_minus_02)
df_02_m = df_02_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_02)
df_02_m = df_02_m.assign(MINUS_LL_rec_V=LL_minus_02)
df_02_m = df_02_m.assign(MINUS_HL_rec_V_pos=HL_minus_02_pos)
df_02_m = df_02_m.assign(MINUS_HL_rec_V_neg=HL_minus_02_neg)
df_02_m = df_02_m.assign(MINUS_B_Perp=b_perp_minus_02)
df_02_m.to_csv(f'{root_folder}visit_02_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 03
'''
ettime_03 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_03)):
    times = spice.str2et(visit_03[p])
    ettime_03.append(times)

# grab first and last timesteps of the visit
first_03 = ettime_03[0]
last_03 = ettime_03[-1]

# run function to determine times and index values of those times for the visit
times_plus_03, idx_plus_03 = between_idxs(plus_error_et,first_03,last_03)
times_minus_03, idx_minus_03 = between_idxs(minus_error_et,first_03,last_03)


# ---------  + % error ------------

iono_time_plus_03 = plus_error[idx_plus_03]
pressure_plus_03 = pressure_array[idx_plus_03]
clock_plus_03 = clock_angle[idx_plus_03]
clock_err_plus_03 = clock_angle_err[idx_plus_03]
LL_plus_03 = LL_rec[idx_plus_03]
HL_plus_03_pos = HL_rec_pos[idx_plus_03] 
HL_plus_03_neg = HL_rec_neg[idx_plus_03]
b_perp_plus_03 = b_perp[idx_plus_03]



# --------- - % error ------------

iono_time_minus_03 = minus_error[idx_minus_03]
pressure_minus_03 = pressure_array[idx_minus_03]
clock_minus_03 = clock_angle[idx_minus_03]
clock_err_minus_03 = clock_angle_err[idx_minus_03]
LL_minus_03 = LL_rec[idx_minus_03]
HL_minus_03_pos = HL_rec_pos[idx_minus_03]
HL_minus_03_neg = HL_rec_neg[idx_minus_03]
b_perp_minus_03 = b_perp[idx_minus_03]


# -------- dataframe -----------

# save data out for use in other modules
df_03_p = pd.DataFrame()
df_03_m = pd.DataFrame()

df_03_p = df_03_p.assign(PLUS_Big_DF_Index=idx_plus_03)
df_03_p = df_03_p.assign(PLUS_Ionosphere_Time=iono_time_plus_03)
df_03_p = df_03_p.assign(PLUS_Pressure=pressure_plus_03)
df_03_p = df_03_p.assign(PLUS_Clock_Angle=clock_plus_03)
df_03_p = df_03_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_03)
df_03_p = df_03_p.assign(PLUS_LL_rec_V=LL_plus_03)
df_03_p = df_03_p.assign(PLUS_HL_rec_V_pos=HL_plus_03_pos)
df_03_p = df_03_p.assign(PLUS_HL_rec_V_neg=HL_plus_03_neg)
df_03_p = df_03_p.assign(PLUS_B_Perp=b_perp_plus_03)
df_03_p.to_csv(f'{root_folder}visit_03_times_plus_error_{error}_new.csv',index=False)


df_03_m = df_03_m.assign(MINUS_Big_DF_Index=idx_minus_03)
df_03_m = df_03_m.assign(MINUS_Ionosphere_Time=iono_time_minus_03)
df_03_m = df_03_m.assign(MINUS_Pressure=pressure_minus_03)
df_03_m = df_03_m.assign(MINUS_Clock_Angle=clock_minus_03)
df_03_m = df_03_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_03)
df_03_m = df_03_m.assign(MINUS_LL_rec_V=LL_minus_03)
df_03_m = df_03_m.assign(MINUS_HL_rec_V_pos=HL_minus_03_pos)
df_03_m = df_03_m.assign(MINUS_HL_rec_V_neg=HL_minus_03_neg)
df_03_m = df_03_m.assign(MINUS_B_Perp=b_perp_minus_03)
df_03_m.to_csv(f'{root_folder}visit_03_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 04
'''
ettime_04 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_04)):
    times = spice.str2et(visit_04[p])
    ettime_04.append(times)

# grab first and last timesteps of the visit
first_04 = ettime_04[0]
last_04 = ettime_04[-1]

# run function to determine times and index values of those times for the visit
times_plus_04, idx_plus_04 = between_idxs(plus_error_et,first_04,last_04)
times_minus_04, idx_minus_04 = between_idxs(minus_error_et,first_04,last_04)


# ---------  + % error ------------

iono_time_plus_04 = plus_error[idx_plus_04]
pressure_plus_04 = pressure_array[idx_plus_04]
clock_plus_04 = clock_angle[idx_plus_04]
clock_err_plus_04 = clock_angle_err[idx_plus_04]
LL_plus_04 = LL_rec[idx_plus_04]
HL_plus_04_pos = HL_rec_pos[idx_plus_04] 
HL_plus_04_neg = HL_rec_neg[idx_plus_04]
b_perp_plus_04 = b_perp[idx_plus_04]



# --------- - % error ------------

iono_time_minus_04 = minus_error[idx_minus_04]
pressure_minus_04 = pressure_array[idx_minus_04]
clock_minus_04 = clock_angle[idx_minus_04]
clock_err_minus_04 = clock_angle_err[idx_minus_04]
LL_minus_04 = LL_rec[idx_minus_04]
HL_minus_04_pos = HL_rec_pos[idx_minus_04]
HL_minus_04_neg = HL_rec_neg[idx_minus_04]
b_perp_minus_04 = b_perp[idx_minus_04]


# -------- dataframe -----------

# save data out for use in other modules
df_04_p = pd.DataFrame()
df_04_m = pd.DataFrame()

df_04_p = df_04_p.assign(PLUS_Big_DF_Index=idx_plus_04)
df_04_p = df_04_p.assign(PLUS_Ionosphere_Time=iono_time_plus_04)
df_04_p = df_04_p.assign(PLUS_Pressure=pressure_plus_04)
df_04_p = df_04_p.assign(PLUS_Clock_Angle=clock_plus_04)
df_04_p = df_04_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_04)
df_04_p = df_04_p.assign(PLUS_LL_rec_V=LL_plus_04)
df_04_p = df_04_p.assign(PLUS_HL_rec_V_pos=HL_plus_04_pos)
df_04_p = df_04_p.assign(PLUS_HL_rec_V_neg=HL_plus_04_neg)
df_04_p = df_04_p.assign(PLUS_B_Perp=b_perp_plus_04)
df_04_p.to_csv(f'{root_folder}visit_04_times_plus_error_{error}_new.csv',index=False)


df_04_m = df_04_m.assign(MINUS_Big_DF_Index=idx_minus_04)
df_04_m = df_04_m.assign(MINUS_Ionosphere_Time=iono_time_minus_04)
df_04_m = df_04_m.assign(MINUS_Pressure=pressure_minus_04)
df_04_m = df_04_m.assign(MINUS_Clock_Angle=clock_minus_04)
df_04_m = df_04_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_04)
df_04_m = df_04_m.assign(MINUS_LL_rec_V=LL_minus_04)
df_04_m = df_04_m.assign(MINUS_HL_rec_V_pos=HL_minus_04_pos)
df_04_m = df_04_m.assign(MINUS_HL_rec_V_neg=HL_minus_04_neg)
df_04_m = df_04_m.assign(MINUS_B_Perp=b_perp_minus_04)
df_04_m.to_csv(f'{root_folder}visit_04_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 05
'''
ettime_05 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_05)):
    times = spice.str2et(visit_05[p])
    ettime_05.append(times)

# grab first and last timesteps of the visit
first_05 = ettime_05[0]
last_05 = ettime_05[-1]

# run function to determine times and index values of those times for the visit
times_plus_05, idx_plus_05 = between_idxs(plus_error_et,first_05,last_05)
times_minus_05, idx_minus_05 = between_idxs(minus_error_et,first_05,last_05)


# ---------  + % error ------------

iono_time_plus_05 = plus_error[idx_plus_05]
pressure_plus_05 = pressure_array[idx_plus_05]
clock_plus_05 = clock_angle[idx_plus_05]
clock_err_plus_05 = clock_angle_err[idx_plus_05]
LL_plus_05 = LL_rec[idx_plus_05]
HL_plus_05_pos = HL_rec_pos[idx_plus_05] 
HL_plus_05_neg = HL_rec_neg[idx_plus_05]
b_perp_plus_05 = b_perp[idx_plus_05]


# --------- - % error ------------

iono_time_minus_05 = minus_error[idx_minus_05]
pressure_minus_05 = pressure_array[idx_minus_05]
clock_minus_05 = clock_angle[idx_minus_05]
clock_err_minus_05 = clock_angle_err[idx_minus_05]
LL_minus_05 = LL_rec[idx_minus_05]
HL_minus_05_pos = HL_rec_pos[idx_minus_05]
HL_minus_05_neg = HL_rec_neg[idx_minus_05]
b_perp_minus_05 = b_perp[idx_minus_05]



# -------- dataframe -----------

# save data out for use in other modules
df_05_p = pd.DataFrame()
df_05_m = pd.DataFrame()

df_05_p = df_05_p.assign(PLUS_Big_DF_Index=idx_plus_05)
df_05_p = df_05_p.assign(PLUS_Ionosphere_Time=iono_time_plus_05)
df_05_p = df_05_p.assign(PLUS_Pressure=pressure_plus_05)
df_05_p = df_05_p.assign(PLUS_Clock_Angle=clock_plus_05)
df_05_p = df_05_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_05)
df_05_p = df_05_p.assign(PLUS_LL_rec_V=LL_plus_05)
df_05_p = df_05_p.assign(PLUS_HL_rec_V_pos=HL_plus_05_pos)
df_05_p = df_05_p.assign(PLUS_HL_rec_V_neg=HL_plus_05_neg)
df_05_p = df_05_p.assign(PLUS_B_Perp=b_perp_plus_05)
df_05_p.to_csv(f'{root_folder}visit_05_times_plus_error_{error}_new.csv',index=False)


df_05_m = df_05_m.assign(MINUS_Big_DF_Index=idx_minus_05)
df_05_m = df_05_m.assign(MINUS_Ionosphere_Time=iono_time_minus_05)
df_05_m = df_05_m.assign(MINUS_Pressure=pressure_minus_05)
df_05_m = df_05_m.assign(MINUS_Clock_Angle=clock_minus_05)
df_05_m = df_05_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_05)
df_05_m = df_05_m.assign(MINUS_LL_rec_V=LL_minus_05)
df_05_m = df_05_m.assign(MINUS_HL_rec_V_pos=HL_minus_05_pos)
df_05_m = df_05_m.assign(MINUS_HL_rec_V_neg=HL_minus_05_neg)
df_05_m = df_05_m.assign(MINUS_B_Perp=b_perp_minus_05)
df_05_m.to_csv(f'{root_folder}visit_05_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 08
'''
ettime_08 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_08)):
    times = spice.str2et(visit_08[p])
    ettime_08.append(times)

# grab first and last timesteps of the visit
first_08 = ettime_08[0]
last_08 = ettime_08[-1]

# run function to determine times and index values of those times for the visit
times_plus_08, idx_plus_08 = between_idxs(plus_error_et,first_08,last_08)
times_minus_08, idx_minus_08 = between_idxs(minus_error_et,first_08,last_08)


# ---------  + % error ------------

iono_time_plus_08 = plus_error[idx_plus_08]
pressure_plus_08 = pressure_array[idx_plus_08]
clock_plus_08 = clock_angle[idx_plus_08]
clock_err_plus_08 = clock_angle_err[idx_plus_08]
LL_plus_08 = LL_rec[idx_plus_08]
HL_plus_08_pos = HL_rec_pos[idx_plus_08] 
HL_plus_08_neg = HL_rec_neg[idx_plus_08]
b_perp_plus_08 = b_perp[idx_plus_08]



# --------- - % error ------------

iono_time_minus_08 = minus_error[idx_minus_08]
pressure_minus_08 = pressure_array[idx_minus_08]
clock_minus_08 = clock_angle[idx_minus_08]
clock_err_minus_08 = clock_angle_err[idx_minus_08]
LL_minus_08 = LL_rec[idx_minus_08]
HL_minus_08_pos = HL_rec_pos[idx_minus_08]
HL_minus_08_neg = HL_rec_neg[idx_minus_08]
b_perp_minus_08 = b_perp[idx_minus_08]



# -------- dataframe -----------

# save data out for use in other modules
df_08_p = pd.DataFrame()
df_08_m = pd.DataFrame()

df_08_p = df_08_p.assign(PLUS_Big_DF_Index=idx_plus_08)
df_08_p = df_08_p.assign(PLUS_Ionosphere_Time=iono_time_plus_08)
df_08_p = df_08_p.assign(PLUS_Pressure=pressure_plus_08)
df_08_p = df_08_p.assign(PLUS_Clock_Angle=clock_plus_08)
df_08_p = df_08_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_08)
df_08_p = df_08_p.assign(PLUS_LL_rec_V=LL_plus_08)
df_08_p = df_08_p.assign(PLUS_HL_rec_V_pos=HL_plus_08_pos)
df_08_p = df_08_p.assign(PLUS_HL_rec_V_neg=HL_plus_08_neg)
df_08_p = df_08_p.assign(PLUS_B_Perp=b_perp_plus_08)
df_08_p.to_csv(f'{root_folder}visit_08_times_plus_error_{error}_new.csv',index=False)


df_08_m = df_08_m.assign(MINUS_Big_DF_Index=idx_minus_08)
df_08_m = df_08_m.assign(MINUS_Ionosphere_Time=iono_time_minus_08)
df_08_m = df_08_m.assign(MINUS_Pressure=pressure_minus_08)
df_08_m = df_08_m.assign(MINUS_Clock_Angle=clock_minus_08)
df_08_m = df_08_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_08)
df_08_m = df_08_m.assign(MINUS_LL_rec_V=LL_minus_08)
df_08_m = df_08_m.assign(MINUS_HL_rec_V_pos=HL_minus_08_pos)
df_08_m = df_08_m.assign(MINUS_HL_rec_V_neg=HL_minus_08_neg)
df_08_m = df_08_m.assign(MINUS_B_Perp=b_perp_minus_08)
df_08_m.to_csv(f'{root_folder}visit_08_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 09
'''
ettime_09 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_09)):
    times = spice.str2et(visit_09[p])
    ettime_09.append(times)

# grab first and last timesteps of the visit
first_09 = ettime_09[0]
last_09 = ettime_09[-1]

# run function to determine times and index values of those times for the visit
times_plus_09, idx_plus_09 = between_idxs(plus_error_et,first_09,last_09)
times_minus_09, idx_minus_09 = between_idxs(minus_error_et,first_09,last_09)


# ---------  + % error ------------

iono_time_plus_09 = plus_error[idx_plus_09]
pressure_plus_09 = pressure_array[idx_plus_09]
clock_plus_09 = clock_angle[idx_plus_09]
clock_err_plus_09 = clock_angle_err[idx_plus_09]
LL_plus_09 = LL_rec[idx_plus_09]
HL_plus_09_pos = HL_rec_pos[idx_plus_09] 
HL_plus_09_neg = HL_rec_neg[idx_plus_09]
b_perp_plus_09 = b_perp[idx_plus_09]


# --------- - % error ------------

iono_time_minus_09 = minus_error[idx_minus_09]
pressure_minus_09 = pressure_array[idx_minus_09]
clock_minus_09 = clock_angle[idx_minus_09]
clock_err_minus_09 = clock_angle_err[idx_minus_09]
LL_minus_09 = LL_rec[idx_minus_09]
HL_minus_09_pos = HL_rec_pos[idx_minus_09]
HL_minus_09_neg = HL_rec_neg[idx_minus_09]
b_perp_minus_09 = b_perp[idx_minus_09]



# -------- dataframe -----------

# save data out for use in other modules
df_09_p = pd.DataFrame()
df_09_m = pd.DataFrame()

df_09_p = df_09_p.assign(PLUS_Big_DF_Index=idx_plus_09)
df_09_p = df_09_p.assign(PLUS_Ionosphere_Time=iono_time_plus_09)
df_09_p = df_09_p.assign(PLUS_Pressure=pressure_plus_09)
df_09_p = df_09_p.assign(PLUS_Clock_Angle=clock_plus_09)
df_09_p = df_09_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_09)
df_09_p = df_09_p.assign(PLUS_LL_rec_V=LL_plus_09)
df_09_p = df_09_p.assign(PLUS_HL_rec_V_pos=HL_plus_09_pos)
df_09_p = df_09_p.assign(PLUS_HL_rec_V_neg=HL_plus_09_neg)
df_09_p = df_09_p.assign(PLUS_B_Perp=b_perp_plus_09)
df_09_p.to_csv(f'{root_folder}visit_09_times_plus_error_{error}_new.csv',index=False)


df_09_m = df_09_m.assign(MINUS_Big_DF_Index=idx_minus_09)
df_09_m = df_09_m.assign(MINUS_Ionosphere_Time=iono_time_minus_09)
df_09_m = df_09_m.assign(MINUS_Pressure=pressure_minus_09)
df_09_m = df_09_m.assign(MINUS_Clock_Angle=clock_minus_09)
df_09_m = df_09_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_09)
df_09_m = df_09_m.assign(MINUS_LL_rec_V=LL_minus_09)
df_09_m = df_09_m.assign(MINUS_HL_rec_V_pos=HL_minus_09_pos)
df_09_m = df_09_m.assign(MINUS_HL_rec_V_neg=HL_minus_09_neg)
df_09_m = df_09_m.assign(MINUS_B_Perp=b_perp_minus_09)
df_09_m.to_csv(f'{root_folder}visit_09_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 10
'''
ettime_10 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_10)):
    times = spice.str2et(visit_10[p])
    ettime_10.append(times)

# grab first and last timesteps of the visit
first_10 = ettime_10[0]
last_10 = ettime_10[-1]

# run function to determine times and index values of those times for the visit
times_plus_10, idx_plus_10 = between_idxs(plus_error_et,first_10,last_10)
times_minus_10, idx_minus_10 = between_idxs(minus_error_et,first_10,last_10)


# ---------  + % error ------------

iono_time_plus_10 = plus_error[idx_plus_10]
pressure_plus_10 = pressure_array[idx_plus_10]
clock_plus_10 = clock_angle[idx_plus_10]
clock_err_plus_10 = clock_angle_err[idx_plus_10]
LL_plus_10 = LL_rec[idx_plus_10]
HL_plus_10_pos = HL_rec_pos[idx_plus_10] 
HL_plus_10_neg = HL_rec_neg[idx_plus_10]
b_perp_plus_10 = b_perp[idx_plus_10]


# --------- - % error ------------

iono_time_minus_10 = minus_error[idx_minus_10]
pressure_minus_10 = pressure_array[idx_minus_10]
clock_minus_10 = clock_angle[idx_minus_10]
clock_err_minus_10 = clock_angle_err[idx_minus_10]
LL_minus_10 = LL_rec[idx_minus_10]
HL_minus_10_pos = HL_rec_pos[idx_minus_10] 
HL_minus_10_neg = HL_rec_neg[idx_minus_10]
b_perp_minus_10 = b_perp[idx_minus_10]



# -------- dataframe -----------

# save data out for use in other modules
df_10_p = pd.DataFrame()
df_10_m = pd.DataFrame()

df_10_p = df_10_p.assign(PLUS_Big_DF_Index=idx_plus_10)
df_10_p = df_10_p.assign(PLUS_Ionosphere_Time=iono_time_plus_10)
df_10_p = df_10_p.assign(PLUS_Pressure=pressure_plus_10)
df_10_p = df_10_p.assign(PLUS_Clock_Angle=clock_plus_10)
df_10_p = df_10_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_10)
df_10_p = df_10_p.assign(PLUS_LL_rec_V=LL_plus_10)
df_10_p = df_10_p.assign(PLUS_HL_rec_V_pos=HL_plus_10_pos)
df_10_p = df_10_p.assign(PLUS_HL_rec_V_neg=HL_plus_10_neg)
df_10_p = df_10_p.assign(PLUS_B_Perp=b_perp_plus_10)
df_10_p.to_csv(f'{root_folder}visit_10_times_plus_error_{error}_new.csv',index=False)


df_10_m = df_10_m.assign(MINUS_Big_DF_Index=idx_minus_10)
df_10_m = df_10_m.assign(MINUS_Ionosphere_Time=iono_time_minus_10)
df_10_m = df_10_m.assign(MINUS_Pressure=pressure_minus_10)
df_10_m = df_10_m.assign(MINUS_Clock_Angle=clock_minus_10)
df_10_m = df_10_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_10)
df_10_m = df_10_m.assign(MINUS_LL_rec_V=LL_minus_10)
df_10_m = df_10_m.assign(MINUS_HL_rec_V_pos=HL_minus_10_pos)
df_10_m = df_10_m.assign(MINUS_HL_rec_V_neg=HL_minus_10_neg)
df_10_m = df_10_m.assign(MINUS_B_Perp=b_perp_minus_10)
df_10_m.to_csv(f'{root_folder}visit_10_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 11
'''
ettime_11 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_11)):
    times = spice.str2et(visit_11[p])
    ettime_11.append(times)

# grab first and last timesteps of the visit
first_11 = ettime_11[0]
last_11 = ettime_11[-1]

# run function to determine times and index values of those times for the visit
times_plus_11, idx_plus_11 = between_idxs(plus_error_et,first_11,last_11)
times_minus_11, idx_minus_11 = between_idxs(minus_error_et,first_11,last_11)


# ---------  + % error ------------

iono_time_plus_11 = plus_error[idx_plus_11]
pressure_plus_11 = pressure_array[idx_plus_11]
clock_plus_11 = clock_angle[idx_plus_11]
clock_err_plus_11 = clock_angle_err[idx_plus_11]
LL_plus_11 = LL_rec[idx_plus_11]
HL_plus_11_pos = HL_rec_pos[idx_plus_11] 
HL_plus_11_neg = HL_rec_neg[idx_plus_11]
b_perp_plus_11 = b_perp[idx_plus_11]


# --------- - % error ------------

iono_time_minus_11 = minus_error[idx_minus_11]
pressure_minus_11 = pressure_array[idx_minus_11]
clock_minus_11 = clock_angle[idx_minus_11]
clock_err_minus_11 = clock_angle_err[idx_minus_11]
LL_minus_11 = LL_rec[idx_minus_11]
HL_minus_11_pos = HL_rec_pos[idx_minus_11]
HL_minus_11_neg = HL_rec_neg[idx_minus_11]
b_perp_minus_11 = b_perp[idx_minus_11]



# -------- dataframe -----------

# save data out for use in other modules
df_11_p = pd.DataFrame()
df_11_m = pd.DataFrame()

df_11_p = df_11_p.assign(PLUS_Big_DF_Index=idx_plus_11)
df_11_p = df_11_p.assign(PLUS_Ionosphere_Time=iono_time_plus_11)
df_11_p = df_11_p.assign(PLUS_Pressure=pressure_plus_11)
df_11_p = df_11_p.assign(PLUS_Clock_Angle=clock_plus_11)
df_11_p = df_11_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_11)
df_11_p = df_11_p.assign(PLUS_LL_rec_V=LL_plus_11)
df_11_p = df_11_p.assign(PLUS_HL_rec_V_pos=HL_plus_11_pos)
df_11_p = df_11_p.assign(PLUS_HL_rec_V_neg=HL_plus_11_neg)
df_11_p = df_11_p.assign(PLUS_B_Perp=b_perp_plus_11)
df_11_p.to_csv(f'{root_folder}visit_11_times_plus_error_{error}_new.csv',index=False)


df_11_m = df_11_m.assign(MINUS_Big_DF_Index=idx_minus_11)
df_11_m = df_11_m.assign(MINUS_Ionosphere_Time=iono_time_minus_11)
df_11_m = df_11_m.assign(MINUS_Pressure=pressure_minus_11)
df_11_m = df_11_m.assign(MINUS_Clock_Angle=clock_minus_11)
df_11_m = df_11_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_11)
df_11_m = df_11_m.assign(MINUS_LL_rec_V=LL_minus_11)
df_11_m = df_11_m.assign(MINUS_HL_rec_V_pos=HL_minus_11_pos)
df_11_m = df_11_m.assign(MINUS_HL_rec_V_neg=HL_minus_11_neg)
df_11_m = df_11_m.assign(MINUS_B_Perp=b_perp_minus_11)
df_11_m.to_csv(f'{root_folder}visit_11_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 12
'''
ettime_12 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_12)):
    times = spice.str2et(visit_12[p])
    ettime_12.append(times)

# grab first and last timesteps of the visit
first_12 = ettime_12[0]
last_12 = ettime_12[-1]

# run function to determine times and index values of those times for the visit
times_plus_12, idx_plus_12 = between_idxs(plus_error_et,first_12,last_12)
times_minus_12, idx_minus_12 = between_idxs(minus_error_et,first_12,last_12)


# ---------  + % error ------------

iono_time_plus_12 = plus_error[idx_plus_12]
pressure_plus_12 = pressure_array[idx_plus_12]
clock_plus_12 = clock_angle[idx_plus_12]
clock_err_plus_12 = clock_angle_err[idx_plus_12]
LL_plus_12 = LL_rec[idx_plus_12]
HL_plus_12_pos = HL_rec_pos[idx_plus_12] 
HL_plus_12_neg = HL_rec_neg[idx_plus_12]
b_perp_plus_12 = b_perp[idx_plus_12]


# --------- - % error ------------

iono_time_minus_12 = minus_error[idx_minus_12]
pressure_minus_12 = pressure_array[idx_minus_12]
clock_minus_12 = clock_angle[idx_minus_12]
clock_err_minus_12 = clock_angle_err[idx_minus_12]
LL_minus_12 = LL_rec[idx_minus_12]
HL_minus_12_pos = HL_rec_pos[idx_minus_12]
HL_minus_12_neg = HL_rec_neg[idx_minus_12]
b_perp_minus_12 = b_perp[idx_minus_12]



# -------- dataframe -----------

# save data out for use in other modules
df_12_p = pd.DataFrame()
df_12_m = pd.DataFrame()

df_12_p = df_12_p.assign(PLUS_Big_DF_Index=idx_plus_12)
df_12_p = df_12_p.assign(PLUS_Ionosphere_Time=iono_time_plus_12)
df_12_p = df_12_p.assign(PLUS_Pressure=pressure_plus_12)
df_12_p = df_12_p.assign(PLUS_Clock_Angle=clock_plus_12)
df_12_p = df_12_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_12)
df_12_p = df_12_p.assign(PLUS_LL_rec_V=LL_plus_12)
df_12_p = df_12_p.assign(PLUS_HL_rec_V_pos=HL_plus_12_pos)
df_12_p = df_12_p.assign(PLUS_HL_rec_V_neg=HL_plus_12_neg)
df_12_p = df_12_p.assign(PLUS_B_Perp=b_perp_plus_12)
df_12_p.to_csv(f'{root_folder}visit_12_times_plus_error_{error}_new.csv',index=False)


df_12_m = df_12_m.assign(MINUS_Big_DF_Index=idx_minus_12)
df_12_m = df_12_m.assign(MINUS_Ionosphere_Time=iono_time_minus_12)
df_12_m = df_12_m.assign(MINUS_Pressure=pressure_minus_12)
df_12_m = df_12_m.assign(MINUS_Clock_Angle=clock_minus_12)
df_12_m = df_12_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_12)
df_12_m = df_12_m.assign(MINUS_LL_rec_V=LL_minus_12)
df_12_m = df_12_m.assign(MINUS_HL_rec_V_pos=HL_minus_12_pos)
df_12_m = df_12_m.assign(MINUS_HL_rec_V_neg=HL_minus_12_neg)
df_12_m = df_12_m.assign(MINUS_B_Perp=b_perp_minus_12)
df_12_m.to_csv(f'{root_folder}visit_12_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 15
'''
ettime_15 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_15)):
    times = spice.str2et(visit_15[p])
    ettime_15.append(times)

# grab first and last timesteps of the visit
first_15 = ettime_15[0]
last_15 = ettime_15[-1]

# run function to determine times and index values of those times for the visit
times_plus_15, idx_plus_15 = between_idxs(plus_error_et,first_15,last_15)
times_minus_15, idx_minus_15 = between_idxs(minus_error_et,first_15,last_15)


# ---------  + % error ------------

iono_time_plus_15 = plus_error[idx_plus_15]
pressure_plus_15 = pressure_array[idx_plus_15]
clock_plus_15 = clock_angle[idx_plus_15]
clock_err_plus_15 = clock_angle_err[idx_plus_15]
LL_plus_15 = LL_rec[idx_plus_15]
HL_plus_15_pos = HL_rec_pos[idx_plus_15] 
HL_plus_15_neg = HL_rec_neg[idx_plus_15]
b_perp_plus_15 = b_perp[idx_plus_15]


# --------- - % error ------------

iono_time_minus_15 = minus_error[idx_minus_15]
pressure_minus_15 = pressure_array[idx_minus_15]
clock_minus_15 = clock_angle[idx_minus_15]
clock_err_minus_15 = clock_angle_err[idx_minus_15]
LL_minus_15 = LL_rec[idx_minus_15]
HL_minus_15_pos = HL_rec_pos[idx_minus_15]
HL_minus_15_neg = HL_rec_neg[idx_minus_15]
b_perp_minus_15 = b_perp[idx_minus_15]


# -------- dataframe -----------

# save data out for use in other modules
df_15_p = pd.DataFrame()
df_15_m = pd.DataFrame()

df_15_p = df_15_p.assign(PLUS_Big_DF_Index=idx_plus_15)
df_15_p = df_15_p.assign(PLUS_Ionosphere_Time=iono_time_plus_15)
df_15_p = df_15_p.assign(PLUS_Pressure=pressure_plus_15)
df_15_p = df_15_p.assign(PLUS_Clock_Angle=clock_plus_15)
df_15_p = df_15_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_15)
df_15_p = df_15_p.assign(PLUS_LL_rec_V=LL_plus_15)
df_15_p = df_15_p.assign(PLUS_HL_rec_V_pos=HL_plus_15_pos)
df_15_p = df_15_p.assign(PLUS_HL_rec_V_neg=HL_plus_15_neg)
df_15_p = df_15_p.assign(PLUS_B_Perp=b_perp_plus_15)
df_15_p.to_csv(f'{root_folder}visit_15_times_plus_error_{error}_new.csv',index=False)


df_15_m = df_15_m.assign(MINUS_Big_DF_Index=idx_minus_15)
df_15_m = df_15_m.assign(MINUS_Ionosphere_Time=iono_time_minus_15)
df_15_m = df_15_m.assign(MINUS_Pressure=pressure_minus_15)
df_15_m = df_15_m.assign(MINUS_Clock_Angle=clock_minus_15)
df_15_m = df_15_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_15)
df_15_m = df_15_m.assign(MINUS_LL_rec_V=LL_minus_15)
df_15_m = df_15_m.assign(MINUS_HL_rec_V_pos=HL_minus_15_pos)
df_15_m = df_15_m.assign(MINUS_HL_rec_V_neg=HL_minus_15_neg)
df_15_m = df_15_m.assign(MINUS_B_Perp=b_perp_minus_15)
df_15_m.to_csv(f'{root_folder}visit_15_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 16
'''
ettime_16 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_16)):
    times = spice.str2et(visit_16[p])
    ettime_16.append(times)

# grab first and last timesteps of the visit
first_16 = ettime_16[0]
last_16 = ettime_16[-1]

# run function to determine times and index values of those times for the visit
times_plus_16, idx_plus_16 = between_idxs(plus_error_et,first_16,last_16)
times_minus_16, idx_minus_16 = between_idxs(minus_error_et,first_16,last_16)


# ---------  + % error ------------

iono_time_plus_16 = plus_error[idx_plus_16]
pressure_plus_16 = pressure_array[idx_plus_16]
clock_plus_16 = clock_angle[idx_plus_16]
clock_err_plus_16 = clock_angle_err[idx_plus_16]
LL_plus_16 = LL_rec[idx_plus_16]
HL_plus_16_pos = HL_rec_pos[idx_plus_16] 
HL_plus_16_neg = HL_rec_neg[idx_plus_16]
b_perp_plus_16 = b_perp[idx_plus_16]


# --------- - % error ------------

iono_time_minus_16 = minus_error[idx_minus_16]
pressure_minus_16  = pressure_array[idx_minus_16]
clock_minus_16 = clock_angle[idx_minus_16]
clock_err_minus_16 = clock_angle_err[idx_minus_16]
LL_minus_16 = LL_rec[idx_minus_16]
HL_minus_16_pos = HL_rec_pos[idx_minus_16]
HL_minus_16_neg = HL_rec_neg[idx_minus_16]
b_perp_minus_16 = b_perp[idx_minus_16]



# -------- dataframe -----------

# save data out for use in other modules
df_16_p = pd.DataFrame()
df_16_m = pd.DataFrame()

df_16_p = df_16_p.assign(PLUS_Big_DF_Index=idx_plus_16)
df_16_p = df_16_p.assign(PLUS_Ionosphere_Time=iono_time_plus_16)
df_16_p = df_16_p.assign(PLUS_Pressure=pressure_plus_16)
df_16_p = df_16_p.assign(PLUS_Clock_Angle=clock_plus_16)
df_16_p = df_16_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_16)
df_16_p = df_16_p.assign(PLUS_LL_rec_V=LL_plus_16)
df_16_p = df_16_p.assign(PLUS_HL_rec_V_pos=HL_plus_16_pos)
df_16_p = df_16_p.assign(PLUS_HL_rec_V_neg=HL_plus_16_neg)
df_16_p = df_16_p.assign(PLUS_B_Perp=b_perp_plus_16)
df_16_p.to_csv(f'{root_folder}visit_16_times_plus_error_{error}_new.csv',index=False)


df_16_m = df_16_m.assign(MINUS_Big_DF_Index=idx_minus_16)
df_16_m = df_16_m.assign(MINUS_Ionosphere_Time=iono_time_minus_16)
df_16_m = df_16_m.assign(MINUS_Pressure=pressure_minus_16)
df_16_m = df_16_m.assign(MINUS_Clock_Angle=clock_minus_16)
df_16_m = df_16_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_16)
df_16_m = df_16_m.assign(MINUS_LL_rec_V=LL_minus_16)
df_16_m = df_16_m.assign(MINUS_HL_rec_V_pos=HL_minus_16_pos)
df_16_m = df_16_m.assign(MINUS_HL_rec_V_neg=HL_minus_16_neg)
df_16_m = df_16_m.assign(MINUS_B_Perp=b_perp_minus_16)
df_16_m.to_csv(f'{root_folder}visit_16_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 17
'''
ettime_17 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_17)):
    times = spice.str2et(visit_17[p])
    ettime_17.append(times)

# grab first and last timesteps of the visit
first_17 = ettime_17[0]
last_17 = ettime_17[-1]

# run function to determine times and index values of those times for the visit
times_plus_17, idx_plus_17 = between_idxs(plus_error_et,first_17,last_17)
times_minus_17, idx_minus_17 = between_idxs(minus_error_et,first_17,last_17)


# ---------  + % error ------------

iono_time_plus_17 = plus_error[idx_plus_17]
pressure_plus_17 = pressure_array[idx_plus_17]
clock_plus_17 = clock_angle[idx_plus_17]
clock_err_plus_17 = clock_angle_err[idx_plus_17]
LL_plus_17 = LL_rec[idx_plus_17]
HL_plus_17_pos = HL_rec_pos[idx_plus_17] 
HL_plus_17_neg = HL_rec_neg[idx_plus_17]
b_perp_plus_17 = b_perp[idx_plus_17]


# --------- - % error ------------

iono_time_minus_17 = minus_error[idx_minus_17]
pressure_minus_17  = pressure_array[idx_minus_17]
clock_minus_17 = clock_angle[idx_minus_17]
clock_err_minus_17 = clock_angle_err[idx_minus_17]
LL_minus_17 = LL_rec[idx_minus_17]
HL_minus_17_pos = HL_rec_pos[idx_minus_17]
HL_minus_17_neg = HL_rec_neg[idx_minus_17]
b_perp_minus_17 = b_perp[idx_minus_17]



# -------- dataframe -----------

# save data out for use in other modules
df_17_p = pd.DataFrame()
df_17_m = pd.DataFrame()

df_17_p = df_17_p.assign(PLUS_Big_DF_Index=idx_plus_17)
df_17_p = df_17_p.assign(PLUS_Ionosphere_Time=iono_time_plus_17)
df_17_p = df_17_p.assign(PLUS_Pressure=pressure_plus_17)
df_17_p = df_17_p.assign(PLUS_Clock_Angle=clock_plus_17)
df_17_p = df_17_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_17)
df_17_p = df_17_p.assign(PLUS_LL_rec_V=LL_plus_17)
df_17_p = df_17_p.assign(PLUS_HL_rec_V_pos=HL_plus_17_pos)
df_17_p = df_17_p.assign(PLUS_HL_rec_V_neg=HL_plus_17_neg)
df_17_p = df_17_p.assign(PLUS_B_Perp=b_perp_plus_17)
df_17_p.to_csv(f'{root_folder}visit_17_times_plus_error_{error}_new.csv',index=False)


df_17_m = df_17_m.assign(MINUS_Big_DF_Index=idx_minus_17)
df_17_m = df_17_m.assign(MINUS_Ionosphere_Time=iono_time_minus_17)
df_17_m = df_17_m.assign(MINUS_Pressure=pressure_minus_17)
df_17_m = df_17_m.assign(MINUS_Clock_Angle=clock_minus_17)
df_17_m = df_17_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_17)
df_17_m = df_17_m.assign(MINUS_LL_rec_V=LL_minus_17)
df_17_m = df_17_m.assign(MINUS_HL_rec_V_pos=HL_minus_17_pos)
df_17_m = df_17_m.assign(MINUS_HL_rec_V_neg=HL_minus_17_neg)
df_17_m = df_17_m.assign(MINUS_B_Perp=b_perp_minus_17)
df_17_m.to_csv(f'{root_folder}visit_17_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 18
'''
ettime_18 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_18)):
    times = spice.str2et(visit_18[p])
    ettime_18.append(times)

# grab first and last timesteps of the visit
first_18 = ettime_18[0]
last_18 = ettime_18[-1]

# run function to determine times and index values of those times for the visit
times_plus_18, idx_plus_18 = between_idxs(plus_error_et,first_18,last_18)
times_minus_18, idx_minus_18 = between_idxs(minus_error_et,first_18,last_18)


# ---------  + % error ------------

iono_time_plus_18 = plus_error[idx_plus_18]
pressure_plus_18 = pressure_array[idx_plus_18]
clock_plus_18 = clock_angle[idx_plus_18]
clock_err_plus_18 = clock_angle_err[idx_plus_18]
LL_plus_18 = LL_rec[idx_plus_18]
HL_plus_18_pos = HL_rec_pos[idx_plus_18] 
HL_plus_18_neg = HL_rec_neg[idx_plus_18]
b_perp_plus_18 = b_perp[idx_plus_18]


# --------- - % error ------------

iono_time_minus_18 = minus_error[idx_minus_18]
pressure_minus_18  = pressure_array[idx_minus_18]
clock_minus_18 = clock_angle[idx_minus_18]
clock_err_minus_18 = clock_angle_err[idx_minus_18]
LL_minus_18 = LL_rec[idx_minus_18]
HL_minus_18_pos = HL_rec_pos[idx_minus_18]
HL_minus_18_neg = HL_rec_neg[idx_minus_18]
b_perp_minus_18 = b_perp[idx_minus_18]



# -------- dataframe -----------

# save data out for use in other modules
df_18_p = pd.DataFrame()
df_18_m = pd.DataFrame()

df_18_p = df_18_p.assign(PLUS_Big_DF_Index=idx_plus_18)
df_18_p = df_18_p.assign(PLUS_Ionosphere_Time=iono_time_plus_18)
df_18_p = df_18_p.assign(PLUS_Pressure=pressure_plus_18)
df_18_p = df_18_p.assign(PLUS_Clock_Angle=clock_plus_18)
df_18_p = df_18_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_18)
df_18_p = df_18_p.assign(PLUS_LL_rec_V=LL_plus_18)
df_18_p = df_18_p.assign(PLUS_HL_rec_V_pos=HL_plus_18_pos)
df_18_p = df_18_p.assign(PLUS_HL_rec_V_neg=HL_plus_18_neg)
df_18_p = df_18_p.assign(PLUS_B_Perp=b_perp_plus_18)
df_18_p.to_csv(f'{root_folder}visit_18_times_plus_error_{error}_new.csv',index=False)


df_18_m = df_18_m.assign(MINUS_Big_DF_Index=idx_minus_18)
df_18_m = df_18_m.assign(MINUS_Ionosphere_Time=iono_time_minus_18)
df_18_m = df_18_m.assign(MINUS_Pressure=pressure_minus_18)
df_18_m = df_18_m.assign(MINUS_Clock_Angle=clock_minus_18)
df_18_m = df_18_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_18)
df_18_m = df_18_m.assign(MINUS_LL_rec_V=LL_minus_18)
df_18_m = df_18_m.assign(MINUS_HL_rec_V_pos=HL_minus_18_pos)
df_18_m = df_18_m.assign(MINUS_HL_rec_V_neg=HL_minus_18_neg)
df_18_m = df_18_m.assign(MINUS_B_Perp=b_perp_minus_18)
df_18_m.to_csv(f'{root_folder}visit_18_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 19
'''
ettime_19 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_19)):
    times = spice.str2et(visit_19[p])
    ettime_19.append(times)

# grab first and last timesteps of the visit
first_19 = ettime_19[0]
last_19 = ettime_19[-1]

# run function to determine times and index values of those times for the visit
times_plus_19, idx_plus_19 = between_idxs(plus_error_et,first_19,last_19)
times_minus_19, idx_minus_19 = between_idxs(minus_error_et,first_19,last_19)


# ---------  + % error ------------

iono_time_plus_19 = plus_error[idx_plus_19]
pressure_plus_19 = pressure_array[idx_plus_19]
clock_plus_19 = clock_angle[idx_plus_19]
clock_err_plus_19 = clock_angle_err[idx_plus_19]
LL_plus_19 = LL_rec[idx_plus_19]
HL_plus_19_pos = HL_rec_pos[idx_plus_19] 
HL_plus_19_neg = HL_rec_neg[idx_plus_19]
b_perp_plus_19 = b_perp[idx_plus_19]


# --------- - % error ------------

iono_time_minus_19 = minus_error[idx_minus_19]
pressure_minus_19  = pressure_array[idx_minus_19]
clock_minus_19 = clock_angle[idx_minus_19]
clock_err_minus_19 = clock_angle_err[idx_minus_19]
LL_minus_19 = LL_rec[idx_minus_19]
HL_minus_19_pos = HL_rec_pos[idx_minus_19]
HL_minus_19_neg = HL_rec_neg[idx_minus_19]
b_perp_minus_19 = b_perp[idx_minus_19]



# -------- dataframe -----------

# save data out for use in other modules
df_19_p = pd.DataFrame()
df_19_m = pd.DataFrame()

df_19_p = df_19_p.assign(PLUS_Big_DF_Index=idx_plus_19)
df_19_p = df_19_p.assign(PLUS_Ionosphere_Time=iono_time_plus_19)
df_19_p = df_19_p.assign(PLUS_Pressure=pressure_plus_19)
df_19_p = df_19_p.assign(PLUS_Clock_Angle=clock_plus_19)
df_19_p = df_19_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_19)
df_19_p = df_19_p.assign(PLUS_LL_rec_V=LL_plus_19)
df_19_p = df_19_p.assign(PLUS_HL_rec_V_pos=HL_plus_19_pos)
df_19_p = df_19_p.assign(PLUS_HL_rec_V_neg=HL_plus_19_neg)
df_19_p = df_19_p.assign(PLUS_B_Perp=b_perp_plus_19)
df_19_p.to_csv(f'{root_folder}visit_19_times_plus_error_{error}_new.csv',index=False)


df_19_m = df_19_m.assign(MINUS_Big_DF_Index=idx_minus_19)
df_19_m = df_19_m.assign(MINUS_Ionosphere_Time=iono_time_minus_19)
df_19_m = df_19_m.assign(MINUS_Pressure=pressure_minus_19)
df_19_m = df_19_m.assign(MINUS_Clock_Angle=clock_minus_19)
df_19_m = df_19_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_19)
df_19_m = df_19_m.assign(MINUS_LL_rec_V=LL_minus_19)
df_19_m = df_19_m.assign(MINUS_HL_rec_V_pos=HL_minus_19_pos)
df_19_m = df_19_m.assign(MINUS_HL_rec_V_neg=HL_minus_19_neg)
df_19_m = df_19_m.assign(MINUS_B_Perp=b_perp_minus_19)
df_19_m.to_csv(f'{root_folder}visit_19_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 20
'''
ettime_20 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_20)):
    times = spice.str2et(visit_20[p])
    ettime_20.append(times)

# grab first and last timesteps of the visit
first_20 = ettime_20[0]
last_20 = ettime_20[-1]

# run function to determine times and index values of those times for the visit
times_plus_20, idx_plus_20 = between_idxs(plus_error_et,first_20,last_20)
times_minus_20, idx_minus_20 = between_idxs(minus_error_et,first_20,last_20)


# ---------  + % error ------------

iono_time_plus_20 = plus_error[idx_plus_20]
pressure_plus_20 = pressure_array[idx_plus_20]
clock_plus_20 = clock_angle[idx_plus_20]
clock_err_plus_20 = clock_angle_err[idx_plus_20]
LL_plus_20 = LL_rec[idx_plus_20]
HL_plus_20_pos = HL_rec_pos[idx_plus_20] 
HL_plus_20_neg = HL_rec_neg[idx_plus_20]
b_perp_plus_20 = b_perp[idx_plus_20]


# --------- - % error ------------

iono_time_minus_20 = minus_error[idx_minus_20]
pressure_minus_20 = pressure_array[idx_minus_20]
clock_minus_20 = clock_angle[idx_minus_20]
clock_err_minus_20 = clock_angle_err[idx_minus_20]
LL_minus_20 = LL_rec[idx_minus_20]
HL_minus_20_pos = HL_rec_pos[idx_minus_20]
HL_minus_20_neg = HL_rec_neg[idx_minus_20]
b_perp_minus_20 = b_perp[idx_minus_20]



# -------- dataframe -----------

# save data out for use in other modules
df_20_p = pd.DataFrame()
df_20_m = pd.DataFrame()

df_20_p = df_20_p.assign(PLUS_Big_DF_Index=idx_plus_20)
df_20_p = df_20_p.assign(PLUS_Ionosphere_Time=iono_time_plus_20)
df_20_p = df_20_p.assign(PLUS_Pressure=pressure_plus_20)
df_20_p = df_20_p.assign(PLUS_Clock_Angle=clock_plus_20)
df_20_p = df_20_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_20)
df_20_p = df_20_p.assign(PLUS_LL_rec_V=LL_plus_20)
df_20_p = df_20_p.assign(PLUS_HL_rec_V_pos=HL_plus_20_pos)
df_20_p = df_20_p.assign(PLUS_HL_rec_V_neg=HL_plus_20_neg)
df_20_p = df_20_p.assign(PLUS_B_Perp=b_perp_plus_20)
df_20_p.to_csv(f'{root_folder}visit_20_times_plus_error_{error}_new.csv',index=False)


df_20_m = df_20_m.assign(MINUS_Big_DF_Index=idx_minus_20)
df_20_m = df_20_m.assign(MINUS_Ionosphere_Time=iono_time_minus_20)
df_20_m = df_20_m.assign(MINUS_Pressure=pressure_minus_20)
df_20_m = df_20_m.assign(MINUS_Clock_Angle=clock_minus_20)
df_20_m = df_20_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_20)
df_20_m = df_20_m.assign(MINUS_LL_rec_V=LL_minus_20)
df_20_m = df_20_m.assign(MINUS_HL_rec_V_pos=HL_minus_20_pos)
df_20_m = df_20_m.assign(MINUS_HL_rec_V_neg=HL_minus_20_neg)
df_20_m = df_20_m.assign(MINUS_B_Perp=b_perp_minus_20)
df_20_m.to_csv(f'{root_folder}visit_20_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 21
'''
ettime_21 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_21)):
    times = spice.str2et(visit_21[p])
    ettime_21.append(times)

# grab first and last timesteps of the visit
first_21 = ettime_21[0]
last_21 = ettime_21[-1]

# run function to determine times and index values of those times for the visit
times_plus_21, idx_plus_21 = between_idxs(plus_error_et,first_21,last_21)
times_minus_21, idx_minus_21 = between_idxs(minus_error_et,first_21,last_21)


# ---------  + % error ------------

iono_time_plus_21 = plus_error[idx_plus_21]
pressure_plus_21 = pressure_array[idx_plus_21]
clock_plus_21 = clock_angle[idx_plus_21]
clock_err_plus_21 = clock_angle_err[idx_plus_21]
LL_plus_21 = LL_rec[idx_plus_21]
HL_plus_21_pos = HL_rec_pos[idx_plus_21] 
HL_plus_21_neg = HL_rec_neg[idx_plus_21]
b_perp_plus_21 = b_perp[idx_plus_21]


# --------- - % error ------------

iono_time_minus_21 = minus_error[idx_minus_21]
pressure_minus_21 = pressure_array[idx_minus_21]
clock_minus_21 = clock_angle[idx_minus_21]
clock_err_minus_21 = clock_angle_err[idx_minus_21]
LL_minus_21 = LL_rec[idx_minus_21]
HL_minus_21_pos = HL_rec_pos[idx_minus_21]
HL_minus_21_neg = HL_rec_neg[idx_minus_21]
b_perp_minus_21 = b_perp[idx_minus_21]


# -------- dataframe -----------

# save data out for use in other modules
df_21_p = pd.DataFrame()
df_21_m = pd.DataFrame()

df_21_p = df_21_p.assign(PLUS_Big_DF_Index=idx_plus_21)
df_21_p = df_21_p.assign(PLUS_Ionosphere_Time=iono_time_plus_21)
df_21_p = df_21_p.assign(PLUS_Pressure=pressure_plus_21)
df_21_p = df_21_p.assign(PLUS_Clock_Angle=clock_plus_21)
df_21_p = df_21_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_21)
df_21_p = df_21_p.assign(PLUS_LL_rec_V=LL_plus_21)
df_21_p = df_21_p.assign(PLUS_HL_rec_V_pos=HL_plus_21_pos)
df_21_p = df_21_p.assign(PLUS_HL_rec_V_neg=HL_plus_21_neg)
df_21_p = df_21_p.assign(PLUS_B_Perp=b_perp_plus_21)
df_21_p.to_csv(f'{root_folder}visit_21_times_plus_error_{error}_new.csv',index=False)


df_21_m = df_21_m.assign(MINUS_Big_DF_Index=idx_minus_21)
df_21_m = df_21_m.assign(MINUS_Ionosphere_Time=iono_time_minus_21)
df_21_m = df_21_m.assign(MINUS_Pressure=pressure_minus_21)
df_21_m = df_21_m.assign(MINUS_Clock_Angle=clock_minus_21)
df_21_m = df_21_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_21)
df_21_m = df_21_m.assign(MINUS_LL_rec_V=LL_minus_21)
df_21_m = df_21_m.assign(MINUS_HL_rec_V_pos=HL_minus_21_pos)
df_21_m = df_21_m.assign(MINUS_HL_rec_V_neg=HL_minus_21_neg)
df_21_m = df_21_m.assign(MINUS_B_Perp=b_perp_minus_21)
df_21_m.to_csv(f'{root_folder}visit_21_times_minus_error_{error}_new.csv',index=False)



'''
VISIT 24
'''
ettime_24 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_24)):
    times = spice.str2et(visit_24[p])
    ettime_24.append(times)

# grab first and last timesteps of the visit
first_24 = ettime_24[0]
last_24 = ettime_24[-1]

# run function to determine times and index values of those times for the visit
times_plus_24, idx_plus_24 = between_idxs(plus_error_et,first_24,last_24)
times_minus_24, idx_minus_24 = between_idxs(minus_error_et,first_24,last_24)


# ---------  + % error ------------

iono_time_plus_24 = plus_error[idx_plus_24]
pressure_plus_24 = pressure_array[idx_plus_24]
clock_plus_24 = clock_angle[idx_plus_24]
clock_err_plus_24 = clock_angle_err[idx_plus_24]
LL_plus_24 = LL_rec[idx_plus_24]
HL_plus_24_pos = HL_rec_pos[idx_plus_24] 
HL_plus_24_neg = HL_rec_neg[idx_plus_24]
b_perp_plus_24 = b_perp[idx_plus_24]



# --------- - % error ------------

iono_time_minus_24 = minus_error[idx_minus_24]
pressure_minus_24 = pressure_array[idx_minus_24]
clock_minus_24 = clock_angle[idx_minus_24]
clock_err_minus_24 = clock_angle_err[idx_minus_24]
LL_minus_24 = LL_rec[idx_minus_24]
HL_minus_24_pos = HL_rec_pos[idx_minus_24]
HL_minus_24_neg = HL_rec_neg[idx_minus_24]
b_perp_minus_24 = b_perp[idx_minus_24]


# -------- dataframe -----------

# save data out for use in other modules
df_24_p = pd.DataFrame()
df_24_m = pd.DataFrame()

df_24_p = df_24_p.assign(PLUS_Big_DF_Index=idx_plus_24)
df_24_p = df_24_p.assign(PLUS_Ionosphere_Time=iono_time_plus_24)
df_24_p = df_24_p.assign(PLUS_Pressure=pressure_plus_24)
df_24_p = df_24_p.assign(PLUS_Clock_Angle=clock_plus_24)
df_24_p = df_24_p.assign(PLUS_Clock_Angle_Error=clock_err_plus_24)
df_24_p = df_24_p.assign(PLUS_LL_rec_V=LL_plus_24)
df_24_p = df_24_p.assign(PLUS_HL_rec_V_pos=HL_plus_24_pos)
df_24_p = df_24_p.assign(PLUS_HL_rec_V_neg=HL_plus_24_neg)
df_24_p = df_24_p.assign(PLUS_B_Perp=b_perp_plus_24)
df_24_p.to_csv(f'{root_folder}visit_24_times_plus_error_{error}_new.csv',index=False)


df_24_m = df_24_m.assign(MINUS_Big_DF_Index=idx_minus_24)
df_24_m = df_24_m.assign(MINUS_Ionosphere_Time=iono_time_minus_24)
df_24_m = df_24_m.assign(MINUS_Pressure=pressure_minus_24)
df_24_m = df_24_m.assign(MINUS_Clock_Angle=clock_minus_24)
df_24_m = df_24_m.assign(MINUS_Clock_Angle_Error=clock_err_minus_24)
df_24_m = df_24_m.assign(MINUS_LL_rec_V=LL_minus_24)
df_24_m = df_24_m.assign(MINUS_HL_rec_V_pos=HL_minus_24_pos)
df_24_m = df_24_m.assign(MINUS_HL_rec_V_neg=HL_minus_24_neg)
df_24_m = df_24_m.assign(MINUS_B_Perp=b_perp_minus_24)
df_24_m.to_csv(f'{root_folder}visit_24_times_minus_error_{error}_new.csv',index=False)

# --------------- doy 155 -----------------

'''
VISIT 25
'''
ettime_25 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_25)):
    times = spice.str2et(visit_25[p])
    ettime_25.append(times)

# grab first and last timesteps of the visit
first_25 = ettime_25[0]
last_25 = ettime_25[-1]

# run function to determine times and index values of those times for the visit
times_plus_25, idx_plus_25 = between_idxs(plus_error_et_155,first_25,last_25)
times_minus_25, idx_minus_25 = between_idxs(minus_error_et_155,first_25,last_25)


# ---------  + % error ------------

iono_time_plus_25 = plus_error_155[idx_plus_25]
pressure_plus_25 = pressure_array_155[idx_plus_25]
clock_plus_25 = clock_angle_155[idx_plus_25]
LL_plus_25 = LL_rec_155[idx_plus_25]
HL_plus_25_pos = HL_rec_pos_155[idx_plus_25] 
HL_plus_25_neg = HL_rec_neg_155[idx_plus_25]
b_perp_plus_25 = b_perp_155[idx_plus_25]


# --------- - % error ------------

iono_time_minus_25 = minus_error_155[idx_minus_25]
pressure_minus_25 = pressure_array_155[idx_minus_25]
clock_minus_25 = clock_angle_155[idx_minus_25]
LL_minus_25 = LL_rec_155[idx_minus_25]
HL_minus_25_pos = HL_rec_pos_155[idx_minus_25]
HL_minus_25_neg = HL_rec_neg_155[idx_minus_25]
b_perp_minus_25 = b_perp_155[idx_minus_25]


# -------- dataframe -----------

# save data out for use in other modules
df_25_p = pd.DataFrame()
df_25_m = pd.DataFrame()

df_25_p = df_25_p.assign(PLUS_Big_DF_Index=idx_plus_25)
df_25_p = df_25_p.assign(PLUS_Ionosphere_Time=iono_time_plus_25)
df_25_p = df_25_p.assign(PLUS_Pressure=pressure_plus_25)
df_25_p = df_25_p.assign(PLUS_Clock_Angle=clock_plus_25)
df_25_p = df_25_p.assign(PLUS_LL_rec_V=LL_plus_25)
df_25_p = df_25_p.assign(PLUS_HL_rec_V_pos=HL_plus_25_pos)
df_25_p = df_25_p.assign(PLUS_HL_rec_V_neg=HL_plus_25_neg)
df_25_p = df_25_p.assign(PLUS_B_Perp=b_perp_plus_25)
df_25_p.to_csv(f'{root_folder}visit_25_times_plus_error_{error}_new.csv',index=False)


df_25_m = df_25_m.assign(MINUS_Big_DF_Index=idx_minus_25)
df_25_m = df_25_m.assign(MINUS_Ionosphere_Time=iono_time_minus_25)
df_25_m = df_25_m.assign(MINUS_Pressure=pressure_minus_25)
df_25_m = df_25_m.assign(MINUS_Clock_Angle=clock_minus_25)
df_25_m = df_25_m.assign(MINUS_LL_rec_V=LL_minus_25)
df_25_m = df_25_m.assign(MINUS_HL_rec_V_pos=HL_minus_25_pos)
df_25_m = df_25_m.assign(MINUS_HL_rec_V_neg=HL_minus_25_neg)
df_25_m = df_25_m.assign(MINUS_B_Perp=b_perp_minus_25)
df_25_m.to_csv(f'{root_folder}visit_25_times_minus_error_{error}_new.csv',index=False)


# ------------ doys 156 + -----------

'''
VISIT 26
'''
ettime_26 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_26)):
    times = spice.str2et(visit_26[p])
    ettime_26.append(times)

# grab first and last timesteps of the visit
first_26 = ettime_26[0]
last_26 = ettime_26[-1]

# run function to determine times and index values of those times for the visit
times_plus_26, idx_plus_26 = between_idxs(plus_error_et_156_plus,first_26,last_26)
times_minus_26, idx_minus_26 = between_idxs(minus_error_et_156_plus,first_26,last_26)


# ---------  + % error ------------

iono_time_plus_26 = plus_error_156_plus[idx_plus_26]
pressure_plus_26 = pressure_array_156_plus[idx_plus_26]
clock_plus_26 = clock_angle_156_plus[idx_plus_26]
LL_plus_26 = LL_rec_156_plus[idx_plus_26]
HL_plus_26_pos = HL_rec_pos_156_plus[idx_plus_26] 
HL_plus_26_neg = HL_rec_neg_156_plus[idx_plus_26]
b_perp_plus_26 = b_perp_156_plus[idx_plus_26]


# --------- - % error ------------

iono_time_minus_26 = minus_error_156_plus[idx_minus_26]
pressure_minus_26 = pressure_array_156_plus[idx_minus_26]
clock_minus_26 = clock_angle_156_plus[idx_minus_26]
LL_minus_26 = LL_rec_156_plus[idx_minus_26]
HL_minus_26_pos = HL_rec_pos_156_plus[idx_minus_26]
HL_minus_26_neg = HL_rec_neg_156_plus[idx_minus_26]
b_perp_minus_26 = b_perp_156_plus[idx_minus_26]


# -------- dataframe -----------

# save data out for use in other modules
df_26_p = pd.DataFrame()
df_26_m = pd.DataFrame()

df_26_p = df_26_p.assign(PLUS_Big_DF_Index=idx_plus_26)
df_26_p = df_26_p.assign(PLUS_Ionosphere_Time=iono_time_plus_26)
df_26_p = df_26_p.assign(PLUS_Pressure=pressure_plus_26)
df_26_p = df_26_p.assign(PLUS_Clock_Angle=clock_plus_26)
df_26_p = df_26_p.assign(PLUS_LL_rec_V=LL_plus_26)
df_26_p = df_26_p.assign(PLUS_HL_rec_V_pos=HL_plus_26_pos)
df_26_p = df_26_p.assign(PLUS_HL_rec_V_neg=HL_plus_26_neg)
df_26_p = df_26_p.assign(PLUS_B_Perp=b_perp_plus_26)
df_26_p.to_csv(f'{root_folder}visit_26_times_plus_error_{error}_new.csv',index=False)


df_26_m = df_26_m.assign(MINUS_Big_DF_Index=idx_minus_26)
df_26_m = df_26_m.assign(MINUS_Ionosphere_Time=iono_time_minus_26)
df_26_m = df_26_m.assign(MINUS_Pressure=pressure_minus_26)
df_26_m = df_26_m.assign(MINUS_Clock_Angle=clock_minus_26)
df_26_m = df_26_m.assign(MINUS_LL_rec_V=LL_minus_26)
df_26_m = df_26_m.assign(MINUS_HL_rec_V_pos=HL_minus_26_pos)
df_26_m = df_26_m.assign(MINUS_HL_rec_V_neg=HL_minus_26_neg)
df_26_m = df_26_m.assign(MINUS_B_Perp=b_perp_minus_26)
df_26_m.to_csv(f'{root_folder}visit_26_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 27
'''
ettime_27 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_27)):
    times = spice.str2et(visit_27[p])
    ettime_27.append(times)

# grab first and last timesteps of the visit
first_27 = ettime_27[0]
last_27 = ettime_27[-1]

# run function to determine times and index values of those times for the visit
times_plus_27, idx_plus_27 = between_idxs(plus_error_et_156_plus,first_27,last_27)
times_minus_27, idx_minus_27 = between_idxs(minus_error_et_156_plus,first_27,last_27)


# ---------  + % error ------------

iono_time_plus_27 = plus_error_156_plus[idx_plus_27]
pressure_plus_27 = pressure_array_156_plus[idx_plus_27]
clock_plus_27 = clock_angle_156_plus[idx_plus_27]
LL_plus_27 = LL_rec_156_plus[idx_plus_27]
HL_plus_27_pos = HL_rec_pos_156_plus[idx_plus_27] 
HL_plus_27_neg = HL_rec_neg_156_plus[idx_plus_27]
b_perp_plus_27 = b_perp_156_plus[idx_plus_27]


# --------- - % error ------------

iono_time_minus_27 = minus_error_156_plus[idx_minus_27]
pressure_minus_27 = pressure_array_156_plus[idx_minus_27]
clock_minus_27 = clock_angle_156_plus[idx_minus_27]
LL_minus_27 = LL_rec_156_plus[idx_minus_27]
HL_minus_27_pos = HL_rec_pos_156_plus[idx_minus_27]
HL_minus_27_neg = HL_rec_neg_156_plus[idx_minus_27]
b_perp_minus_27 = b_perp_156_plus[idx_minus_27]


# -------- dataframe -----------

# save data out for use in other modules
df_27_p = pd.DataFrame()
df_27_m = pd.DataFrame()

df_27_p = df_27_p.assign(PLUS_Big_DF_Index=idx_plus_27)
df_27_p = df_27_p.assign(PLUS_Ionosphere_Time=iono_time_plus_27)
df_27_p = df_27_p.assign(PLUS_Pressure=pressure_plus_27)
df_27_p = df_27_p.assign(PLUS_Clock_Angle=clock_plus_27)
df_27_p = df_27_p.assign(PLUS_LL_rec_V=LL_plus_27)
df_27_p = df_27_p.assign(PLUS_HL_rec_V_pos=HL_plus_27_pos)
df_27_p = df_27_p.assign(PLUS_HL_rec_V_neg=HL_plus_27_neg)
df_27_p = df_27_p.assign(PLUS_B_Perp=b_perp_plus_27)
df_27_p.to_csv(f'{root_folder}visit_27_times_plus_error_{error}_new.csv',index=False)


df_27_m = df_27_m.assign(MINUS_Big_DF_Index=idx_minus_27)
df_27_m = df_27_m.assign(MINUS_Ionosphere_Time=iono_time_minus_27)
df_27_m = df_27_m.assign(MINUS_Pressure=pressure_minus_27)
df_27_m = df_27_m.assign(MINUS_Clock_Angle=clock_minus_27)
df_27_m = df_27_m.assign(MINUS_LL_rec_V=LL_minus_27)
df_27_m = df_27_m.assign(MINUS_HL_rec_V_pos=HL_minus_27_pos)
df_27_m = df_27_m.assign(MINUS_HL_rec_V_neg=HL_minus_27_neg)
df_27_m = df_27_m.assign(MINUS_B_Perp=b_perp_minus_27)
df_27_m.to_csv(f'{root_folder}visit_27_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 28
'''
ettime_28 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_28)):
    times = spice.str2et(visit_28[p])
    ettime_28.append(times)

# grab first and last timesteps of the visit
first_28 = ettime_28[0]
last_28 = ettime_28[-1]

# run function to determine times and index values of those times for the visit
times_plus_28, idx_plus_28 = between_idxs(plus_error_et_156_plus,first_28,last_28)
times_minus_28, idx_minus_28 = between_idxs(minus_error_et_156_plus,first_28,last_28)


# ---------  + % error ------------

iono_time_plus_28 = plus_error_156_plus[idx_plus_28]
pressure_plus_28 = pressure_array_156_plus[idx_plus_28]
clock_plus_28 = clock_angle_156_plus[idx_plus_28]
LL_plus_28 = LL_rec_156_plus[idx_plus_28]
HL_plus_28_pos = HL_rec_pos_156_plus[idx_plus_28] 
HL_plus_28_neg = HL_rec_neg_156_plus[idx_plus_28]
b_perp_plus_28 = b_perp_156_plus[idx_plus_28]


# --------- - % error ------------

iono_time_minus_28 = minus_error_156_plus[idx_minus_28]
pressure_minus_28 = pressure_array_156_plus[idx_minus_28]
clock_minus_28 = clock_angle_156_plus[idx_minus_28]
LL_minus_28 = LL_rec_156_plus[idx_minus_28]
HL_minus_28_pos = HL_rec_pos_156_plus[idx_minus_28]
HL_minus_28_neg = HL_rec_neg_156_plus[idx_minus_28]
b_perp_minus_28 = b_perp_156_plus[idx_minus_28]


# -------- dataframe -----------

# save data out for use in other modules
df_28_p = pd.DataFrame()
df_28_m = pd.DataFrame()

df_28_p = df_28_p.assign(PLUS_Big_DF_Index=idx_plus_28)
df_28_p = df_28_p.assign(PLUS_Ionosphere_Time=iono_time_plus_28)
df_28_p = df_28_p.assign(PLUS_Pressure=pressure_plus_28)
df_28_p = df_28_p.assign(PLUS_Clock_Angle=clock_plus_28)
df_28_p = df_28_p.assign(PLUS_LL_rec_V=LL_plus_28)
df_28_p = df_28_p.assign(PLUS_HL_rec_V_pos=HL_plus_28_pos)
df_28_p = df_28_p.assign(PLUS_HL_rec_V_neg=HL_plus_28_neg)
df_28_p = df_28_p.assign(PLUS_B_Perp=b_perp_plus_28)
df_28_p.to_csv(f'{root_folder}visit_28_times_plus_error_{error}_new.csv',index=False)


df_28_m = df_28_m.assign(MINUS_Big_DF_Index=idx_minus_28)
df_28_m = df_28_m.assign(MINUS_Ionosphere_Time=iono_time_minus_28)
df_28_m = df_28_m.assign(MINUS_Pressure=pressure_minus_28)
df_28_m = df_28_m.assign(MINUS_Clock_Angle=clock_minus_28)
df_28_m = df_28_m.assign(MINUS_LL_rec_V=LL_minus_28)
df_28_m = df_28_m.assign(MINUS_HL_rec_V_pos=HL_minus_28_pos)
df_28_m = df_28_m.assign(MINUS_HL_rec_V_neg=HL_minus_28_neg)
df_28_m = df_28_m.assign(MINUS_B_Perp=b_perp_minus_28)
df_28_m.to_csv(f'{root_folder}visit_28_times_minus_error_{error}_new.csv',index=False)


# ------------- doy 175 and 176 ---------

'''
VISIT 34
'''
ettime_34 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_34)):
    times = spice.str2et(visit_34[p])
    ettime_34.append(times)

# grab first and last timesteps of the visit
first_34 = ettime_34[0]
last_34 = ettime_34[-1]

# run function to determine times and index values of those times for the visit
times_plus_34, idx_plus_34 = between_idxs(plus_error_et_175,first_34,last_34)
times_minus_34, idx_minus_34 = between_idxs(minus_error_et_175,first_34,last_34)


# ---------  + % error ------------

iono_time_plus_34 = plus_error_175[idx_plus_34]
pressure_plus_34 = pressure_array_175[idx_plus_34]
clock_plus_34 = clock_angle_175[idx_plus_34]
LL_plus_34 = LL_rec_175[idx_plus_34]
HL_plus_34_pos = HL_rec_pos_175[idx_plus_34] 
HL_plus_34_neg = HL_rec_neg_175[idx_plus_34]
b_perp_plus_34 = b_perp_175[idx_plus_34]


# --------- - % error ------------

iono_time_minus_34 = minus_error_175[idx_minus_34]
pressure_minus_34 = pressure_array_175[idx_minus_34]
clock_minus_34 = clock_angle_175[idx_minus_34]
LL_minus_34 = LL_rec_175[idx_minus_34]
HL_minus_34_pos = HL_rec_pos_175[idx_minus_34]
HL_minus_34_neg = HL_rec_neg_175[idx_minus_34]
b_perp_minus_34 = b_perp_175[idx_minus_34]


# -------- dataframe -----------

# save data out for use in other modules
df_34_p = pd.DataFrame()
df_34_m = pd.DataFrame()

df_34_p = df_34_p.assign(PLUS_Big_DF_Index=idx_plus_34)
df_34_p = df_34_p.assign(PLUS_Ionosphere_Time=iono_time_plus_34)
df_34_p = df_34_p.assign(PLUS_Pressure=pressure_plus_34)
df_34_p = df_34_p.assign(PLUS_Clock_Angle=clock_plus_34)
df_34_p = df_34_p.assign(PLUS_LL_rec_V=LL_plus_34)
df_34_p = df_34_p.assign(PLUS_HL_rec_V_pos=HL_plus_34_pos)
df_34_p = df_34_p.assign(PLUS_HL_rec_V_neg=HL_plus_34_neg)
df_34_p = df_34_p.assign(PLUS_B_Perp=b_perp_plus_34)
df_34_p.to_csv(f'{root_folder}visit_34_times_plus_error_{error}_new.csv',index=False)


df_34_m = df_34_m.assign(MINUS_Big_DF_Index=idx_minus_34)
df_34_m = df_34_m.assign(MINUS_Ionosphere_Time=iono_time_minus_34)
df_34_m = df_34_m.assign(MINUS_Pressure=pressure_minus_34)
df_34_m = df_34_m.assign(MINUS_Clock_Angle=clock_minus_34)
df_34_m = df_34_m.assign(MINUS_LL_rec_V=LL_minus_34)
df_34_m = df_34_m.assign(MINUS_HL_rec_V_pos=HL_minus_34_pos)
df_34_m = df_34_m.assign(MINUS_HL_rec_V_neg=HL_minus_34_neg)
df_34_m = df_34_m.assign(MINUS_B_Perp=b_perp_minus_34)
df_34_m.to_csv(f'{root_folder}visit_34_times_minus_error_{error}_new.csv',index=False)


'''
VISIT 35
'''
ettime_35 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_35)):
    times = spice.str2et(visit_35[p])
    ettime_35.append(times)

# grab first and last timesteps of the visit
first_35 = ettime_35[0]
last_35 = ettime_35[-1]

# run function to determine times and index values of those times for the visit
times_plus_35, idx_plus_35 = between_idxs(plus_error_et_175,first_35,last_35)
times_minus_35, idx_minus_35 = between_idxs(minus_error_et_175,first_35,last_35)


# ---------  + % error ------------

iono_time_plus_35 = plus_error_175[idx_plus_35]
pressure_plus_35 = pressure_array_175[idx_plus_35]
clock_plus_35 = clock_angle_175[idx_plus_35]
LL_plus_35 = LL_rec_175[idx_plus_35]
HL_plus_35_pos = HL_rec_pos_175[idx_plus_35] 
HL_plus_35_neg = HL_rec_neg_175[idx_plus_35]
b_perp_plus_35 = b_perp_175[idx_plus_35]


# --------- - % error ------------

iono_time_minus_35 = minus_error_175[idx_minus_35]
pressure_minus_35 = pressure_array_175[idx_minus_35]
clock_minus_35 = clock_angle_175[idx_minus_35]
LL_minus_35 = LL_rec_175[idx_minus_35]
HL_minus_35_pos = HL_rec_pos_175[idx_minus_35]
HL_minus_35_neg = HL_rec_neg_175[idx_minus_35]
b_perp_minus_35 = b_perp_175[idx_minus_35]


# -------- dataframe -----------

# save data out for use in other modules
df_35_p = pd.DataFrame()
df_35_m = pd.DataFrame()

df_35_p = df_35_p.assign(PLUS_Big_DF_Index=idx_plus_35)
df_35_p = df_35_p.assign(PLUS_Ionosphere_Time=iono_time_plus_35)
df_35_p = df_35_p.assign(PLUS_Pressure=pressure_plus_35)
df_35_p = df_35_p.assign(PLUS_Clock_Angle=clock_plus_35)
df_35_p = df_35_p.assign(PLUS_LL_rec_V=LL_plus_35)
df_35_p = df_35_p.assign(PLUS_HL_rec_V_pos=HL_plus_35_pos)
df_35_p = df_35_p.assign(PLUS_HL_rec_V_neg=HL_plus_35_neg)
df_35_p = df_35_p.assign(PLUS_B_Perp=b_perp_plus_35)
df_35_p.to_csv(f'{root_folder}visit_35_times_plus_error_{error}_new.csv',index=False)


df_35_m = df_35_m.assign(MINUS_Big_DF_Index=idx_minus_35)
df_35_m = df_35_m.assign(MINUS_Ionosphere_Time=iono_time_minus_35)
df_35_m = df_35_m.assign(MINUS_Pressure=pressure_minus_35)
df_35_m = df_35_m.assign(MINUS_Clock_Angle=clock_minus_35)
df_35_m = df_35_m.assign(MINUS_LL_rec_V=LL_minus_35)
df_35_m = df_35_m.assign(MINUS_HL_rec_V_pos=HL_minus_35_pos)
df_35_m = df_35_m.assign(MINUS_HL_rec_V_neg=HL_minus_35_neg)
df_35_m = df_35_m.assign(MINUS_B_Perp=b_perp_minus_35)
df_35_m.to_csv(f'{root_folder}visit_35_times_minus_error_{error}_new.csv',index=False)
