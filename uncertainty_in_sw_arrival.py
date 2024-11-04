"""
Created on Mon Sep  2 13:26:22 2024

@author: hannah
"""

# load in relevant modules
import pandas as pd
import spiceypy as spice
import numpy as np
import datetime as dt

# leap seconds kernal - need this for ephemerous time
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

#  ------------------- dataframes -----------------
visit_times = pd.read_csv(root_folder+'visit_times.csv')
df = pd.read_csv(root_folder+'juno_data_big_df_realtime_new.csv')
df_155 = pd.read_csv(root_folder+'juno_data_doy_155.csv')
df_156_plus = pd.read_csv(root_folder+'juno_data_doy_156_plus.csv')
df_175 = pd.read_csv(root_folder+'juno_data_doy_175.csv')

travel_times_main = df['Total_Travel_Time'].to_numpy()
travel_times_155 = df_155['Total_Travel_Time'].to_numpy()
travel_times_156_plus = df_156_plus['Total_Travel_Time'].to_numpy()
travel_times_175 = df_175['Total_Travel_Time'].to_numpy()


# -------------- error calculation -------------

error='10'
error_div = int(error)/100


times = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(df['Juno_Detection_Time'])):
    timetest = dt.datetime.strptime(df['Juno_Detection_Time'][p],'%Y-%m-%d %H:%M:%S.%f')
    times.append(timetest)
    
times_array = np.array(times)

plus_tt = []
minus_tt = []
errors = []
for i in range(len(travel_times_main)):
    calc = travel_times_main[i] * error_div
    plus = travel_times_main[i] + calc
    minus = travel_times_main[i] - calc
    
    errors.append(calc)
    plus_tt.append(plus)
    minus_tt.append(minus)

plus_times  = []
minus_times = []
for q in range(len(travel_times_main)):
    plus_calc = times[q] + dt.timedelta(seconds=plus_tt[q])
    minus_calc = times[q] + dt.timedelta(seconds=minus_tt[q])
    plus_times.append(plus_calc)
    minus_times.append(minus_calc)
    
df_extra = pd.DataFrame()

df_extra = df_extra.assign(Error_in_Travel_Time=errors)
df_extra = df_extra.assign(Ionosphere_Times_Plus_Error=plus_times)
df_extra = df_extra.assign(Ionosphere_Times_Minus_Error=minus_times)

df_new = pd.concat([df, df_extra],axis=1)


df_new.to_csv(f'{root_folder}juno_data_big_df_realtime_new_errors_{error}.csv',index=False)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
df_iono =  df_new

df_iono.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
df_iono.to_csv(f'{root_folder}juno_data_big_df_ionotime_new_errors_{error}.csv',index=False)


# -------- visit 25 -------------

times_155 = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(df_155['UTC'])):
    timetest = dt.datetime.strptime(df_155['UTC'][p],'%Y-%m-%d %H:%M:%S.%f')
    times_155.append(timetest)
    
times_155 = np.array(times_155)

plus_tt_155 = []
minus_tt_155 = []
errors_155 = []
for i in range(len(travel_times_155)):
    calc_155 = travel_times_155[i] * error_div
    plus_155 = travel_times_155[i] + calc_155
    minus_155 = travel_times_155[i] - calc_155
    
    errors_155.append(calc_155)
    plus_tt_155.append(plus_155)
    minus_tt_155.append(minus_155)

    
plus_times_155  = []
minus_times_155 = []
for q in range(len(travel_times_155)):
    plus_calc_155 = times_155[q] + dt.timedelta(seconds=plus_tt_155[q])
    minus_calc_155 = times_155[q] + dt.timedelta(seconds=minus_tt_155[q])
    plus_times_155.append(plus_calc_155)
    minus_times_155.append(minus_calc_155)
    

df_extra_155 = pd.DataFrame()

df_extra_155 = df_extra_155.assign(Error_in_Travel_Time=errors_155)
df_extra_155 = df_extra_155.assign(Ionosphere_Times_Plus_Error=plus_times_155)
df_extra_155 = df_extra_155.assign(Ionosphere_Times_Minus_Error=minus_times_155)


df_new_155 = pd.concat([df_155, df_extra_155],axis=1)


df_new_155.to_csv(f'{root_folder}juno_data_doy_155_errors_{error}.csv',index=False)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
df_iono_155 =  df_new_155

df_iono_155.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
df_iono_155.to_csv(f'{root_folder}juno_data_doy_155_ionotime_errors_{error}.csv',index=False)



# -------------- visits 26, 27, 28 ----------------


times_156_plus = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(df_156_plus['UTC'])):
    timetest = dt.datetime.strptime(df_156_plus['UTC'][p],'%Y-%m-%d %H:%M:%S.%f')
    times_156_plus.append(timetest)
    
times_156_plus = np.array(times_156_plus)


plus_tt_156_plus = []
minus_tt_156_plus = []
errors_156_plus = []
for i in range(len(travel_times_156_plus)):
    calc_156_plus = travel_times_156_plus[i] * error_div
    plus_156_plus = travel_times_156_plus[i] + calc_156_plus
    minus_156_plus = travel_times_156_plus[i] - calc_156_plus
    
    errors_156_plus.append(calc_156_plus)
    plus_tt_156_plus.append(plus_156_plus)
    minus_tt_156_plus.append(minus_156_plus)

    
plus_times_156_plus  = []
minus_times_156_plus = []
for q in range(len(travel_times_156_plus)):
    plus_calc_156_plus = times_156_plus[q] + dt.timedelta(seconds=plus_tt_156_plus[q])
    minus_calc_156_plus = times_156_plus[q] + dt.timedelta(seconds=minus_tt_156_plus[q])
    plus_times_156_plus.append(plus_calc_156_plus)
    minus_times_156_plus.append(minus_calc_156_plus)


df_extra_156_plus = pd.DataFrame()

df_extra_156_plus = df_extra_156_plus.assign(Error_in_Travel_Time=errors_156_plus)
df_extra_156_plus = df_extra_156_plus.assign(Ionosphere_Times_Plus_Error=plus_times_156_plus)
df_extra_156_plus = df_extra_156_plus.assign(Ionosphere_Times_Minus_Error=minus_times_156_plus)


df_new_156_plus = pd.concat([df_156_plus, df_extra_156_plus],axis=1)


df_new_156_plus.to_csv(f'{root_folder}juno_data_doy_156_plus_errors_{error}.csv',index=False)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
df_iono_156_plus =  df_new_156_plus

df_iono_156_plus.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
df_iono_156_plus.to_csv(f'{root_folder}juno_data_doy_156_plus_ionotime_errors_{error}.csv',index=False)

    
# ----------------- visits 34 and 35 ---------------------

times_175 = []
# need to make this into for loop as can only convert one thing into date time at a time
for p in range(len(df_175['UTC'])):
    timetest = dt.datetime.strptime(df_175['UTC'][p],'%Y-%m-%d %H:%M:%S.%f')
    times_175.append(timetest)
    
times_175 = np.array(times_175)
    

plus_tt_175 = []
minus_tt_175 = []
errors_175 = []
for i in range(len(travel_times_175)):
    calc_175 = travel_times_175[i] * error_div
    plus_175 = travel_times_175[i] + calc_175
    minus_175 = travel_times_175[i] - calc_175
    
    errors_175.append(calc_175)
    plus_tt_175.append(plus_175)
    minus_tt_175.append(minus_175)

    
plus_times_175  = []
minus_times_175 = []
for q in range(len(travel_times_175)):
    plus_calc_175 = times_175[q] + dt.timedelta(seconds=plus_tt_175[q])
    minus_calc_175 = times_175[q] + dt.timedelta(seconds=minus_tt_175[q])
    plus_times_175.append(plus_calc_175)
    minus_times_175.append(minus_calc_175)
    
    

df_extra_175 = pd.DataFrame()

df_extra_175 = df_extra_175.assign(Error_in_Travel_Time=errors_175)
df_extra_175 = df_extra_175.assign(Ionosphere_Times_Plus_Error=plus_times_175)
df_extra_175 = df_extra_175.assign(Ionosphere_Times_Minus_Error=minus_times_175)


df_new_175 = pd.concat([df_175, df_extra_175],axis=1)


df_new_175.to_csv(f'{root_folder}juno_data_doy_175_errors_{error}.csv',index=False)

# rearrange dataframe to be sorted by time sw/mag data paramters effect ionosphere
df_iono_175 =  df_new_175

df_iono_175.sort_values(by='Time_Impacts_Ionosphere', inplace=True)

# # export SORTED dataframe to read into other files
df_iono_175.to_csv(f'{root_folder}juno_data_doy_175_ionotime_errors_{error}.csv',index=False)
