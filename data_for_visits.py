"""
Created on Mon Feb  5 11:55:41 2024

@author: hannah

organises and saves out time data (ionospheric time and travel time to ionosphere) for each visit
also saves out index location for big dataframe

precursor to per_visit_compiler
(could combine with?)
"""

# load in relevant modules
import pandas as pd
import spiceypy as spice
import numpy as np

# leap seconds kernal - need this for ephemerous time
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

# dataframes
visit_times = pd.read_csv(root_folder+'visit_times.csv')
pre_155_sw = pd.read_csv(root_folder+'juno_data_big_df_ionotime_updated.csv')

sw_155 = pd.read_csv(root_folder+'juno_data_doy_155_ionotime.csv')
sw_156_plus = pd.read_csv(root_folder+'juno_data_doy_156_plus_ionotime.csv')
sw_175 = pd.read_csv(root_folder+'juno_data_doy_175_ionotime.csv')


# grab all visit tiem data
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

# grab solar wind time reaching ionosphere data in ephemerous time
# need ephemerous time for less than and more than calculations
iono_time = pre_155_sw['ET_Iono_Time'].to_numpy()
iono_time_date = pre_155_sw['Time_Impacts_Ionosphere'].to_numpy()
travel_time = pre_155_sw['Total_Travel_Time'].to_numpy()


iono_time_date_155 = sw_155['Time_Impacts_Ionosphere'].to_numpy()
travel_time_155 = sw_155['Total_Travel_Time'].to_numpy()

iono_time_date_156_plus = sw_156_plus['Time_Impacts_Ionosphere'].to_numpy()
travel_time_156_plus = sw_156_plus['Total_Travel_Time'].to_numpy()

iono_time_date_175 = sw_175['Time_Impacts_Ionosphere'].to_numpy()
travel_time_175 = sw_175['Total_Travel_Time'].to_numpy()


def et_time(time_date):
    iono_time = []
    for i in range(len(time_date)):
        time = spice.str2et(time_date[i])
        iono_time.append(time)
    return iono_time

iono_time_155 = et_time(iono_time_date_155)
iono_time_156_plus = et_time(iono_time_date_156_plus)
iono_time_175 = et_time(iono_time_date_175)


# this function records tne index and value of any data with dates between (and equal to) first and last value 
# need the dataset (ie the ionosphere impact time data), first timestep of the visit and last timestep of the visit
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


def travel_times_grabber(dataset,idx_start,idx_end):
    t_time = []
    for i in range(len(dataset)):
        if (dataset[i] >= dataset[idx_start] and dataset[i] <= dataset[idx_end]):
            t_time.append(dataset)
    return t_time

            

# ----------- main dataframe -----------

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
times_01, idx_01 = between_idxs(iono_time,first_01,last_01)

travel_time_01 = travel_time[idx_01]
iono_time_01 = iono_time_date[idx_01]

# save data out for use in other modules
df_01 = pd.DataFrame()
df_01 = df_01.assign(Big_DF_Index=idx_01)
df_01 = df_01.assign(Ionosphere_Time=iono_time_01)
df_01 = df_01.assign(SW_Travel_Time=travel_time_01)

df_01.to_csv(root_folder+'visit_01_times.csv',index=False)


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
times_02, idx_02 = between_idxs(iono_time,first_02,last_02)

travel_time_02 = travel_time[idx_02]
iono_time_02 = iono_time_date[idx_02]

# save data out for use in other modules
df_02 = pd.DataFrame()
df_02 = df_02.assign(Big_DF_Index=idx_02)
df_02 = df_02.assign(Ionosphere_Time=iono_time_02)
df_02 = df_02.assign(SW_Travel_Time=travel_time_02)

df_02.to_csv(root_folder+'visit_02_times.csv',index=False)


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
times_03, idx_03 = between_idxs(iono_time,first_03,last_03)

travel_time_03 = travel_time[idx_03]
iono_time_03 = iono_time_date[idx_03]

# save data out for use in other modules
df_03 = pd.DataFrame()
df_03 = df_03.assign(Big_DF_Index=idx_03)
df_03 = df_03.assign(Ionosphere_Time=iono_time_03)
df_03 = df_03.assign(SW_Travel_Time=travel_time_03)

df_03.to_csv(root_folder+'visit_03_times.csv',index=False)


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
times_04, idx_04 = between_idxs(iono_time,first_04,last_04)

travel_time_04 = travel_time[idx_04]
iono_time_04 = iono_time_date[idx_04]

# save data out for use in other modules
df_04 = pd.DataFrame()
df_04 = df_04.assign(Big_DF_Index=idx_04)
df_04 = df_04.assign(Ionosphere_Time=iono_time_04)
df_04 = df_04.assign(SW_Travel_Time=travel_time_04)

df_04.to_csv(root_folder+'visit_04_times.csv',index=False)



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
times_05, idx_05 = between_idxs(iono_time,first_05,last_05)

travel_time_05 = travel_time[idx_05]
iono_time_05 = iono_time_date[idx_05]

# save data out for use in other modules
df_05 = pd.DataFrame()
df_05 = df_05.assign(Big_DF_Index=idx_05)
df_05 = df_05.assign(Ionosphere_Time=iono_time_05)
df_05 = df_05.assign(SW_Travel_Time=travel_time_05)

df_05.to_csv(root_folder+'visit_05_times.csv',index=False)
    


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
times_08, idx_08 = between_idxs(iono_time,first_08,last_08)

travel_time_08 = travel_time[idx_08]
iono_time_08 = iono_time_date[idx_08]

# save data out for use in other modules
df_08 = pd.DataFrame()
df_08 = df_08.assign(Big_DF_Index=idx_08)
df_08 = df_08.assign(Ionosphere_Time=iono_time_08)
df_08 = df_08.assign(SW_Travel_Time=travel_time_08)

df_08.to_csv(root_folder+'visit_08_times.csv',index=False)

    

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
times_09, idx_09 = between_idxs(iono_time,first_09,last_09)

travel_time_09 = travel_time[idx_09]
iono_time_09 = iono_time_date[idx_09]

# save data out for use in other modules
df_09 = pd.DataFrame()
df_09 = df_09.assign(Big_DF_Index=idx_09)
df_09 = df_09.assign(Ionosphere_Time=iono_time_09)
df_09 = df_09.assign(SW_Travel_Time=travel_time_09)

df_09.to_csv(root_folder+'visit_09_times.csv',index=False)



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
times_10, idx_10 = between_idxs(iono_time,first_10,last_10)

travel_time_10 = travel_time[idx_10]
iono_time_10 = iono_time_date[idx_10]

# save data out for use in other modules
df_10 = pd.DataFrame()
df_10 = df_10.assign(Big_DF_Index=idx_10)
df_10 = df_10.assign(Ionosphere_Time=iono_time_10)
df_10 = df_10.assign(SW_Travel_Time=travel_time_10)

df_10.to_csv(root_folder+'isit_10_times.csv',index=False)



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
times_11, idx_11 = between_idxs(iono_time,first_11,last_11)

travel_time_11 = travel_time[idx_11]
iono_time_11 = iono_time_date[idx_11]

# save data out for use in other modules
df_11 = pd.DataFrame()
df_11 = df_11.assign(Big_DF_Index=idx_11)
df_11 = df_11.assign(Ionosphere_Time=iono_time_11)
df_11 = df_11.assign(SW_Travel_Time=travel_time_11)

df_11.to_csv(root_folder+'visit_11_times.csv',index=False)



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
times_12, idx_12 = between_idxs(iono_time,first_12,last_12)

travel_time_12 = travel_time[idx_12]
iono_time_12 = iono_time_date[idx_12]

# save data out for use in other modules
df_12 = pd.DataFrame()
df_12 = df_12.assign(Big_DF_Index=idx_12)
df_12 = df_12.assign(Ionosphere_Time=iono_time_12)
df_12 = df_12.assign(SW_Travel_Time=travel_time_12)

df_12.to_csv(root_folder+'visit_12_times.csv',index=False)



'''
VISIT 13
'''
ettime_13 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_13)):
    times = spice.str2et(visit_13[p])
    ettime_13.append(times)

# grab first and last timesteps of the visit
first_13 = ettime_13[0]
last_13 = ettime_13[-1]

# run function to determine times and index values of those times for the visit
times_13, idx_13 = between_idxs(iono_time,first_13,last_13)

travel_time_13 = travel_time[idx_13]
iono_time_13 = iono_time_date[idx_13]

# save data out for use in other modules
df_13 = pd.DataFrame()
df_13 = df_13.assign(Big_DF_Index=idx_13)
df_13 = df_13.assign(Ionosphere_Time=iono_time_13)
df_13 = df_13.assign(SW_Travel_Time=travel_time_13)

df_13.to_csv(root_folder+'visit_13_times.csv',index=False)



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
times_15, idx_15 = between_idxs(iono_time,first_15,last_15)

travel_time_15 = travel_time[idx_15]
iono_time_15 = iono_time_date[idx_15]

# save data out for use in other modules
df_15 = pd.DataFrame()
df_15 = df_15.assign(Big_DF_Index=idx_15)
df_15 = df_15.assign(Ionosphere_Time=iono_time_15)
df_15 = df_15.assign(SW_Travel_Time=travel_time_15)

df_15.to_csv(root_folder+'visit_15_times.csv',index=False)
    


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
times_16, idx_16 = between_idxs(iono_time,first_16,last_16)

travel_time_16 = travel_time[idx_16]
iono_time_16 = iono_time_date[idx_16]

# save data out for use in other modules
df_16 = pd.DataFrame()
df_16 = df_16.assign(Big_DF_Index=idx_16)
df_16 = df_16.assign(Ionosphere_Time=iono_time_16)
df_16 = df_16.assign(SW_Travel_Time=travel_time_16)

df_16.to_csv(root_folder+'visit_16_times.csv',index=False)


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
times_17, idx_17 = between_idxs(iono_time,first_17,last_17)

travel_time_17 = travel_time[idx_17]
iono_time_17 = iono_time_date[idx_17]

# save data out for use in other modules
df_17 = pd.DataFrame()
df_17 = df_17.assign(Big_DF_Index=idx_17)
df_17 = df_17.assign(Ionosphere_Time=iono_time_17)
df_17 = df_17.assign(SW_Travel_Time=travel_time_17)

df_17.to_csv(root_folder+'visit_17_times.csv',index=False)




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
times_18, idx_18 = between_idxs(iono_time,first_18,last_18)

travel_time_18 = travel_time[idx_18]
iono_time_18 = iono_time_date[idx_18]

# save data out for use in other modules
df_18 = pd.DataFrame()
df_18 = df_18.assign(Big_DF_Index=idx_18)
df_18 = df_18.assign(Ionosphere_Time=iono_time_18)
df_18 = df_18.assign(SW_Travel_Time=travel_time_18)

df_18.to_csv(root_folder+'visit_18_times.csv',index=False)



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
times_19, idx_19 = between_idxs(iono_time,first_19,last_19)

travel_time_19 = travel_time[idx_19]
iono_time_19 = iono_time_date[idx_19]

# save data out for use in other modules
df_19 = pd.DataFrame()
df_19 = df_19.assign(Big_DF_Index=idx_19)
df_19 = df_19.assign(Ionosphere_Time=iono_time_19)
df_19 = df_19.assign(SW_Travel_Time=travel_time_19)

df_19.to_csv(root_folder+'visit_19_times.csv',index=False)



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
times_20, idx_20 = between_idxs(iono_time,first_20,last_20)

travel_time_20 = travel_time[idx_20]
iono_time_20 = iono_time_date[idx_20]

# save data out for use in other modules
df_20 = pd.DataFrame()
df_20 = df_20.assign(Big_DF_Index=idx_20)
df_20 = df_20.assign(Ionosphere_Time=iono_time_20)
df_20 = df_20.assign(SW_Travel_Time=travel_time_20)

df_20.to_csv(root_folder+'visit_20_times.csv',index=False)




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
times_21, idx_21 = between_idxs(iono_time,first_21,last_21)

travel_time_21 = travel_time[idx_21]
iono_time_21 = iono_time_date[idx_21]

# save data out for use in other modules
df_21 = pd.DataFrame()
df_21 = df_21.assign(Big_DF_Index=idx_21)
df_21 = df_21.assign(Ionosphere_Time=iono_time_21)
df_21 = df_21.assign(SW_Travel_Time=travel_time_21)

df_21.to_csv(root_folder+'visit_21_times.csv',index=False)


'''
VISIT 23
'''
ettime_23 = []
# convert visit into ephermous time to be able to compare with iono time
for p in range(len(visit_23)):
    times = spice.str2et(visit_23[p])
    ettime_23.append(times)

# grab first and last timesteps of the visit
first_23 = ettime_23[0]
last_23 = ettime_23[-1]

# run function to determine times and index values of those times for the visit
times_23, idx_23 = between_idxs(iono_time,first_23,last_23)

travel_time_23 = travel_time[idx_23]
iono_time_23 = iono_time_date[idx_23]

# save data out for use in other modules
df_23 = pd.DataFrame()
df_23 = df_23.assign(Big_DF_Index=idx_23)
df_23 = df_23.assign(Ionosphere_Time=iono_time_23)
df_23 = df_23.assign(SW_Travel_Time=travel_time_23)

df_23.to_csv(root_folder+'visit_23_times.csv',index=False)


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
times_24, idx_24 = between_idxs(iono_time,first_24,last_24)

travel_time_24 = travel_time[idx_24]
iono_time_24 = iono_time_date[idx_24]

# save data out for use in other modules
df_24 = pd.DataFrame()
df_24 = df_24.assign(Big_DF_Index=idx_24)
df_24 = df_24.assign(Ionosphere_Time=iono_time_24)
df_24 = df_24.assign(SW_Travel_Time=travel_time_24)

df_24.to_csv(root_folder+'visit_24_times.csv',index=False)


# ---------------- 25 -------------

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
times_25, idx_25 = between_idxs(iono_time_155,first_25,last_25)

travel_time_25 = travel_time_155[idx_25]
iono_time_25 = iono_time_date_155[idx_25]

# save data out for use in other modules
df_25 = pd.DataFrame()
df_25 = df_25.assign(Big_DF_Index=idx_25)
df_25 = df_25.assign(Ionosphere_Time=iono_time_25)
df_25 = df_25.assign(SW_Travel_Time=travel_time_25)

df_25.to_csv(root_folder+'visit_25_times.csv',index=False)


# -------------- 26, 27, 28 ------------

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
times_26, idx_26 = between_idxs(iono_time_156_plus,first_26,last_26)

travel_time_26 = travel_time_156_plus[idx_26]
iono_time_26 = iono_time_date_156_plus[idx_26]

# save data out for use in other modules
df_26 = pd.DataFrame()
df_26 = df_26.assign(Big_DF_Index=idx_26)
df_26 = df_26.assign(Ionosphere_Time=iono_time_26)
df_26 = df_26.assign(SW_Travel_Time=travel_time_26)

df_26.to_csv(root_folder+'visit_26_times.csv',index=False)


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
times_27, idx_27 = between_idxs(iono_time_156_plus,first_27,last_27)

travel_time_27 = travel_time_156_plus[idx_27]
iono_time_27 = iono_time_date_156_plus[idx_27]

# save data out for use in other modules
df_27 = pd.DataFrame()
df_27 = df_27.assign(Big_DF_Index=idx_27)
df_27 = df_27.assign(Ionosphere_Time=iono_time_27)
df_27 = df_27.assign(SW_Travel_Time=travel_time_27)

df_27.to_csv(root_folder+'visit_27_times.csv',index=False)


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
times_28, idx_28 = between_idxs(iono_time_156_plus,first_28,last_28)

travel_time_28 = travel_time_156_plus[idx_28]
iono_time_28 = iono_time_date_156_plus[idx_28]

# save data out for use in other modules
df_28 = pd.DataFrame()
df_28 = df_28.assign(Big_DF_Index=idx_28)
df_28 = df_28.assign(Ionosphere_Time=iono_time_28)
df_28 = df_28.assign(SW_Travel_Time=travel_time_28)

df_28.to_csv(root_folder+'visit_28_times.csv',index=False)


# ------------------------- 34 & 35 -------------------------


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
times_34, idx_34 = between_idxs(iono_time_175,first_34,last_34)

travel_time_34 = travel_time_175[idx_34]
iono_time_34 = iono_time_date_175[idx_34]

# save data out for use in other modules
df_34 = pd.DataFrame()
df_34 = df_34.assign(Big_DF_Index=idx_34)
df_34 = df_34.assign(Ionosphere_Time=iono_time_34)
df_34 = df_34.assign(SW_Travel_Time=travel_time_34)

df_34.to_csv(root_folder+'visit_34_times.csv',index=False)


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
times_35, idx_35 = between_idxs(iono_time_175,first_35,last_35)

travel_time_35 = travel_time_175[idx_35]
iono_time_35 = iono_time_date_175[idx_35]

# save data out for use in other modules
df_35 = pd.DataFrame()
df_35 = df_35.assign(Big_DF_Index=idx_35)
df_35 = df_35.assign(Ionosphere_Time=iono_time_35)
df_35 = df_35.assign(SW_Travel_Time=travel_time_35)

df_35.to_csv(root_folder+'visit_35_times.csv',index=False)
