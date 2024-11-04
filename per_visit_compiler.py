"""
Created on Fri May 24 10:44:51 2024

@author: hannah

this script reads in all solar wind related dataframes and saves out new dataframes breaking down
storing the datasets the previous dataframes into smaller parts to organise by visit

saves out as visit_{visit}_data
"""
import pandas as pd

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

df = pd.read_csv(root_folder+'juno_data_big_df_ionotime_new.csv')
df_155 = pd.read_csv(root_folder+'juno_data_doy_155_ionotime.csv')
df_156_plus = pd.read_csv(root_folder+'juno_data_doy_156_plus_ionotime.csv')
df_175 = pd.read_csv(root_folder+'juno_data_doy_175_ionotime.csv')

visit_01 = pd.read_csv(root_folder+'visit_01_times.csv',delimiter=',')
visit_02 = pd.read_csv(root_folder+'visit_02_times.csv',delimiter=',')
visit_03 = pd.read_csv(root_folder+'visit_03_times.csv',delimiter=',')
visit_04 = pd.read_csv(root_folder+'visit_04_times.csv',delimiter=',')
visit_05 = pd.read_csv(root_folder+'visit_05_times.csv',delimiter=',')
visit_08 = pd.read_csv(root_folder+'visit_08_times.csv',delimiter=',')
visit_09 = pd.read_csv(root_folder+'visit_09_times.csv',delimiter=',')
visit_10 = pd.read_csv(root_folder+'visit_10_times.csv',delimiter=',')
visit_11 = pd.read_csv(root_folder+'visit_11_times.csv',delimiter=',')
visit_12 = pd.read_csv(root_folder+'visit_12_times.csv',delimiter=',')
visit_13 = pd.read_csv(root_folder+'visit_13_times.csv',delimiter=',')
visit_15 = pd.read_csv(root_folder+'visit_15_times.csv',delimiter=',')
visit_16 = pd.read_csv(root_folder+'visit_16_times.csv',delimiter=',')
visit_17 = pd.read_csv(root_folder+'visit_17_times.csv',delimiter=',')
visit_18 = pd.read_csv(root_folder+'visit_18_times.csv',delimiter=',')
visit_19 = pd.read_csv(root_folder+'visit_19_times.csv',delimiter=',')
visit_20 = pd.read_csv(root_folder+'visit_20_times.csv',delimiter=',')
visit_21 = pd.read_csv(root_folder+'visit_21_times.csv',delimiter=',')
visit_23 = pd.read_csv(root_folder+'visit_23_times.csv',delimiter=',')
visit_24 = pd.read_csv(root_folder+'visit_24_times.csv',delimiter=',')
visit_25 = pd.read_csv(root_folder+'visit_25_times.csv',delimiter=',')
visit_26 = pd.read_csv(root_folder+'visit_26_times.csv',delimiter=',')
visit_27 = pd.read_csv(root_folder+'visit_27_times.csv',delimiter=',')
visit_28 = pd.read_csv(root_folder+'visit_28_times.csv',delimiter=',')
visit_34 = pd.read_csv(root_folder+'visit_34_times.csv',delimiter=',')
visit_35 = pd.read_csv(root_folder+'visit_35_times.csv',delimiter=',')


index_01 = visit_01['Big_DF_Index'].to_numpy()
index_02 = visit_02['Big_DF_Index'].to_numpy()
index_03 = visit_03['Big_DF_Index'].to_numpy()
index_04 = visit_04['Big_DF_Index'].to_numpy()
index_05 = visit_05['Big_DF_Index'].to_numpy()
index_08 = visit_08['Big_DF_Index'].to_numpy()
index_09 = visit_09['Big_DF_Index'].to_numpy()
index_10 = visit_10['Big_DF_Index'].to_numpy()
index_11 = visit_11['Big_DF_Index'].to_numpy()
index_12 = visit_12['Big_DF_Index'].to_numpy()
index_13 = visit_13['Big_DF_Index'].to_numpy()
index_15 = visit_15['Big_DF_Index'].to_numpy()
index_16 = visit_16['Big_DF_Index'].to_numpy()
index_17 = visit_17['Big_DF_Index'].to_numpy()
index_18 = visit_18['Big_DF_Index'].to_numpy()
index_19 = visit_19['Big_DF_Index'].to_numpy()
index_20 = visit_20['Big_DF_Index'].to_numpy()
index_21 = visit_21['Big_DF_Index'].to_numpy()
index_23 = visit_23['Big_DF_Index'].to_numpy()
index_24 = visit_24['Big_DF_Index'].to_numpy()
index_25 = visit_25['Big_DF_Index'].to_numpy()
index_26 = visit_26['Big_DF_Index'].to_numpy()
index_27 = visit_27['Big_DF_Index'].to_numpy()
index_28 = visit_28['Big_DF_Index'].to_numpy()
index_34 = visit_34['Big_DF_Index'].to_numpy()
index_35 = visit_35['Big_DF_Index'].to_numpy()


iono_time_01 = visit_01['Ionosphere_Time'].to_numpy()
iono_time_02 = visit_02['Ionosphere_Time'].to_numpy()
iono_time_03 = visit_03['Ionosphere_Time'].to_numpy()
iono_time_04 = visit_04['Ionosphere_Time'].to_numpy()
iono_time_05 = visit_05['Ionosphere_Time'].to_numpy()
iono_time_08 = visit_08['Ionosphere_Time'].to_numpy()
iono_time_09 = visit_09['Ionosphere_Time'].to_numpy()
iono_time_10 = visit_10['Ionosphere_Time'].to_numpy()
iono_time_11 = visit_11['Ionosphere_Time'].to_numpy()
iono_time_12 = visit_12['Ionosphere_Time'].to_numpy()
iono_time_13 = visit_13['Ionosphere_Time'].to_numpy()
iono_time_15 = visit_15['Ionosphere_Time'].to_numpy()
iono_time_16 = visit_16['Ionosphere_Time'].to_numpy()
iono_time_17 = visit_17['Ionosphere_Time'].to_numpy()
iono_time_18 = visit_18['Ionosphere_Time'].to_numpy()
iono_time_19 = visit_19['Ionosphere_Time'].to_numpy()
iono_time_20 = visit_20['Ionosphere_Time'].to_numpy()
iono_time_21 = visit_21['Ionosphere_Time'].to_numpy()
iono_time_23 = visit_23['Ionosphere_Time'].to_numpy()
iono_time_24 = visit_24['Ionosphere_Time'].to_numpy()
iono_time_25 = visit_25['Ionosphere_Time'].to_numpy()
iono_time_26 = visit_26['Ionosphere_Time'].to_numpy()
iono_time_27 = visit_27['Ionosphere_Time'].to_numpy()
iono_time_28 = visit_28['Ionosphere_Time'].to_numpy()
iono_time_34 = visit_34['Ionosphere_Time'].to_numpy()
iono_time_35 = visit_35['Ionosphere_Time'].to_numpy()


tt_01 = visit_01['SW_Travel_Time'].to_numpy()
tt_02 = visit_02['SW_Travel_Time'].to_numpy()
tt_03 = visit_03['SW_Travel_Time'].to_numpy()
tt_04 = visit_04['SW_Travel_Time'].to_numpy()
tt_05 = visit_05['SW_Travel_Time'].to_numpy()
tt_08 = visit_08['SW_Travel_Time'].to_numpy()
tt_09 = visit_09['SW_Travel_Time'].to_numpy()
tt_10 = visit_10['SW_Travel_Time'].to_numpy()
tt_11 = visit_11['SW_Travel_Time'].to_numpy()
tt_12 = visit_12['SW_Travel_Time'].to_numpy()
tt_13 = visit_13['SW_Travel_Time'].to_numpy()
tt_15 = visit_15['SW_Travel_Time'].to_numpy()
tt_16 = visit_16['SW_Travel_Time'].to_numpy()
tt_17 = visit_17['SW_Travel_Time'].to_numpy()
tt_18 = visit_18['SW_Travel_Time'].to_numpy()
tt_19 = visit_19['SW_Travel_Time'].to_numpy()
tt_20 = visit_20['SW_Travel_Time'].to_numpy()
tt_21 = visit_21['SW_Travel_Time'].to_numpy()
tt_23 = visit_23['SW_Travel_Time'].to_numpy()
tt_24 = visit_24['SW_Travel_Time'].to_numpy()
tt_25 = visit_25['SW_Travel_Time'].to_numpy()
tt_26 = visit_26['SW_Travel_Time'].to_numpy()
tt_27 = visit_27['SW_Travel_Time'].to_numpy()
tt_28 = visit_28['SW_Travel_Time'].to_numpy()
tt_34 = visit_34['SW_Travel_Time'].to_numpy()
tt_35 = visit_35['SW_Travel_Time'].to_numpy()


pressure_array = df['SW_Pressure'].to_numpy()
velocity_array = df['SW_Velocity'].to_numpy()
clock_angle = df['Clock_Angle'].to_numpy()
#clock_angle_err = df['Clock_Angle_Error'].to_numpy()
b_perp = df['B_Perp'].to_numpy()
Br = df['Br'].to_numpy()
Bt = df['Bt'].to_numpy()
Bn = df['Bn'].to_numpy()
juno_time = df['Juno_Detection_Time'].to_numpy()
LL_rec = df['Low_Latitude_Reconnection_Voltage'].to_numpy()
HL_rec_neg = df['High_Latitude_Reconnection_Voltage_BY_NEG'].to_numpy()
HL_rec_pos = df['High_Latitude_Reconnection_Voltage_BY_POS'].to_numpy()
bs_loc = df['Bow_Shock_Stand_Off_RJ'].to_numpy()
mp_loc = df['Magnetopause_Stand_Off_RJ'].to_numpy()

pressure_array_155 = df_155['Pressure'].to_numpy()
velocity_array_155 = df_155['Velocity'].to_numpy()
clock_angle_155 = df_155['Clock_Angle'].to_numpy()
b_perp_155 = df_155['B_perp'].to_numpy()
Br_155 = df_155['Br'].to_numpy()
Bt_155 = df_155['Bt'].to_numpy()
Bn_155 = df_155['Bn'].to_numpy()
juno_time_155 = df_155['UTC'].to_numpy()
LL_rec_155 = df_155['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_155 = df_155['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_155 = df_155['High_Latitude_Reconnection_By_Pos'].to_numpy()
bs_loc_155 = df_155['Bow_Shock_Location'].to_numpy()
mp_loc_155 = df_155['Magnetopause_Location'].to_numpy()

pressure_array_156_plus = df_156_plus['Pressure'].to_numpy()
velocity_array_156_plus = df_156_plus['Velocity'].to_numpy()
clock_angle_156_plus = df_156_plus['Clock_Angle'].to_numpy()
b_perp_156_plus = df_156_plus['B_perp'].to_numpy()
Br_156_plus = df_156_plus['Br'].to_numpy()
Bt_156_plus = df_156_plus['Bt'].to_numpy()
Bn_156_plus = df_156_plus['Bn'].to_numpy()
juno_time_156_plus = df_156_plus['UTC'].to_numpy()
LL_rec_156_plus = df_156_plus['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_156_plus = df_156_plus['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_156_plus = df_156_plus['High_Latitude_Reconnection_By_Pos'].to_numpy()
bs_loc_156_plus = df_156_plus['Bow_Shock_Location'].to_numpy()
mp_loc_156_plus = df_156_plus['Magnetopause_Location'].to_numpy()

pressure_array_175 = df_175['Pressure'].to_numpy()
velocity_array_175 = df_175['Velocity'].to_numpy()
clock_angle_175 = df_175['Clock_Angle'].to_numpy()
b_perp_175 = df_175['B_perp'].to_numpy()
Br_175 = df_175['Br'].to_numpy()
Bt_175 = df_175['Bt'].to_numpy()
Bn_175 = df_175['Bn'].to_numpy()
juno_time_175 = df_175['UTC'].to_numpy()
LL_rec_175 = df_175['Low_Latitude_Reconnection'].to_numpy()
HL_rec_neg_175 = df_175['High_Latitude_Reconnection_By_Neg'].to_numpy()
HL_rec_pos_175 = df_175['High_Latitude_Reconnection_By_Pos'].to_numpy()
bs_loc_175 = df_175['Bow_Shock_Location'].to_numpy()
mp_loc_175 = df_175['Magnetopause_Location'].to_numpy()


def index_compare(index, variable):
    index_values = []
    for i  in index:
        index_values.append(variable[i])
    return index_values


bs_01 = index_compare(index_01, bs_loc)
bs_02 = index_compare(index_02, bs_loc)
bs_03 = index_compare(index_03, bs_loc)
bs_04 = index_compare(index_04, bs_loc)
bs_05 = index_compare(index_05, bs_loc)
bs_08 = index_compare(index_08, bs_loc)
bs_09 = index_compare(index_09, bs_loc)
bs_10 = index_compare(index_10, bs_loc)
bs_11 = index_compare(index_11, bs_loc)
bs_12 = index_compare(index_12, bs_loc)
bs_13 = index_compare(index_13, bs_loc)
bs_15 = index_compare(index_15, bs_loc)
bs_16 = index_compare(index_16, bs_loc)
bs_17 = index_compare(index_17, bs_loc)
bs_18 = index_compare(index_18, bs_loc)
bs_19 = index_compare(index_19, bs_loc)
bs_20 = index_compare(index_20, bs_loc)
bs_21 = index_compare(index_21, bs_loc)
bs_23 = index_compare(index_23, bs_loc)
bs_24 = index_compare(index_24, bs_loc)

bs_25 = index_compare(index_25, bs_loc_155)

bs_26 = index_compare(index_26, bs_loc_156_plus)
bs_27 = index_compare(index_27, bs_loc_156_plus)
bs_28 = index_compare(index_28, bs_loc_156_plus)

bs_34 = index_compare(index_34, bs_loc_175)
bs_35 = index_compare(index_35, bs_loc_175)


mp_01 = index_compare(index_01, mp_loc)
mp_02 = index_compare(index_02, mp_loc)
mp_03 = index_compare(index_03, mp_loc)
mp_04 = index_compare(index_04, mp_loc)
mp_05 = index_compare(index_05, mp_loc)
mp_08 = index_compare(index_08, mp_loc)
mp_09 = index_compare(index_09, mp_loc)
mp_10 = index_compare(index_10, mp_loc)
mp_11 = index_compare(index_11, mp_loc)
mp_12 = index_compare(index_12, mp_loc)
mp_13 = index_compare(index_13, mp_loc)
mp_15 = index_compare(index_15, mp_loc)
mp_16 = index_compare(index_16, mp_loc)
mp_17 = index_compare(index_17, mp_loc)
mp_18 = index_compare(index_18, mp_loc)
mp_19 = index_compare(index_19, mp_loc)
mp_20 = index_compare(index_20, mp_loc)
mp_21 = index_compare(index_21, mp_loc)
mp_23 = index_compare(index_23, mp_loc)
mp_24 = index_compare(index_24, mp_loc)

mp_25 = index_compare(index_25, mp_loc_155)

mp_26 = index_compare(index_26, mp_loc_156_plus)
mp_27 = index_compare(index_27, mp_loc_156_plus)
mp_28 = index_compare(index_28, mp_loc_156_plus)

mp_34 = index_compare(index_34, mp_loc_175)
mp_35 = index_compare(index_35, mp_loc_175)


br_01 = index_compare(index_01, Br)
br_02 = index_compare(index_02, Br)
br_03 = index_compare(index_03, Br)
br_04 = index_compare(index_04, Br)
br_05 = index_compare(index_05, Br)
br_08 = index_compare(index_08, Br)
br_09 = index_compare(index_09, Br)
br_10 = index_compare(index_10, Br)
br_11 = index_compare(index_11, Br)
br_12 = index_compare(index_12, Br)
br_13 = index_compare(index_13, Br)
br_15 = index_compare(index_15, Br)
br_16 = index_compare(index_16, Br)
br_17 = index_compare(index_17, Br)
br_18 = index_compare(index_18, Br)
br_19 = index_compare(index_19, Br)
br_20 = index_compare(index_20, Br)
br_21 = index_compare(index_21, Br)
br_23 = index_compare(index_23, Br)
br_24 = index_compare(index_24, Br)

br_25 = index_compare(index_25, Br_155)

br_26 = index_compare(index_26, Br_156_plus)
br_27 = index_compare(index_27, Br_156_plus)
br_28 = index_compare(index_28, Br_156_plus)

br_34 = index_compare(index_34, Br_175)
br_35 = index_compare(index_35, Br_175)


bt_01 = index_compare(index_01, Bt)
bt_02 = index_compare(index_02, Bt)
bt_03 = index_compare(index_03, Bt)
bt_04 = index_compare(index_04, Bt)
bt_05 = index_compare(index_05, Bt)
bt_08 = index_compare(index_08, Bt)
bt_09 = index_compare(index_09, Bt)
bt_10 = index_compare(index_10, Bt)
bt_11 = index_compare(index_11, Bt)
bt_12 = index_compare(index_12, Bt)
bt_13 = index_compare(index_13, Bt)
bt_15 = index_compare(index_15, Bt)
bt_16 = index_compare(index_16, Bt)
bt_17 = index_compare(index_17, Bt)
bt_18 = index_compare(index_18, Bt)
bt_19 = index_compare(index_19, Bt)
bt_20 = index_compare(index_20, Bt)
bt_21 = index_compare(index_21, Bt)
bt_23 = index_compare(index_23, Bt)
bt_24 = index_compare(index_24, Bt)

bt_25 = index_compare(index_25, Bt_155)

bt_26 = index_compare(index_26, Bt_156_plus)
bt_27 = index_compare(index_27, Bt_156_plus)
bt_28 = index_compare(index_28, Bt_156_plus)

bt_34 = index_compare(index_34, Bt_175)
bt_35 = index_compare(index_35, Bt_175)


bn_01 = index_compare(index_01, Bn)
bn_02 = index_compare(index_02, Bn)
bn_03 = index_compare(index_03, Bn)
bn_04 = index_compare(index_04, Bn)
bn_05 = index_compare(index_05, Bn)
bn_08 = index_compare(index_08, Bn)
bn_09 = index_compare(index_09, Bn)
bn_10 = index_compare(index_10, Bn)
bn_11 = index_compare(index_11, Bn)
bn_12 = index_compare(index_12, Bn)
bn_13 = index_compare(index_13, Bn)
bn_15 = index_compare(index_15, Bn)
bn_16 = index_compare(index_16, Bn)
bn_17 = index_compare(index_17, Bn)
bn_18 = index_compare(index_18, Bn)
bn_19 = index_compare(index_19, Bn)
bn_20 = index_compare(index_20, Bn)
bn_21 = index_compare(index_21, Bn)
bn_23 = index_compare(index_23, Bn)
bn_24 = index_compare(index_24, Bn)

bn_25 = index_compare(index_25, Bn_155)

bn_26 = index_compare(index_26, Bn_156_plus)
bn_27 = index_compare(index_27, Bn_156_plus)
bn_28 = index_compare(index_28, Bn_156_plus)

bn_34 = index_compare(index_34, Bn_175)
bn_35 = index_compare(index_35, Bn_175)


clock_01 = index_compare(index_01, clock_angle)
clock_02 = index_compare(index_02, clock_angle)
clock_03 = index_compare(index_03, clock_angle)
clock_04 = index_compare(index_04, clock_angle)
clock_05 = index_compare(index_05, clock_angle)
clock_08 = index_compare(index_08, clock_angle)
clock_09 = index_compare(index_09, clock_angle)
clock_10 = index_compare(index_10, clock_angle)
clock_11 = index_compare(index_11, clock_angle)
clock_12 = index_compare(index_12, clock_angle)
clock_13 = index_compare(index_13, clock_angle)
clock_15 = index_compare(index_15, clock_angle)
clock_16 = index_compare(index_16, clock_angle)
clock_17 = index_compare(index_17, clock_angle)
clock_18 = index_compare(index_18, clock_angle)
clock_19 = index_compare(index_19, clock_angle)
clock_20 = index_compare(index_20, clock_angle)
clock_21 = index_compare(index_21, clock_angle)
clock_23 = index_compare(index_23, clock_angle)
clock_24 = index_compare(index_24, clock_angle)

clock_25 = index_compare(index_25, clock_angle_155)

clock_26 = index_compare(index_26, clock_angle_156_plus)
clock_27 = index_compare(index_27, clock_angle_156_plus)
clock_28 = index_compare(index_28, clock_angle_156_plus)

clock_34 = index_compare(index_34, clock_angle_175)
clock_35 = index_compare(index_35, clock_angle_175)

b_perp_01 = index_compare(index_01, b_perp)
b_perp_02 = index_compare(index_02, b_perp)
b_perp_03 = index_compare(index_03, b_perp)
b_perp_04 = index_compare(index_04, b_perp)
b_perp_05 = index_compare(index_05, b_perp)
b_perp_08 = index_compare(index_08, b_perp)
b_perp_09 = index_compare(index_09, b_perp)
b_perp_10 = index_compare(index_10, b_perp)
b_perp_11 = index_compare(index_11, b_perp)
b_perp_12 = index_compare(index_12, b_perp)
b_perp_13 = index_compare(index_13, b_perp)
b_perp_15 = index_compare(index_15, b_perp)
b_perp_16 = index_compare(index_16, b_perp)
b_perp_17 = index_compare(index_17, b_perp)
b_perp_18 = index_compare(index_18, b_perp)
b_perp_19 = index_compare(index_19, b_perp)
b_perp_20 = index_compare(index_20, b_perp)
b_perp_21 = index_compare(index_21, b_perp)
b_perp_23 = index_compare(index_23, b_perp)
b_perp_24 = index_compare(index_24, b_perp)

b_perp_25 = index_compare(index_25, b_perp_155)

b_perp_26 = index_compare(index_26, b_perp_156_plus)
b_perp_27 = index_compare(index_27, b_perp_156_plus)
b_perp_28 = index_compare(index_28, b_perp_156_plus)

b_perp_34 = index_compare(index_34, b_perp_175)
b_perp_35 = index_compare(index_35, b_perp_175)

velocity_01 = index_compare(index_01, velocity_array)
velocity_02 = index_compare(index_02, velocity_array)
velocity_03 = index_compare(index_03, velocity_array)
velocity_04 = index_compare(index_04, velocity_array)
velocity_05 = index_compare(index_05, velocity_array)
velocity_08 = index_compare(index_08, velocity_array)
velocity_09 = index_compare(index_09, velocity_array)
velocity_10 = index_compare(index_10, velocity_array)
velocity_11 = index_compare(index_11, velocity_array)
velocity_12 = index_compare(index_12, velocity_array)
velocity_13 = index_compare(index_13, velocity_array)
velocity_15 = index_compare(index_15, velocity_array)
velocity_16 = index_compare(index_16, velocity_array)
velocity_17 = index_compare(index_17, velocity_array)
velocity_18 = index_compare(index_18, velocity_array)
velocity_19 = index_compare(index_19, velocity_array)
velocity_20 = index_compare(index_20, velocity_array)
velocity_21 = index_compare(index_21, velocity_array)
velocity_23 = index_compare(index_23, velocity_array)
velocity_24 = index_compare(index_24, velocity_array)

velocity_25 = index_compare(index_25, velocity_array_155)

velocity_26 = index_compare(index_26, velocity_array_156_plus)
velocity_27 = index_compare(index_27, velocity_array_156_plus)
velocity_28 = index_compare(index_28, velocity_array_156_plus)

velocity_34 = index_compare(index_34, velocity_array_175)
velocity_35 = index_compare(index_35, velocity_array_175)

pressure_01 = index_compare(index_01, pressure_array)
pressure_02 = index_compare(index_02, pressure_array)
pressure_03 = index_compare(index_03, pressure_array)
pressure_04 = index_compare(index_04, pressure_array)
pressure_05 = index_compare(index_05, pressure_array)
pressure_08 = index_compare(index_08, pressure_array)
pressure_09 = index_compare(index_09, pressure_array)
pressure_10 = index_compare(index_10, pressure_array)
pressure_11 = index_compare(index_11, pressure_array)
pressure_12 = index_compare(index_12, pressure_array)
pressure_13 = index_compare(index_13, pressure_array)
pressure_15 = index_compare(index_15, pressure_array)
pressure_16 = index_compare(index_16, pressure_array)
pressure_17 = index_compare(index_17, pressure_array)
pressure_18 = index_compare(index_18, pressure_array)
pressure_19 = index_compare(index_19, pressure_array)
pressure_20 = index_compare(index_20, pressure_array)
pressure_21 = index_compare(index_21, pressure_array)
pressure_23 = index_compare(index_23, pressure_array)
pressure_24 = index_compare(index_24, pressure_array)

pressure_25 = index_compare(index_25, pressure_array_155)

pressure_26 = index_compare(index_26, pressure_array_156_plus)
pressure_27 = index_compare(index_27, pressure_array_156_plus)
pressure_28 = index_compare(index_28, pressure_array_156_plus)

pressure_34 = index_compare(index_34, pressure_array_175)
pressure_35 = index_compare(index_35, pressure_array_175)

juno_time_01 = index_compare(index_01, juno_time)
juno_time_02 = index_compare(index_02, juno_time)
juno_time_03 = index_compare(index_03, juno_time)
juno_time_04 = index_compare(index_04, juno_time)
juno_time_05 = index_compare(index_05, juno_time)
juno_time_08 = index_compare(index_08, juno_time)
juno_time_09 = index_compare(index_09, juno_time)
juno_time_10 = index_compare(index_10, juno_time)
juno_time_11 = index_compare(index_11, juno_time)
juno_time_12 = index_compare(index_12, juno_time)
juno_time_13 = index_compare(index_13, juno_time)
juno_time_15 = index_compare(index_15, juno_time)
juno_time_16 = index_compare(index_16, juno_time)
juno_time_17 = index_compare(index_17, juno_time)
juno_time_18 = index_compare(index_18, juno_time)
juno_time_19 = index_compare(index_19, juno_time)
juno_time_20 = index_compare(index_20, juno_time)
juno_time_21 = index_compare(index_21, juno_time)
juno_time_23 = index_compare(index_23, juno_time)
juno_time_24 = index_compare(index_24, juno_time)

juno_time_25 = index_compare(index_25, juno_time_155)

juno_time_26 = index_compare(index_26, juno_time_156_plus)
juno_time_27 = index_compare(index_27, juno_time_156_plus)
juno_time_28 = index_compare(index_28, juno_time_156_plus)

juno_time_34 = index_compare(index_34, juno_time_175)
juno_time_35 = index_compare(index_35, juno_time_175)


LL_01 = index_compare(index_01, LL_rec)
LL_02 = index_compare(index_02, LL_rec)
LL_03 = index_compare(index_03, LL_rec)
LL_04 = index_compare(index_04, LL_rec)
LL_05 = index_compare(index_05, LL_rec)
LL_08 = index_compare(index_08, LL_rec)
LL_09 = index_compare(index_09, LL_rec)
LL_10 = index_compare(index_10, LL_rec)
LL_11 = index_compare(index_11, LL_rec)
LL_12 = index_compare(index_12, LL_rec)
LL_13 = index_compare(index_13, LL_rec)
LL_15 = index_compare(index_15, LL_rec)
LL_16 = index_compare(index_16, LL_rec)
LL_17 = index_compare(index_17, LL_rec)
LL_18 = index_compare(index_18, LL_rec)
LL_19 = index_compare(index_19, LL_rec)
LL_20 = index_compare(index_20, LL_rec)
LL_21 = index_compare(index_21, LL_rec)
LL_23 = index_compare(index_23, LL_rec)
LL_24 = index_compare(index_24, LL_rec)

LL_25 = index_compare(index_25, LL_rec_155)

LL_26 = index_compare(index_26, LL_rec_156_plus)
LL_27 = index_compare(index_27, LL_rec_156_plus)
LL_28 = index_compare(index_28, LL_rec_156_plus)

LL_34 = index_compare(index_34, LL_rec_175)
LL_35 = index_compare(index_35, LL_rec_175)


HL_01_pos = index_compare(index_01, HL_rec_pos)
HL_02_pos = index_compare(index_02, HL_rec_pos)
HL_03_pos = index_compare(index_03, HL_rec_pos)
HL_04_pos = index_compare(index_04, HL_rec_pos)
HL_05_pos = index_compare(index_05, HL_rec_pos)
HL_08_pos = index_compare(index_08, HL_rec_pos)
HL_09_pos = index_compare(index_09, HL_rec_pos)
HL_10_pos = index_compare(index_10, HL_rec_pos)
HL_11_pos = index_compare(index_11, HL_rec_pos)
HL_12_pos = index_compare(index_12, HL_rec_pos)
HL_13_pos = index_compare(index_13, HL_rec_pos)
HL_15_pos = index_compare(index_15, HL_rec_pos)
HL_16_pos = index_compare(index_16, HL_rec_pos)
HL_17_pos = index_compare(index_17, HL_rec_pos)
HL_18_pos = index_compare(index_18, HL_rec_pos)
HL_19_pos = index_compare(index_19, HL_rec_pos)
HL_20_pos = index_compare(index_20, HL_rec_pos)
HL_21_pos = index_compare(index_21, HL_rec_pos)
HL_23_pos = index_compare(index_23, HL_rec_pos)
HL_24_pos = index_compare(index_24, HL_rec_pos)

HL_25_pos = index_compare(index_25, HL_rec_pos_155)

HL_26_pos = index_compare(index_26, HL_rec_pos_156_plus)
HL_27_pos = index_compare(index_27, HL_rec_pos_156_plus)
HL_28_pos = index_compare(index_28, HL_rec_pos_156_plus)

HL_34_pos = index_compare(index_34, HL_rec_pos_175)
HL_35_pos = index_compare(index_35, HL_rec_pos_175)

HL_01_neg = index_compare(index_01, HL_rec_neg)
HL_02_neg = index_compare(index_02, HL_rec_neg)
HL_03_neg = index_compare(index_03, HL_rec_neg)
HL_04_neg = index_compare(index_04, HL_rec_neg)
HL_05_neg = index_compare(index_05, HL_rec_neg)
HL_08_neg = index_compare(index_08, HL_rec_neg)
HL_09_neg = index_compare(index_09, HL_rec_neg)
HL_10_neg = index_compare(index_10, HL_rec_neg)
HL_11_neg = index_compare(index_11, HL_rec_neg)
HL_12_neg = index_compare(index_12, HL_rec_neg)
HL_13_neg = index_compare(index_13, HL_rec_neg)
HL_15_neg = index_compare(index_15, HL_rec_neg)
HL_16_neg = index_compare(index_16, HL_rec_neg)
HL_17_neg = index_compare(index_17, HL_rec_neg)
HL_18_neg = index_compare(index_18, HL_rec_neg)
HL_19_neg = index_compare(index_19, HL_rec_neg)
HL_20_neg = index_compare(index_20, HL_rec_neg)
HL_21_neg = index_compare(index_21, HL_rec_neg)
HL_23_neg = index_compare(index_23, HL_rec_neg)
HL_24_neg = index_compare(index_24, HL_rec_neg)

HL_25_neg = index_compare(index_25, HL_rec_neg_155)

HL_26_neg = index_compare(index_26, HL_rec_neg_156_plus)
HL_27_neg = index_compare(index_27, HL_rec_neg_156_plus)
HL_28_neg = index_compare(index_28, HL_rec_neg_156_plus)

HL_34_neg = index_compare(index_34, HL_rec_neg_175)
HL_35_neg = index_compare(index_35, HL_rec_neg_175)



'''
# VISIT 01
'''

visit_01_df = pd.DataFrame()

visit_01_df = visit_01_df.assign(Index=index_01)
visit_01_df = visit_01_df.assign(Ionospheric_Time=iono_time_01)
visit_01_df = visit_01_df.assign(SW_Travel_Time=tt_01)
visit_01_df = visit_01_df.assign(Juno_Detection_Time=juno_time_01)
visit_01_df = visit_01_df.assign(Sw_Pressure=pressure_01)
visit_01_df = visit_01_df.assign(SW_Velocity=velocity_01)
visit_01_df = visit_01_df.assign(Clock_Angle=clock_01)
visit_01_df = visit_01_df.assign(B_Perp=b_perp_01)
visit_01_df = visit_01_df.assign(Br=br_01)
visit_01_df = visit_01_df.assign(Bt=bt_01)
visit_01_df = visit_01_df.assign(Bn=bn_01)
visit_01_df = visit_01_df.assign(BS_Location=bs_01)
visit_01_df = visit_01_df.assign(MP_Location=mp_01)
visit_01_df = visit_01_df.assign(LL_Recon_V=LL_01)
visit_01_df = visit_01_df.assign(HL_Recon_V_By_Plus=HL_01_pos)
visit_01_df = visit_01_df.assign(HL_Recon_V_By_Neg=HL_01_neg)

visit_01_df.to_csv(root_folder+'visit_01_data.csv',index=False)


'''
# VISIT 02
'''

visit_02_df = pd.DataFrame()

visit_02_df = visit_02_df.assign(Index=index_02)
visit_02_df = visit_02_df.assign(Ionospheric_Time=iono_time_02)
visit_02_df = visit_02_df.assign(SW_Travel_Time=tt_02)
visit_02_df = visit_02_df.assign(Juno_Detection_Time=juno_time_02)
visit_02_df = visit_02_df.assign(Sw_Pressure=pressure_02)
visit_02_df = visit_02_df.assign(SW_Velocity=velocity_02)
visit_02_df = visit_02_df.assign(Clock_Angle=clock_02)
visit_02_df = visit_02_df.assign(B_Perp=b_perp_02)
visit_02_df = visit_02_df.assign(Br=br_02)
visit_02_df = visit_02_df.assign(Bt=bt_02)
visit_02_df = visit_02_df.assign(Bn=bn_02)
visit_02_df = visit_02_df.assign(BS_Location=bs_02)
visit_02_df = visit_02_df.assign(MP_Location=mp_02)
visit_02_df = visit_02_df.assign(LL_Recon_V=LL_02)
visit_02_df = visit_02_df.assign(HL_Recon_V_By_Plus=HL_02_pos)
visit_02_df = visit_02_df.assign(HL_Recon_V_By_Neg=HL_02_neg)

visit_02_df.to_csv(root_folder+'visit_02_data.csv',index=False)


'''
# VISIT 03
'''

visit_03_df = pd.DataFrame()

visit_03_df = visit_03_df.assign(Index=index_03)
visit_03_df = visit_03_df.assign(Ionospheric_Time=iono_time_03)
visit_03_df = visit_03_df.assign(SW_Travel_Time=tt_03)
visit_03_df = visit_03_df.assign(Juno_Detection_Time=juno_time_03)
visit_03_df = visit_03_df.assign(Sw_Pressure=pressure_03)
visit_03_df = visit_03_df.assign(SW_Velocity=velocity_03)
visit_03_df = visit_03_df.assign(Clock_Angle=clock_03)
visit_03_df = visit_03_df.assign(B_Perp=b_perp_03)
visit_03_df = visit_03_df.assign(Br=br_03)
visit_03_df = visit_03_df.assign(Bt=bt_03)
visit_03_df = visit_03_df.assign(Bn=bn_03)
visit_03_df = visit_03_df.assign(BS_Location=bs_03)
visit_03_df = visit_03_df.assign(MP_Location=mp_03)
visit_03_df = visit_03_df.assign(LL_Recon_V=LL_03)
visit_03_df = visit_03_df.assign(HL_Recon_V_By_Plus=HL_03_pos)
visit_03_df = visit_03_df.assign(HL_Recon_V_By_Neg=HL_03_neg)

visit_03_df.to_csv(root_folder+'visit_03_data.csv',index=False)


'''
# VISIT 04
'''

visit_04_df = pd.DataFrame()

visit_04_df = visit_04_df.assign(Index=index_04)
visit_04_df = visit_04_df.assign(Ionospheric_Time=iono_time_04)
visit_04_df = visit_04_df.assign(SW_Travel_Time=tt_04)
visit_04_df = visit_04_df.assign(Juno_Detection_Time=juno_time_04)
visit_04_df = visit_04_df.assign(Sw_Pressure=pressure_04)
visit_04_df = visit_04_df.assign(SW_Velocity=velocity_04)
visit_04_df = visit_04_df.assign(Clock_Angle=clock_04)
visit_04_df = visit_04_df.assign(B_Perp=b_perp_04)
visit_04_df = visit_04_df.assign(Br=br_04)
visit_04_df = visit_04_df.assign(Bt=bt_04)
visit_04_df = visit_04_df.assign(Bn=bn_04)
visit_04_df = visit_04_df.assign(BS_Location=bs_04)
visit_04_df = visit_04_df.assign(MP_Location=mp_04)
visit_04_df = visit_04_df.assign(LL_Recon_V=LL_04)
visit_04_df = visit_04_df.assign(HL_Recon_V_By_Plus=HL_04_pos)
visit_04_df = visit_04_df.assign(HL_Recon_V_By_Neg=HL_04_neg)

visit_04_df.to_csv(root_folder+'visit_04_data.csv',index=False)



'''
# VISIT 05
'''

visit_05_df = pd.DataFrame()

visit_05_df = visit_05_df.assign(Index=index_05)
visit_05_df = visit_05_df.assign(Ionospheric_Time=iono_time_05)
visit_05_df = visit_05_df.assign(SW_Travel_Time=tt_05)
visit_05_df = visit_05_df.assign(Juno_Detection_Time=juno_time_05)
visit_05_df = visit_05_df.assign(Sw_Pressure=pressure_05)
visit_05_df = visit_05_df.assign(SW_Velocity=velocity_05)
visit_05_df = visit_05_df.assign(Clock_Angle=clock_05)
visit_05_df = visit_05_df.assign(B_Perp=b_perp_05)
visit_05_df = visit_05_df.assign(Br=br_05)
visit_05_df = visit_05_df.assign(Bt=bt_05)
visit_05_df = visit_05_df.assign(Bn=bn_05)
visit_05_df = visit_05_df.assign(BS_Location=bs_05)
visit_05_df = visit_05_df.assign(MP_Location=mp_05)
visit_05_df = visit_05_df.assign(LL_Recon_V=LL_05)
visit_05_df = visit_05_df.assign(HL_Recon_V_By_Plus=HL_05_pos)
visit_05_df = visit_05_df.assign(HL_Recon_V_By_Neg=HL_05_neg)

visit_05_df.to_csv(root_folder+'isit_05_data.csv',index=False)




'''
# VISIT 08
'''

visit_08_df = pd.DataFrame()

visit_08_df = visit_08_df.assign(Index=index_08)
visit_08_df = visit_08_df.assign(Ionospheric_Time=iono_time_08)
visit_08_df = visit_08_df.assign(SW_Travel_Time=tt_08)
visit_08_df = visit_08_df.assign(Juno_Detection_Time=juno_time_08)
visit_08_df = visit_08_df.assign(Sw_Pressure=pressure_08)
visit_08_df = visit_08_df.assign(SW_Velocity=velocity_08)
visit_08_df = visit_08_df.assign(Clock_Angle=clock_08)
visit_08_df = visit_08_df.assign(B_Perp=b_perp_08)
visit_08_df = visit_08_df.assign(Br=br_08)
visit_08_df = visit_08_df.assign(Bt=bt_08)
visit_08_df = visit_08_df.assign(Bn=bn_08)
visit_08_df = visit_08_df.assign(BS_Location=bs_08)
visit_08_df = visit_08_df.assign(MP_Location=mp_08)
visit_08_df = visit_08_df.assign(LL_Recon_V=LL_08)
visit_08_df = visit_08_df.assign(HL_Recon_V_By_Plus=HL_08_pos)
visit_08_df = visit_08_df.assign(HL_Recon_V_By_Neg=HL_08_neg)

visit_08_df.to_csv(root_folder+'visit_08_data.csv',index=False)


'''
# VISIT 09
'''

visit_09_df = pd.DataFrame()

visit_09_df = visit_09_df.assign(Index=index_09)
visit_09_df = visit_09_df.assign(Ionospheric_Time=iono_time_09)
visit_09_df = visit_09_df.assign(SW_Travel_Time=tt_09)
visit_09_df = visit_09_df.assign(Juno_Detection_Time=juno_time_09)
visit_09_df = visit_09_df.assign(Sw_Pressure=pressure_09)
visit_09_df = visit_09_df.assign(SW_Velocity=velocity_09)
visit_09_df = visit_09_df.assign(Clock_Angle=clock_09)
visit_09_df = visit_09_df.assign(B_Perp=b_perp_09)
visit_09_df = visit_09_df.assign(Br=br_09)
visit_09_df = visit_09_df.assign(Bt=bt_09)
visit_09_df = visit_09_df.assign(Bn=bn_09)
visit_09_df = visit_09_df.assign(BS_Location=bs_09)
visit_09_df = visit_09_df.assign(MP_Location=mp_09)
visit_09_df = visit_09_df.assign(LL_Recon_V=LL_09)
visit_09_df = visit_09_df.assign(HL_Recon_V_By_Plus=HL_09_pos)
visit_09_df = visit_09_df.assign(HL_Recon_V_By_Neg=HL_09_neg)

visit_09_df.to_csv(root_folder+'visit_09_data.csv',index=False)


'''
# VISIT 10
'''

visit_10_df = pd.DataFrame()

visit_10_df = visit_10_df.assign(Index=index_10)
visit_10_df = visit_10_df.assign(Ionospheric_Time=iono_time_10)
visit_10_df = visit_10_df.assign(SW_Travel_Time=tt_10)
visit_10_df = visit_10_df.assign(Juno_Detection_Time=juno_time_10)
visit_10_df = visit_10_df.assign(Sw_Pressure=pressure_10)
visit_10_df = visit_10_df.assign(SW_Velocity=velocity_10)
visit_10_df = visit_10_df.assign(Clock_Angle=clock_10)
visit_10_df = visit_10_df.assign(B_Perp=b_perp_10)
visit_10_df = visit_10_df.assign(Br=br_10)
visit_10_df = visit_10_df.assign(Bt=bt_10)
visit_10_df = visit_10_df.assign(Bn=bn_10)
visit_10_df = visit_10_df.assign(BS_Location=bs_10)
visit_10_df = visit_10_df.assign(MP_Location=mp_10)
visit_10_df = visit_10_df.assign(LL_Recon_V=LL_10)
visit_10_df = visit_10_df.assign(HL_Recon_V_By_Plus=HL_10_pos)
visit_10_df = visit_10_df.assign(HL_Recon_V_By_Neg=HL_10_neg)

visit_10_df.to_csv(root_folder+'visit_10_data.csv',index=False)


'''
# VISIT 11
'''

visit_11_df = pd.DataFrame()

visit_11_df = visit_11_df.assign(Index=index_11)
visit_11_df = visit_11_df.assign(Ionospheric_Time=iono_time_11)
visit_11_df = visit_11_df.assign(SW_Travel_Time=tt_11)
visit_11_df = visit_11_df.assign(Juno_Detection_Time=juno_time_11)
visit_11_df = visit_11_df.assign(Sw_Pressure=pressure_11)
visit_11_df = visit_11_df.assign(SW_Velocity=velocity_11)
visit_11_df = visit_11_df.assign(Clock_Angle=clock_11)
visit_11_df = visit_11_df.assign(B_Perp=b_perp_11)
visit_11_df = visit_11_df.assign(Br=br_11)
visit_11_df = visit_11_df.assign(Bt=bt_11)
visit_11_df = visit_11_df.assign(Bn=bn_11)
visit_11_df = visit_11_df.assign(BS_Location=bs_11)
visit_11_df = visit_11_df.assign(MP_Location=mp_11)
visit_11_df = visit_11_df.assign(LL_Recon_V=LL_11)
visit_11_df = visit_11_df.assign(HL_Recon_V_By_Plus=HL_11_pos)
visit_11_df = visit_11_df.assign(HL_Recon_V_By_Neg=HL_11_neg)

visit_11_df.to_csv(root_folder+'visit_11_data.csv',index=False)


'''
# VISIT 12
'''

visit_12_df = pd.DataFrame()

visit_12_df = visit_12_df.assign(Index=index_12)
visit_12_df = visit_12_df.assign(Ionospheric_Time=iono_time_12)
visit_12_df = visit_12_df.assign(SW_Travel_Time=tt_12)
visit_12_df = visit_12_df.assign(Juno_Detection_Time=juno_time_12)
visit_12_df = visit_12_df.assign(Sw_Pressure=pressure_12)
visit_12_df = visit_12_df.assign(SW_Velocity=velocity_12)
visit_12_df = visit_12_df.assign(Clock_Angle=clock_12)
visit_12_df = visit_12_df.assign(B_Perp=b_perp_12)
visit_12_df = visit_12_df.assign(Br=br_12)
visit_12_df = visit_12_df.assign(Bt=bt_12)
visit_12_df = visit_12_df.assign(Bn=bn_12)
visit_12_df = visit_12_df.assign(BS_Location=bs_12)
visit_12_df = visit_12_df.assign(MP_Location=mp_12)
visit_12_df = visit_12_df.assign(LL_Recon_V=LL_12)
visit_12_df = visit_12_df.assign(HL_Recon_V_By_Plus=HL_12_pos)
visit_12_df = visit_12_df.assign(HL_Recon_V_By_Neg=HL_12_neg)

visit_12_df.to_csv(root_folder+'visit_12_data.csv',index=False)


'''
# VISIT 13
'''

visit_13_df = pd.DataFrame()

visit_13_df = visit_13_df.assign(Index=index_13)
visit_13_df = visit_13_df.assign(Ionospheric_Time=iono_time_13)
visit_13_df = visit_13_df.assign(SW_Travel_Time=tt_13)
visit_13_df = visit_13_df.assign(Juno_Detection_Time=juno_time_13)
visit_13_df = visit_13_df.assign(Sw_Pressure=pressure_13)
visit_13_df = visit_13_df.assign(SW_Velocity=velocity_13)
visit_13_df = visit_13_df.assign(Clock_Angle=clock_13)
visit_13_df = visit_13_df.assign(B_Perp=b_perp_13)
visit_13_df = visit_13_df.assign(Br=br_13)
visit_13_df = visit_13_df.assign(Bt=bt_13)
visit_13_df = visit_13_df.assign(Bn=bn_13)
visit_13_df = visit_13_df.assign(BS_Location=bs_13)
visit_13_df = visit_13_df.assign(MP_Location=mp_13)
visit_13_df = visit_13_df.assign(LL_Recon_V=LL_13)
visit_13_df = visit_13_df.assign(HL_Recon_V_By_Plus=HL_13_pos)
visit_13_df = visit_13_df.assign(HL_Recon_V_By_Neg=HL_13_neg)

visit_13_df.to_csv(root_folder+'visit_13_data.csv',index=False)


'''
# VISIT 15
'''

visit_15_df = pd.DataFrame()

visit_15_df = visit_15_df.assign(Index=index_15)
visit_15_df = visit_15_df.assign(Ionospheric_Time=iono_time_15)
visit_15_df = visit_15_df.assign(SW_Travel_Time=tt_15)
visit_15_df = visit_15_df.assign(Juno_Detection_Time=juno_time_15)
visit_15_df = visit_15_df.assign(Sw_Pressure=pressure_15)
visit_15_df = visit_15_df.assign(SW_Velocity=velocity_15)
visit_15_df = visit_15_df.assign(Clock_Angle=clock_15)
visit_15_df = visit_15_df.assign(B_Perp=b_perp_15)
visit_15_df = visit_15_df.assign(Br=br_15)
visit_15_df = visit_15_df.assign(Bt=bt_15)
visit_15_df = visit_15_df.assign(Bn=bn_15)
visit_15_df = visit_15_df.assign(BS_Location=bs_15)
visit_15_df = visit_15_df.assign(MP_Location=mp_15)
visit_15_df = visit_15_df.assign(LL_Recon_V=LL_15)
visit_15_df = visit_15_df.assign(HL_Recon_V_By_Plus=HL_15_pos)
visit_15_df = visit_15_df.assign(HL_Recon_V_By_Neg=HL_15_neg)

visit_15_df.to_csv(root_folder+'visit_15_data.csv',index=False)


'''
# VISIT 16
'''

visit_16_df = pd.DataFrame()

visit_16_df = visit_16_df.assign(Index=index_16)
visit_16_df = visit_16_df.assign(Ionospheric_Time=iono_time_16)
visit_16_df = visit_16_df.assign(SW_Travel_Time=tt_16)
visit_16_df = visit_16_df.assign(Juno_Detection_Time=juno_time_16)
visit_16_df = visit_16_df.assign(Sw_Pressure=pressure_16)
visit_16_df = visit_16_df.assign(SW_Velocity=velocity_16)
visit_16_df = visit_16_df.assign(Clock_Angle=clock_16)
visit_16_df = visit_16_df.assign(B_Perp=b_perp_16)
visit_16_df = visit_16_df.assign(Br=br_16)
visit_16_df = visit_16_df.assign(Bt=bt_16)
visit_16_df = visit_16_df.assign(Bn=bn_16)
visit_16_df = visit_16_df.assign(BS_Location=bs_16)
visit_16_df = visit_16_df.assign(MP_Location=mp_16)
visit_16_df = visit_16_df.assign(LL_Recon_V=LL_16)
visit_16_df = visit_16_df.assign(HL_Recon_V_By_Plus=HL_16_pos)
visit_16_df = visit_16_df.assign(HL_Recon_V_By_Neg=HL_16_neg)

visit_16_df.to_csv(root_folder+'visit_16_data.csv',index=False)


'''
# VISIT 17
'''

visit_17_df = pd.DataFrame()

visit_17_df = visit_17_df.assign(Index=index_17)
visit_17_df = visit_17_df.assign(Ionospheric_Time=iono_time_17)
visit_17_df = visit_17_df.assign(SW_Travel_Time=tt_17)
visit_17_df = visit_17_df.assign(Juno_Detection_Time=juno_time_17)
visit_17_df = visit_17_df.assign(Sw_Pressure=pressure_17)
visit_17_df = visit_17_df.assign(SW_Velocity=velocity_17)
visit_17_df = visit_17_df.assign(Clock_Angle=clock_17)
visit_17_df = visit_17_df.assign(B_Perp=b_perp_17)
visit_17_df = visit_17_df.assign(Br=br_17)
visit_17_df = visit_17_df.assign(Bt=bt_17)
visit_17_df = visit_17_df.assign(Bn=bn_17)
visit_17_df = visit_17_df.assign(BS_Location=bs_17)
visit_17_df = visit_17_df.assign(MP_Location=mp_17)
visit_17_df = visit_17_df.assign(LL_Recon_V=LL_17)
visit_17_df = visit_17_df.assign(HL_Recon_V_By_Plus=HL_17_pos)
visit_17_df = visit_17_df.assign(HL_Recon_V_By_Neg=HL_17_neg)

visit_17_df.to_csv(root_folder+'visit_17_data.csv',index=False)


'''
# VISIT 18
'''

visit_18_df = pd.DataFrame()

visit_18_df = visit_18_df.assign(Index=index_18)
visit_18_df = visit_18_df.assign(Ionospheric_Time=iono_time_18)
visit_18_df = visit_18_df.assign(SW_Travel_Time=tt_18)
visit_18_df = visit_18_df.assign(Juno_Detection_Time=juno_time_18)
visit_18_df = visit_18_df.assign(Sw_Pressure=pressure_18)
visit_18_df = visit_18_df.assign(SW_Velocity=velocity_18)
visit_18_df = visit_18_df.assign(Clock_Angle=clock_18)
visit_18_df = visit_18_df.assign(B_Perp=b_perp_18)
visit_18_df = visit_18_df.assign(Br=br_18)
visit_18_df = visit_18_df.assign(Bt=bt_18)
visit_18_df = visit_18_df.assign(Bn=bn_18)
visit_18_df = visit_18_df.assign(BS_Location=bs_18)
visit_18_df = visit_18_df.assign(MP_Location=mp_18)
visit_18_df = visit_18_df.assign(LL_Recon_V=LL_18)
visit_18_df = visit_18_df.assign(HL_Recon_V_By_Plus=HL_18_pos)
visit_18_df = visit_18_df.assign(HL_Recon_V_By_Neg=HL_18_neg)

visit_18_df.to_csv(root_folder+'visit_18_data.csv',index=False)


'''
#VISIT 19
'''

visit_19_df = pd.DataFrame()

visit_19_df = visit_19_df.assign(Index=index_19)
visit_19_df = visit_19_df.assign(Ionospheric_Time=iono_time_19)
visit_19_df = visit_19_df.assign(SW_Travel_Time=tt_19)
visit_19_df = visit_19_df.assign(Juno_Detection_Time=juno_time_19)
visit_19_df = visit_19_df.assign(Sw_Pressure=pressure_19)
visit_19_df = visit_19_df.assign(SW_Velocity=velocity_19)
visit_19_df = visit_19_df.assign(Clock_Angle=clock_19)
visit_19_df = visit_19_df.assign(B_Perp=b_perp_19)
visit_19_df = visit_19_df.assign(Br=br_19)
visit_19_df = visit_19_df.assign(Bt=bt_19)
visit_19_df = visit_19_df.assign(Bn=bn_19)
visit_19_df = visit_19_df.assign(BS_Location=bs_19)
visit_19_df = visit_19_df.assign(MP_Location=mp_19)
visit_19_df = visit_19_df.assign(LL_Recon_V=LL_19)
visit_19_df = visit_19_df.assign(HL_Recon_V_By_Plus=HL_19_pos)
visit_19_df = visit_19_df.assign(HL_Recon_V_By_Neg=HL_19_neg)


visit_19_df.to_csv(root_folder+'visit_19_data.csv',index=False)


'''
#VISIT 20
'''

visit_20_df = pd.DataFrame()

visit_20_df = visit_20_df.assign(Index=index_20)
visit_20_df = visit_20_df.assign(Ionospheric_Time=iono_time_20)
visit_20_df = visit_20_df.assign(SW_Travel_Time=tt_20)
visit_20_df = visit_20_df.assign(Juno_Detection_Time=juno_time_20)
visit_20_df = visit_20_df.assign(Sw_Pressure=pressure_20)
visit_20_df = visit_20_df.assign(SW_Velocity=velocity_20)
visit_20_df = visit_20_df.assign(Clock_Angle=clock_20)
visit_20_df = visit_20_df.assign(B_Perp=b_perp_20)
visit_20_df = visit_20_df.assign(Br=br_20)
visit_20_df = visit_20_df.assign(Bt=bt_20)
visit_20_df = visit_20_df.assign(Bn=bn_20)
visit_20_df = visit_20_df.assign(BS_Location=bs_20)
visit_20_df = visit_20_df.assign(MP_Location=mp_20)
visit_20_df = visit_20_df.assign(LL_Recon_V=LL_20)
visit_20_df = visit_20_df.assign(HL_Recon_V_By_Plus=HL_20_pos)
visit_20_df = visit_20_df.assign(HL_Recon_V_By_Neg=HL_20_neg)

visit_20_df.to_csv(root_folder+'visit_20_data.csv',index=False)


'''
#VISIT 21
'''

visit_21_df = pd.DataFrame()

visit_21_df = visit_21_df.assign(Index=index_21)
visit_21_df = visit_21_df.assign(Ionospheric_Time=iono_time_21)
visit_21_df = visit_21_df.assign(SW_Travel_Time=tt_21)
visit_21_df = visit_21_df.assign(Juno_Detection_Time=juno_time_21)
visit_21_df = visit_21_df.assign(Sw_Pressure=pressure_21)
visit_21_df = visit_21_df.assign(SW_Velocity=velocity_21)
visit_21_df = visit_21_df.assign(Clock_Angle=clock_21)
visit_21_df = visit_21_df.assign(B_Perp=b_perp_21)
visit_21_df = visit_21_df.assign(Br=br_21)
visit_21_df = visit_21_df.assign(Bt=bt_21)
visit_21_df = visit_21_df.assign(Bn=bn_21)
visit_21_df = visit_21_df.assign(BS_Location=bs_21)
visit_21_df = visit_21_df.assign(MP_Location=mp_21)
visit_21_df = visit_21_df.assign(LL_Recon_V=LL_21)
visit_21_df = visit_21_df.assign(HL_Recon_V_By_Plus=HL_21_pos)
visit_21_df = visit_21_df.assign(HL_Recon_V_By_Neg=HL_21_neg)

visit_21_df.to_csv(root_folder+'visit_21_data.csv',index=False)


'''
# VISIT 23
'''

visit_23_df = pd.DataFrame()

visit_23_df = visit_23_df.assign(Index=index_23)
visit_23_df = visit_23_df.assign(Ionospheric_Time=iono_time_23)
visit_23_df = visit_23_df.assign(SW_Travel_Time=tt_23)
visit_23_df = visit_23_df.assign(Juno_Detection_Time=juno_time_23)
visit_23_df = visit_23_df.assign(Sw_Pressure=pressure_23)
visit_23_df = visit_23_df.assign(SW_Velocity=velocity_23)
visit_23_df = visit_23_df.assign(Clock_Angle=clock_23)
visit_23_df = visit_23_df.assign(B_Perp=b_perp_23)
visit_23_df = visit_23_df.assign(Br=br_23)
visit_23_df = visit_23_df.assign(Bt=bt_23)
visit_23_df = visit_23_df.assign(Bn=bn_23)
visit_23_df = visit_23_df.assign(BS_Location=bs_23)
visit_23_df = visit_23_df.assign(MP_Location=mp_23)
visit_23_df = visit_23_df.assign(LL_Recon_V=LL_23)
visit_23_df = visit_23_df.assign(HL_Recon_V_By_Plus=HL_23_pos)
visit_23_df = visit_23_df.assign(HL_Recon_V_By_Neg=HL_23_neg)

visit_23_df.to_csv(root_folder+'visit_23_data.csv',index=False)


'''
# VISIT 24
'''

visit_24_df = pd.DataFrame()

visit_24_df = visit_24_df.assign(Index=index_24)
visit_24_df = visit_24_df.assign(Ionospheric_Time=iono_time_24)
visit_24_df = visit_24_df.assign(SW_Travel_Time=tt_24)
visit_24_df = visit_24_df.assign(Juno_Detection_Time=juno_time_24)
visit_24_df = visit_24_df.assign(Sw_Pressure=pressure_24)
visit_24_df = visit_24_df.assign(SW_Velocity=velocity_24)
visit_24_df = visit_24_df.assign(Clock_Angle=clock_24)
visit_24_df = visit_24_df.assign(B_Perp=b_perp_24)
visit_24_df = visit_24_df.assign(Br=br_24)
visit_24_df = visit_24_df.assign(Bt=bt_24)
visit_24_df = visit_24_df.assign(Bn=bn_24)
visit_24_df = visit_24_df.assign(BS_Location=bs_24)
visit_24_df = visit_24_df.assign(MP_Location=mp_24)
visit_24_df = visit_24_df.assign(LL_Recon_V=LL_24)
visit_24_df = visit_24_df.assign(HL_Recon_V_By_Plus=HL_24_pos)
visit_24_df = visit_24_df.assign(HL_Recon_V_By_Neg=HL_24_neg)

visit_24_df.to_csv(root_folder+'visit_24_data.csv',index=False)


'''
# VISIT 25
'''

visit_25_df = pd.DataFrame()

visit_25_df = visit_25_df.assign(Index=index_25)
visit_25_df = visit_25_df.assign(Ionospheric_Time=iono_time_25)
visit_25_df = visit_25_df.assign(SW_Travel_Time=tt_25)
visit_25_df = visit_25_df.assign(Juno_Detection_Time=juno_time_25)
visit_25_df = visit_25_df.assign(Sw_Pressure=pressure_25)
visit_25_df = visit_25_df.assign(SW_Velocity=velocity_25)
visit_25_df = visit_25_df.assign(Clock_Angle=clock_25)
visit_25_df = visit_25_df.assign(B_Perp=b_perp_25)
visit_25_df = visit_25_df.assign(Br=br_25)
visit_25_df = visit_25_df.assign(Bt=bt_25)
visit_25_df = visit_25_df.assign(Bn=bn_25)
visit_25_df = visit_25_df.assign(BS_Location=bs_25)
visit_25_df = visit_25_df.assign(MP_Location=mp_25)
visit_25_df = visit_25_df.assign(LL_Recon_V=LL_25)
visit_25_df = visit_25_df.assign(HL_Recon_V_By_Plus=HL_25_pos)
visit_25_df = visit_25_df.assign(HL_Recon_V_By_Neg=HL_25_neg)

visit_25_df.to_csv(root_folder+'visit_25_data.csv',index=False)


'''
# VISIT 26
'''

visit_26_df = pd.DataFrame()

visit_26_df = visit_26_df.assign(Index=index_26)
visit_26_df = visit_26_df.assign(Ionospheric_Time=iono_time_26)
visit_26_df = visit_26_df.assign(SW_Travel_Time=tt_26)
visit_26_df = visit_26_df.assign(Juno_Detection_Time=juno_time_26)
visit_26_df = visit_26_df.assign(Sw_Pressure=pressure_26)
visit_26_df = visit_26_df.assign(SW_Velocity=velocity_26)
visit_26_df = visit_26_df.assign(Clock_Angle=clock_26)
visit_26_df = visit_26_df.assign(B_Perp=b_perp_26)
visit_26_df = visit_26_df.assign(Br=br_26)
visit_26_df = visit_26_df.assign(Bt=bt_26)
visit_26_df = visit_26_df.assign(Bn=bn_26)
visit_26_df = visit_26_df.assign(BS_Location=bs_26)
visit_26_df = visit_26_df.assign(MP_Location=mp_26)
visit_26_df = visit_26_df.assign(LL_Recon_V=LL_26)
visit_26_df = visit_26_df.assign(HL_Recon_V_By_Plus=HL_26_pos)
visit_26_df = visit_26_df.assign(HL_Recon_V_By_Neg=HL_26_neg)

visit_26_df.to_csv(root_folder+'visit_26_data.csv',index=False)


'''
# VISIT 27
'''

visit_27_df = pd.DataFrame()

visit_27_df = visit_27_df.assign(Index=index_27)
visit_27_df = visit_27_df.assign(Ionospheric_Time=iono_time_27)
visit_27_df = visit_27_df.assign(SW_Travel_Time=tt_27)
visit_27_df = visit_27_df.assign(Juno_Detection_Time=juno_time_27)
visit_27_df = visit_27_df.assign(Sw_Pressure=pressure_27)
visit_27_df = visit_27_df.assign(SW_Velocity=velocity_27)
visit_27_df = visit_27_df.assign(Clock_Angle=clock_27)
visit_27_df = visit_27_df.assign(B_Perp=b_perp_27)
visit_27_df = visit_27_df.assign(Br=br_27)
visit_27_df = visit_27_df.assign(Bt=bt_27)
visit_27_df = visit_27_df.assign(Bn=bn_27)
visit_27_df = visit_27_df.assign(BS_Location=bs_27)
visit_27_df = visit_27_df.assign(MP_Location=mp_27)
visit_27_df = visit_27_df.assign(LL_Recon_V=LL_27)
visit_27_df = visit_27_df.assign(HL_Recon_V_By_Plus=HL_27_pos)
visit_27_df = visit_27_df.assign(HL_Recon_V_By_Neg=HL_27_neg)

visit_27_df.to_csv(root_folder+'visit_27_data.csv',index=False)


'''
# VISIT 28
'''

visit_28_df = pd.DataFrame()

visit_28_df = visit_28_df.assign(Index=index_28)
visit_28_df = visit_28_df.assign(Ionospheric_Time=iono_time_28)
visit_28_df = visit_28_df.assign(SW_Travel_Time=tt_28)
visit_28_df = visit_28_df.assign(Juno_Detection_Time=juno_time_28)
visit_28_df = visit_28_df.assign(Sw_Pressure=pressure_28)
visit_28_df = visit_28_df.assign(SW_Velocity=velocity_28)
visit_28_df = visit_28_df.assign(Clock_Angle=clock_28)
visit_28_df = visit_28_df.assign(B_Perp=b_perp_28)
visit_28_df = visit_28_df.assign(Br=br_28)
visit_28_df = visit_28_df.assign(Bt=bt_28)
visit_28_df = visit_28_df.assign(Bn=bn_28)
visit_28_df = visit_28_df.assign(BS_Location=bs_28)
visit_28_df = visit_28_df.assign(MP_Location=mp_28)
visit_28_df = visit_28_df.assign(LL_Recon_V=LL_28)
visit_28_df = visit_28_df.assign(HL_Recon_V_By_Plus=HL_28_pos)
visit_28_df = visit_28_df.assign(HL_Recon_V_By_Neg=HL_28_neg)

visit_28_df.to_csv(root_folder+'visit_28_data.csv',index=False)


'''
# VISIT 34
'''

visit_34_df = pd.DataFrame()

visit_34_df = visit_24_df.assign(Index=index_24)
visit_34_df = visit_24_df.assign(Ionospheric_Time=iono_time_24)
visit_34_df = visit_24_df.assign(SW_Travel_Time=tt_24)
visit_34_df = visit_24_df.assign(Juno_Detection_Time=juno_time_24)
visit_34_df = visit_24_df.assign(Sw_Pressure=pressure_24)
visit_34_df = visit_24_df.assign(SW_Velocity=velocity_24)
visit_34_df = visit_24_df.assign(Clock_Angle=clock_24)
visit_34_df = visit_24_df.assign(B_Perp=b_perp_24)
visit_34_df = visit_24_df.assign(Br=br_24)
visit_34_df = visit_24_df.assign(Bt=bt_24)
visit_34_df = visit_24_df.assign(Bn=bn_24)
visit_34_df = visit_24_df.assign(BS_Location=bs_24)
visit_34_df = visit_24_df.assign(MP_Location=mp_24)
visit_34_df = visit_24_df.assign(LL_Recon_V=LL_24)
visit_34_df = visit_24_df.assign(HL_Recon_V_By_Plus=HL_24_pos)
visit_34_df = visit_24_df.assign(HL_Recon_V_By_Neg=HL_24_neg)

visit_34_df.to_csv(root_folder+'visit_34_data.csv',index=False)


'''
# VISIT 35
'''

visit_35_df = pd.DataFrame()

visit_35_df = visit_35_df.assign(Index=index_35)
visit_35_df = visit_35_df.assign(Ionospheric_Time=iono_time_35)
visit_35_df = visit_35_df.assign(SW_Travel_Time=tt_35)
visit_35_df = visit_35_df.assign(Juno_Detection_Time=juno_time_35)
visit_35_df = visit_35_df.assign(Sw_Pressure=pressure_35)
visit_35_df = visit_35_df.assign(SW_Velocity=velocity_35)
visit_35_df = visit_35_df.assign(Clock_Angle=clock_35)
visit_35_df = visit_35_df.assign(B_Perp=b_perp_35)
visit_35_df = visit_35_df.assign(Br=br_35)
visit_35_df = visit_35_df.assign(Bt=bt_35)
visit_35_df = visit_35_df.assign(Bn=bn_35)
visit_35_df = visit_35_df.assign(BS_Location=bs_35)
visit_35_df = visit_35_df.assign(MP_Location=mp_35)
visit_35_df = visit_35_df.assign(LL_Recon_V=LL_35)
visit_35_df = visit_35_df.assign(HL_Recon_V_By_Plus=HL_35_pos)
visit_35_df = visit_35_df.assign(HL_Recon_V_By_Neg=HL_35_neg)

visit_35_df.to_csv(root_folder+'visit_35_data.csv',index=False)
