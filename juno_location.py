"""
Created on Mon Oct  9 12:55:01 2023

@author: hannah

this script calculates the location of juno relative to jupiter using solar wind data
"""
import matplotlib.pyplot as plt
import numpy as np
import spiceypy as sp

import pandas as pd

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'

# read in solar wind data
sw_df = pd.read_csv(root_folder+'solar_wind_data.csv',delimiter=',')
# replace T in dataframe w/ :: to make compatible w/ spiceypy
sw_df[['UTC']] = sw_df[['UTC']].replace('T','::', regex=True)# help(sp.str2et)

'''
kernals section
'''
# leap seconds kernal
sp.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")
# this one is for iau - planetary constants
sp.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has co-oridinate systems in it
sp.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")
# juno location
sp.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160522_160729_160909.bsp")
sp.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/spk_rec_160312_160522_160614.bsp")

#et = np.linspace(sp.str2et("2016-180::00:00:00"),sp.str2et("2016-181::00:00:00"))

# convert to ephermeous time for main bulk dataframe
et = np.linspace(sp.str2et("2016-139::00:02:30"),sp.str2et("2016-155::00:59:30"),54785)

# et time for doy 155 - 60s avg
et_155 = np.linspace(sp.str2et("2016-06-03 01:00:08.94500"),sp.str2et("2016-06-03 23:59:09.91000"),1380)

# et time for doy 166 - 158
et_156_plus = np.linspace(sp.str2et("2016-06-04 00:07:38.91300"),sp.str2et("2016-06-06 23:57:41.93000"),432)

# et time for doy 175/176
et_175 = np.linspace(sp.str2et("2016-06-23 00:00:27.082"),sp.str2et("2016-06-24 08:15:28.441"),3871)

# also convert to UTC
UTC = sp.et2utc(et, 'ISOD', 6)
UTC_155 = sp.et2utc(et_155,'ISOD',6)
UTC_156_plus = sp.et2utc(et_156_plus,'ISOD',6)
UTC_175 = sp.et2utc(et_175,'ISOD',6)

# spk - position
# 599 = jupiter, 5th planet, 99 = centre of jupes
# owlt = one way light time
pos, owlt = sp.spkpos("Juno", et, "JUNO_JSO","none","599")

pos_155, owlt_155 = sp.spkpos("Juno", et_155, "JUNO_JSO","none","599")
pos_156_plus, owlt_155_plus = sp.spkpos("Juno", et_156_plus, "JUNO_JSO","none","599")
pos_175, owlt_175 = sp.spkpos("Juno", et_175, "JUNO_JSO","none","599")

# ,0 = x, 1 = y, 2 = zfor pos
marker_size=1
plt.scatter(pos[:,0]/71492, pos[:,1]/71492,marker='.',s=marker_size)
plt.show()

plt.scatter(pos_155[:,0]/71492, pos_155[:,1]/71492,marker='.',s=marker_size)
plt.show()

marker_size=1
plt.scatter(pos_156_plus[:,0]/71492, pos_156_plus[:,1]/71492,marker='.',s=marker_size)
plt.show()

marker_size=1
plt.scatter(pos_175[:,0]/71492, pos_175[:,1]/71492,marker='.',s=marker_size)
plt.show()


# want to write the distances in x, y and z along with time into a dataframe
juno_pos = pd.DataFrame()
juno_pos = juno_pos.assign(UTC=UTC)
juno_pos = juno_pos.assign(ET=et)
juno_pos = juno_pos.assign(X=pos[:,0])
juno_pos = juno_pos.assign(XRJ=(pos[:,0]/71492))
juno_pos = juno_pos.assign(Y=pos[:,1])
juno_pos = juno_pos.assign(YRJ=(pos[:,1]/71492))
juno_pos = juno_pos.assign(Z=pos[:,2])
juno_pos = juno_pos.assign(ZRJ=(pos[:,2]/71492))

juno_pos.to_csv(root_folder+'juno_position_df.csv',index=False)


# want to write the distances in x, y and z along with time into a dataframe
juno_pos_155 = pd.DataFrame()
juno_pos_155 = juno_pos_155.assign(UTC=UTC_155)
juno_pos_155 = juno_pos_155.assign(ET=et_155)
juno_pos_155 = juno_pos_155.assign(X=pos_155[:,0])
juno_pos_155 = juno_pos_155.assign(XRJ=(pos_155[:,0]/71492))
juno_pos_155 = juno_pos_155.assign(Y=pos_155[:,1])
juno_pos_155 = juno_pos_155.assign(YRJ=(pos_155[:,1]/71492))
juno_pos_155 = juno_pos_155.assign(Z=pos_155[:,2])
juno_pos_155 = juno_pos_155.assign(ZRJ=(pos_155[:,2]/71492))

juno_pos_155.to_csv(root_folder+'juno_position_155_df.csv',index=False)


# want to write the distances in x, y and z along with time into a dataframe
juno_pos_156_plus = pd.DataFrame()
juno_pos_156_plus = juno_pos_156_plus.assign(UTC=UTC_156_plus)
juno_pos_156_plus = juno_pos_156_plus.assign(ET=et_156_plus)
juno_pos_156_plus = juno_pos_156_plus.assign(X=pos_156_plus[:,0])
juno_pos_156_plus = juno_pos_156_plus.assign(XRJ=(pos_156_plus[:,0]/71492))
juno_pos_156_plus = juno_pos_156_plus.assign(Y=pos_156_plus[:,1])
juno_pos_156_plus = juno_pos_156_plus.assign(YRJ=(pos_156_plus[:,1]/71492))
juno_pos_156_plus = juno_pos_156_plus.assign(Z=pos_156_plus[:,2])
juno_pos_156_plus = juno_pos_156_plus.assign(ZRJ=(pos_156_plus[:,2]/71492))

juno_pos_156_plus.to_csv(root_folder+'juno_position_156_plus_df.csv',index=False)


# want to write the distances in x, y and z along with time into a dataframe
juno_pos_175 = pd.DataFrame()
juno_pos_175 = juno_pos_175.assign(UTC=UTC_175)
juno_pos_175 = juno_pos_175.assign(ET=et_175)
juno_pos_175 = juno_pos_175.assign(X=pos_175[:,0])
juno_pos_175 = juno_pos_175.assign(XRJ=(pos_175[:,0]/71492))
juno_pos_175 = juno_pos_175.assign(Y=pos_175[:,1])
juno_pos_175 = juno_pos_175.assign(YRJ=(pos_175[:,1]/71492))
juno_pos_175 = juno_pos_175.assign(Z=pos_175[:,2])
juno_pos_175 = juno_pos_175.assign(ZRJ=(pos_175[:,2]/71492))

juno_pos_175.to_csv(root_folder+'juno_position_175_df.csv',index=False)
juno_pos = juno_pos.assign(ZRJ=(pos[:,2]/71492))
