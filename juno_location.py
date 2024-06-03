"""
Created on Mon Oct  9 12:55:01 2023

@author: hannah
"""
import matplotlib.pyplot as plt
import numpy as np
import spiceypy as sp

import pandas as pd

# read in solar wind data
sw_df = pd.read_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/solar_wind_data.csv',delimiter=',')
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

# convert to ephermeous time
et = np.linspace(sp.str2et("2016-139::00:02:30"),sp.str2et("2016-155::00:59:30"),54785)
# also convert to UTC
UTC = sp.et2utc(et, 'ISOD', 6)

# spk - position
# 599 = jupiter, 5th planet, 99 = centre of jupes
# owlt = one way light time
pos, owlt = sp.spkpos("Juno", et, "JUNO_JSO","none","599")

# ,0 = x, 1 = y, 2 = zfor pos
marker_size=1
plt.scatter(pos[:,0]/71492, pos[:,1]/71492,marker='.',s=marker_size)
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

juno_pos.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/juno_position_df.csv',index=False)
