#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 10:22:41 2023

@author: hannah
"""

# code to calculate clock angle and perpendicular componenent of solar wind

# file to process pds data
def clock_angle_calculator(Bn,Bt,check,plotting):
    # import relevant modules
    import numpy as np
    import matplotlib.pyplot as plt

    # empty grid for clock angle
    clock_angle = []
    clock_angle_range = []
    # calculate clock angle for each moving avg
    for j in range(len(Bn)):
    # use arctan2 electric boogaloo
    # 0 is Bz+, +/-180 is Bz-, +90 is By+, -90 is By-
        theta_c = np.arctan2(Bt[j],Bn[j])
        theta_c_deg = np.degrees(theta_c)
        
        c_angle_min = theta_c_deg - 10.25
        c_angle_max = theta_c_deg + 10.25
        
        '''
        need to finish sorting this but implement it as an error bar?
        '''
        
        c_angle_range = c_angle_max - c_angle_min
        
        
        # can 'alter' to a 0 - 360 clock by +360 to -ve data points
        # ie 0 is +Bz, +90 is By+, 180 is -Bz and 270 is -By
        # if theta_c_deg < 0.0:
        #     theta_c_deg += 360
        # save clock angle to array of clock angles
        clock_angle.append(theta_c_deg)
        clock_angle_range.append(c_angle_range)
    
    return clock_angle, clock_angle_range
        
    
    # check clock angle values are within bounds   
    if check == 'yes':
        print(np.max(clock_angle))
        print(np.min(clock_angle))
        
    if plotting == 'yes':
        plt.hist(clock_angle)
        plt.xlabel('Clock Angle ($^o$)')
        plt.ylabel('Counts')
        
    return clock_angle

def B_perp_calculator(Bn,Bt,check):
    import numpy as np
    
    B_perp = []
    
    for k in range(len(Bn)):
        square = (Bt[k]**2) + (Bn[k]**2)
        root = np.sqrt(square)
        B_perp.append(root)
        
    if check == 'yes':
        print(np.max(B_perp))
        print(np.min(B_perp))
    
    return B_perp
