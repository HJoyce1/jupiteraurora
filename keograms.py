#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 18:12:46 2026

@author: hannah

makes auroral keograms to compare with solar wind data 
keogram code adapted from work by Leah Claire, Lancaster University

"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import pandas as pd
from dateutil.parser import parse
import datetime as dt
import matplotlib.patheffects as pe
from matplotlib.path import Path
import spiceypy as spice
import collections
from matplotlib import path
import scipy.constants as c
# import collections.abc as collections
from scipy import signal
from tqdm import tqdm
import glob
import datetime
import matplotlib.dates as mdates
import os
from datetime import timedelta
import pandas as pd

from mpl_toolkits.axes_grid1 import make_axes_locatable

visit_list = ['34'] # start with 18 in morning
visit='34' # 01 produces error, remember to go back to it - fixed!

root = '/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/'


def total_mag(Bx,By,Bz):
    Btot = np.sqrt((Bx**2)+(By**2)+(Bz**2))
    return Btot

# data for each visit
visit_02_data = pd.read_csv(root+'visit_02_data_aug.csv', delimiter=',')
visit_03_data = pd.read_csv(root+'visit_03_data_aug.csv', delimiter=',')
visit_04_data = pd.read_csv(root+'visit_04_data_aug.csv', delimiter=',')
visit_05_data = pd.read_csv(root+'visit_05_data_aug.csv', delimiter=',')
visit_08_data = pd.read_csv(root+'visit_08_data_aug.csv', delimiter=',')
visit_09_data = pd.read_csv(root+'visit_09_data_aug.csv', delimiter=',')
visit_10_data = pd.read_csv(root+'visit_10_data_aug.csv', delimiter=',')
visit_11_data = pd.read_csv(root+'visit_11_data_aug.csv', delimiter=',')
visit_16_data = pd.read_csv(root+'visit_16_data_aug.csv', delimiter=',')
visit_17_data = pd.read_csv(root+'visit_17_data_aug.csv', delimiter=',')
visit_18_data = pd.read_csv(root+'visit_18_data_aug.csv', delimiter=',')
visit_19_data = pd.read_csv(root+'visit_19_data_aug.csv', delimiter=',')
visit_20_data = pd.read_csv(root+'visit_20_data_aug.csv', delimiter=',')
visit_21_data = pd.read_csv(root+'visit_21_data_aug.csv', delimiter=',')
visit_24_data = pd.read_csv(root+'visit_24_data_aug.csv', delimiter=',')
visit_25_data = pd.read_csv(root+'visit_25_data_aug.csv', delimiter=',')
visit_27_data = pd.read_csv(root+'visit_27_data_aug.csv', delimiter=',')
visit_28_data = pd.read_csv(root+'visit_28_data_aug.csv', delimiter=',')
visit_34_data = pd.read_csv(root+'visit_34_data_aug.csv', delimiter=',')
visit_35_data = pd.read_csv(root+'visit_35_data_aug.csv', delimiter=',')


# # high cml - 01, 12, 15, 26
visit_01_data = pd.read_csv(root+'visit_01_data_aug.csv', delimiter=',')
visit_12_data = pd.read_csv(root+'visit_12_data_aug.csv', delimiter=',')
visit_15_data = pd.read_csv(root+'visit_15_data_aug.csv', delimiter=',')
visit_26_data = pd.read_csv(root+'visit_26_data_aug.csv', delimiter=',')

bx_01 = visit_01_data['Bx'].to_numpy()
by_01 = visit_01_data['By'].to_numpy()
bz_01 = visit_01_data['Bz'].to_numpy()

btot_01 = total_mag(bx_01,by_01,bz_01)

bx_12 = visit_12_data['Bx'].to_numpy()
by_12 = visit_12_data['By'].to_numpy()
bz_12 = visit_12_data['Bz'].to_numpy()

btot_12 = total_mag(bx_12,by_12,bz_12)

bx_15 = visit_15_data['Bx'].to_numpy()
by_15 = visit_15_data['By'].to_numpy()
bz_15 = visit_15_data['Bz'].to_numpy()

btot_15 = total_mag(bx_15,by_15,bz_15)

bx_26 = visit_26_data['Bx'].to_numpy()
by_26 = visit_26_data['By'].to_numpy()
bz_26 = visit_26_data['Bz'].to_numpy()

btot_26 = total_mag(bx_26,by_26,bz_26)



bx_02 = visit_02_data['Bx'].to_numpy()
by_02 = visit_02_data['By'].to_numpy()
bz_02 = visit_02_data['Bz'].to_numpy()

btot_02 = total_mag(bx_02,by_02,bz_02)


bx_03 = visit_03_data['Bx'].to_numpy()
by_03 = visit_03_data['By'].to_numpy()
bz_03 = visit_03_data['Bz'].to_numpy()

btot_03 = total_mag(bx_03,by_03,bz_03)


bx_04 = visit_04_data['Bx'].to_numpy()
by_04 = visit_04_data['By'].to_numpy()
bz_04 = visit_04_data['Bz'].to_numpy()

btot_04 = total_mag(bx_04,by_04,bz_04)


bx_05 = visit_05_data['Bx'].to_numpy()
by_05 = visit_05_data['By'].to_numpy()
bz_05 = visit_05_data['Bz'].to_numpy()

btot_05 = total_mag(bx_05,by_05,bz_05)


bx_08 = visit_08_data['Bx'].to_numpy()
by_08 = visit_08_data['By'].to_numpy()
bz_08 = visit_08_data['Bz'].to_numpy()

btot_08 = total_mag(bx_08,by_08,bz_08)


bx_09 = visit_09_data['Bx'].to_numpy()
by_09 = visit_09_data['By'].to_numpy()
bz_09 = visit_09_data['Bz'].to_numpy()

btot_09 = total_mag(bx_09,by_09,bz_09)


bx_10 = visit_10_data['Bx'].to_numpy()
by_10 = visit_10_data['By'].to_numpy()
bz_10 = visit_10_data['Bz'].to_numpy()

btot_10 = total_mag(bx_10,by_10,bz_10)


bx_11 = visit_11_data['Bx'].to_numpy()
by_11 = visit_11_data['By'].to_numpy()
bz_11 = visit_11_data['Bz'].to_numpy()

btot_11 = total_mag(bx_11,by_11,bz_11)


bx_16 = visit_16_data['Bx'].to_numpy()
by_16 = visit_16_data['By'].to_numpy()
bz_16 = visit_16_data['Bz'].to_numpy()

btot_16 = total_mag(bx_16,by_16,bz_16)


bx_17 = visit_17_data['Bx'].to_numpy()
by_17 = visit_17_data['By'].to_numpy()
bz_17 = visit_17_data['Bz'].to_numpy()

btot_17 = total_mag(bx_17,by_17,bz_17)


bx_18 = visit_18_data['Bx'].to_numpy()
by_18 = visit_18_data['By'].to_numpy()
bz_18 = visit_18_data['Bz'].to_numpy()

btot_18 = total_mag(bx_18,by_18,bz_18)


bx_19 = visit_19_data['Bx'].to_numpy()
by_19 = visit_19_data['By'].to_numpy()
bz_19 = visit_19_data['Bz'].to_numpy()

btot_19 = total_mag(bx_19,by_19,bz_19)


bx_20 = visit_20_data['Bx'].to_numpy()
by_20 = visit_20_data['By'].to_numpy()
bz_20 = visit_20_data['Bz'].to_numpy()

btot_20 = total_mag(bx_20,by_20,bz_20)


bx_21 = visit_21_data['Bx'].to_numpy()
by_21 = visit_21_data['By'].to_numpy()
bz_21 = visit_21_data['Bz'].to_numpy()

btot_21 = total_mag(bx_21,by_21,bz_21)


bx_24 = visit_24_data['Bx'].to_numpy()
by_24 = visit_24_data['By'].to_numpy()
bz_24 = visit_24_data['Bz'].to_numpy()

btot_24 = total_mag(bx_24,by_24,bz_24)


bx_25 = visit_25_data['Bx'].to_numpy()
by_25 = visit_25_data['By'].to_numpy()
bz_25 = visit_25_data['Bz'].to_numpy()

btot_25 = total_mag(bx_25,by_25,bz_25)


bx_27 = visit_27_data['Bx'].to_numpy()
by_27 = visit_27_data['By'].to_numpy()
bz_27 = visit_27_data['Bz'].to_numpy()

btot_27 = total_mag(bx_27,by_27,bz_27)


bx_28 = visit_28_data['Bx'].to_numpy()
by_28 = visit_28_data['By'].to_numpy()
bz_28 = visit_28_data['Bz'].to_numpy()

btot_28 = total_mag(bx_28,by_28,bz_28)


bx_34 = visit_34_data['Bx'].to_numpy()
by_34 = visit_34_data['By'].to_numpy()
bz_34 = visit_34_data['Bz'].to_numpy()

btot_34 = total_mag(bx_34,by_34,bz_34)


bx_35 = visit_35_data['Bx'].to_numpy()
by_35 = visit_35_data['By'].to_numpy()
bz_35 = visit_35_data['Bz'].to_numpy()

btot_35 = total_mag(bx_35,by_35,bz_35)


# extract clock angle for each visit
clock_02 = visit_02_data['Clock_Angle'].to_numpy()
clock_03 = visit_03_data['Clock_Angle'].to_numpy()
clock_04 = visit_04_data['Clock_Angle'].to_numpy()
clock_05 = visit_05_data['Clock_Angle'].to_numpy()
clock_08 = visit_08_data['Clock_Angle'].to_numpy()
clock_09 = visit_09_data['Clock_Angle'].to_numpy()
clock_10 = visit_10_data['Clock_Angle'].to_numpy()
clock_11 = visit_11_data['Clock_Angle'].to_numpy()
clock_16 = visit_16_data['Clock_Angle'].to_numpy()
clock_17 = visit_17_data['Clock_Angle'].to_numpy()
clock_18 = visit_18_data['Clock_Angle'].to_numpy()
clock_19 = visit_19_data['Clock_Angle'].to_numpy()
clock_20 = visit_20_data['Clock_Angle'].to_numpy()
clock_21 = visit_21_data['Clock_Angle'].to_numpy()
clock_24 = visit_24_data['Clock_Angle'].to_numpy()
clock_25 = visit_25_data['Clock_Angle'].to_numpy()
clock_27 = visit_27_data['Clock_Angle'].to_numpy()
clock_28 = visit_28_data['Clock_Angle'].to_numpy()
clock_34 = visit_34_data['Clock_Angle'].to_numpy()
clock_35 = visit_35_data['Clock_Angle'].to_numpy()

# clock_13 = visit_13_data['Clock_Angle'].to_numpy()

# high cml
clock_01 = visit_01_data['Clock_Angle'].to_numpy()
clock_12 = visit_12_data['Clock_Angle'].to_numpy()
clock_15 = visit_15_data['Clock_Angle'].to_numpy()
clock_26 = visit_26_data['Clock_Angle'].to_numpy()

# ------------- pressures ----------------

# extract pressure for each visit
pressure_02 = visit_02_data['Sw_Pressure'].to_numpy()
pressure_03 = visit_03_data['Sw_Pressure'].to_numpy()
pressure_04 = visit_04_data['Sw_Pressure'].to_numpy()
pressure_05 = visit_05_data['Sw_Pressure'].to_numpy()
pressure_08 = visit_08_data['Sw_Pressure'].to_numpy()
pressure_09 = visit_09_data['Sw_Pressure'].to_numpy()
pressure_10 = visit_10_data['Sw_Pressure'].to_numpy()
pressure_11 = visit_11_data['Sw_Pressure'].to_numpy()
pressure_16 = visit_16_data['Sw_Pressure'].to_numpy()
pressure_17 = visit_17_data['Sw_Pressure'].to_numpy()
pressure_18 = visit_18_data['Sw_Pressure'].to_numpy()
pressure_19 = visit_19_data['Sw_Pressure'].to_numpy()
pressure_20 = visit_20_data['Sw_Pressure'].to_numpy()
pressure_21 = visit_21_data['Sw_Pressure'].to_numpy()
pressure_24 = visit_24_data['Sw_Pressure'].to_numpy()
pressure_25 = visit_25_data['Sw_Pressure'].to_numpy()
pressure_27 = visit_27_data['Sw_Pressure'].to_numpy()
pressure_28 = visit_28_data['Sw_Pressure'].to_numpy()
pressure_34 = visit_34_data['Sw_Pressure'].to_numpy()
pressure_35 = visit_35_data['Sw_Pressure'].to_numpy()

#pressure_13 = visit_13_data['Sw_Pressure'].to_numpy()

# high cml
pressure_01 = visit_01_data['Sw_Pressure'].to_numpy()
pressure_12 = visit_12_data['Sw_Pressure'].to_numpy()
pressure_15 = visit_15_data['Sw_Pressure'].to_numpy()
pressure_26 = visit_26_data['Sw_Pressure'].to_numpy()

# ------------ perpendicular magnetic field -------------

b_perp_02 = visit_02_data['B_Perp'].to_numpy()
b_perp_03 = visit_03_data['B_Perp'].to_numpy()
b_perp_04 = visit_04_data['B_Perp'].to_numpy()
b_perp_05 = visit_05_data['B_Perp'].to_numpy()
b_perp_08 = visit_08_data['B_Perp'].to_numpy()
b_perp_09 = visit_09_data['B_Perp'].to_numpy()
b_perp_10 = visit_10_data['B_Perp'].to_numpy()
b_perp_11 = visit_11_data['B_Perp'].to_numpy()
b_perp_16 = visit_16_data['B_Perp'].to_numpy()
b_perp_17 = visit_17_data['B_Perp'].to_numpy()
b_perp_18 = visit_18_data['B_Perp'].to_numpy()
b_perp_19 = visit_19_data['B_Perp'].to_numpy()
b_perp_20 = visit_20_data['B_Perp'].to_numpy()
b_perp_21 = visit_21_data['B_Perp'].to_numpy()
b_perp_24 = visit_24_data['B_Perp'].to_numpy()
b_perp_25 = visit_25_data['B_Perp'].to_numpy()
b_perp_27 = visit_27_data['B_Perp'].to_numpy()
b_perp_28 = visit_28_data['B_Perp'].to_numpy()
b_perp_34 = visit_34_data['B_Perp'].to_numpy()
b_perp_35 = visit_35_data['B_Perp'].to_numpy()

# high cml
b_perp_01 = visit_01_data['B_Perp'].to_numpy()
b_perp_12 = visit_12_data['B_Perp'].to_numpy()
b_perp_15 = visit_15_data['B_Perp'].to_numpy()
b_perp_26 = visit_26_data['B_Perp'].to_numpy()

# leap seconds kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/naif0012.tls")
# this one is for iau - planetary constants
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/pck00010.tpc")
# this one has coordinate systems in it
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/juno_v12.tf")
# kernal for planet locations
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/de436s.bsp")
# hst kernal
spice.furnsh("/Users/hannah/OneDrive - Lancaster University/aurora/hst.bsp")

# ------------------------------------------------------------------------------
# Nichols constants:
au_to_km = 1.495978707e8
dpr = 180. / np.pi
autokm = 1.4959787066E8

# Nichols spice stuff:
planets = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn',
           'Uranus', 'Neptune', 'Pluto', 'Vulcan']
inx = planets.index('Jupiter')
naifobj = 99 + (inx + 1) * 100
frame = 'IAU_' + 'Jupiter'
radii = spice.bodvcd(naifobj, 'RADII', 3)
#print(radii)
rpeqkm = radii[1][0]
rpplkm = radii[1][2]
oblt = 1. - rpplkm / rpeqkm
obt = oblt

# Hardwire epoch
epoch = 'J2000'
corr = 'LT'

# -----------------------------------------------------------------------------

root_folder = '/Users/hannah/OneDrive - Lancaster University/aurora/'
root_saves = '/Users/hannah/OneDrive - Lancaster University/aurora/keograms/'

#visit_list = ['02','03','04','05','08','09','10','11','16','17','18','19','20','21','24','25','27','28','34','35']#['04','05','08','09','10','11','12','15','16','17','18','19','20','21']#'04','05','08','09','10','11','12','15','16','17','18','19','20','21'] #'04','05','08','09','10','11','12','13','15','16','17','18','19','20','21'
#visit_list = ['01','12','15','26']

def ignore_nan_counter(data):
    count = 0
    for i in data:
        if not np.isnan(i):
            count+=1
    return count

# This function turns python datetime into spice/ET time - From time_conversions.py
# Input: datetime 1-D array/list or single value
def datetime2et(pytimes):
    isscalar = False
    if not isinstance(pytimes, collections.Iterable):
        isscalar = True
        pytimes = [pytimes]
    utctimes = np.array([dt.strftime(iii, '%Y-%m-%d, %H:%M:%S.%f') for iii in pytimes])
    ettimes = np.array([spice.utc2et(iii) for iii in utctimes])
    if isscalar:
        return np.ndarray.item(ettimes)
    else:
        return ettimes


class shape():
    
    def __init__(self, verticies):
        
        """
        Creates a shape bounded by input verticies
        
        Params:
        --------------------
        
        verticies: list
            list of verticies of shape e.g: [[lat1, lon1], [lat2, lon2], 
                                              [lat3, lon3]]
        """
        
        self.verticies = verticies
        
        # seperate lat and lon from verticies (useful for plotting verticies)
        self.lat_verts = [vert[0] for vert in self.verticies]
        self.lon_verts = [vert[1] for vert in self.verticies]
        # re-add first point so lines will connect when plotting
        self.lat_verts.append(self.lat_verts[0])
        self.lon_verts.append(self.lon_verts[0])
        
        return
    
    def insert_points(self, lats, lons, data):
        
        """
        Saves pixels whoes latitudes/longitudes are contained withing
        this shapes verticies
        
        Params
        ------
        
        lats: list or np.ndarray
            latitude location of each intensity point
        
        lons: list or np.ndarray
            longitude location of each intensity point
            
        intensities: list or np.ndarray
            intensities of each point
            
        shape of lats, lons and intensities must be equal
        """

        # check shapes are equal
        if lats.shape != data.shape:
            raise ValueError("lats and intensities must be same shape")
        if lons.shape != data.shape:
            raise ValueError("lons and intensities must be same shape")
        if lats.shape != lons.shape:
            raise ValueError("lats and lons must be same shape")
        
        # arange coordinates for use in path function
        coords = np.vstack([lats.flatten(), lons.flatten()]).T
        #breakpoint()

        # check if points within polygon
        poly = path.Path(self.verticies)
        inside_mask = poly.contains_points(coords) # this is already making the mask 
        
        # store intensities to self - just in case
        self.data = data.flatten()[inside_mask]
        self.lats = lats.flatten()[inside_mask]
        self.lons = lons.flatten()[inside_mask]
        
        return inside_mask.reshape(720,1440)
    
def apply_mask(region,image,visit,plotting): # needs region, im_flip and visit, 
    # this bit establishes what type of visit we have and what coordinates to use  

    if  visit == '01' or visit == '12' or visit == '15' or visit == '26':
        print('Using High CML Regions')
        # dusk_active_region = [[20,192.25],[30,200],[20,220],[15,230],[15,220]]
        # swirl_region = [[6.75,111],[17,155],[18,185],[10,190]]
        # noon_active_region = [[18,154],[24,154],[28,192],[22,192]]
        
        dusk_active_region = [[20,192.5],[30,200],[20,220],[15,230],[15,220]]
        swirl_region = [[7,112],[17,155],[18,185],[10,190]] #[7,112],[17,155],[18,185],[10,190]
        noon_active_region = [[18,154],[24,154],[28,192],[22,192]]
        
        test_region = [[15,210],[15,230],[25,230],[25,210]]#10,179.75],[10,200],[20,200],[20,180]
        
        test = shape(test_region)
        
        # dusk_active_region = [[20,167.25],[30,160],[20,140],[15,130],[15,140]]
        # swirl_region = [[2.75,247],[17,205],[18,175],[4,170]]
        # noon_active_region = [[18,206],[24,206],[28,168],[22,168]]


    else:
        print('Using Standard CML Regions')
        # noon_active_region = [[23,174.75],[29,174.75],[29,202],[23,202]]
        # dusk_active_region = [[3,184.75],[22,184.75],[22,205],[3,205]]
        # swirl_region = [[6.75,99],[17,143],[18,173],[10,178]]
        
        noon_active_region = [[23,175],[29,175],[29,202],[23,202]]
        dusk_active_region = [[2,185],[22,185],[22,205],[3,205]]
        swirl_region = [[7,100],[17,143],[18,173],[10,178]] #[7,100],[17,143],[18,173],[10,178]
        
        test_region = [[15,210],[15,230],[25,230],[25,210]]
        
        test = shape(test_region)
    
        # swirl_region = [[2.75,261],[17,217],[18,175],[4,182]]
        # noon_active_region = [[23,184.75],[29,184.75],[29,158],[23,158]]
        # dusk_active_region = [[3,174.75],[22,174.75],[22,155],[3,155]]
        
    
    noon_active = shape(noon_active_region)
    dusk_active = shape(dusk_active_region)
    swirl = shape(swirl_region)

    lats1 = np.arange(0,180,0.25)
    lons1 = np.arange(0,360,0.25)
    
    # make 2D grid of lat/lons
    llons, llats = np.meshgrid(lons1, lats1)
    
    image_extract = image

    lons = np.arange(0,1440,1)
    lats = np.arange(0,160,1)   # 0-40 degrees colat in image pixel res
    
    if region == 'test' or region == 'Test':
        print('Masking Test Region')
        mask_test = test.insert_points(llats, llons, image)
        #mask_noon160 = mask_noon[0:160,:]
        
        #breakpoint()
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=30.)
            
            coord1 = [40,720]
            coord2 = [40,800]
            coord3 = [80,800]
            coord4 = [80,720]
            
            # Overplot the region of interest, e.g. a lat-lon box here:
            plt.plot([coord1[1],coord2[1],coord3[1],coord4[1],coord1[1]],  # corner A repeated to
                     [coord1[0],coord2[0],coord3[0],coord4[0],coord1[0]],  # close the box.
                     color='red',linewidth=3.)
            plt.show()
            
        if plotting == 'yes':
            plt.figure()
            plt.imshow(mask_test, origin='lower')  # zoom to see tiny squares! #roi mask
            plt.title('Mask in image space')
            plt.show()
            
        image_extract[mask_test==False] = np.nan
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=30.)
            plt.title('Image masked off by the ROI')
            plt.show()
        
        im_full          = np.zeros((720,1440))        # new array full of nans
        im_full[0:160,:] = image_extract[0:160,:] #roi_im_full[0:160,:]
        
        # -------------
        if plotting == 'yes':
            plt.figure()
            plt.imshow(im_full, origin='lower')
            plt.title('ROI intensities in full image space')
            plt.show()
    
    elif region == 'noon' or region == 'Noon':
        print('Masking Noon Active Region')
        mask_noon = noon_active.insert_points(llats, llons, image)
        #mask_noon160 = mask_noon[0:160,:]
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons1,lats1[0:160],image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=30.)
            
            coord1 = noon_active_region[0]
            coord2 = noon_active_region[1]
            coord3 = noon_active_region[2]
            coord4 = noon_active_region[3]
            
            # Overplot the region of interest, e.g. a lat-lon box here:
            plt.plot([coord1[1],coord2[1],coord3[1],coord4[1],coord1[1]],  # corner A repeated to
                     [coord1[0],coord2[0],coord3[0],coord4[0],coord1[0]],  # close the box.
                     color='red',linewidth=3.)
            plt.show()
            
        if plotting == 'yes':
            plt.figure()
            plt.imshow(mask_noon, origin='lower')  # zoom to see tiny squares! #roi mask
            plt.title('Mask in image space')
            plt.show()
            
        image_extract[mask_noon==False] = np.nan
        
        noon_px = noon_active.data
        
        tot_size = len(noon_px)
        num_px = ignore_nan_counter(noon_px)
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=1000.)
            plt.title('Image masked off by the ROI')
            plt.show()
        
        im_full          = np.zeros((720,1440))        # new array full of nans
        im_full[0:160,:] = image_extract[0:160,:] #roi_im_full[0:160,:]
        
        # -------------
        if plotting == 'yes':
            plt.figure()
            plt.imshow(im_full, origin='lower')
            plt.title('ROI intensities in full image space')
            plt.show()
            

    elif region == 'dusk' or region == 'Dusk':
        print('Masking Dusk Active Region')
        mask_dusk = dusk_active.insert_points(llats, llons, image)
        #mask_dusk160 = mask_dusk[0:160,:]
            
        if plotting == 'yes':
            if visit == '01' or visit == '12' or visit == '15' or visit == '26':
                plt.figure(figsize=(8,6))
                plt.pcolormesh(lons1,lats1[0:160],image_extract[0:160,:],cmap='cubehelix',
                               vmin=0.,vmax=1000.)
                plt.xlabel('SIII longitude')
                plt.ylabel('co-latitude')
                
                coord1 = dusk_active_region[0]
                coord2 = dusk_active_region[1]
                coord3 = dusk_active_region[2]
                coord4 = dusk_active_region[3]
                
                # Overplot the region of interest, e.g. a lat-lon box here:
                plt.plot([coord1[1],coord2[1],coord3[1],coord4[1],coord1[1]],  # corner A repeated to
                         [coord1[0],coord2[0],coord3[0],coord4[0],coord1[0]],  # close the box.
                         color='red',linewidth=3.)
                plt.show()
            else:
                plt.figure(figsize=(8,6))
                plt.pcolormesh(lons1,lats1,image_extract,cmap='cubehelix',
                               vmin=0.,vmax=30.)
                
                coord1 = dusk_active_region[0]
                coord2 = dusk_active_region[1]
                coord3 = dusk_active_region[2]
                coord4 = dusk_active_region[3]
                #coord5 = dusk_active_region[4]
                
                # Overplot the region of interest, e.g. a lat-lon box here:
                plt.plot([coord1[1],coord2[1],coord3[1],coord4[1],coord1[1]],  # corner A repeated to
                         [coord1[0],coord2[0],coord3[0],coord4[0],coord1[0]],  # close the box.
                         color='red',linewidth=3.)
                plt.show()
            
        if plotting == 'yes':
            plt.figure()
            plt.imshow(mask_dusk, origin='lower')  # zoom to see tiny squares! #roi mask
            plt.title('Mask in image space')
            plt.show()
            
        image_extract[mask_dusk==False] = np.nan
        
        dusk_px = dusk_active.data
        
        tot_size = len(dusk_px)
        num_px = ignore_nan_counter(dusk_px)
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=1000.)
            plt.title('Image masked off by the ROI')
            plt.show()
        
        im_full          = np.zeros((720,1440))        # new array full of nans
        im_full[0:160,:] = image_extract[0:160,:] #roi_im_full[0:160,:]
        
        if plotting == 'yes':
            plt.figure()
            plt.imshow(im_full, origin='lower')
            plt.title('ROI intensities in full image space')
            plt.show()
            
        # dar_boundary = path.Path([(20,192.25), (30,200), (20,220), (15,230), (15,220), (20,192.25)]) 
        # testlons = np.arange(0,360,0.25)
        # testcolats = np.arange(0,40,0.25)
        # llons, llats = np.meshgrid(testlons, testcolats) # checked correct
        # coords = np.vstack([llats.flatten(), llons.flatten()]).T
        # dar_mask = dar_boundary.contains_points(coords) # ([[llats], [llons]])
        # dar_mask_2D = dar_mask.reshape(160,1440)

        # roi_mask_full          = np.zeros((720,1440))
        # roi_mask_full[0:160,:] = dar_mask_2D
        

    elif region == 'swirl' or region == 'Swirl':
        print('Masking Swirl Region')
        mask_swirl = swirl.insert_points(llats, llons, image)
        #mask_swirl160 = mask_swirl[0:160,:]
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons1,lats1,image_extract,cmap='cubehelix',
                           vmin=0.,vmax=30.)
            
            coord1 = swirl_region[0]
            coord2 = swirl_region[1]
            coord3 = swirl_region[2]
            coord4 = swirl_region[3]
            
            # Overplot the region of interest, e.g. a lat-lon box here:
            plt.plot([coord1[1],coord2[1],coord3[1],coord4[1],coord1[1]],  # corner A repeated to
                     [coord1[0],coord2[0],coord3[0],coord4[0],coord1[0]],  # close the box.
                     color='red',linewidth=3.)
            plt.show()
            
        if plotting == 'yes':
            plt.figure()
            plt.imshow(mask_swirl, origin='lower')  # zoom to see tiny squares! #roi mask
            plt.title('Mask in image space')
            plt.show()
            
        image_extract[mask_swirl==False] = np.nan
        
        swirl_px = swirl.data
        
        tot_size = len(swirl_px)
        num_px = ignore_nan_counter(swirl_px)
        
        #breakpoint()
        
        if plotting == 'yes':
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,image_extract[0:160,:],cmap='cubehelix',
                           vmin=0.,vmax=1000.)
            plt.title('Image masked off by the ROI')
            plt.show()
        
        im_full          = np.zeros((720,1440))        # new array full of nans
        im_full[0:160,:] = image_extract[0:160,:] #roi_im_full[0:160,:]
        
        if plotting == 'yes':
            plt.figure()
            plt.imshow(im_full, origin='lower')
            plt.title('ROI intensities in full image space')
            plt.show()
        
    else:
        print('Not a Valid Region')
    
    #breakpoint()
        
    return im_full, tot_size, num_px #roi_mask_full
       
# actual function to call to cut out regions
def image_processing(file,n,visit,visit_name,plotting):
    if __name__ == "__main__":

        # open intensities file
        infolist = fits.open(file)
        header = infolist[1].header
        image_data = infolist[1].data
        # hdu_list = fits.open(file)  # opens the FITS files, accessing data plus header info
        # hdu_list.info()                      # print file information

        # accessing specific header info entries:
        exp_time  = header['EXPT']
        aperture  = header['APERTURE']   # filter type - important as it determines the Gustin conversion factors for intensity/counts/powers etc.
        cml       = header['CML']
        dece      = header['DECE']
        hem       = header['HEMISPH']
        # obt     = hdu_list[0].header['OBLT']       # can't see oblateness in the fits header?
        dist_org  = header['DIST_ORG']   # Earth-planet distance in AU before reduction
        pcx       = header['PCX']        # Planet centre pixel
        pcy       = header['PCY']        # Planet centre pixel
        nppa_org  = header['NPPA_ORG']   # North pole position angle before reduction
        nppa      = header['NPPA']       # North pole position angle
        pixsize   = header['PXSEC']    # Pixel size in arc seconds
        pxsec     = pixsize                          # just matching a variable name here that's used in the broject function
        dist      = header['DIST']       # Standard (scaled) Earth-planet distance in AU
        dmeq_orig = header['DMEQ_ORG']   # Original diameter of planet equator in arcsecond
        # ------------------------------------------------------------------------------
        cts2kr    = header['CTS2KR']     # reciprocal of conversion factor in counts/sec/kR
        
        #breakpoint()
        
        # If 1 / conversion factor is ~3994, this implies a colour ratio of 1.10
        # for Saturn with a STIS SrF2 image (see Gustin+ 2012 Table 1):
        colour_ratio = 2.5 #1.10 for Saturn

        # And this in turn means that the counts-per-second to total emitted power (Watts)
        # conversion factor is 9.04e-10 (Gustin+2012 Table 2), for STIS SrF2:
        gustin_conv_factor = 1.02e-9 #9.04e-10 for Saturn
        gustin_conv_factor_swirl = 1.16e-09#(1.16 * 10**-9) 1.16 × 10−9
        # "Conversion factor to be multiplied by the squared HST-planet distance (km)
        # to determine the total emitted power (Watts) from observed counts per second."

        # ------------------------------------------------------------------------------

        # In some fits files (Jupiter), these 'delta' values are listed as DELRPKM in the header.
        # If not (Saturn), it's hard-wired in here depending on the target planet (probably Saturn!).

        deltas = {'Mars': 0., 'Jupiter': 240., 'Saturn': 1100., 'Uranus': 0.}   # auroral emission altitudes at homopause in km
        delrpkm = deltas['Jupiter']
        rpkm = rpeqkm                                # just matching a variable name here that's used in the broject function

        # if find delrpkm fine, otherwise set some default 'deltas' as above. **********
        # ------------------------------------------------------------------------------

        # convert HST timestamp to time at Jupiter using light travel time:
        start_time = parse(header['UDATE'])     # create datetime object
        try:
            dist_org = header['DIST_ORG']
            ltime = dist_org*c.au/c.c
            lighttime = dt.timedelta(seconds=ltime)
        except KeyError:
            ltime = dt.timedelta(seconds=2524.42) 
            
        exposure = dt.timedelta(seconds=exp_time)
        start_time_jup = start_time - lighttime          # correct for light travel time
        end_time_jup = start_time_jup + exposure      # end of exposure time
        mid_ex_jup = start_time_jup + (exposure/2.)   # mid-exposure time at Jupiter
        mid_ex = start_time + (exposure/2.)                 # mid-exposure time at HST
        
        # --------------- make polar projection plot -----------------
        
        if plotting == 'yes':
            image_data = infolist[1].data
            # make a quick-look plot to check image array content:
            plt.figure()
            plt.title('Image array in fits file.')
            plt.imshow(image_data, cmap='cubehelix',origin='lower',vmin=1, vmax=1000)
            plt.xlabel('longitude pixels')
            plt.ylabel('co-latitude pixels')
            # plt.colorbar()
            cbar = plt.colorbar(pad=0.05)
            cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
    
            plt.show()
            
        # ------------------------------------------------------------------------------
        
        # perform limb trimming based on angle of surface vector normal to the sun
        lonm = np.radians(np.linspace(0.,360.,num=1440))
        latm = np.radians(np.linspace(-90.,90.,num=720))

        limb_mask = np.zeros((720,1440))   # rows by columns
        cmlr = np.radians(cml)             # convert CML to radians
        dec  = np.radians(dece)            # convert declination angle to radians
        for i in range(0,720):
            limb_mask[i,:] = np.sin(latm[i])*np.sin(dec) + np.cos(latm[i])*np.cos(dec)*np.cos(lonm-cmlr)

        limb_mask    = np.flip(limb_mask,axis=1)     # flip the mask horizontally, not sure why this is needed
        cliplim = np.cos(np.radians(88.))            # set a minimum required vector normal surface-sun angle
        clipind = np.squeeze([limb_mask >= cliplim]) # False = out of bounds (over the limb)
        image_data[clipind  == False] = np.nan  # set image array values outside clip mask to nans
        
            #infolist.close()                # close the file once you're done with it

        im_clean = np.flip(image_data,0)
        lons = np.arange(0,1440,1)
        lats = np.arange(0,720,1)   # 0-40 degrees colat in image pixel res.
        
        im_4broject = im_clean.copy()
        
        if plotting == 'yes':
            # Quick plot check of the centred, limb-trimmed image:
            plt.figure(figsize=(8,6))
            plt.pcolormesh(lons,lats,np.flip(np.flip(im_clean, axis=1)),cmap='cubehelix',
                            vmin=0.,vmax=1000.)
            plt.title('Image centred, limb-trimmed.')
            plt.xlabel('longitude pixels')
            plt.ylabel('co-latitude pixels')
    
            plt.show()

        # flip image vertically if required (ease of indexing) and extract auroral region:
        if  hem == 'north':
            im_flip = np.flip(image_data,0)
            image_extract = im_flip[0:160,:] # extract image in colat range 0-40 deg (4*40 = 160 pixels in image lat space):
        elif hem == 'south':
            image_extract = image_data[0:160,:]
            
        
        im_flip[im_flip < -100] = np.nan
        

        if plotting == 'yes':
            plt.figure()
            plt.imshow(image_extract,cmap='cubehelix',vmin=0.,vmax=1000.)
            plt.title('Auroral region extracted')
            plt.xlabel('longitude pixels')
            plt.ylabel('co-latitude pixels')
            cbar = plt.colorbar(pad=0.05)
            cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
            plt.show()
        
        # ----------------
        
            # set up polar coords
            rho   = np.linspace(0,40,     num=160 ) # colat vector with image pixel resolution steps
            theta = np.linspace(0,2*np.pi,num=1440) # longitude vector in radian space and image pixel resolution steps

            plt.figure(figsize=(8,6))
            fs = 12
            ax = plt.subplot(projection='polar')           # initialize polar projection
            ax.set_title('Polar projection.')
            plt.fill_between(theta, 0, 40, alpha=0.2,hatch="/",color='gray')
            ax.set_theta_zero_location("N")                # set angle 0.0 to top of plot
            ax.set_xticklabels(['0','45','90','135','180','225','270','315'],
                                fontweight='bold',fontsize=fs)
            ax.tick_params(axis='x',pad=-1.)               # shift position of LT labels
            ax.set_yticklabels(['','','','',''])           # turn off auto lat labels
            ax.set_yticks([0,10,20,30,40])                    # but set grid spacing
            ax.set_ylim([0,40])                            # max colat range
    
            # plot image data in log-colour scale:
            # plt.pcolormesh(theta,rho,image_extract,cmap='cubehelix',
            #                norm=LogNorm(vmin=.1,vmax=100.))
    
            # plot image data in linear colour-scale:
            plt.pcolormesh(theta,rho,image_extract,cmap='cubehelix',
                            vmin=0.,vmax=1000.)
    
            # Add colourbar: 
            cbar = plt.colorbar(ticks=[0.,100.,500.,900.],pad=0.05)
            cbar.ax.set_yticklabels(['0','100','500','900'])
            cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
    
            plt.show()

        return im_flip, cml, dece, dist, pcx, pcy, pxsec, nppa, rpkm, delrpkm, oblt, cts2kr, gustin_conv_factor, gustin_conv_factor_swirl, im_4broject, mid_ex_jup, dmeq_orig # may need other bits

# Nicked from Jonny's pypline.py file:
def _cylbroject(pimage, cml, dece, dmeq, xcen, ycen, psize, nppa, req, obt, ndiv=2, correct=True):
    
    # print(cml)
    # print(dece)
    # print(dmeq)
    # print(xcen)
    # print(ycen)
    # print(psize)
    # print(nppa)
    # print(req)
    # print(obt)
    # print(ndiv)

    if nppa == 999:
        nppa = 0

    ny, nx = pimage.shape
    xsize, ysize = 1400, 1400
    rimage = np.zeros((ysize, xsize))
    cimage = np.zeros((ysize, xsize))

    rad2deg = 180. / np.pi
    deg2rad = np.pi / 180.

   # # # third, calculate the pixel size. 1 au = 1.49598e8 km */
    plen = (psize / 3600.) * deg2rad * dmeq * 1.49598e8
    #print(plen) # in km

    # # # first initialize global variables */
    sn = np.sin(nppa * deg2rad)
    cn = np.cos(nppa * deg2rad)
    sd = np.sin(dece * deg2rad)
    cd = np.cos(dece * deg2rad)
    td = np.tan((dece + 90.0) * deg2rad)
    a = req / plen  # /* rquatorial planet radius in  pixel scale */
    o2 = 2.0 * obt
    oo = obt * obt
    a1o = a * (1.0 - obt)
    px = 0
    py = 0
    ii = 0
    jj = 0

    lo = np.zeros(nx * ndiv)
    sb = np.zeros(nx * ndiv)
    cb = np.zeros(nx * ndiv)
    ll = np.zeros(nx * ndiv)
    la = np.zeros(ny * ndiv)
    sa = np.zeros(ny * ndiv)
    ca = np.zeros(ny * ndiv)
    caca = np.zeros(ny * ndiv)
    r = np.zeros(ny * ndiv)

    d = 360.0 / nx / ndiv
    d2 = 360.0 / nx / ndiv / 2.0
    temp = (1.0 - obt) * (1.0 - obt) * td
    for i in range(nx * ndiv):  # / * longitude * /
        lo[i] = i * d + d2 + cml
        sb[i] = np.sin(lo[i] * deg2rad)
        cb[i] = np.cos(lo[i] * deg2rad)
        ll[i] = np.arctan(temp * cb[i]) * rad2deg

    d = 180.0 / ny / ndiv
    d2 = 180.0 / ny / ndiv / 2.0
    for i in range(ny * ndiv):  # / * latitude * /
        la[i] = i * d + d2 - 90.0
        sa[i] = np.sin(la[i] * deg2rad)
        ca[i] = np.cos(la[i] * deg2rad)
        caca[i] = ca[i] * ca[i]
        r[i] = a1o / (1.0 - (o2 - oo) * caca[i])**0.5

#  / * here's the big loop. start from the longitude corresponding to the left edge of the brojected image.
#    that is, start from 360-cml+90 and follow down (leftward) on the pimage while the apparent longitude on the
#    brojected image increases to the right. * /
    for i in range(nx * ndiv):
        if dece < 0.0:
            start = 0
            end = int(((ll[i] + 90.0) / 180.0 * ndiv * ny) + 1)
        else:
            start = int(((ll[i] + 90.0) / 180.0 * ndiv * ny))
            end = ndiv * ny
        for j in range(start, end):
            x = r[j] * ca[j] * sb[i]
            y = r[j] * sa[j]
            z = r[j] * ca[j] * cb[i]
            px = x
            py = y * cd - z * sd
            temp = px
            px = int(px * cn - py * sn + xcen)
            py = int(temp * sn + py * cn + ycen)
            if (px >= 0) & (px < xsize) & (py >= 0) & (py < ysize):
                ii = int(i / ndiv)
                jj = int(j / ndiv)
                value = pimage[jj, ii] / ndiv / ndiv
                # print(i, j, ii, jj, py, px, value)
                # return
                rimage[py, px] += value


#  /  *  Here, correct for area effect. Added by Juwhan Kim, 03 / 01 / 2005.
    if correct is True:
        value = 1.0 / ndiv / ndiv
        for i in range(nx * ndiv):
            if dece < 0.0:
                start = 0
                end = int(((ll[i] + 90.0) / 180.0 * ndiv * ny) + 1)
            else:
                start = int(((ll[i] + 90.0) / 180.0 * ndiv * ny))
                end = ndiv * ny
            for j in range(start, end):
                x = r[j] * ca[j] * sb[i]
                y = r[j] * sa[j]
                z = r[j] * ca[j] * cb[i]
                px = x
                py = y * cd - z * sd
                temp = px
                px = int(px * cn - py * sn + xcen)
                py = int(temp * sn + py * cn + ycen)
                if (px >= 0) & (px < xsize) & (py >= 0) & (py < ysize):
                    cimage[py, px] += value

        for i in range(xsize):
            for j in range(ysize):
                cimval = cimage[j, i]
                if cimval != 0:
                    rimage[j, i] = rimage[j, i] / cimval

    return rimage


def power_calculator(visit_list, year, prefix, extra, time, region, plotting):
    # Initialize lists at the START of the function
    cmls, deces, all_times = [], [], []
    all_keo_strips, all_lat_strips = [], []
    au_to_km = 1.496e8
    
    # We will return these at the very end
    final_powers, final_total_region, final_region_pixels = [], [], []

    for i in visit_list:
        # Reset per-visit powers so they don't accumulate across different visits
        powers, total_region, region_pixels = [], [], []
        
        arch = '*_v'+ i
        ti = str('/*0'+time+'*')
        visit_name = prefix+arch[-2:]
        ab = glob.glob(f'{root_folder}data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti) 
        ab.sort() 
        
        #ab = ab[3:]
        
        # Initialize a reference CML for this specific visit
        visit_start_cml = None #temp

        for n, f_path in tqdm(enumerate(ab), desc=f"Visit {i}"):
            res = image_processing(f_path, n, i, visit_name, plotting)
            image, cml, dece, dist, pcx, pcy, pxsec, nppa, rpkm, delrpkm, oblt, cts2kr, \
            gustin_conv_factor, gustin_conv_factor_swirl, im_4broject, mid_ex_jup, dmeq_orig = res
            
            # masks the RAW image before projection
            im_full_masked, tot_size, num_px = apply_mask(region, image, i, plotting)
            
            # DEFINE CYLBROJECT
            def cylbroject(img, ndiv=2):
                return _cylbroject(img, cml, dece, dist, pcx, pcy, pxsec, nppa, rpkm + delrpkm, oblt, ndiv, True)
            
            # PROJECT THE MASKED IMAGE
            bimage = cylbroject(np.flip(np.flip(im_full_masked, axis=1)), ndiv=2)
            
            # --- grid ---
            ny, nx = bimage.shape
            rows, cols = np.where(~np.isnan(bimage) & (bimage > 0))
            active_pixels = bimage[rows, cols]
            
            
# =========================this currently does NOT work =======================
#           idea is to try and regrid the lons after the power conversion as 
#           cylbroject remaps to a non 2D linear grid? but this doesn't work 
#           and I can't figure out why
#
#             # Capture the very first CML to use as a "grounding" point
#             if n == 0:
#                 visit_start_cml = cml
#             
#             lons_relative = (cols / nx) * 360
#             # abs_lons = (lons_relative)
#             rotation_correction = (cml - visit_start_cml)
# 
#             abs_lons = (lons_relative) % 360
#             #abs_lons = (lons_relative - (cml - visit_start_cml)) % 360
# =============================================================================

            # ny = 1400 (The full height of the cylindrical projection)
            # This scale maps the 1400 pixels to 180 degrees
            deg_per_pixel = 180.0 / 720
            print(f"DEBUG: pcy={pcy}, ny={ny}, min_row={np.min(rows)}")
            
            # Print this for each region (Swirl, Noon, Dusk)
            print(f"--- REGION DEBUG ---")
            print(f"Region: {region}") # (e.g., 'swirl')
            print(f"Row Range: {np.min(rows)} to {np.max(rows)}")
            print(f"Pixel Rows: {np.min(rows)} to {np.max(rows)}")
            print(f"Width in Pixels: {np.max(rows) - np.min(rows)}") 

            reso = 0.458
            pole_offset = 210.7

            #abs_lats = (rows-pole_offset) * reso
            # The Rigorous Linear Fit based on your Region Debug data
           # Select the scale based on the region to ensure perfect mapping
            if  visit == '01' or visit == '12' or visit == '15' or visit == '26':# or visit == '0111':
                if region == 'noon':
                    slope = -0.278
                    intercept = 74.30
                elif region == 'dusk':
                    slope = -0.306
                    intercept = 93.65
                elif region == 'swirl':
                    slope = -0.458
                    intercept = 121.5
            else:
                if region == 'noon':
                    # Maps Row 203 to 23 and Row 167 to 29
                    slope = -0.1666
                    intercept = 56.8
                elif region == 'dusk':
                    # Maps Row 257 to 2 and Row 208 to 22
                    slope = -0.408
                    intercept = 106.9
                elif region == 'swirl':
                    # Maps Row 250 to 7 and Row 226 to 18
                    slope = -0.458
                    intercept = 121.5
                    
            # Apply directly to pixels
            abs_lats = intercept + (slope * rows)
            
        
            # power calculation
            distance_squared = (dist * au_to_km)**2
            # np.nansum(bimage) now only sees what apply_mask allowed through
            current_power_gw = np.nansum(bimage) * (1/5120 if region == 'swirl' else cts2kr) * \
                               (gustin_conv_factor_swirl if region == 'swirl' else gustin_conv_factor) * \
                               distance_squared / 1e9

            print(f"{mid_ex_jup.strftime('%H:%M:%S')} | Power: {current_power_gw:.2f} GW")

            # store stats
            powers.append(current_power_gw)
            total_region.append(tot_size)
            region_pixels.append(num_px)
            cmls.append(cml); deces.append(dece); all_times.append(mid_ex_jup)
            
            
            # --- BINNING ---
            total_raw_intensity = np.nansum(active_pixels)
            strip_scale = current_power_gw / total_raw_intensity if total_raw_intensity > 0 else 0
            
            print(f"DEBUG: Lat Range = {np.nanmin(abs_lats):.1f} to {np.nanmax(abs_lats):.1f}")
            # Latitude Power Sum
            lat_p, _ = np.histogram(abs_lats, bins=np.arange(0, 40.5, 0.5), 
                         weights=active_pixels * strip_scale)
            
            # Apply the exact same shift and scale to the mask
            m_rows, m_cols = np.where(~np.isnan(bimage))
            #m_lats = (m_rows-pole_offset) * reso
            #m_lats = intercept + (slope * m_rows)
            # m_lats = lat_grid[m_rows, m_cols]
            # m_intercept, m_slope, m_rowss = get_calibrated_lats(region, m_rows, np.min(rows), np.max(rows))
            m_lats = m_lats = intercept + (slope * m_rows)
            
            lat_vis_counts, _ = np.histogram(m_lats, bins=np.arange(0, 40.5, 0.5))

            # 1. Dynamic Expected Pixels: Use the actual width of the mask for this region
            # This ensures the ratio is 1.0 when the region is fully in view
            expected_lat_pixels = np.max(lat_vis_counts)
            
            # 2. Correction with safety floor
            if expected_lat_pixels > 0:
                lat_visibility_ratio = lat_vis_counts / expected_lat_pixels
                
                # Only correct bins that are at least 5% visible to avoid noise spikes
                visible_mask = lat_visibility_ratio > 0.05
                
                lat_prof_corrected = np.zeros_like(lat_p)
                # Power / Ratio: If 100% visible, power remains unchanged.
                lat_prof_corrected[visible_mask] = lat_p[visible_mask] / lat_visibility_ratio[visible_mask]
            else:
                lat_prof_corrected = lat_p

            all_lat_strips.append(lat_prof_corrected)
            
            '''
            all_keo_strips.append(lon_prof)
            all_lat_strips.append(lat_prof)
            '''

        # using [-len(ab):] ensures that if you run multiple visits, the keogram only shows the current one
        # plot_stacked_keograms(all_times[-len(ab):], all_lat_strips[-len(ab):], region, i, np.mean(deces[-len(ab):]))
        
        # Accumulate the visit-specific stats into the master return lists
        final_powers.extend(powers)
        final_total_region.extend(total_region)
        final_region_pixels.extend(region_pixels)

    return final_powers, cmls, deces, final_total_region, final_region_pixels, all_keo_strips, all_times, all_lat_strips


def plot_stacked_keograms(times, lon_strips, lat_strips, region_name, visit_id, dece_val):
    # lon_matrix = np.array(lon_strips, dtype=float).T
    lat_matrix = np.array(lat_strips, dtype=float).T
    
    # 1. Bleach noise to white
    # lon_matrix[lon_matrix <= 0.05] = np.nan
    lat_matrix[lat_matrix <= 0.05] = np.nan

    
    # automatically find the data bounds
    has_lat = ~np.all(np.isnan(lat_matrix), axis=1)
    idx_lat = np.where(has_lat)[0]
    if len(idx_lat) == 0: return
    
    '''
    #fiinding longitude bounds
    has_lon = ~np.all(np.isnan(lon_matrix), axis=1)
    idx_lon = np.where(has_lon)[0]
    '''
    # buffer (2 degrees for lat, 5 for lon)
    lat_min, lat_max = max(0, idx_lat[0] - 2), min(180, idx_lat[-1] + 2)
    # lon_min, lon_max = max(0, idx_lon[0] - 5), min(360, idx_lon[-1] + 5)
    
    lat_edges = np.arange(0, 181, 1)
    # lon_edges = np.arange(0, 361, 1)
    
    # convert to numpy array
    times_arr = np.array(times)
    
    # time edges to match flat shading
    if len(times_arr) > 1:
        diff = times_arr[1] - times_arr[0]
        t_edges = np.append(times_arr - diff/2, times_arr[-1] + diff/2)
    else:
        t_edges = np.array([times_arr[0] - datetime.timedelta(minutes=1), 
                            times_arr[0] + datetime.timedelta(minutes=1)])
        
    print(f"Latitude Matrix Shape: {lat_matrix.shape}")
    print(f"Latitude Edges Range: {lat_edges[lat_min]} to {lat_edges[lat_max]}")

    fig, (ax_lon, ax_lat) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    current_cmap = plt.colormaps['magma'].copy()
    current_cmap.set_bad(color='white') 

    # --- LATITUDE ---
    im2 = ax_lat.pcolormesh(t_edges, lat_edges, lat_matrix, 
                            cmap=current_cmap, shading='flat', vmin=0.1)
    ax_lat.set_ylim(lat_max, lat_min) 
    ax_lat.set_ylabel('Colatitude ($^\circ$)')
    plt.colorbar(im2, ax=ax_lat, label='GW / deg')
    '''
    # --- LONGITUDE ---
    im1 = ax_lon.pcolormesh(t_edges, lon_edges, lon_matrix, 
                            cmap=current_cmap, shading='flat', vmin=0.1)
    ax_lon.set_ylim(lon_min, lon_max)
    ax_lon.set_ylabel('SIII Longitude ($^\circ$)')
    plt.colorbar(im1, ax=ax_lon, label='GW / deg')

    ax_lat.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.xticks(rotation=45)
    '''
    # FIND THE MAX INDEX IN  DATA
    max_row, max_col = np.unravel_index(np.nanargmax(lat_matrix), lat_matrix.shape)
    peak_time = times[max_col]
    print("--- TIMING VERIFICATION ---")
    print(f"Brightest pixel found at index: {max_col}")
    print(f"Timestamp for that index: {peak_time}")
    print("---------------------------")

    plt.savefig(f'{root_folder}keograms/keogram_{region_name}_{visit_id}.jpg', bbox_inches='tight', dpi=400)
    plt.close(fig)    
    

def plot_full_auroral_stack(keo_data, keo_times, phys_times, clock, press, bperp, btot, visit_id):
    # 6 Axes: 3 Keograms, 1 Clock, 1 Pressure, 1 Combined Mag
    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(12, 18), sharex=True)
    plt.subplots_adjust(hspace=0.1)

    cmap = plt.colormaps['magma'].copy()
    cmap.set_bad(color='white')
    # lat_edges = np.arange(0, 181, 1)
    lat_edges = np.arange(0, 40.5, 0.5)
    
    # --- HELPER: Process Matrix & Find Bounds for a Single Region ---
    def get_region_plot_params(strips):
        matrix = np.array(strips, dtype=float).T
        matrix[matrix <= 0.05] = np.nan # Bleach noise
        
        # Find bounds for this specific matrix
        has_lat = ~np.all(np.isnan(matrix), axis=1)
        idx_lat = np.where(has_lat)[0]
        
        if len(idx_lat) == 0:
            return matrix, 0, 40 # Default if empty
            
        # l_min, l_max = max(0, idx_lat[0] - 2), min(180, idx_lat[-1] + 2)
        # return matrix, l_min, l_max
        # Convert index back to degrees (0.25 deg per index)
        l_min = idx_lat[0] * 0.5
        l_max = (idx_lat[-1] + 1) * 0.5 # +1 to get the far edge of the last pixel
        # change 0.5 if lat edges changes to 0.5
        return matrix, l_min, l_max
    
    def process_keo(strips):
        m = np.array(strips, dtype=float).T
        m[m <= 0.05] = np.nan
        return m

    def get_edges(t):
        t_arr = np.array(t)
        d = t_arr[1] - t_arr[0] if len(t_arr) > 1 else datetime.timedelta(minutes=1)
        return np.append(t_arr - d/2, t_arr[-1] + d/2)

    k_edges = get_edges(keo_times)

    
    # --- AX1: SWIRL ---
    s_mat, s_min, s_max = get_region_plot_params(keo_data['swirl'][1])
    k1=ax1.pcolormesh(k_edges, lat_edges, s_mat, cmap=cmap, vmin=0.1)
    ax1.set_ylim(s_max+0.5, s_min-0.5)
    #ax1.set_ylim(6,18) # Zoom in on the Swirl region
    ax1.set_ylabel('Swirl\nColat ($^\circ$)',fontsize=14)
    ax1.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax1.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)
    # ax1.invert_yaxis()

    # --- AX2: NOON ---
    n_mat, n_min, n_max = get_region_plot_params(keo_data['noon'][1])
    k2=ax2.pcolormesh(k_edges, lat_edges, n_mat, cmap=cmap, vmin=0.1)
    ax2.set_ylim(n_max+0.5, n_min-0.5)
    #ax2.set_ylim(23, 29) # Zoom in on the Noon region
    ax2.set_ylabel('Noon\nColat ($^\circ$)',fontsize=14)
    ax2.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax2.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)
    # ax2.invert_yaxis()

    # --- AX3: DUSK ---
    d_mat, d_min, d_max = get_region_plot_params(keo_data['dusk'][1])
    k3=ax3.pcolormesh(k_edges, lat_edges, d_mat, cmap=cmap, vmin=0.1)
    ax3.set_ylim(d_max+0.5, d_min-0.5)
    #ax3.set_ylim(2, 22) # Zoom in on the Dusk region
    ax3.set_ylabel('Dusk\nColat ($^\circ$)',fontsize=14)
    ax3.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax3.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)
    # ax3.invert_yaxis()

    # --- CLOCK ANGLE (AX4) ---
    ax4.scatter(phys_times, clock, color='tab:blue', label='Clock Angle')
    ax4.plot(phys_times, clock,color='tab:blue')
    ax4.set_ylabel('Clock ($^\circ$)',fontsize=14)
    # ax4.axhline(0, color='black', lw=0.5, ls='--')
    ax4.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax4.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)

    # --- PRESSURE (AX5) ---
    ax5.scatter(phys_times, press, color='tab:red')
    ax5.plot(phys_times,press,color='tab:red')
    ax5.set_ylabel('Pressure\n(nPa)',fontsize=14)
    ax5.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax5.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)

    # --- COMBINED MAGNETIC FIELD (AX6) ---
    ax6.scatter(phys_times, bperp, color='blue', label='$B_{\perp}$')
    ax6.scatter(phys_times, btot, color='orange', label='$B_{tot}$', linestyle='--')
    ax6.plot(phys_times, btot,color='orange')
    ax6.plot(phys_times,bperp,color='blue')
    ax6.set_ylabel('Magnetic Field \nStrength (nT)',fontsize=14)
    ax6.legend(loc='upper right', fontsize='small', frameon=False)
    ax6.tick_params(which='major',direction='in',bottom=True, top=True, left=True, right=True, labelsize=14, length=5)
    ax6.tick_params(which='minor', direction='in', bottom=True, top=True, left=True, right=True, length=2)
    ax6.legend(framealpha=0.5,fontsize=14, loc='best')#,loc ="upper left")
    
    # Final Formatting
    ax6.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    ax6.set_xlabel('Ionospheric Time', fontsize=14)
    
    # Replace your colorbar calls with this pattern for ax1, ax2, and ax3:
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="3%", pad=0.1) # 3% width of the main plot
    cb1 = fig.colorbar(k1, cax=cax1)
    cb1.set_label('GW / deg', fontsize=14)
    cb1.ax.tick_params(labelsize=14)

    
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="3%", pad=0.1)
    cb2 = fig.colorbar(k2, cax=cax2, label='GW / deg')
    cb2.set_label('GW / deg', fontsize=14)
    cb2.ax.tick_params(labelsize=14)
    
    divider3 = make_axes_locatable(ax3)
    cax3 = divider3.append_axes("right", size="3%", pad=0.1)
    cb3 = fig.colorbar(k3, cax=cax3, label='GW / deg')
    cb3.set_label('GW / deg', fontsize=14)
    cb3.ax.tick_params(labelsize=14)
    
    for ax in [ax4, ax5, ax6]:
        divider = make_axes_locatable(ax)
        _ = divider.append_axes("right", size="3%", pad=0.1).set_visible(False)
    '''
    plt.colorbar(k1, ax=ax1, label='GW / deg')
    plt.colorbar(k2, ax=ax2, label='GW / deg')
    plt.colorbar(k3, ax=ax3, label='GW / deg')
    '''
    return fig

    
# user input to access right file
year = input("Year of the visit:  \n")
if year == '2016':
    pre = input('Campaign from Jonny or Denis? (1/2)  ')
    if pre == '1' or pre == 'jonny' or pre == 'j' or pre == 'J' or pre == 'Jonny':
        prefix = 'ocx8'
        extra = 'nichols/' # this may be have to be removed if your directory system did not differentiate between Jonny's and Dennis's campaigns
    else:
        prefix = 'od8k'
        extra = 'grodent/' # this may be have to be removed if your directory system did not differentiate between Jonny's and Dennis's campaigns
elif year == '2019':
    prefix = 'odxc'
    extra = ''
elif year == '2021':
    prefix = 'oef4'
    extra = ''
elif year == '2017' or  year == '2018':
    prefix = 'od8k'
    extra = ''
elif year == '2022':
    prefix = 'oeow'
    extra = ''
time = str(input('Exposure time (in seconds: 10, 30, 100...): \n'))

'''
# --------- timestamps ------------
times_02 = visit_02_data['Ionospheric_Times'].
'''

swirl_p, cmls, deces, s_tot, s_px, swirl_strips, swirl_times, swirl_lat_strips = power_calculator(visit_list, year, prefix, extra, time, 'swirl','no')
noon_p, cmls, deces, n_tot, n_px, noon_strips, noon_times, noon_lat_strips = power_calculator(visit_list, year, prefix, extra, time, 'noon','no')
dusk_p, cmls, deces, d_tot, d_px, dusk_strips, dusk_times, dusk_lat_strips = power_calculator(visit_list, year, prefix, extra, time, 'dusk','no')

# Now you can actually plot it: (all_times[-len(ab):], all_keo_strips[-len(ab):], all_lat_strips[-len(ab):], region, i, np.mean(deces[-len(ab):]))
# plot_stacked_keograms(swirl_times, swirl_lat_strips, 'Swirl', visit, np.mean(deces))
# plot_stacked_keograms(noon_times, noon_lat_strips, 'Noon', visit, np.mean(deces))
# plot_stacked_keograms(dusk_times, dusk_lat_strips, 'Dusk', visit, np.mean(deces))

# --- EXECUTION EXAMPLE ---
k_dict = {
    'swirl': (swirl_times, swirl_lat_strips),
    'noon': (noon_times, noon_lat_strips),
    'dusk': (dusk_times, dusk_lat_strips)
}

# Fetch the column using globals
raw_times = globals()[f"visit_{visit}_data"]['Ionospheric_Time']

# Convert to actual datetime objects
phys_times = pd.to_datetime(raw_times)

# Fetch your data arrays dynamically
clock = globals()[f"clock_{visit}"]
pressure = globals()[f"pressure_{visit}"]
b_perp = globals()[f"b_perp_{visit}"]
btot = globals()[f"btot_{visit}"]

# call the plotter
# use 'noon_times' as the 'phys_times' argument as all times should be same
fig = plot_full_auroral_stack(k_dict,noon_times,phys_times,clock, pressure, b_perp, btot, visit)
saveloc = (f'{root_saves}keogram_stack_{visit}_plot.jpg') 
plt.savefig(saveloc,bbox_inches='tight',dpi=400)
