#=== originally: HST_emission_power.py ========================================
# Script to read in projected/unprojected HST STIS .fits images, define a
# spatial region, and calculate the UV emission power.
#
# Translation between unprojected/projected images is partly handled by IDL
# pipeline legacy code from Boston University, and partly from python code
# provided by Jonny Nichols at Leicester. (broject function)

# Emission power computed as per Gustin+ 2012.

# Should allow e.g., calculation of total auroral oval UV emission power using
# statistical auroral boundaries, planetary auroral comparisons, or
# application to Voronoi image segmentations, etc.

# === Base code: Joe Kinrade - 8 January 2025 =================================
# Modified by Hannah Joyce to load in and process on a per visit basis as well
# as to mask out specific auroral regions to extract power from

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
#spice.furnsh("/Users/joe/data/SPICE/cassini_generic.mk") # SPICE kernels
#import time_conversions as tc
import collections
from matplotlib import path
import scipy.constants as c
# import collections.abc as collections # If using Python version 3.10 and above
from scipy import signal
from tqdm import tqdm
import glob

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

#visit_list = ['02','03','04','05','08','09','10','11','16','17','18','19','20','21','24','25','27','28','34','35']#['04','05','08','09','10','11','12','15','16','17','18','19','20','21']#'04','05','08','09','10','11','12','15','16','17','18','19','20','21'] #'04','05','08','09','10','11','12','13','15','16','17','18','19','20','21'
#visit_list = ['01','12','15','26']
visit_list = ['01'] # start with 18 in morning
visit='01'

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
            if visit == '01' or visit == '12' or visit == '15' or 'visit' == '26':
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

        # convert HST timestamp to time at Saturn using light travel time:
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

        return im_flip, cml, dece, dist, pcx, pcy, pxsec, nppa, rpkm, delrpkm, oblt, cts2kr, gustin_conv_factor, gustin_conv_factor_swirl, im_4broject # may need other bits
    
# Load fits file and image/header ----------------------------------------------

#file = '/Users/hannah/OneDrive - Lancaster University/aurora/data/2016/extract/nichols/151_v19/nopolar100/jup_16-151-13-39-57_0100_v19_stis_f25srf2_proj.fits'

# ==============================================================================
# Now at the point we can try back-projecting this projected image mask
# (or masked-out image) using the back-project function:
# ==============================================================================

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

# this cylbroject function definition is feeding inputs into _cylbroject - ? JK
# def cylbroject(image, ndiv=6):
#     # self._check_image_loaded(proj=True)      # commented out JK
#     print('Brojecting with ndiv = ', ndiv)
#     bimage = _cylbroject(image,
#                                   cml, dece, dist,
#                                   pcx, pcy, pxsec,    # ***** need to define pxsec
#                                   nppa, rpkm + delrpkm,
#                                   oblt, ndiv, True)
#     return bimage
#
#

def power_calculator(visit_list,year,prefix,extra,time,region,plotting): # region
    cmls=[]
    deces=[]
    # loop for visits
    for i in visit_list:
        
        powers = []
        total_region = []
        region_pixels = []
        print(f'VISIT {i} \n \n')
        
        arch = '*_v'+ i
        ti = str('/*0'+time+'*')
        visit_name = prefix+arch[-2:]
        ab = glob.glob(f'{root_folder}data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti) 
        ab.sort() 
    
        # loop through each file in the visit
        for  n,i in tqdm(enumerate(ab)):
            
            # if n == 1:
            #     breakpoint()
    
            print(n)
            print(i)
            
            filename = str(i)[-51:-5]
            visit = filename[-20:-18]
            
            # use test region
            #i = '/Users/hannah/OneDrive - Lancaster University/aurora/data/2016/extract/nichols/137_v01/nopolar100/jup_16-137-23-43-30_0100_v01_stis_f25srf2_proj.fits'
            
            image, cml, dece, dist, pcx, pcy, pxsec, nppa, rpkm, delrpkm, oblt, cts2kr, gustin_conv_factor, gustin_conv_factor_swirl, im_4broject = image_processing(i, n, visit, visit_name, plotting)
            
            cmls.append(cml)
            deces.append(dece)

            distance_squared = (dist * au_to_km)**2          # AU in km

            im_full, tot_size, num_px = apply_mask(region, image, visit, plotting)
            
            total_region.append(tot_size)
            region_pixels.append(num_px)
            
            def cylbroject(image, ndiv=2):
                # self._check_image_loaded(proj=True)      # commented out JK
                print('Brojecting with ndiv = ', ndiv)
                bimage = _cylbroject(image,
                                              cml, dece, dist,
                                              pcx, pcy, pxsec,    # ***** need to define pxsec
                                              nppa, rpkm + delrpkm,
                                              oblt, ndiv, True)
                return bimage
            
            #breakpoint()
            bimage = cylbroject(np.flip(np.flip(im_full, axis=1)), ndiv=2)   # flips required to get bimage looking right?
            #full_image = cylbroject(np.flip(np.flip(im_4broject, axis=1)), ndiv=2)
            #bimage = cylbroject(np.flip(im_full),ndiv=2)
            #breakpoint()

            if plotting == 'yes':
            
                plt.figure()
                plt.title('Image array in fits file.')
                plt.imshow(bimage, cmap='cubehelix',origin='lower')
                plt.show()
                
                # plt.figure()
                # plt.title('Brojected full image.')
                # plt.imshow(full_image, cmap='cubehelix',origin='lower',vmin=0,vmax=1000)
                # plt.xlabel('pixels')
                # plt.ylabel('pixels')
                # cbar = plt.colorbar(pad=0.05)
                # cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)

                # # plt.colorbar()
                # # cbar = plt.colorbar(pad=0.05)
                # # cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
                # fignamei = 'broject_full.pdf'
                # plt.savefig(fignamei, dpi=350) #, bbox_inches='tight')  

                # plt.show()
                
            print(np.nansum(bimage))
            if region == 'swirl':
                print('Using CR 12')
                # colour ratio 12
                total_power_emitted_from_roi = np.nansum(bimage) * (1/5120) * distance_squared * gustin_conv_factor_swirl / 1e9
            else:
                print('Using CR 2.5')
                # colour ratio 2.5
                # calculate emitted power from ROI in GW (exposure time not required here as kR intensities are per second):
                total_power_emitted_from_roi = np.nansum(bimage) * cts2kr * distance_squared * gustin_conv_factor / 1e9

            print('Total power emitted from ROI in GW:')
            print(total_power_emitted_from_roi)
            
            powers.append(total_power_emitted_from_roi)
            
            # # swirl region
            # swirl_full = apply_mask('swirl', image, visit, plotting)

            # # Try and see  what happens! Not sure if image input needs to be full [1440,720] , centred projection?
            # simage = cylbroject(np.flip(swirl_full),ndiv=2)
            # #bimage = cylbroject(image_centred,ndiv=2)
            # print('Swirl Successfully Projected')

            # plt.figure()
            # plt.title('Image array in fits file.')
            # plt.imshow(simage, cmap='cubehelix',origin='lower')
            # # plt.colorbar()
            # # cbar = plt.colorbar(pad=0.05)
            # # cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
            # plt.show()

            # # calculate emitted power from ROI in GW (exposure time not required here as kR intensities are per second):
            # total_power_swirl = np.nansum(simage) * cts2kr * distance_squared * gustin_conv_factor / 1e9
            
            # swirl_powers.append(total_power_swirl)
            
            # # --------------------            

            # # dusk active region
            # dusk_full = apply_mask('dusk', image, visit, plotting)

            # # Try and see  what happens! Not sure if image input needs to be full [1440,720] , centred projection?
            # dimage = cylbroject(np.flip(dusk_full),ndiv=2)
            # #bimage = cylbroject(image_centred,ndiv=2)
            # print('Dusk Successfully Projected')

            # plt.figure()
            # plt.title('Image array in fits file.')
            # plt.imshow(dimage, cmap='cubehelix',origin='lower')
            # # plt.colorbar()
            # # cbar = plt.colorbar(pad=0.05)
            # # cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
            # plt.show()

            # # calculate emitted power from ROI in GW (exposure time not required here as kR intensities are per second):
            # total_power_dusk = np.nansum(dimage) * cts2kr * distance_squared * gustin_conv_factor / 1e9
            
            # dusk_powers.append(total_power_dusk)
            
            # # --------------------    
            
            # # noon active region
            # noon_full = apply_mask('noon', image, visit, plotting)

            # # Try and see  what happens! Not sure if image input needs to be full [1440,720] , centred projection?
            # nimage = cylbroject(np.flip(noon_full),ndiv=2)
            # #bimage = cylbroject(image_centred,ndiv=2)
            # print('Noon Successfully Projected')

            # plt.figure()
            # plt.title('Image array in fits file.')
            # plt.imshow(nimage, cmap='cubehelix',origin='lower')
            # # plt.colorbar()
            # # cbar = plt.colorbar(pad=0.05)
            # # cbar.ax.set_ylabel('Intensity [kR]',fontsize=12)
            # plt.show()
    
            # # calculate emitted power from ROI in GW (exposure time not required here as kR intensities are per second):
            # total_power_noon = np.nansum(nimage) * cts2kr * distance_squared * gustin_conv_factor / 1e9
            
            # noon_powers.append(total_power_noon)

    return powers, cmls, deces, total_region, region_pixels #swirl_powers, dusk_powers, noon_powers, 
    

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
#region = str(input("Region of Interest:  \n"))


# # ==============================================================================
# # Once the back-projected image looks OK, we can proceed with the emission power
# # calculation here.
# # ==============================================================================
#
# # ISOLATE THE ROI INTENSITIES IN A FULL 1440*720 PROJECTED IMAGE (all other pixels set to nans)

swirl_powers, cmls, deces, total_swirl, swirl_pixels = power_calculator(visit_list,year,prefix,extra,time,'swirl','no') # region
noon_powers, cmls, deces, total_noon, noon_pixels = power_calculator(visit_list,year,prefix,extra,time,'noon','no')
dusk_powers, cmls, deces, total_dusk, dusk_pixels = power_calculator(visit_list,year,prefix,extra,time,'dusk','no')
'''
#test_powers, cmls, deces = power_calculator(visit_list,year,prefix,extra,time,'test','yes')

#dataframe for total power per image
total_power_regions_df = pd.DataFrame()

total_power_regions_df = total_power_regions_df.assign(CML=cmls)
total_power_regions_df = total_power_regions_df.assign(DECE=deces)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Swirl=swirl_powers)
total_power_regions_df = total_power_regions_df.assign(Total_Region_Swirl=total_swirl)
total_power_regions_df = total_power_regions_df.assign(Total_Pixels_Swirl=swirl_pixels)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Noon_Active=noon_powers)
total_power_regions_df = total_power_regions_df.assign(Total_Region_Noon=total_noon)
total_power_regions_df = total_power_regions_df.assign(Total_Pixels_Noon=noon_pixels)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Dusk_Active=dusk_powers)
total_power_regions_df = total_power_regions_df.assign(Total_Region_Dusk=total_dusk)
total_power_regions_df = total_power_regions_df.assign(Total_Pixels_Dusk=dusk_pixels)

total_power_regions_df.to_csv(f'{root_folder}python_scripts/dataframes/new_total_powers_regions_{visit}_CR_adj_26update.csv',index=False)
'''
