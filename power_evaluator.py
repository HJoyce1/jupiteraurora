"""
Created on Wed Mar  6 13:39:13 2024

@author: hannah

script that converts intensities to power values and can write out to several dataframes

this needs re-looking at as currently requires uncommenting and commenting out sections to work
it is also unclear exactly what dataframes are being saved to what
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
import numpy as np
from pathlib import Path
from dateutil.parser import parse
import datetime as dt
import pandas as pd
from matplotlib.colors import LogNorm
from reading_mfp import moonfploc
from scipy import stats
import glob
import imageio
from astropy.io import fits
from tqdm import tqdm
from matplotlib import path
import scipy.constants as c
import rayleigh_to_power as rtp


visit_list = ['03']#['02','03','04','05','08','09','10','16','17','18','19','20','21','24','25','27','34','35']#['04','05','08','09','10','11','12','15','16','17','18','19','20','21']#'04','05','08','09','10','11','12','15','16','17','18','19','20','21'] #'04','05','08','09','10','11','12','13','15','16','17','18','19','20','21'
visit='03'


# create class to handle sorting pixels into regions
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
    
    def insert_points(self, lats, lons, power):
        
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
        if lats.shape != power.shape:
            raise ValueError("lats and intensities must be same shape")
        if lons.shape != power.shape:
            raise ValueError("lons and intensities must be same shape")
        if lats.shape != lons.shape:
            raise ValueError("lats and lons must be same shape")
        
        # arange coordinates for use in path function
        coords = np.vstack([lats.flatten(), lons.flatten()]).T
        
        # check if points within polygon
        poly = path.Path(self.verticies)
        inside_mask = poly.contains_points(coords) # this is already making the mask 
        
        # mask = np.ma.masked_where(inside_mask)
        # print(np.where(inside_mask)) 
        
        # false_array = np.zeros((720,1440), dtype=bool)
        
        '''
        set inside mask to true? + then multiply image by the mask??????
        I honestly don't know need to look into this with a better frame of mind
        '''
    
        
        '''
        need to make a self.mask function where all values False in inside_mask are set to np.nan
        or just export the mask and do this outside?
        
        ie can get self.mask = inside_mask 
        
        then on outside it'll be swirl.mask to call the mask?
        and then we can be like if false set to np.nan
        not sure this will directly work as other values will then be true but then we can overlay maybe?
        '''
        
        # store intensities to self
        self.power = power.flatten()[inside_mask]
        self.lats = lats.flatten()[inside_mask]
        self.lons = lons.flatten()[inside_mask]
        
        return
    


# actual function to call to cut out regions
def box_cutter(file,n,visit,visit_name, plotting,histogram):
    if __name__ == "__main__":
     
        '''
        doy = filename[-39:-36]
        tinti = int(filename[-25:-22])
        tint = str(tinti)
        hora = filename[11:19]
        year = filename[4:6]
        visit = 'v' + str(v) + '_20' + str(year)

        # one of the two titles for every plot
        plt.suptitle(f'Visit {prefix}{v} (DOY: {doy}/20{year}, {hora[0:2]}:{hora[3:5]}:{hora[6:]})', y=0.99, fontsize=14)
        
        we want DOYTHrMinSec
         filename is DOY-Hour-Min-Sec 
        '''
        
        # open intensities file
        infolist = fits.open(file)
        header = infolist[1].header
        image = infolist[1].data
        
        cml = header['CML']
        
        #breakpoint()
        
        # jupiter times
        start_time = parse(header['UDATE'])     # create datetime object
        try:
            dist_org = header['DIST_ORG']
            ltime = dist_org*c.au/c.c
            lighttime = dt.timedelta(seconds=ltime)
        except KeyError:
            lighttime = dt.timedelta(seconds=2524.42) 
        
        time_date = start_time - lighttime  # correct for light travel time
        
        image_flip = np.flip(image,0)
        image_final = image_flip[:(int((image.shape[0])/1)),:]
                        # if using power conversion
        powers, app_area, abs_area = rtp.rayleigh_to_power(image_final, n, visit)
        #intensity = np.where(powers <= 0, np.nan, powers)
        #print(powers)
      
        lats = np.arange(0,180,0.25)
        lons = np.arange(0,360,0.25)
        
        # make 2D grid of lat/lons
        llons, llats = np.meshgrid(lons, lats)
        
        lon_flat = llons.flatten()
        lat_flat = llats.flatten()
        
        lon_rad = np.radians(lon_flat)
        lat_rad = np.radians(lat_flat)
    
        # define each region
        # assign verticies to confine regions
        # swirl = shape([[4,100],[6,100],[22,170],[22,185],[16,185]])
        # #swirl = shape([[100,4],[100,6],[170,22],[185,22],[185,16]])
        # noon_active = shape([[24,170],[28,170],[32,190],[28,205],[24,205]])
        # #noon_active = shape([[170,24],[170,28],[190,32],[205,28],[205,24]])
        # dusk_active = shape([[0,100],[3,100],[16,188],[24,188],[24,205],[10,205],[0,160]])
        # #dusk_active = shape([[100,0],[100,3],[188,16],[188,24],[205,24],[205,10],[160,0]])
        
        if visit == '12' or visit == '15':
            # print(visit)
            # swirl = shape([[4,100],[6,100],[24,170],[24,188],[16,188]])
            # #swirl = shape([[100,4],[100,6],[170,22],[185,22],[185,16]])
            # noon_active = shape([[24,170],[28,170],[32,190],[28,205],[24,205]])
            # #noon_active = shape([[170,24],[170,28],[190,32],[205,28],[205,24]])
            # dusk_active = shape([[0,100],[4,100],[16,188],[24,188],[24,205],[10,205],[0,160]])
            # #dusk_active = shape([[100,0],[100,3],[188,16],[188,24],[205,24],[205,10],[160,0]])
            
            # # deep rarefaction
            # noon_12 = [[22,197],[28,197],[33,167],[30,150],[25,150]]
            # dusk_12 = [[8,120],[22,171],[24,150],[16,120]]
            # swirl_12 = [[0,270],[22,200],[22,171],[8,120]]
            
            # noon_15 = [[22,200],[26,200],[32,170],[31,168],[24,140]]
            # dusk_15 = [[14,140],[22,173],[24,150],[16,130]]
            # swirl_15 = [[0,270],[22,200],[22,173],[14,140]]
            
            # noon_avg = np.mean(np.array([noon_12, noon_15]),axis=0)
            # dusk_avg = np.mean(np.array([dusk_12, dusk_15]),axis=0)
            # swirl_avg = np.mean(np.array([swirl_12, swirl_15]),axis=0)
        
            # swirl = shape(swirl_avg)
            # noon_active = shape(noon_avg)  
            # dusk_active = shape(dusk_avg)
            print('correct place')
            
            dusk_active = shape([[8,270],[11,230],[23,180],[24,215],[18,230],[15,270]])
            swirl = shape([[0,90],[2,90],[15,140],[22,158],[23,180],[11,230],[8,270],[0,270]])
            noon_active = shape([[23,180],[22,158],[27,162.5],[32.5,191.5],[30.5,201],[24,215]])
            
            #breakpoint()
            
        else:
            # assign verticies to confine regions
            swirl = shape([[4,100],[6,100],[24,170],[24,188],[16,188]])
            #swirl = shape([[100,4],[100,6],[170,22],[185,22],[185,16]])
            noon_active = shape([[24,170],[28,170],[32,190],[28,205],[24,205]])
            #noon_active = shape([[170,24],[170,28],[190,32],[205,28],[205,24]])
            dusk_active = shape([[0,100],[4,100],[16,188],[24,188],[24,205],[10,205],[0,160]])
            #dusk_active = shape([[100,0],[100,3],[188,16],[188,24],[205,24],[205,10],[160,0]])
            
            
        #breakpoint()
        inside_io_footprint = shape([[0,359.9],[8,359.8],[8.4,282.3],[13.4,264.1],[19.5,251.7],[25.1,242.0],[28.5,234.7],[31.3,228.9],[34.0,224.1],[36.8,218.8],[38.7,214.3],[40.1,210.1],[40.9,205.8],
                              [41.4,201.3],[41.2,196.4],[40.9,191.5],[40.5,187.1],[39.5,182.1],[38.0,176.9],[36.7,171.1],[34.9,165.5],[33.3,160.0],[31.7,154.0],[30.0,147.4],
                              [28.6,140.4],[26.5,132.1],[24.6,122.9],[23.0,113.0],[21.7,102.5],[20.5,91.5],[19.8,79.7],[19.0,67.1],[17.7,53.1],[15.9,37.8],[8,0.1],[0,0]])
                             #shape([[81.6,77.7],[76.6,95.9],[70.5,108.3],[64.9,118.0],[61.5,125.3],[58.7,131.1],[56.0,135.9],[53.2,141.2],[51.3,145.7],[49.9,149.9],[49.1,154.2],[48.6,158.7],[48.8,163.6],[49.1,168.5],[49.5,172.9],[50.5,177.9],[52.0,183.1],[53.3,188.9],[55.1,194.5],[56.7,200.0],[58.3,206.0],[60.0,212.5],[61.4,219.6],[63.5,227.9],[65.4,237.1],[67.0,247.0],[68.3,257.5],[69.5,268.5],[70.2,280.3],[71.0,292.9],[72.3,306.9],[74.1,322.2]])

        # insert data into regions
        swirl.insert_points(llats, llons, powers)
        noon_active.insert_points(llats, llons, powers)
        dusk_active.insert_points(llats, llons, powers)
        
        inside_io_footprint.insert_points(llats,llons, powers)
        
        '''
        now you can get data within each shape with:
        intensities within the shape = shapeA.intensities
        latitudes of those intensities = shapeA.latitudes
        longitudes of those intensities = shapeA.longitudes
        '''


        #print(dusk_active.intensities)
        
        # if plotting is on this will plot the location of all points in the
        # grid in black, and then will highlight in different colours which
        # of those points belong to the specified regions
        if plotting == 'yes' or plotting == 'y' or plotting == 'Yes' or plotting == 'Y':
            rlim=40
            
            radials = np.linspace(0,rlim,6,dtype='int') # get 6 evenly spaced values from 0 to rlim
            radials = np.arange(0,rlim,10,dtype='int')
            
            # convert values to radians ffor plotting
            lats_rad_swirl = (swirl.lats)
            lat_verts_rad_swirl = np.radians(swirl.lat_verts)
            lons_rad_swirl = np.radians(swirl.lons)
            lon_verts_rad_swirl = np.radians(swirl.lon_verts)
            
            lats_rad_noon_active = (noon_active.lats)
            lat_verts_rad_noon_active = np.radians(noon_active.lat_verts)
            lons_rad_noon_active = np.radians(noon_active.lons)
            lon_verts_rad_noon_active = np.radians(noon_active.lon_verts)
            
            lats_rad_dusk_active = (dusk_active.lats)
            lat_verts_rad_dusk_active = np.radians(dusk_active.lat_verts)
            lons_rad_dusk_active = np.radians(dusk_active.lons)
            lon_verts_rad_dusk_active = np.radians(dusk_active.lon_verts)
            
        
            # creature figure
            fig = plt.figure()
            ax = fig.add_subplot(projection='polar')
            ax.set_theta_offset(np.pi / 2.0)
            
            # plot location of all grid points
            plt.scatter(lon_rad, lat_rad, color="white", s=1, label="All Powers")
            
            # swirl region
            # plt.scatter(lat_verts_rad_swirl, lon_verts_rad_swirl, color="r")
            plt.scatter(lons_rad_swirl, lats_rad_swirl, color="r", s=1, label="Swirl Region")
            
            # noon active region
            #plt.scatter(lat_verts_rad_noon_active, lon_verts_rad_noon_active, color="b")
            plt.scatter(lons_rad_noon_active, lats_rad_noon_active, color="b", s=1, label="Noon Active Region")
            
            # dusk active region
            #plt.scatter(lat_verts_rad_dusk_active, lon_verts_rad_dusk_active, color="g")
            plt.scatter(lons_rad_dusk_active, lats_rad_dusk_active, color="g", s=1, label="Dusk Active Region")
            
            plt.scatter()
            
            # confine plot to auroral region
            #ax.set_rlim(0,np.radians(40))
            ax.set_rticks(np.arange(radials[1],rlim,10,dtype='int'))           
            ax.set_rlabel_position(0) 
        
            '''
            atm the vertices don't plot properly so I've commented them out
            issue: they don't follow the polar grid when plotting
            '''
        
        
            plt.legend()
            plt.show()
            
            
        '''
        plotting for histogram of intensities of each region
        '''
            
            
        if histogram == 'yes' or plotting == 'y' or plotting == 'Yes' or plotting == 'Y':
            
            
            fig = plt.figure(figsize=(5,10))
            plt.subplots_adjust(wspace=0, hspace=0.15) 
            ax1 = fig.add_subplot(3,1,1)
            ax1.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
            ax1.hist(swirl.intensities,bins=(0,250,500,1000,1250,1500,1750,2000,2250,2500,2700,3000),color='red',label='Swirl Region')
            ax1.set_ylabel('Counts')
            ax1.legend(bbox_to_anchor=[0.95,0.95])
            
            ax2 = fig.add_subplot(3,1,2)
            ax2.hist(noon_active.intensities,bins=(0,250,500,1000,1250,1500,1750,2000,2250,2500,2700,3000),color='blue',label='Noon Active Region')
            ax2.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
            ax2.set_ylabel('Counts')
            ax2.legend(bbox_to_anchor=[0.95,0.95])
            
            ax3 = fig.add_subplot(3,1,3)
            ax3.hist(dusk_active.intensities,bins=(0,250,500,1000,1250,1500,1750,2000,2250,2500,2700,3000),color='green',label='Dusk Active Region')
            ax3.tick_params(which='both',direction='in',bottom=True, top=True, left=True, right=True)
            ax3.set_xlabel('Intensity (kR)')
            ax3.set_ylabel('Counts')
            ax3.legend(bbox_to_anchor=[0.95,0.95])
            #plt.show()
            
            
            if not os.path.exists('/Users/hannah/OneDrive - Lancaster University/aurora/histrograms/'+visit_name+'/'):
                os.makedirs('/Users/hannah/OneDrive - Lancaster University/aurora/histograms/'+visit_name+'/',exist_ok=True)
                
            filename = str(file)[-51:-5]   
            saveloc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/histograms/{visit_name}/'+filename+'.jpg')
            plt.savefig(saveloc,bbox_inches='tight',dpi=400)
            
            
        #if viewing_percent == 'yes' or viewing_percent == 'y':

            
         
        # swirl region
        # mean_swirl = np.nanmean(swirl.power)
        # median_swirl = np.nanmedian(swirl.power)
        # total_swirl = np.nansum(swirl.power)
        # min_swirl = np.nanmin(swirl.power)
        # max_swirl = np.nanmax(swirl.power)
        # range_swirl = max_swirl - min_swirl 

        # # noon active
        # mean_noon_active = np.nanmean(noon_active.power)
        # median_noon_active = np.nanmedian(noon_active.power)
        # total_noon_active = np.nansum(noon_active.power)
        # min_noon_active = np.nanmin(noon_active.power)
        # max_noon_active = np.nanmax(noon_active.power)
        # range_noon_active = max_noon_active - min_noon_active

        # # dusk active region
        # mean_dusk_active = np.nanmean(dusk_active.power)
        # median_dusk_active = np.nanmedian(dusk_active.power)
        # total_dusk_active = np.nansum(dusk_active.power)   
        # min_dusk_active = np.nanmin(dusk_active.power)
        # max_dusk_active = np.nanmax(dusk_active.power)
        # range_dusk_active = max_dusk_active - min_dusk_active
        
        swirl_power = swirl.power
        noon_active_power = noon_active.power
        dusk_active_power = dusk_active.power
        
        inside_io_footprint_power = inside_io_footprint.power
        
        # swirl_size = len(swirl_power)
        # dusk_active_size = len(dusk_active_power)
        # noon_active_size = len(noon_active_power)
        
        swirl_tot = len(swirl_power)
        noon_active_tot = len(noon_active_power)
        dusk_active_tot = len(dusk_active_power)
        
        swirl_size = ignore_nan_counter(swirl_power)
        dusk_active_size = ignore_nan_counter(dusk_active_power)
        noon_active_size = ignore_nan_counter(noon_active_power)
         
        #return median_swirl, range_swirl, median_noon_active, range_noon_active, median_dusk_active, range_dusk_active, time_date
        return powers, swirl_power, noon_active_power, dusk_active_power, time_date, cml, swirl_size, noon_active_size, dusk_active_size, swirl_tot, noon_active_tot, dusk_active_tot, inside_io_footprint_power
    
def ignore_nan_counter(data):
    count = 0
    for i in data:
        if not np.isnan(i):
            count+=1
    return count
    
# function to glob the files for each visit
def power_evaluator(visit_list,year,prefix,extra,time):
    
    # define arrays will need to store data
    dates = []
    cmls=[]
    
    tot_pixel_swirl=[]
    tot_pixel_noon_active=[]
    tot_pixel_dusk_active=[]
    
    tot_all_pix_swirl = []
    tot_all_pix_noon_active=[]
    tot_all_pix_dusk_active=[]
    
    inside_maxs = []
    inside_mins = []
    inside_medians = []
    
    swirl_all_powers = []
    noon_all_powers = []
    dusk_all_powers = []
    
    medians_total_swirl = []
    medians_total_dusk_active = []
    medians_total_noon_active = []
    
    means_total_swirl = []
    means_total_dusk_active = []
    means_total_noon_active = []
    
    # # used if want to compile all powers across all visits
    # swirl_powers_all = []
    # noon_active_powers_all = []
    # dusk_active_powers_all = []
    
    # loop for visits
    for i in visit_list:
        i = str(i)
        print(f'VISIT {i} \n \n')
        
        arch = '*_v'+ i
        ti = str('/*0'+time+'*')
        visit_name = prefix+arch[-2:]
        ab = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti) # ab = glob.glob(f'C:/Users/moralpom/phd/data/HST/Jupiter/{year}/extract/{extra}'+arch+'/nopolar'+time+ti)
        ab.sort() 
        
        # arrays to store data for image loop
        time_dates = []

        inside_power = []
        
        swirl_powers = []
        dusk_active_powers = []
        noon_active_powers = []
        
        
        # loop through each file in the visit
        for  n,i in tqdm(enumerate(ab)):
            print(n)
            print(i)
            
            filename = str(i)[-51:-5]
            visit = filename[-20:-18]
            #breakpoint()
            
            powers, swirl_power, noon_active_power, dusk_active_power, time_date, cml, pixels_swirl, pixels_noon_active, pixels_dusk_active, swirl_tot, noon_active_tot, dusk_active_tot, inside_io_footprint_power = box_cutter(i,n,visit,visit_name,'no','no')
            #powers = np.where(powers <= 0, np.nan, powers)
            print(np.where(np.isinf(powers)))
            #powers =  np.where(powers > 75*10**6,np.nan, powers)
            
            
            #### POWERS RETURNED AT THIS POINT ARE FROM A SINGLE IMAGE ####
            
            # store cml of image
            cmls.append(cml)
            # save out times on each image
            time_dates.append(time_date)
            
            
            '''
            SECTION FOR COUNTING PIXELS
            '''
            # how many pixels seen
            tot_pixel_swirl.append(pixels_swirl)
            tot_pixel_noon_active.append(pixels_noon_active)
            tot_pixel_dusk_active.append(pixels_dusk_active)
            
            # how many pixels in each region
            tot_all_pix_swirl.append(swirl_tot)
            tot_all_pix_noon_active.append(noon_active_tot)
            tot_all_pix_dusk_active.append(dusk_active_tot)
            
            
            '''
            SECTION FOR LOOKING AT WHOLE AURORAL REGION
            '''
            
            # looking at whole auroral region
            inside_io_power = np.ma.masked_invalid(inside_io_footprint_power).sum()
            inside_power.append(inside_io_power)
            
            
            '''
            SECTION FOR LOOKING AT REGIONS
            '''
            
            total_swirls = np.ma.masked_invalid(swirl_power).sum()
            swirl_powers.append(total_swirls)
            
            total_dusk = np.ma.masked_invalid(dusk_active_power).sum()
            dusk_active_powers.append(total_dusk)
            
            total_noon = np.ma.masked_invalid(noon_active_power).sum()
            noon_active_powers.append(total_noon)
            
            swirl_all_powers.append(swirl_power)
            noon_all_powers.append(noon_active_power)
            dusk_all_powers.append(dusk_active_powers)
            
            
            # medians
            
            
            
 
            
        
        ##### OUTER LOOP FOR EACH VISIT ###### - ONLY RUNS ONCE IF LOOKING AT ONE VISIT
        
        # get midpoint for date for each visit
        midpoint_date_time = time_dates[len(time_dates)//2]
        dates.append(midpoint_date_time) 
        
        
        '''
        SECTION FOR LOOKING AT WHOLE AURORAL REGION
        '''
        
        # calculate averages for all images in visit
        median_inside_io = np.nanmedian(inside_power)
        max_inside_io = np.nanmax(inside_power)
        min_inside_io = np.nanmin(inside_power)
        
        # append averages to big grid for all visits
        inside_medians.append(median_inside_io)
        inside_maxs.append(max_inside_io)
        inside_mins.append(min_inside_io)
        
        
        '''
        SECTION FOR LOOKING AT REGIONS - CURRENTLY NOT USED / RETURNED
        '''
        
        # # median - median in this case is the median for each visit from the total per image
        
        # median_total_swirl = np.nanmedian(swirl_powers)
        # medians_total_swirl.append(median_total_swirl)
        
        # median_total_dusk = np.nanmedian(dusk_active_powers)
        # medians_total_dusk_active.append(median_total_dusk)
        
        # median_total_noon = np.nanmedian(noon_active_powers)
        # medians_total_noon_active.append(median_total_noon)
        
        # mean_total_swirl = np.nanmean(swirl_powers)
        # means_total_swirl.append(mean_total_swirl)
        
        # mean_total_dusk = np.nanmean(dusk_active_powers)
        # means_total_dusk_active.append(mean_total_dusk)
        
        # mean_total_noon = np.nanmean(noon_active_powers)
        # means_total_noon_active.append(mean_total_noon)
    
        # #  append total power values from each visit into a grid for all visits
        # total_swirl_active_all.append(total_swirl)
        # total_noon_active_all.append(total_noon_active)
        # total_dusk_active_all.append(total_dusk_active)
        
        # # append power values to make list of list - list of powers per image for each visit (probably a bad idea)
        # swirl_powers_all.append(swirl_powers)
        # dusk_active_powers_all.append(dusk_active_powers)
        # noon_active_powers_all.append(noon_active_powers)
    
            
        
    return powers, cmls, dates, swirl_powers, dusk_active_powers, noon_active_powers, inside_power, inside_maxs, inside_mins, inside_medians, tot_pixel_swirl, tot_pixel_noon_active, tot_pixel_dusk_active, tot_all_pix_swirl, tot_all_pix_noon_active, tot_all_pix_dusk_active, swirl_all_powers, noon_all_powers, dusk_all_powers#, medians_total_swirl, medians_total_dusk_active, medians_total_noon_active, means_total_swirl, means_total_dusk_active, means_total_noon_active


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

    
powers, cmls, dates, swirl_powers, dusk_active_powers, noon_active_powers, inside_power, inside_maxs, inside_mins, inside_medians, tot_pixel_swirl, tot_pixel_noon_active, tot_pixel_dusk_active, tot_all_pix_swirl, tot_all_pix_noon_active, tot_all_pix_dusk_active, swirl_all_powers, noon_all_powers, dusk_all_powers = power_evaluator(visit_list,year,prefix,extra,time) #, medians_total_swirl, medians_total_dusk_active, medians_total_noon_active, means_total_swirl, means_total_dusk_active, means_total_noon_active

##### 

'''
SECTION FOR COUNTING PIXELS
'''

# change to array format to calculate percentages
tot_pixel_swirl = np.array(tot_pixel_swirl)
tot_all_pix_swirl = np.array(tot_all_pix_swirl)

tot_pixel_noon_active = np.array(tot_pixel_noon_active)
tot_all_pix_noon_active = np.array(tot_all_pix_noon_active)

tot_pixel_dusk_active = np.array(tot_pixel_dusk_active)
tot_all_pix_dusk_active = np.array(tot_all_pix_dusk_active)

# percentage calculation for pixels seen
percentage_swirl = (tot_pixel_swirl/tot_all_pix_swirl)*100
percentage_noon_active = (tot_pixel_noon_active/tot_all_pix_noon_active)*100
percentage_dusk_active = (tot_pixel_dusk_active/tot_all_pix_dusk_active)*100




###### DATAFRAMES ########


'''
TOTAL POWER AVERAGES ACROSS AURORAL REGION
'''

# # put in new dataframe
# average_power_data_df = pd.DataFrame()

# average_power_data_df = average_power_data_df.assign(Date=dates)
# average_power_data_df = average_power_data_df.assign(Medians=inside_medians)
# average_power_data_df = average_power_data_df.assign(Mins=inside_mins)
# average_power_data_df = average_power_data_df.assign(Maxs=inside_maxs)

# # export dataframe to read into other files
# average_power_data_df.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/average_power_data_intensity_inside_io.csv',index=False)


'''
TOTAL POWER ACROSS AURORAL REGION PER IMAGE
'''

# total_power_data_df = pd.DataFrame()

# total_power_data_df = total_power_data_df.assign(CML=cmls)
# total_power_data_df = total_power_data_df.assign(Total_Power=inside_power)

# total_power_data_df.to_csv(f'/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/total_power_data_intensity_inside_io_{visit}.csv',index=False)



'''
PIXEL COUNT PER REGION
'''

# # data frame for total pixel percentage/ratio
# total_pixels_df = pd.DataFrame()

# total_pixels_df = total_pixels_df.assign(CML=cmls)
# total_pixels_df = total_pixels_df.assign(Pixels_Seen_Swirl_Region=tot_pixel_swirl)
# total_pixels_df = total_pixels_df.assign(Pixels_Tot_Swirl_Region=tot_all_pix_swirl)
# total_pixels_df = total_pixels_df.assign(Percentage_Seen_Swirl_Region=percentage_swirl)
# total_pixels_df = total_pixels_df.assign(Pixels_Seen_Noon_Active=tot_pixel_noon_active)
# total_pixels_df = total_pixels_df.assign(Pixels_Tot_Noon_active_Region=tot_all_pix_noon_active)
# total_pixels_df = total_pixels_df.assign(Percentage_Seen_Noon_Active_Region=percentage_noon_active)
# total_pixels_df = total_pixels_df.assign(Pixels_Seen_Dusk_Active=tot_pixel_dusk_active)
# total_pixels_df = total_pixels_df.assign(Pixels_Tot_Dusk_Active_Region=tot_all_pix_dusk_active)
# total_pixels_df = total_pixels_df.assign(Percentage_Seen_Dusk_Active_Region=percentage_dusk_active)

# # export
# total_pixels_df.to_csv('/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/main_regions_pixels_seen.csv',index=False)


'''
ALL POWER VALUES PER REGION
'''

# # dataframe for actual power values
# power_values_df = pd.DataFrame()

# power_values_df = power_values_df.assign(CML=cmls)
# power_values_df = power_values_df.assign(Swirl_Powers=swirl_all_powers)
# power_values_df = power_values_df.assign(Dusk_Active_Powers=dusk_all_powers)
# power_values_df = power_values_df.assign(Noon_Active_Powers=noon_all_powers)

# power_values_df.to_csv(f'/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/power_values_{visit}_df.csv',index=False)


'''
TOTAL POWER PER REGION PER IMAGE
'''

#dataframe for total power per image
total_power_regions_df = pd.DataFrame()

total_power_regions_df = total_power_regions_df.assign(CML=cmls)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Swirl=swirl_powers)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Noon_Active=noon_active_powers)
total_power_regions_df = total_power_regions_df.assign(Total_Power_Dusk_Active=dusk_active_powers)

total_power_regions_df.to_csv(f'/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/total_power_values_regions_test{visit}_df.csv',index=False)
#(f'/Users/hannah/OneDrive - Lancaster University/aurora/python_scripts/dataframes/total_power_values_regions_{visit}_df.csv',index=False)
