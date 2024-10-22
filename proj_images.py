"""
Created on Tue May 23 13:01:48 2023
updated version
@author: dmoral, edited by @hjoyce

this script generates easy to visualise images of the processed auroral images for HST visits in either intensity
or power depending on user input - these images are then turned into looping gifs, displaying time of visit and also
showing location of the CML line. This script is fully automated to user inputs and also provides the option to 
draw boxes on the different polar auroral regions to visualise their location as well as 'boxing off' to whole auroral 
region if desired 

- main emission boxes are currently commented out as they need to be improved upon and are not filled in yet
- southern hemisphere currently has no options to draw on any regional boxes
- local time view of images also has no current options to draw on any regional boxes

older versions of this file as proj_images and proj_images_boxes (v2).
"""
# import all the necessary modules
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

datapath = Path.cwd() # a function from the path module - identifies current working directory

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
     
     def insert_points(self, lats, lons, intensities):
         
         #whole_grid_mask = np.full((int(intensities.shape[0],intensities.shape[1])), False)
         
         """
         Saves intensities whoes latitudes/longitudes are contained withing
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
         if lats.shape != intensities.shape:
             raise ValueError("lats and intensities must be same shape")
         if lons.shape != intensities.shape:
             raise ValueError("lons and intensities must be same shape")
         if lats.shape != lons.shape:
             raise ValueError("lats and lons must be same shape")
         
         # arange coordinates for use in path function
         coords = np.vstack([lats.flatten(), lons.flatten()]).T
         '''
         if path.Path(self.verticies) == True:
             poly = path.Path(self.verticies)
             inside_mask = poly.contains_points(coords)
        '''
             
             
         # check if points within polygon
         poly = path.Path(self.verticies)
         inside_mask = poly.contains_points(coords) # this is already making the mask 
         
         
         # store intensities to self
         self.intensities = intensities.flatten()[inside_mask]
         self.lats = lats.flatten()[inside_mask]
         self.lons = lons.flatten()[inside_mask]
         
         
         
         return
     
 # actual function to call to cut out regions
    
#this is the funcion used for plotting the images
def moind(aa, visit, n, header, filename, prefix, dpi = 300, crop = 1, rlim = 30, fixed = 'lon',
          hemis = 'North', full = True, moonfp = True, output = 'w', regions = True, fill_in = True, whole_region_box = True, photo=0):
    # dpi? crop? rlim? 
    # fixed - lon (long) or lat, hemis is North/N or South/S, moonfp is what moon footprint or not
    
    # scatter plot:
    cml = header['CML'] # store CML vcentral meridian longitude, (viewing angle from Earth) value
    dece = header['DECE']
    exp_time = header['EXPT']
    iolon = header['IOLON'] # io longitudes - start
    iolon1 = header['IOLON1'] # io middle
    iolon2 = header['IOLON2'] # io end
    eulon = header['EULON'] # europa longitudes
    eulon1 = header['EULON1']
    eulon2 = header['EULON2']
    galon = header['GALON'] # ganymede longitudes
    galon1 = header['GALON1']
    galon2 = header['GALON2']
    # need these defined
    
    # jupiter times
    start_time = parse(header['UDATE'])     # create datetime object
    try:
        dist_org = header['DIST_ORG']
        ltime = dist_org*c.au/c.c
        lighttime = dt.timedelta(seconds=ltime)
    except KeyError:
        lighttime = dt.timedelta(seconds=2524.42) 
    exposure = dt.timedelta(seconds=exp_time)
    
    start_time_jup = start_time - lighttime       # correct for light travel time
    end_time_jup = start_time_jup + exposure      # end of exposure time
    mid_ex_jup = start_time_jup + (exposure/2.)   # mid-exposure
    
    # plot
    latbins = np.radians(np.linspace(-90, 90, num=aa.shape[0])) # latitude, radians conversion
    lonbins = np.radians(np.linspace(0, 360, num=aa.shape[1])) # longitude, radians conversion
    
    # polar variables
    rho   = np.linspace(0,180, num=int(aa.shape[0])) # colatitude vector with image pixel resolution steps
    theta = np.linspace(0,2*np.pi,num=aa.shape[1])
    if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
        rho = rho[::-1]      # if south, changes the orientation 

    # creating the mask
    mask = np.zeros((int(aa.shape[0]),aa.shape[1]))   # rows by columns, fill with zeroes for now
    cmlr = np.radians(cml)                  # convert CML to radians
    dec = np.radians(dece)        # convert declination angle to radians
    
    # filling the mask
    for i in range(0,mask.shape[0]):
        mask[i,:] = np.sin(latbins[i])*np.sin(dec) + np.cos(latbins[i])*np.cos(dec)*np.cos(lonbins-cmlr)
    
    mask = np.flip(mask,axis=1) # flip the mask horizontally, not sure why this is needed?
    cliplim = np.cos(np.radians(89))       # set a minimum required vector normal surface-sun angle
    clipind = np.squeeze([mask >= cliplim])
    
    # applying the mask
    aa[clipind == False] = np.nan
    if __name__ == "__main__":
        
        # open intensities file
        image2 = aa
        image_flip2 = np.flip(image2,0)
      
        lats = np.arange(0,180,0.25)
        lons = np.arange(0,360,0.25)
        
        intensities = image_flip2[:(int((aa.shape[0])/crop)),:]
        #breakpoint()
        # intensities = np.roll(intensities,180*4,axis=1) - if south
        #intensities = intensities.T #image_flip2.T
        
        # make 2D grid of lat/lons
        llons, llats = np.meshgrid(lons, lats)
        
        # lon_flat = llons.flatten()
        # lat_flat = llats.flatten()
        
        #lon_rad = np.radians(lon_flat)
        #lat_rad = np.radians(lat_flat)
    
        # define each region
        v_visit = visit
        
        print(visit)
        if visit == '12' or visit == '15':        
            # # deep rarefaction
            # noon_12 = [[22,163],[28,163],[33,193],[30,210],[25,210]] #[[22,197],[28,197],[33,167],[30,150],[25,150]]
            # dusk_12 = [[8,240],[22,189],[24,210],[16,240]]#[[8,120],[22,171],[24,150],[16,120]]
            # swirl_12 = [[0,90],[22,160],[22,189],[8,240]]#[[0,270],[22,200],[22,171],[8,120]]
            
            # noon_15 = [[22,160],[26,160],[32,190],[31,192],[24,220]]#[[22,200],[26,200],[32,170],[31,168],[24,140]]
            # dusk_15 = [[14,220],[22,187],[24,210],[16,230]]#[14,140],[22,173],[24,150],[16,130]]
            # swirl_15 = [[0,90],[22,160],[22,187],[14,220]]#[[0,270],[22,200],[22,173],[14,140]]
            
            # noon_avg = np.mean(np.array([noon_12, noon_15]),axis=0)
            # dusk_avg = np.mean(np.array([dusk_12, dusk_15]),axis=0)
            # swirl_avg = np.mean(np.array([swirl_12, swirl_15]),axis=0)
            # #breakpoint()    
        
            # swirl = shape(swirl_avg)
            # noon_active = shape(noon_avg)  
            # dusk_active = shape(dusk_avg)
            
            dusk_active = shape([[8,270],[11,230],[23,180],[24,215],[18,230],[15,270]])
            swirl = shape([[0,90],[2,90],[15,140],[22,158],[23,180],[11,230],[8,270],[0,270]])
            noon_active = shape([[23,180],[22,158],[27,162.5],[32.5,191.5],[30.5,201],[24,215]])
            
        else:
            #assign verticies to confine regions
            swirl = shape([[4,100],[6,100],[24,170],[24,188],[16,188]])
            #swirl = shape([[100,4],[100,6],[170,22],[185,22],[185,16]])
            noon_active = shape([[24,170],[28,170],[32,190],[28,205],[24,205]])
            #noon_active = shape([[170,24],[170,28],[190,32],[205,28],[205,24]])
            dusk_active = shape([[0,100],[4,100],[16,188],[24,188],[24,205],[10,205],[0,160]])
            #dusk_active = shape([[100,0],[100,3],[188,16],[188,24],[205,24],[205,10],[160,0]])
        
        #breakpoint()
        if whole_region_box == True:
            inside_io = shape([[0,359.9],[8,359.8],[8.4,282.3],[13.4,264.1],[19.5,251.7],[25.1,242.0],[28.5,234.7],[31.3,228.9],[34.0,224.1],[36.8,218.8],[38.7,214.3],[40.1,210.1],[40.9,205.8],
                              [41.4,201.3],[41.2,196.4],[40.9,191.5],[40.5,187.1],[39.5,182.1],[38.0,176.9],[36.7,171.1],[34.9,165.5],[33.3,160.0],[31.7,154.0],[30.0,147.4],
                              [28.6,140.4],[26.5,132.1],[24.6,122.9],[23.0,113.0],[21.7,102.5],[20.5,91.5],[19.8,79.7],[19.0,67.1],[17.7,53.1],[15.9,37.8],[8,0.1],[0,0]])

            inside_io.insert_points(llats,llons,intensities)
            
            lons_inside_io = np.radians(inside_io.lons)
            lats_inside_io = (inside_io.lats)
        
        # insert data into regions
        swirl.insert_points(llats, llons,intensities)
        noon_active.insert_points(llats, llons, intensities)
        dusk_active.insert_points(llats, llons, intensities)
        
        
        
        '''
        now you can get data within each shape with:
        intensities within the shape = shapeA.intensities
        latitudes of those intensities = shapeA.latitudes
        longitudes of those intensities = shapeA.longitudes
        '''
        
        # if plotting is on this will plot the location of all points in the
        # grid in black, and then will highlight in different colours which
        # of those points belong to the specified regions
            
        
        # convert values to radians ffor plotting
        lats_rad_swirl = (swirl.lats)
        #lat_verts_rad_swirl = np.radians(swirl.lat_verts)
        lons_rad_swirl = np.radians(swirl.lons)
        #lon_verts_rad_swirl = np.radians(swirl.lon_verts)
        
        lats_rad_noon_active = (noon_active.lats)
        #lat_verts_rad_noon_active = np.radians(noon_active.lat_verts)
        lons_rad_noon_active = np.radians(noon_active.lons)
        #lon_verts_rad_noon_active = np.radians(noon_active.lon_verts)
        
        lats_rad_dusk_active = (dusk_active.lats)
        #lat_verts_rad_dusk_active = np.radians(dusk_active.lat_verts)
        lons_rad_dusk_active = np.radians(dusk_active.lons)
        ##lon_verts_rad_dusk_active = np.radians(dusk_active.lon_verts)
    

#   KR_MIN = cliplim
    #aa[aa < KR_MIN] = cliplim
##########################################################################

    # plotting the polar projection of the image
    plt.figure(figsize=(7,6))
        
    ax = plt.subplot(projection='polar')
    radials = np.linspace(0,rlim,6,dtype='int') # get 6 evenly spaced values from 0 to rlim
    radials = np.arange(0,rlim,10,dtype='int') # rearange values w/ step size of 10

    # shifting the image to have CML pointing southwards in the image
    if fixed == 'lon':
        image_centred = aa
        im_flip = np.flip(image_centred,0) # reverse the image along the longitudinal (x, theta) axis
    # cropping the image
        if output == 'kr':
            if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
                corte = im_flip[:(int((aa.shape[0])/crop)),:] # if crop=1, nothing changes # NO ??
                corte = np.roll(corte,180*4,axis=1)
                plt.plot(np.roll([np.radians(180-cml),np.radians(180-cml)],180*4),[0, 180], 'r--', lw=1.2) # cml
                plt.text(np.radians(180-cml), 3+rlim, 'CML', fontsize=11, color='r',
                     horizontalalignment='center', verticalalignment='center', fontweight='bold') 
            else: # if North
                corte = im_flip[:(int((aa.shape[0])/crop)),:]
                plt.plot(np.roll([np.radians(360-cml),np.radians(360-cml)],180*4),[0, 180], 'r--', lw=1.2) # cml
                plt.text(np.radians(360-cml), 3+rlim, 'CML', fontsize=11, color='r',
                     horizontalalignment='center', verticalalignment='center', fontweight='bold') 
        
        elif output == 'w':
            if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
                corte = im_flip[:(int((aa.shape[0])/crop)),:] # if crop=1, nothing changes # NO ??
                corte = np.roll(corte,180*4,axis=1)
                corte, app_area, abs_area = rtp.rayleigh_to_power(corte, n, visit)
                plt.plot(np.roll([np.radians(180-cml),np.radians(180-cml)],180*4),[0, 180], 'r--', lw=1.2) # cml
                plt.text(np.radians(180-cml), 3+rlim, 'CML', fontsize=11, color='r',
                        horizontalalignment='center', verticalalignment='center', fontweight='bold') 
            else: # if North
                corte = im_flip[:(int((aa.shape[0])/crop)),:]
                corte, app_area, abs_area = rtp.rayleigh_to_power(corte, n, visit)
                plt.plot(np.roll([np.radians(360-cml),np.radians(360-cml)],180*4),[0, 180], 'r--', lw=1.2) # cml
                plt.text(np.radians(360-cml), 3+rlim, 'CML', fontsize=11, color='r',
                        horizontalalignment='center', verticalalignment='center', fontweight='bold') 
            
        # label image
        if full == True:
            if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
                ax.set_xticklabels(['180°','','','','','','','','','90°','','','','','','','','',
                                    '0°','','','','','','','','','270°','','','','','','','',''], fontweight='bold')   
                ax.set_xticks(np.linspace(0,2*np.pi,37))                      # but set grid spacing#
                poshem = np.radians(45) # position of the "S" marker
    
            else: # if North
                ax.set_xticklabels(['0°','','','','','','','','','270°','','','','','','','','',
                                    '180°','','','','','','','','','90°','','','','','','','',''], fontweight='bold')   
                ax.set_xticks(np.linspace(0,2*np.pi,37))                      # but set grid spacing
                poshem = np.radians(45) # position of the "N" marker
                #ax.set_xticklabels(['0°','','270°','','180°','','90°'], fontweight='bold')
            
            ytickl = []
            for i in radials:
                ytickl.append(str(i)+'°')

            ax.set_yticklabels(ytickl,color='w',fontsize=10) #, weight='bold') # turn off auto lat labels
            ax.set_rticks(np.arange(radials[1],rlim,10,dtype='int'))           
            ax.set_rlabel_position(0)   # position of the radial labels  
            shrink = 1. # size of the colorbar
            possub = 1.05 # position in the y axis of the subtitle
                   
        else:
            if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
                ax.set_xticklabels(['90°','45°','0°','315°','270°'], fontweight='bold')
                ax.set_xticks(np.linspace(np.pi/2,3*np.pi/2,5))                       # but set grid spacing
                poshem = np.radians(-135) # position of the "S" marker
    
            else:  # if North
                ax.set_xticklabels(['270°','225°','180°','135°','90°'], fontweight='bold')
                ax.set_xticks(np.linspace(np.pi/2,3*np.pi/2,5))  
                poshem = np.radians(135) # position of the "N" marker
            
            shrink = 0.75 # size of the colorbar
            possub = 1.03 # position in the y axis of the subtitle
            ax.set_thetalim([np.pi/2,3*np.pi/2])
            ax.set_yticklabels(['',str(radials[1])+'°',str(radials[2])+'°',str(radials[3])+'°',
                                str(radials[4])+'°',str(radials[5])+'°'],color='w',fontsize=10) #, weight='bold') # turn off auto lat labels
            ax.set_rticks(np.linspace(radials[1],rlim,5,dtype='int'))                       # but set grid spacing
        ax.set_theta_zero_location("N")                   # set angle 0.0 to top of plot

    elif fixed == 'lt':
        image_centred = np.roll(aa,int(cml-180.)*4,axis=1)
        im_flip = np.flip(image_centred,0)
        corte = im_flip[:(int((aa.shape[0])/crop)),:] # if crop=1, nothing changes
        if output == 'w':
            corte, app_area, abs_area = rtp.rayleigh_to_power(corte, n, visit)
            
        if full == True:
            # IndexError: index 4 is out of bounds for axis 0 with size 4 - investigate this
            ax.set_yticklabels(['',str(radials[1])+'°',str(radials[2])+'°',str(radials[3])+'°'], color='w', fontsize=10) # turn off auto lat labels #,str(radials[4])+'°' ,str(radials[5])+'°'
            ax.set_rticks(np.linspace(radials[1],rlim,5,dtype='int'))                       # but set grid spacing
            ax.set_xticklabels(['00','','','03','','','06','','','09','','','12',
                '','','15','','','18','','','21','',''], fontweight='bold')
            ax.set_xticks(np.linspace(0,2*np.pi,25))                      # but set grid spacing
            ax.set_rlabel_position(0)         
            shrink = 1. # size of the colorbar
            possub = 1.05 # position in the y axis of the subtitle

            if hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
                ax.set_theta_zero_location("N")                   # set angle 0.0 to bottom of plot
                poshem = np.radians(45) # position of the "S" marker
            else:  # if North
                ax.set_theta_zero_location("N")                   # set angle 0.0 to bottom of plot
                poshem = np.radians(45) # position of the "N" marker
        else: # if not full, but half circle
            ax.set_xticklabels(['06','9','12','15','18'], fontweight='bold')
            ax.set_xticks(np.linspace(np.pi/2,3*np.pi/2,5))                       # but set grid spacing
            ax.set_thetalim([np.pi/2,3*np.pi/2])
            ax.set_yticklabels(['',str(radials[1]),'',str(radials[3]),'',str(radials[5])],color='w', fontsize=10) # turn off auto lat labels
            ax.set_rticks(np.linspace(radials[1],rlim,5,dtype='int'))                       # but set grid spacing
            shrink = 0.75 # size of the colorbar
            possub = 1.02 # position in the y axis of the subtitle
            poshem = np.radians(135) # position of the "S" marker

    ax.set_facecolor('k') # black background
    ax.set_rlim([0,rlim]) # max colat range
    ax.tick_params(axis='both',pad=2.)    # shift position of LT labels
    plt.rgrids(radials) #, color='white')

    # naming variables
    v = filename[-20:-18]
    #print(v)
    visita = filename[-51:-5]
    if moonfp == True:
        visita = 'mfp_'+visita
    doy = filename[-39:-36]
    tinti = int(filename[-25:-22])
    tint = str(tinti)
    hora = filename[11:19]
    year = filename[4:6]
    visit = 'v' + str(v) + '_20' + str(year)
    time_jup = start_time_jup.strftime('%H:%M:%S') 

    # title text displayed on gif
    plt.suptitle(f'Visit {prefix}{v} - DOY: {doy}/20{year}, {time_jup}',y=0.99, fontsize=14)

    cmlround = np.round(cml, decimals=1)
    
    # the other title + the 0° longitudinal meridian for the LT fixed case
    if fixed == 'lon':
        plt.title(f'Integration time={tint} seconds. CML: {cmlround}°', y=possub, fontsize=12)
    elif fixed == 'lt':
        plt.title(f'Fixed LT. Integration time={tint} s. CML: {cmlround}°', fontsize=12)
        plt.text(np.radians(cml)+np.pi, 4+rlim, '0°', color='coral', fontsize=12,
                 horizontalalignment='center', verticalalignment='bottom', fontweight='bold')
        plt.plot([np.radians(cml)+np.pi,np.radians(cml)+np.pi],[0, 180], color='coral', \
                  path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                  linestyle='-.', lw=1) # prime meridian (longitude 0)

    max_watt = 5*10**14
    # actual plot and colorbar (change the vmin and vmax to play with the limits
    # of the colorbars, recommended to enhance/saturate certain features)
    if int(tint) < 30:
        plt.pcolormesh(theta,rho[:(int((aa.shape[0])/crop))],corte,norm=LogNorm(vmin=10., vmax=1500.), cmap='inferno')
        cbar = plt.colorbar(ticks=[10.,40.,100.,200.,400.,800.,1500.], shrink=shrink, pad=0.06)
        cbar.ax.set_yticklabels(['10','40','100','200','400','800','1500'])
    else:
        if output == 'kr':
            plt.pcolormesh(theta,rho[:(int((aa.shape[0])/crop))],corte,norm=LogNorm(vmin=10., vmax=3000.))
            cbar = plt.colorbar(ticks=[10.,40.,100.,200.,400.,1000.,3000.], shrink=shrink, pad=0.06)
            cbar.ax.set_yticklabels(['10','40','100','200','400','1000','3000'])
            cbar.ax.set_ylabel('Intensity [kR]', rotation=270.)
        elif output == 'w':
            plt.pcolormesh(theta,rho[:(int((aa.shape[0])/crop))],corte, norm=LogNorm(vmin=10**6, vmax=10**9))
            cbar = plt.colorbar(shrink=shrink, pad=0.06)
            cbar.ax.set_ylabel('Auroral Power [W]', rotation=270.)


#####################################################
    
    # grids (major and minor)
    plt.grid(True, which='major', color='w', alpha=0.6, linestyle='-')
    plt.minorticks_on()
    plt.grid(True, which='minor', color='w', alpha=0.2, linestyle='--')
    
    # stronger meridional lines for the 0, 90, 180, 270 degrees:
    plt.plot([np.radians(0),np.radians(0)],[0, 180], 'w', lw=0.9) 
    plt.plot([np.radians(90),np.radians(90)],[0, 180], 'w', lw=0.9) 
    plt.plot([np.radians(180),np.radians(180)],[0, 180], 'w', lw=0.9) 
    plt.plot([np.radians(270),np.radians(270)],[0, 180], 'w', lw=0.9) 

    # deprecated variable (maybe useful for the South/LT fixed definition?)
    shift = 0# cml-180.

    # print which hemisphere are we in:
    plt.text(poshem, 1.3*rlim, str(hemis).capitalize(), fontsize=21, color='k', 
             horizontalalignment='center', verticalalignment='center', fontweight='bold') 
    
    # io footprint cut-off drawing
    if whole_region_box == True:
        io_lon_1 = np.linspace(np.radians(282.3), np.radians(264.1),100)
        io_lat_1 = np.linspace(8.4,13.4,100)

        io_lon_2 = np.linspace(np.radians(264.1), np.radians(251.7),100)
        io_lat_2 = np.linspace(13.4,19.5,100)

        io_lon_3 = np.linspace(np.radians(251.7), np.radians(242),100)
        io_lat_3 = np.linspace(19.5,25.1,100)

        io_lon_4 = np.linspace(np.radians(242), np.radians(234.7),100)
        io_lat_4 = np.linspace(25.1,28.5,100)

        io_lon_5 = np.linspace(np.radians(234.7), np.radians(228.9),100)
        io_lat_5 = np.linspace(28.5,31.3,100)

        io_lon_6 = np.linspace(np.radians(228.9), np.radians(224.1),100)
        io_lat_6 = np.linspace(31.3,34,100)

        io_lon_7 = np.linspace(np.radians(224.1), np.radians(218.8),100)
        io_lat_7 = np.linspace(34,36.8,100)

        io_lon_8 = np.linspace(np.radians(218.8), np.radians(214.3),100)
        io_lat_8 = np.linspace(36.8,38.7,100)

        io_lon_9 = np.linspace(np.radians(214.3), np.radians(210.1),100)
        io_lat_9 = np.linspace(38.7,40.1,100)

        io_lon_10 = np.linspace(np.radians(210.1), np.radians(205.8),100)
        io_lat_10 = np.linspace(40.1,40.9,100)

        io_lon_11 = np.linspace(np.radians(205.8), np.radians(201.3),100)
        io_lat_11 = np.linspace(40.9,41.4,100)

        io_lon_12 = np.linspace(np.radians(201.3), np.radians(196.4),100)
        io_lat_12 = np.linspace(41.4,41.2,100)

        io_lon_13 = np.linspace(np.radians(196.4), np.radians(191.5),100)
        io_lat_13 = np.linspace(41.2,40.9,100)

        io_lon_14 = np.linspace(np.radians(191.5), np.radians(187.1),100)
        io_lat_14 = np.linspace(40.9,40.5,100)

        io_lon_15 = np.linspace(np.radians(187.1), np.radians(182.1),100)
        io_lat_15 = np.linspace(40.5,39.5,100)

        io_lon_16 = np.linspace(np.radians(182.1), np.radians(176.9),100)
        io_lat_16 = np.linspace(39.5,38.0,100)

        io_lon_17 = np.linspace(np.radians(176.9), np.radians(171.1),100)
        io_lat_17 = np.linspace(38.0,36.7,100)

        io_lon_18 = np.linspace(np.radians(171.1), np.radians(165.5),100)
        io_lat_18 = np.linspace(36.7,34.9,100)

        io_lon_19 = np.linspace(np.radians(165.5), np.radians(160.0),100)
        io_lat_19 = np.linspace(34.9,33.3,100)

        io_lon_20 = np.linspace(np.radians(160.0), np.radians(154.0),100)
        io_lat_20 = np.linspace(33.3,31.7,100)

        io_lon_21 = np.linspace(np.radians(154.0), np.radians(147.4),100)
        io_lat_21 = np.linspace(31.7,30.0,100)

        io_lon_22 = np.linspace(np.radians(147.4), np.radians(140.4),100)
        io_lat_22 = np.linspace(30.0,28.6,100)
    
        io_lon_23 = np.linspace(np.radians(140.4), np.radians(132.1),100)
        io_lat_23 = np.linspace(28.6,26.5,100)

        io_lon_24 = np.linspace(np.radians(132.1), np.radians(122.9),100)
        io_lat_24 = np.linspace(26.5,24.6,100)

        io_lon_25 = np.linspace(np.radians(122.9), np.radians(113.0),100)
        io_lat_25 = np.linspace(24.6,23.0,100)

        io_lon_26 = np.linspace(np.radians(113.0), np.radians(102.5),100)
        io_lat_26 = np.linspace(23.0,21.7,100)

        io_lon_27 = np.linspace(np.radians(102.5), np.radians(91.5),100)
        io_lat_27 = np.linspace(21.7,20.5,100)

        io_lon_28 = np.linspace(np.radians(91.5), np.radians(79.7),100)
        io_lat_28 = np.linspace(20.5,19.8,100)

        io_lon_29 = np.linspace(np.radians(79.7), np.radians(67.1),100)
        io_lat_29 = np.linspace(19.8,19.0,100)

        io_lon_30 = np.linspace(np.radians(67.1), np.radians(53.1),100)
        io_lat_30 = np.linspace(19.0,17.7,100)

        io_lon_31 = np.linspace(np.radians(53.1), np.radians(37.8),100)
        io_lat_31 = np.linspace(17.7, 15.9,100)

        io_lon_32 = np.linspace(np.radians(0.1), np.radians(37.8),100)
        io_lat_32 = np.linspace(8, 15.9,100)
    
        io_lon_33 = np.linspace(np.radians(282.3), np.radians(359.8),100)
        io_lat_33 = np.linspace(8.4,8,100)
    
        io_lon_34 = np.linspace(np.radians(359.8),np.radians(359.9),100)
        io_lat_34 = np.linspace(8,0,100)
    
        io_lon_35 = np.linspace(np.radians(0.1), np.radians(0),100)
        io_lat_35 = np.linspace(8,0,100)
        #breakpoint()
        
        #breakpoint()
        plt.plot(io_lon_1, io_lat_1, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_2, io_lat_2, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_3, io_lat_3, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_4, io_lat_4, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_5, io_lat_5, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_6, io_lat_6, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_7, io_lat_7, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_8, io_lat_8, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_9, io_lat_9, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_10, io_lat_10, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_11, io_lat_11, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_12, io_lat_12, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_13, io_lat_13, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_14, io_lat_14, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_15, io_lat_15, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_16, io_lat_16, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_17, io_lat_17, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_18, io_lat_18, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_19, io_lat_19, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_20, io_lat_20, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_21, io_lat_21, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_22, io_lat_22, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_23, io_lat_23, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_24, io_lat_24, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_25, io_lat_25, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_26, io_lat_26, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_27, io_lat_27, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_28, io_lat_28, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_29, io_lat_29, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_30, io_lat_30, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_31, io_lat_31, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_32, io_lat_32, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_33, io_lat_33, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_34, io_lat_34, color='peachpuff', linestyle='dotted', lw=2.5)
        plt.plot(io_lon_35, io_lat_35, color='peachpuff', linestyle='dotted', lw=2.5)
        
        plt.scatter(lons_inside_io, lats_inside_io, color='peachpuff',alpha=0.03, s=1, zorder=6, label='Auroral Emission')
        
    
    # drawing the regions (if marked; only available for the North so far)
    if regions == True:
        
        # ------------ for visits 12 & 15 --------------
        # swirl
        swirls = np.linspace(np.radians(90+shift),np.radians(140+shift),200)
        swirl_line = np.linspace(2,15,200)
        swirls2 = np.linspace(np.radians(140+shift),np.radians(158+shift),200)
        swirl_line2 = np.linspace(15,22,200)
        swirls3 = np.linspace(np.radians(158+shift),np.radians(180+shift),200)
        swirl_line3 = np.linspace(22,23,200)
        swirls4 = np.linspace(np.radians(180+shift),np.radians(230+shift),200)
        swirl_line4 = np.linspace(23,11,200)
        swirls5 = np.linspace(np.radians(230+shift),np.radians(270+shift),200)
        swirl_line5 = np.linspace(11,8,200)
        
        border1 = np.linspace(0,2,200)
        border2 = np.linspace(8,0,200)
        
        # need two lines plotted 0 - 2 at 90 and 8 - 0 at 270
        
        # dusk
        dusks = np.linspace(np.radians(270+shift),np.radians(230+shift),200)
        dusk_line = np.linspace(8,11,200)
        dusks2 = np.linspace(np.radians(230+shift),np.radians(180+shift),200)
        dusk_line2 = np.linspace(11,23,200)
        dusks3 = np.linspace(np.radians(180+shift),np.radians(215+shift),200)
        dusk_line3 = np.linspace(23,24,200)
        dusks4 = np.linspace(np.radians(215+shift),np.radians(230+shift),200)
        dusk_line4 = np.linspace(24,18,200)
        dusks5 = np.linspace(np.radians(230+shift),np.radians(270+shift),200)
        dusk_line5 = np.linspace(18,15,200)
        
        border3 = np.linspace(8,15,200)
        
        # noon
        noons = np.linspace(np.radians(180+shift),np.radians(158+shift),200)
        noon_line = np.linspace(23,22,200)
        noons2 = np.linspace(np.radians(158+shift),np.radians(162.5+shift),200)
        noon_line2 = np.linspace(22,27,200)
        noons3 = np.linspace(np.radians(162.5+shift),np.radians(191.5+shift),200)
        noon_line3 = np.linspace(27,32.5,200)
        noons4 = np.linspace(np.radians(191.5+shift),np.radians(201+shift),200)
        noon_line4 = np.linspace(32.5,30.5,200)
        noons5 = np.linspace(np.radians(201+shift),np.radians(215+shift),200)
        noon_line5 = np.linspace(30.5,24,200)
        noons6 = np.linspace(np.radians(215+shift),np.radians(180+shift),200)
        noon_line6 = np.linspace(24,23,200)
        
        
        # ----------- for other visits ---------
        # # dawn arc
        # lon_dawn = np.linspace(np.radians(180+shift),np.radians(100+shift),200)
        # uplat_dawn = np.linspace(33, 11, 200)
        # downlat_dawn = np.linspace(39, 20, 200)
        
        # swirl region
        swirl = np.linspace(np.radians(170+shift),np.radians(100+shift),200)
        swirl2 = np.linspace(np.radians(188+shift),np.radians(100+shift),200)
        swirl_linee = np.linspace(24, 6, 200)
        swirl_linee2 = np.linspace(16, 4, 200)
        upnoon = np.linspace(np.radians(188+shift),np.radians(170+shift),200)
        
        # ------------- testing -----------
        # # avg swirl
        # # [14, 237], [25, 190], [23, 170], [7, 190]
        # # [14, 123], [25, 170], [23, 190], [7, 170]
        
        # swirl_1 = np.linspace(np.radians(123), np.radians(170), 200)
        # swirl_2 = np.linspace(np.radians(170), np.radians(190), 200)
        # swirl_3 = np.linspace(np.radians(190), np.radians(170), 200)
        # swirl_4 = np.linspace(np.radians(170), np.radians(123), 200)
        
        # swirl_lat1 = np.linspace(14, 25, 200)
        # swirl_lat2 = np.linspace(25, 23, 200)
        # swirl_lat3 = np.linspace(23, 7, 200)
        # swirl_lat4 = np.linspace(7, 14, 200)
        
        # med swirl
        # [17, 230], [25, 190], [23, 171], [8, 208]
        # [17, 130], [25, 170], [23, 189], [8, 152]
        
        # swirl_1 = np.linspace(np.radians(90), np.radians(170), 200)
        # swirl_2 = np.linspace(np.radians(170), np.radians(189), 200)
        # swirl_3 = np.linspace(np.radians(189), np.radians(130), 200) # 90 = 152
        # swirl_4 = np.linspace(np.radians(130), np.radians(90), 200) # 90 = 152, 160
        
        # swirl_lat1 = np.linspace(8, 25, 200)
        # swirl_lat2 = np.linspace(25, 23, 200)
        # swirl_lat3 = np.linspace(23, 0, 200)
        # swirl_lat4 = np.linspace(0, 8, 200) 
        
        # # avg dusk
        # # [7, 189], [23, 170], [23, 150], [6, 134]
        # # [7, 171], [23, 190], [23, 210], [6, 226]
        
        # dusk_1 = np.linspace(np.radians(171), np.radians(190), 200)
        # dusk_2 = np.linspace(np.radians(190), np.radians(210), 200)
        # dusk_3 = np.linspace(np.radians(210), np.radians(226), 200)
        # dusk_4 = np.linspace(np.radians(226), np.radians(171), 200)
        
        # dusk_lat1 = np.linspace(7, 23, 200)
        # dusk_lat2 = np.linspace(23, 23, 200)
        # dusk_lat3 = np.linspace(23, 6, 200)
        # dusk_lat4 = np.linspace(6, 7, 200)
        
        # med dusk
        # [8, 208], [23, 171], [23, 150], [4, 140]
        # [8, 152], [23, 189], [23, 210], [4, 220]
        # [8, 152], [23, 189], [23, 210], [10, 220], [0, 90]
        
        # dusk_1 = np.linspace(np.radians(152), np.radians(189), 200)
        # dusk_2 = np.linspace(np.radians(189), np.radians(210), 200)
        # dusk_3 = np.linspace(np.radians(210), np.radians(210), 200)
        # #dusk_4 = np.linspace(np.radians(210), np.radians(210), 200)
        # dusk_5 = np.linspace(np.radians(210), np.radians(210), 200)
        # dusk_6 = np.linspace(np.radians(210), np.radians(130), 200)
        # dusk_7 = np.linspace(np.radians(130), np.radians(152), 200)
        
        # dusk_lat1 = np.linspace(8, 23, 200)
        # dusk_lat2 = np.linspace(23, 23, 200)
        # dusk_lat3 = np.linspace(23, 12, 200)
        # #dusk_lat4 = np.linspace(12, 8, 200)
        # dusk_lat5 = np.linspace(12, 0, 200)
        # dusk_lat6 = np.linspace(0, 0, 200)
        # dusk_lat7 = np.linspace(0, 8, 200)
        
        # # avg noon
        # # [25, 187], [29, 183], [32, 166], [28, 153], [23, 150]
        # # [25, 173], [29, 177], [32, 194], [28, 207], [23, 210]
        
        # noon_1 = np.linspace(np.radians(173), np.radians(177), 200)
        # noon_2 = np.linspace(np.radians(177), np.radians(194), 200)
        # noon_3 = np.linspace(np.radians(194), np.radians(207), 200)
        # noon_4 = np.linspace(np.radians(207), np.radians(210), 200)
        # noon_5 = np.linspace(np.radians(210), np.radians(173), 200)
        
        # noon_lat1 = np.linspace(25, 29, 200)
        # noon_lat2 = np.linspace(29, 32, 200)
        # noon_lat3 = np.linspace(32, 28, 200)
        # noon_lat4 = np.linspace(28, 23, 200)
        # noon_lat5 = np.linspace(23, 25, 200)
        
        # med noon
        # [25, 188], [30, 180], [32, 167], [28, 152], [23, 150]
        # [25, 172], [30, 180], [32, 193], [28, 208], [23, 210]
        
        # noon_1 = np.linspace(np.radians(172), np.radians(180), 200)
        # noon_2 = np.linspace(np.radians(180), np.radians(193), 200)
        # noon_3 = np.linspace(np.radians(193), np.radians(208), 200)
        # noon_4 = np.linspace(np.radians(208), np.radians(210), 200)
        # noon_5 = np.linspace(np.radians(210), np.radians(172), 200)
        
        # noon_lat1 = np.linspace(25, 30, 200)
        # noon_lat2 = np.linspace(30, 32, 200)
        # noon_lat3 = np.linspace(32, 28, 200)
        # noon_lat4 = np.linspace(28, 23, 200)
        # noon_lat5 = np.linspace(23, 25, 200)
        
        
        # --------- testing --------
    
        
        # boundary between active regions
        boundary = np.linspace(np.radians(205+shift),np.radians(188+shift),200)
       
        # noon active region
        updusk = np.linspace(np.radians(205+shift),np.radians(170+shift),200) #top of region line
        lon_noon_a = np.linspace(np.radians(205+shift),np.radians(190+shift),100)
        lon_noon_b = np.linspace(np.radians(190+shift),np.radians(170+shift),100)
        downlat_noon_a = np.linspace(28, 32, 100)
        downlat_noon_b = np.linspace(32, 28, 100)
        # dusk active region
        active = np.linspace(np.radians(205+shift),np.radians(160+shift),200) 
        active_line = np.linspace(10, 0, 200)
        active_2 = np.linspace(np.radians(188+shift),np.radians(100+shift),200)
        active_line2 = np.linspace(16, 4, 200)
    
        
        # dusk main emission - kidney & bean - don't use atm
        # kidney1 = np.linspace(np.radians(180+shift),np.radians(190+shift),200)
        # kidney2 = np.linspace(np.radians(190+shift),np.radians(210+shift),200)
        # kidney3 = np.linspace(np.radians(210+shift),np.radians(208+shift),200)
        # kidney4 = np.linspace(np.radians(208+shift),np.radians(205+shift),200)
        # kidney5 = np.linspace(np.radians(270+shift),np.radians(245+shift),200)
        
        # kidney6 = np.linspace(np.radians(230+shift),np.radians(225+shift),200)
        # kidney7 = np.linspace(np.radians(225+shift),np.radians(210+shift),200)
        # kidney8 = np.linspace(np.radians(210+shift),np.radians(182+shift),200)
        # kidney9 = np.linspace(np.radians(205+shift),np.radians(182+shift),200)
        # kidney10 = np.linspace(np.radians(245+shift),np.radians(230+shift),200)
        
        # bean_line1 = np.linspace(0, 4, 200)
        # bean_line2 = np.linspace(4, 10, 200)
        # bean_line3 = np.linspace(10, 20, 200)
        # bean_line4 = np.linspace(20, 30, 200)
        # bean_line5 = np.linspace(9, 10, 200)
        
        # bean_line6 = np.linspace(12,20,200)
        # bean_line7 = np.linspace(20,33,200)
        # bean_line8 = np.linspace(33,38,200)
        # bean_line9 = np.linspace(30, 34, 200)
        # bean_line10 = np.linspace(10, 12, 200)            
            
        if fixed == 'lon':
            if hemis == "North" or hemis == "north" or hemis == "N" or hemis == "n":  
                print(v_visit)
                if v_visit == '12' or v_visit == '15':
                    plt.plot(swirls,swirl_line,color='plum',linestyle='--', lw=2)
                    plt.plot(swirls2,swirl_line2,color='plum',linestyle='--', lw=2)
                    plt.plot(swirls3,swirl_line3,color='plum',linestyle='--', lw=2)
                    plt.plot(swirls4,swirl_line4,color='plum',linestyle='--', lw=2)
                    plt.plot(swirls5,swirl_line5,color='plum',linestyle='--', lw=2)
                    plt.plot(200*[np.radians(90+shift)], border1,color='plum',linestyle='--', lw=2)
                    plt.plot(200*[np.radians(270+shift)],border2,color='plum',linestyle='--', lw=2)
                    
                    plt.plot(dusks,dusk_line,color='khaki',linestyle='-', lw=1)
                    plt.plot(dusks2,dusk_line2,color='khaki',linestyle='-', lw=1)
                    plt.plot(dusks3,dusk_line3,color='khaki',linestyle='-', lw=1)
                    plt.plot(dusks4,dusk_line4,color='khaki',linestyle='-', lw=1)
                    plt.plot(dusks5,dusk_line5,color='khaki',linestyle='-', lw=1)
                    plt.plot(200*[np.radians(270+shift)],border3,color='khaki',linestyle='-',lw=1)
                    
                    plt.plot(noons,noon_line,'c-', lw=1)
                    plt.plot(noons2,noon_line2,'c-', lw=1)
                    plt.plot(noons3,noon_line3,'c-', lw=1)
                    plt.plot(noons4,noon_line4,'c-', lw=1)
                    plt.plot(noons5,noon_line5,'c-', lw=1)
                    plt.plot(noons6,noon_line6,'c-', lw=1)
                    
                    plt.legend(facecolor='black', labelcolor= 'linecolor',bbox_to_anchor=[0, 0.9], loc='upper left')
                    
                    if fill_in == True:
                        plt.scatter(lons_rad_swirl, lats_rad_swirl, color="plum", alpha=0.05, s=1, zorder=9, label="Swirl Region")
                    
                    # noon active region
                    #plt.scatter(lat_verts_rad_noon_active, lon_verts_rad_noon_active, color="b")
                        plt.scatter(lons_rad_noon_active, lats_rad_noon_active, color="c", alpha=0.1, s=1, zorder=7, label="Noon Active Region")
                    
                    # dusk active region
                    #plt.scatter(lat_verts_rad_dusk_active, lon_verts_rad_dusk_active, color="g")
                        plt.scatter(lons_rad_dusk_active, lats_rad_dusk_active, color='khaki', alpha=0.03, s=1, zorder=8, label="Dusk Active Region")
                        
                    else:
                        print("regions not filled in")
                    
                else:
                    print('not visit 12 or 15')
                    # dusk active boundary
                    plt.plot([np.radians(205+shift), np.radians(205+shift)], [24, 10], color='khaki', linestyle='-', lw=1.5, label="Dusk Active Region")
                    plt.plot([np.radians(188+shift), np.radians(188+shift)], [24, 16], color='khaki', linestyle='-',lw=1.5)
                    plt.plot([np.radians(100+shift), np.radians(100+shift)], [4, 0], color='khaki', linestyle='-',lw=1.5)
                    #plt.plot(active, 200*[10], 'r-', lw=1.5)
                    plt.plot(boundary, 200*[24], color='khaki', linestyle='-', lw=2.5)
                    plt.plot(active, active_line, color='khaki', linestyle='-', lw=1.5)
                    plt.plot(active_2, active_line2, color='khaki', linestyle='-', lw=1.5)
                    # plt.plot(active_3, active_line3, 'r-', lw=1.5)
                    
                     # #LHS / 'dawn' arc
                     # plt.plot([np.radians(100+shift), np.radians(100+shift)], [20, 11], 'k-', lw=1) #130
                     # plt.plot([np.radians(180+shift), np.radians(180+shift)], [33, 39], 'k-', lw=1) # changed from 54 to 39 to close the polygon
                     # plt.plot(lon_dawn, uplat_dawn, 'k-', lw=1)
                     # plt.plot(lon_dawn, downlat_dawn, 'k-', lw=1)
                    
                    # noon active boundary
                    plt.plot([np.radians(205+shift), np.radians(205+shift)], [24, 28], 'c-', lw=1.5, label="Noon Active Region")    
                    plt.plot([np.radians(170+shift), np.radians(170+shift)], [24, 28], 'c-', lw=1.5)
                    plt.plot(lon_noon_a, downlat_noon_a, 'c-', lw=1.5)
                    plt.plot(lon_noon_b, downlat_noon_b, 'c-', lw=1.5)
                    plt.plot(updusk, 200*[24], 'c-', lw=1.5)
                    
                    #swirl region
                    plt.plot([np.radians(188+shift), np.radians(188+shift)], [24, 16], color='plum',linestyle='--', lw=1, label="Swirl Region")
                    plt.plot([np.radians(100+shift), np.radians(100+shift)], [6, 4], color='plum',linestyle='--', lw=1)
            
                    plt.plot(upnoon, 200*[24],color='plum',linestyle='--',lw=1)
                    plt.plot(swirl,swirl_linee, color='plum',linestyle='--', lw=1)
                    plt.plot(swirl2, swirl_linee2, color='plum',linestyle='--', lw=1)
                    
                    
                    
                    # ----------- testing --------------
                     
                     
                   #  plt.plot(swirl_1, swirl_lat1,color='plum',linestyle='--',lw=1)
                   #  plt.plot(swirl_2, swirl_lat2,color='plum',linestyle='--',lw=1)
                   #  plt.plot(swirl_3, swirl_lat3,color='plum',linestyle='--',lw=1)
                   #  plt.plot(swirl_4, swirl_lat4,color='plum',linestyle='--',lw=1)
                    
                   #  plt.plot(dusk_1, dusk_lat1,color='khaki',linestyle='-',lw=1)
                   #  plt.plot(dusk_2, dusk_lat2,color='khaki',linestyle='-',lw=1)
                   #  plt.plot(dusk_3, dusk_lat3,color='khaki',linestyle='-',lw=1)
                   # # plt.plot(dusk_4, dusk_lat4,color='khaki',linestyle='-',lw=1)
                   #  plt.plot(dusk_5, dusk_lat5,color='khaki',linestyle='-',lw=1)
                   #  plt.plot(dusk_6, dusk_lat6,color='khaki',linestyle='-',lw=1)
                   #  plt.plot(dusk_7, dusk_lat7,color='khaki',linestyle='-',lw=1)
                    
                   #  plt.plot(noon_1, noon_lat1,color='cyan',linestyle='-',lw=1)
                   #  plt.plot(noon_2, noon_lat2,color='cyan',linestyle='-',lw=1)
                   #  plt.plot(noon_3, noon_lat3,color='cyan',linestyle='-',lw=1)
                   #  plt.plot(noon_4, noon_lat4,color='cyan',linestyle='-',lw=1)
                   #  plt.plot(noon_5, noon_lat5,color='cyan',linestyle='-',lw=1)
                
    
                    plt.legend(facecolor='black', labelcolor= 'linecolor',bbox_to_anchor=[0, 0.9], loc='upper left')
                    
                    if fill_in == True:
                        plt.scatter(lons_rad_swirl, lats_rad_swirl, color="plum", alpha=0.05, s=1, zorder=9, label="Swirl Region")
                    
                    # noon active region
                    #plt.scatter(lat_verts_rad_noon_active, lon_verts_rad_noon_active, color="b")
                        plt.scatter(lons_rad_noon_active, lats_rad_noon_active, color="c", alpha=0.1, s=1, zorder=7, label="Noon Active Region")
                    
                    # dusk active region
                    #plt.scatter(lat_verts_rad_dusk_active, lon_verts_rad_dusk_active, color="g")
                        plt.scatter(lons_rad_dusk_active, lats_rad_dusk_active, color='khaki', alpha=0.03, s=1, zorder=8, label="Dusk Active Region")
                        
                    else:
                        print("regions not filled in")
                    
                    
                    # #RHS / 'dusk' ME 
                    # plt.plot([np.radians(270+shift),np.radians(270+shift)],[0,9],'m-',lw=1)
                    # plt.plot([np.radians(182+shift), np.radians(182+shift)], [34, 38], 'm-', lw=1)
    
                    # plt.plot(kidney1, bean_line1, 'm-', lw=1)
                    # plt.plot(kidney2, bean_line2, 'm-', lw=1)
                    # plt.plot(kidney3, bean_line3, 'm-', lw=1)
                    # plt.plot(kidney4, bean_line4, 'm-', lw=1)
                    # plt.plot(kidney5, bean_line5, 'm-', lw=1)
                    # plt.plot(kidney6, bean_line6, 'm-', lw=1)
                    # plt.plot(kidney7, bean_line7, 'm-', lw=1)
                    # plt.plot(kidney8, bean_line8, 'm-',lw=1)
                    # plt.plot(kidney9, bean_line9, 'm-',lw=1)
                    # plt.plot(kidney10, bean_line10, 'm-',lw=1)
                
    
                    # colour in sections - only works for polar region so far
                    #plt.scatter(lon_rad, lat_flat, color="red", s=1, label="All Intensities")
    
                    # swirl region
                    #plt.scatter(lat_verts_rad_swirl, lon_verts_rad_swirl, color="r")
            
                
            # elif hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
            # !!! not defined yet
        
        # regions for the LT fixed case - would this even work?
        elif fixed == 'lt':
            #dusk active boundary
            plt.plot([np.radians(205+shift), np.radians(205+shift)], [24, 10], color='khaki', linestyle='-', lw=1.5, label="Dusk Active Region")
            plt.plot([np.radians(188+shift), np.radians(188+shift)], [24, 16], color='khaki', linestyle='-',lw=1.5)
            plt.plot([np.radians(100+shift), np.radians(100+shift)], [3, 0], color='khaki', linestyle='-',lw=1.5)
            #plt.plot(active, 200*[10], 'r-', lw=1.5)
            plt.plot(boundary, 200*[24], color='khaki', linestyle='-', lw=2.5)
            plt.plot(active, active_line, color='khaki', linestyle='-', lw=1.5)
            plt.plot(active_2, active_line2, color='khaki', linestyle='-', lw=1.5)
            # plt.plot(active_3, active_line3, 'r-', lw=1.5)
            
            # #LHS / 'dawn' arc
            # plt.plot([np.radians(100+shift), np.radians(100+shift)], [20, 11], 'k-', lw=1) #130
            # plt.plot([np.radians(180+shift), np.radians(180+shift)], [33, 39], 'k-', lw=1) # changed from 54 to 39 to close the polygon
            # plt.plot(lon_dawn, uplat_dawn, 'k-', lw=1)
            # plt.plot(lon_dawn, downlat_dawn, 'k-', lw=1)
            
            # noon active boundary
            plt.plot([np.radians(205+shift), np.radians(205+shift)], [24, 28], 'c-', lw=1.5, label="Noon Active Region")    
            plt.plot([np.radians(170+shift), np.radians(170+shift)], [24, 28], 'c-', lw=1.5)
            plt.plot(lon_noon_a, downlat_noon_a, 'c-', lw=1.5)
            plt.plot(lon_noon_b, downlat_noon_b, 'c-', lw=1.5)
            plt.plot(updusk, 200*[24], 'c-', lw=1.5)
            
            #swirl region
            plt.plot([np.radians(185+shift), np.radians(185+shift)], [22, 16], color='plum',linestyle='--', lw=1, label="Swirl Region")
            plt.plot([np.radians(100+shift), np.radians(100+shift)], [6, 4], color='plum',linestyle='--', lw=1)
    
            plt.plot(upnoon, 200*[22],color='plum',linestyle='--',lw=1)
            plt.plot(swirl,swirl_line, color='plum',linestyle='--', lw=1)
            plt.plot(swirl2, swirl_line2, color='plum',linestyle='--', lw=1)
            
            plt.legend(facecolor='black', labelcolor= 'linecolor',bbox_to_anchor=[0, 0.8], loc='upper left')
            
            # #RHS / 'dusk' ME 
            # plt.plot([np.radians(270+shift),np.radians(270+shift)],[0,9],'m-',lw=1)
            # plt.plot([np.radians(182+shift), np.radians(182+shift)], [34, 38], 'm-', lw=1)

            # plt.plot(kidney1, bean_line1, 'm-', lw=1)
            # plt.plot(kidney2, bean_line2, 'm-', lw=1)
            # plt.plot(kidney3, bean_line3, 'm-', lw=1)
            # plt.plot(kidney4, bean_line4, 'm-', lw=1)
            # plt.plot(kidney5, bean_line5, 'm-', lw=1)
            # plt.plot(kidney6, bean_line6, 'm-', lw=1)
            # plt.plot(kidney7, bean_line7, 'm-', lw=1)
            # plt.plot(kidney8, bean_line8, 'm-',lw=1)
            # plt.plot(kidney9, bean_line9, 'm-',lw=1)
            # plt.plot(kidney10, bean_line10, 'm-',lw=1)
            
            #plt.scatter(lon_rad, lat_flat, color="red", s=1, label="All Intensities")
            # swirl region
            #plt.scatter(lat_verts_rad_swirl, lon_verts_rad_swirl, color="r")
            #breakpoint()
            plt.scatter(lons_rad_swirl, lats_rad_swirl, color="plum", alpha=0.05, s=1, zorder=9, label="Swirl Region")
            
            # noon active region
            #plt.scatter(lat_verts_rad_noon_active, lon_verts_rad_noon_active, color="b")
            plt.scatter(lons_rad_noon_active, lats_rad_noon_active, color="c", alpha=0.1, s=1, zorder=7, label="Noon Active Region")
            
            # dusk active region
            #plt.scatter(lat_verts_rad_dusk_active, lon_verts_rad_dusk_active, color="g")
            plt.scatter(lons_rad_dusk_active, lats_rad_dusk_active, color='khaki', alpha=0.03, s=1, zorder=8, label="Dusk Active Region")
            # # delimitation (region)
            # updusk = np.linspace(np.radians(205+shift),np.radians(170+shift),200)
            # lon_dawn = np.linspace(np.radians(180+shift),np.radians(130+shift),200)
            # uplat_dawn = np.linspace(33, 15, 200)
            # downlat_dawn = np.linspace(39, 23, 200)
            
            # lon_noon_a = np.linspace(np.radians(205+shift),np.radians(190+shift),100)
            # lon_noon_b = np.linspace(np.radians(190+shift),np.radians(170+shift),100)
            # downlat_noon_a = np.linspace(28, 32, 100)
            # downlat_noon_b = np.linspace(32, 27, 100)
            
            # if hemis == "North" or hemis == "north" or hemis == "N" or hemis == "n":
            #     # dusk boundary
            #     plt.plot([np.radians(205+shift), np.radians(205+shift)], [20, 10], 'r-', lw=1.5)
            #     plt.plot([np.radians(170+shift), np.radians(170+shift)], [10, 20], 'r-', lw=1.5)
            #     plt.plot(updusk, 200*[10], 'r-', lw=1.5)
            #     plt.plot(updusk, 200*[20], 'r-', lw=1.5)
            #     # dawn boundary
            #     plt.plot([np.radians(130+shift), np.radians(130+shift)], [23, 15], 'b-', lw=1)
            #     plt.plot([np.radians(180+shift), np.radians(180+shift)], [33, 39], 'b-', lw=1) # changed from 54 to 39 to close the polygon
            #     plt.plot(lon_dawn, uplat_dawn, 'b-', lw=1)
            #     plt.plot(lon_dawn, downlat_dawn, 'b-', lw=1)
                
            
            #     # noon boundary
            #     plt.plot([np.radians(205+shift), np.radians(205+shift)], [22, 28], 'y-', lw=1.5)    
            #     plt.plot([np.radians(170+shift), np.radians(170+shift)], [27, 22], 'y-', lw=1.5)
            #     plt.plot(lon_noon_a, downlat_noon_a, 'y-', lw=1.5)
            #     plt.plot(lon_noon_b, downlat_noon_b, 'y-', lw=1.5)
            #     plt.plot(updusk, 200*[22], 'y-', lw=1.5)
            #     # polar boundary
            #     plt.plot([np.radians(205+shift), np.radians(205+shift)], [10, 28], 'w--', lw=1)
            #     plt.plot(lon_noon_a, downlat_noon_a, 'w--', lw=1)
            #     plt.plot(lon_noon_b, downlat_noon_b, 'w--', lw=1)
            #     plt.plot([np.radians(170+shift), np.radians(170+shift)], [27, 10], 'w--', lw=1)
            #     plt.plot(updusk, 200*[10], 'w--', lw=1)
                
            
            # elif hemis == "South" or hemis == "south" or hemis == "S" or hemis == "s":
            # !!! not defined yet
            
            
#################################################################################
                # MUST SHIFT THE MOONFP FOR THE LT FIXED CASE #
################################################################################

    # drawing the moon footprints 
    if moonfp == True:
        #retrieve their expected longitude and latitude (from Hess et al., 2011)
        nlonio, ncolatio, slonio, scolatio, nloneu, ncolateu, sloneu, scolateu, nlonga, ncolatga, slonga, scolatga = moonfploc(iolon,eulon,galon)
        nlonio1, ncolatio1, slonio1, scolatio1, nloneu1, ncolateu1, sloneu1, scolateu1, nlonga1, ncolatga1, slonga1, scolatga1 = moonfploc(iolon1,eulon1,galon1)
        nlonio2, ncolatio2, slonio2, scolatio2, nloneu2, ncolateu2, sloneu2, scolateu2, nlonga2, ncolatga2, slonga2, scolatga2 = moonfploc(iolon2,eulon2,galon2)
       
        # plot a colored mark in their expected location, together with their name
        if fixed == 'lon':
            if hemis == "North" or hemis == "north" or hemis == "N" or hemis == "n":
                # we define some intervals for plotting the moon footprints because if they
                # are supposed to be way inside the "night" hemisphere (only within +-120degrees
                # from CML), if not, we do not plot them
                if abs(cml-nlonio1) < 120 or abs(cml-nlonio1) > 240:
                    plt.plot([2*np.pi-(np.radians(nlonio1)),2*np.pi-(np.radians(nlonio2))],[ncolatio, ncolatio], 'k-', lw=4)
                    plt.plot([2*np.pi-(np.radians(nlonio1)),2*np.pi-(np.radians(nlonio2))],[ncolatio, ncolatio], color='gold', linestyle='-', lw=2.5)
                    plt.text(2*np.pi-(np.radians(nlonio)), 3.5+ncolatio, 'IO', color='gold', fontsize=10,  fontweight='bold',alpha=0.5,\
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')], horizontalalignment='center', verticalalignment='center')
                
                if abs(cml-nloneu1) < 120 or abs(cml-nloneu1) > 240:
                    plt.plot([2*np.pi-(np.radians(nloneu1)),2*np.pi-(np.radians(nloneu2))],[ncolateu, ncolateu], 'k-', lw=4)
                    plt.plot([2*np.pi-(np.radians(nloneu1)),2*np.pi-(np.radians(nloneu2))],[ncolateu, ncolateu], color='aquamarine', linestyle='-', lw=2.5)
                    plt.text(2*np.pi-(np.radians(nloneu)), 3.5+ncolateu, 'EUR', color='aquamarine', fontsize=10, fontweight='bold',alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center')
                
                if abs(cml-nlonga1) < 120 or abs(cml-nlonga1) > 250:
                    plt.plot([2*np.pi-(np.radians(nlonga1)),2*np.pi-(np.radians(nlonga2))],[ncolatga, ncolatga], 'k-', lw=4)
                    plt.plot([2*np.pi-(np.radians(nlonga1)),2*np.pi-(np.radians(nlonga2))],[ncolatga, ncolatga], 'w-', lw=2.5)
                    plt.text(2*np.pi-(np.radians(nlonga)), 3.5+ncolatga, 'GAN', color='w', fontsize=10, fontweight='bold',alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center')
            else: #if we are in the Southern hemisphere
                if abs(cml-slonio1) < 120 or abs(cml-slonio1) > 240:
                    plt.plot([(np.radians(180-slonio1)),(np.radians(180-slonio2))],[scolatio, scolatio], 'k-', lw=4)
                    plt.plot([(np.radians(180-slonio1)),(np.radians(180-slonio2))],[scolatio, scolatio], color='gold', linestyle='-', lw=2.5)
                    plt.text((np.radians(180-slonio)), 3.5+scolatio, 'IO', color='gold', fontsize=10, alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
                
                if abs(cml-sloneu1) < 120 or abs(cml-sloneu1) > 240:
                    plt.plot([(np.radians(180-sloneu1)),(np.radians(180-sloneu2))],[scolateu, scolateu], 'k-', lw=4)
                    plt.plot([(np.radians(180-sloneu1)),(np.radians(180-sloneu2))],[scolateu, scolateu], color='aquamarine', linestyle='-', lw=2.5)
                    plt.text((np.radians(180-sloneu)), 3.5+scolateu, 'EUR', color='aquamarine', fontsize=10, alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
                
                if abs(cml-slonga1) < 120 or abs(cml-slonga1) > 240:
                    plt.plot([(np.radians(180-slonga1)),(np.radians(180-slonga2))],[scolatga, scolatga], 'k-', lw=4)
                    plt.plot([(np.radians(180-slonga1)),(np.radians(180-slonga2))],[scolatga, scolatga], 'w-', lw=2.5)
                    plt.text((np.radians(180-slonga)), 3.5+scolatga, 'GAN', color='w', fontsize=10, alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
    
        elif fixed == 'lt':
            if hemis == "North" or hemis == "north" or hemis == "N" or hemis == "n":
                if abs(cml-nlonio1) < 120 or abs(cml-nlonio1) > 240:
                    plt.plot([(np.radians(180+cml-nlonio1)),(np.radians(180+cml-nlonio2))],[ncolatio, ncolatio], 'k-', lw=4)
                    plt.plot([(np.radians(180+cml-nlonio1)),(np.radians(180+cml-nlonio2))],[ncolatio, ncolatio], color='gold', linestyle='-', lw=2.5)
                    plt.text((np.radians(180+cml-nlonio)), 3.5+ncolatio, 'IO ', color='gold', fontsize=10,  fontweight='bold',alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center')
                
                if abs(cml-nloneu1) < 120 or abs(cml-nloneu1) > 240:
                    plt.plot([2*np.pi-(np.radians(180+cml-nloneu1)),2*np.pi-(np.radians(180+cml-nloneu2))],[ncolateu, ncolateu], 'k-', lw=4)
                    plt.plot([2*np.pi-(np.radians(180+cml-nloneu1)),2*np.pi-(np.radians(180+cml-nloneu2))],[ncolateu, ncolateu], color='aquamarine', linestyle='-', lw=2.5)
                    plt.text((np.radians(180+cml-nloneu)), 3.5+ncolateu, 'EUR', color='aquamarine', fontsize=10, fontweight='bold',alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center')
                
                if abs(cml-nlonga1) < 120 or abs(cml-nlonga1) > 240:
                    plt.plot([(np.radians(180+cml-nlonga1)),(np.radians(180+cml-nlonga2))],[ncolatga, ncolatga], 'k-', lw=4)
                    plt.plot([(np.radians(180+cml-nlonga1)),(np.radians(180+cml-nlonga2))],[ncolatga, ncolatga], 'w-', lw=2.5)
                    plt.text((np.radians(180+cml-nlonga)), 3.5+ncolatga, 'GAN', color='w', fontsize=10, fontweight='bold',alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center')
       
            else: # south hemisphere
                if abs(cml-slonio1) < 120 or abs(cml-slonio1) > 240:
                    plt.plot([(np.radians(180+cml-slonio1)),(np.radians(180+cml-slonio2))],[scolatio, scolatio], 'k-', lw=4)
                    plt.plot([(np.radians(180+cml-slonio1)),(np.radians(180+cml-slonio2))],[scolatio, scolatio], color='gold', linestyle='-', lw=2.5)
                    plt.text((np.radians(180+cml-slonio)), 3.5+scolatio, 'IO ', color='gold', fontsize=10, alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
                
                if abs(cml-sloneu1) < 120 or abs(cml-sloneu1) > 240:
                    plt.plot([(np.radians(180+cml-sloneu1)),(np.radians(180+cml-sloneu2))],[scolateu, scolateu], 'k-', lw=4)
                    plt.plot([(np.radians(180+cml-sloneu1)),(np.radians(180+cml-sloneu2))],[scolateu, scolateu], color='aquamarine', linestyle='-', lw=2.5)
                    plt.text((np.radians(180+cml-sloneu)), 3.5+scolateu, 'EUR', color='aquamarine', fontsize=10,alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
                
                if abs(cml-slonga1) < 120 or abs(cml-slonga1) > 240:
                    plt.plot([(np.radians(180+cml-slonga1)),(np.radians(180+cml-slonga2))],[scolatga, scolatga], 'k-', lw=4)
                    plt.plot([(np.radians(180+cml-slonga1)),(np.radians(180+cml-slonga2))],[scolatga, scolatga], 'w-', lw=2.5)
                    plt.text((np.radians(180+cml-slonga)), 3.5+scolatga, 'GAN', color='w', fontsize=10,alpha=0.5,
                             path_effects=[patheffects.withStroke(linewidth=1, foreground='black')],\
                             horizontalalignment='center', verticalalignment='center', fontweight='bold')
                        

    # defining the final filename, adding sufixes depending n the stuff we are showing
    namesave = str(filename) #filename[5:] for jonny 2016 01 as for some reason files don't follow right convention, neither does 02

    if fixed == 'lon': 
        #print(namesave)
        if whole_region_box == True:
            if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/whole_emission_boxed/'+prefix+''+v+'/'):
                os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/whole_emission_boxed/'+prefix+''+v+'/')
                
            saveloc= (f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/whole_emission_boxed/'+prefix+''+v+'/mo_'+namesave+'_fixlon.jpg')
        
        else:
         if regions == True:
            if output == 'kr':
                 if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_intensity/'+prefix+''+v+'/'):
                     os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_intensity/'+prefix+''+v+'/')
                 saveloc= (f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_intensity/'+prefix+''+v+'/mo_'+namesave+'_fixlon.jpg')
            
            elif output == 'w':
                if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_power/'+prefix+''+v+'/'):
                        os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_power/'+prefix+''+v+'/')
                saveloc= (f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/boxed_power/'+prefix+''+v+'/mo_'+namesave+'_fixlon.jpg')
        
        if regions == False:
            if output == 'kr':
                if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/'+prefix+''+v+'/'):
                    os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/'+prefix+''+v+'/')
                saveloc= (f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/'+prefix+''+v+'/mo_'+namesave+'_fixlon.jpg')
            
            elif output == 'w':
                if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/nobox_power/'+prefix+''+v+'/'):
                    os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/nobox_power/'+prefix+''+v+'/')
                saveloc= (f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/nobox_power/'+prefix+''+v+'/mo_'+namesave+'_fixlon.jpg')
        # if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{prefix}{v}/fin/{tint}s/'):
        #     os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{prefix}{v}/fin/{tint}s/')
            
        print('Name of the saved image is mo_'+str(namesave)+"_fixlon.jpg")
        plt.savefig(saveloc, dpi=dpi) #/{prefix}{v}/fin/{tint}s/mo_{namesave}_fixlon.jpg' 
        print(saveloc)
        
        
    elif fixed == 'lt':
        # first for need f before '/Users/ if have {prefix} in2018
        if not os.path.exists(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/'+prefix+''+v+'/'): #{prefix}{v}/fin/
            os.makedirs(f'/Users/hannah/OneDrive - Lancaster University/aurora/{fixed}/'+prefix+''+v+'/') #/{prefix}{v}/fin/
        # if not os.path.exists('/Users/hannah/OneDrive - Lancaster University/aurora/'):#/{prefix}{v}/fin/{tint}s_fixedlt/
        #     os.makedirs('/Users/hannah/OneDrive - Lancaster University/aurora/')#/{prefix}{v}/fin/{tint}s_fixedlt/
       
        print('Name of the saved image is mo_'+str(filename)+"_fixedlt.jpg")
        saveloc = (f'/Users/hannah/OneDrive - Lancaster University//aurora/{fixed}/'+prefix+''+v+'/mo_'+namesave+'_fixedlt.jpg')
        plt.savefig(saveloc, dpi=300)# save location/{prefix}{v}/fin/{tint}s_fixedlt/ goes between aurora and movisit
        print(saveloc)

    plt.close()
    
# and this chunk is to call the function:
def multigif(lista, year, prefix, extra, time, radius, moonfp, full, fixed, output, regions, whole_region_box, polarf):#, mf, indf, polarf, secondf
    for l in lista:
        l = str(l)
        print(f'VISIT {l} \n \n') #just to check out the visit we are plotting
        
        # we grab all the files we are interested in plotting
        arch = '*_v'+ l
        ti = str('/*0'+time+'*')
        loc = (f'/Users/hannah/OneDrive - Lancaster University/aurora/data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti)
        ab = glob.glob(loc)#('/Users/hannah/OneDrive - Lancaster University/aurora/data/'+arch+'/nopolar'+time+ti)#(f'/Users/hannah/OneDrive - Lancaster University/aurora/data/{year}/extract/{extra}'+arch+'/nopolar'+time+ti)#('/Users/hannah/OneDrive - Lancaster University/aurora/data/'+arch+'/nopolar'+time+ti)#+ti)#/nopolar'+time+ti) #
        print(loc)
        #breakpoint()
        ab.sort()   
        

        # print(f"Length of ab is: {len(ab)}")
        
        # define the image resolution for when there are many images so the 
        # final GIF is not too large
        if int(time) < 100:
            cal = 150
        else:
            cal = 300
            
        # and now we loop to apply the plotting function (moind) to every image
        # in the visit
        for n,i in tqdm(enumerate(ab)):
            hdulist = fits.open(i)
            header = hdulist[1].header
            image = hdulist[1].data
            
            #raise ValueError
            
            try:
                hemis = header['HEMISPH']
            except NameError:
                hemis = str(input('Input hemisphere manually:  ("north" or "south")  '))
            filename = str(i)[-51:-5]
            visit = filename[-20:-18]
            moind(image, visit, n, header, filename, prefix, dpi = cal, crop=1, rlim = radius, fixed = fixed, hemis = hemis, full=full, moonfp=moonfp, output=output, regions=regions, fill_in=fill_in, whole_region_box=whole_region_box, photo=n)#, mf=mf, indf=indf, polarf=polarf, secondf=secondf
            hdulist.close()
            print(f'Image {n} of {str(i)[-51:-5]} created.')
        
        # and now we start creating the GIF
        visita = prefix+arch[-2:]
        imagesgif = []
        
        # durat = input('GIF duration (in seconds): ') #deprecated
        fps = 5
        if int(time) < 50:
            fps = 40
        
        # name of the GIF
        gifname = visita + '_t'+time
        
        print(f'On to the GIF ({gifname})')
        
        nam = str('*0'+time+'*')

        # optional suffixes to the GIF's name
        # + grabbing the correct images recently generated by moind
        if polarf == True:
            if fixed == 'lon':
                gifname += '_flon'
                if whole_region_box == True:
                    gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/whole_emission_boxed/{visita}/*.jpg')
                if whole_region_box == False:
                    if regions == True:
                        if output == 'w':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/boxed_power/{visita}/*.jpg')
                        if output == 'kr':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/boxed_intensity/{visita}/*.jpg')
                    if regions == False:
                        if output == 'w':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/nobox_power/{visita}/*.jpg')
                        if output == 'kr':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/{visita}/*.jpg')
            elif fixed == 'lt':
                gifname += '_flt'
                gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lt/{visita}/*.jpg') # C:/Users/moralpom/phd/Pictures/polar/ #+str(time)+'s_fixedlt/mopf'+nam #

        else:
            if fixed == 'lon':
                gifname += '_flon'
                if whole_region_box == True:
                    gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/whole_emission_boxed/{visita}/*.jpg')
                if whole_region_box == False:
                    if regions == True:
                        if output == 'w':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/boxed_power/{visita}/*.jpg')
                        if output == 'kr':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/boxed_intensity/{visita}/*.jpg')
                    if regions == False:
                        if output == 'w':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/nobox_power/{visita}/*.jpg')
                        if output == 'kr':
                            gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lon/{visita}/*.jpg')
            elif fixed == 'lt':
                gifname += '_flt'
                gifphotos = glob.glob(f'/Users/hannah/OneDrive - Lancaster University/aurora/lt/{visita}/*.jpg')#+visita+'/fin') # C:/Users/moralpom/phd/Pictures/polar/ #+str(time)+'s_fixedlt/'+nam
        
        if moonfp == True:
            gifname += '_mfp'
        
        if full == False:
            gifname += '_half'
        
        # final suffix (termination)
        gifname += '.gif'
        
        # GIF creation:
        gifphotos.sort()
        for file in gifphotos:
            imagesgif.append(imageio.imread(file))

        # saving the GIF - this saves to the right place
        if whole_region_box == True:
            imageio.mimsave('/Users/hannah/OneDrive - Lancaster University/aurora/gifs/whole_emission_boxed/' + gifname, imagesgif, fps=fps)
        if regions == True:
            if output == 'w':
                imageio.mimsave('/Users/hannah/OneDrive - Lancaster University/aurora/gifs/boxed_power/' + gifname, imagesgif, fps=fps)
            if output == 'kr':
                imageio.mimsave('/Users/hannah/OneDrive - Lancaster University/aurora/gifs/boxed_intensity/' + gifname, imagesgif, fps=fps)
        if regions == False:
            if output == 'w':
                imageio.mimsave('/Users/hannah/OneDrive - Lancaster University/aurora/gifs/nobox_power/' + gifname, imagesgif, fps=fps)
            if output == 'kr':
                imageio.mimsave('/Users/hannah/OneDrive - Lancaster University/aurora/gifs/' + gifname, imagesgif, fps=fps) # C:/Users/moralpom/phd/Pictures/gifs/Case studies/fin/

# and this last part is the one that must be run every time:
#os.chdir('/Users/hannah/OneDrive - Lancaster University/aurora/') # C:/Users/moralpom/phd/pypeline/

# you have to input the year of the visits you are plotting (so cannot mix visits
# from different years in the same "run" of the code, for filepathing reasons)
year = input("Year of the visit:  \n")
if year == '2016':
    pre = input('Campaign from Jonny or Denis? (1/2)  ')
    if pre == '1' or pre == 'jonny' or pre == 'j' or pre == 'J' or pre == 'Jonny':
        prefix = 'ocx8'
        extra = 'nichols/' #this may be have to be removed if your directory system did not differentiate between Jonny's and Dennis's campaigns
    else:
        prefix = 'od8k'
        extra = 'grodent/' #this may be have to be removed if your directory system did not differentiate between Jonny's and Dennis's campaigns
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
    
time = str(input('Exposure time (in seconds: 10, 30, 100...): \n')) #usually 100
radius = int(input('Max. radius (in degrees of colatitude): \n'))   #usually 40
moonfp = not bool(input('Moon footprints printed? (Default: Enter for YES)\n')) #usually yes
full = not bool(input('Show the whole hemisphere? (Default: Enter for YES)\n')) #usually yes
fixed = str(input('Fix longitude (lon) or Local Time (lt):\n')) 
output = str(input('Output Values in Intensity (kr) or Power (w):\n')) 
regions = not bool(input('Show Regions (Default: Enter for YES)\n')) 
fill_in = not bool(input('Fill in Regions (Default: Enter for YES)\n'))
whole_region_box = not bool(input('Show Whole Auroral Region Boxed (Default: Enter for YES)\n')) 

# the only part you have to add manually is the particular set of visit numbers you want
# to plot. If you want only one that is perfectly fine but it must be IN a list
lista = ['26','27','28','24','25','34','35']#, '05', '17']#,'10','11','12','13','15','16','18','19','20','21'] #for example, or lista = ['0v'] 'v04','v10','v11'

multigif(lista, year, prefix, extra, time, radius, moonfp, full, fixed, output, regions, whole_region_box, True) # this is what I need to call
