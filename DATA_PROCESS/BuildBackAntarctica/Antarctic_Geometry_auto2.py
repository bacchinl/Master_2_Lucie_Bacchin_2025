#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 10:04:19 2023

@author: cmosbeux

Description: This script reconstructs the ice sheet as it was back in time based on current geometry and dhdt trend. The grounding line position and flotation criterion 
             are based on hydrostatic equilibrium hypothesis. 
             The reconstrution for a given year is based on the hypothesis that bedmachine geometry is representative  of the year 2020.

Requierment: 
        - topography file (containing ice thickness and bedrock geometry): e.g.: BedMachineAntarctica-v03.nc
        - dhdt files for the shelves and the ice sheet, I uses dhdt observations from IceSat1-2 missions over 2003-2019

Potential Requierment:
        - Path to data are hard coded for my system, you should define the filename in the different functions
        - if you want to make some plots, you will need my Antarctica_Background Module and add it to your paths (available on my github)
        - I handle netcdf with a format_reading Module of my own too (also available on my github)
        - You might want to adjust water and ice density depending on the densities you will use in your modelling (ri and rw).
"""


from Antarctica_Background import Plot_Antarctica, basin, scale,plot_GL
from format_reading import netcdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import cm
import scipy
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
import h5py
from FunctionsDataPlot import plot_dhdt_obs
import copy
from scipy.ndimage import gaussian_filter, zoom
import sys
from scipy import interpolate
import os

#-------------------------------------------------------------------------
# USER KEYWORDS
#-------------------------------------------------------------------------
plot = True
ri = 0.917
rw = 1.028
z_sl = 0
years = [2017 ]
resolution = 10 #km
power_rate = 1 # for long period, we can use dhdt rates that decrease with t
## Défine x_fin and y_fin for cropping
#x_fin = np.array(list(range(-1694617 , -1506616, 250)))
#y_fin = np.array(list(range( -371235,-204734, 250)))

#-------------------------------------------------------------------------

#%%Read dHdT from IceSat1-2 (over the shelves)
def Build_IceSat2_ContinuousMap(x_final=None, y_final=None, masked=False):
    print("Building IceSat dhdt...")
    #Ice Sheet
    filename  = '../../../DATA/AIS_mass_change.h5'
    
    f = h5py.File(filename, "r")
    
    a_group_key = list(f.keys())[0]
    # Get the data
    x_ground = f['x'][()].T
    y_ground = f['y'][()].T[::-1,:]
    
    data_ground = f['dHdt'][()].T[::-1,:]
    
    #Ice Shelf
    filename = '../../../DATA/ICE1_ICE2_AnIS_dHdt_2003_2018_R209_05KM_FLOAT_MASS_F2.h5'
    f = h5py.File(filename, "r")
    a_group_key = list(f.keys())[0]
    # Get the data
    x_shelf = f['x'][()]
    y_shelf = f['y'][()]
    data_shelf = f['dhdt'][()]
    data_shelf = data_shelf*9.26
    
    #Interpolation of the shelf data on the ground data
    print("\tInterpolating over same grid...")
    lim = 1
    x_share, y_share = x_ground[lim:-lim,lim:-lim], y_ground[lim:-lim,lim:-lim] 
    # print('xground = %d' % x_share.min())
    points_share = list(zip(x_share.flatten(), y_share.flatten()))
    
    #make the interpolator
    interpolator = RegularGridInterpolator((x_shelf[0,:], y_shelf[::-1,0]), data_shelf[::-1,:].T)
    #interpolation
    data_shelf_share = interpolator(points_share)
    data_shelf = data_shelf_share.reshape(x_share.shape)
   

    #mask 
    print("\tInterpolating over same grid...")
    lim = 1
    x_share, y_share = x_ground[lim:-lim,lim:-lim], y_ground[lim:-lim,lim:-lim] 
    # print('xground = %d' % x_share.min())
    points_share = list(zip(x_share.flatten(), y_share.flatten()))
    
    #make the interpolator
    interpolator = RegularGridInterpolator((x_shelf[0,:], y_shelf[::-1,0]), data_shelf[::-1,:].T)
    #interpolation
    data_shelf_share = interpolator(points_share)
    data_shelf = data_shelf_share.reshape(x_share.shape)
    
    #mask 
    mask = np.logical_and(np.isnan(data_shelf), np.isnan(data_ground[lim:-lim,lim:-lim]))
    n, m = mask.shape
    
    mask2 = copy.copy(mask)
    for i in np.arange(n):
        for j in np.arange(m):
            if mask[i,j]:
                if False in mask[i-1:i+2,j-1:j+2]:
                    mask2[i,j] = False
                    
    data = np.nan_to_num(data_shelf) + np.nan_to_num(data_ground[lim:-lim,lim:-lim])
    data[data==0.0]=np.nan
    
    #mask NaNs
    data = np.ma.masked_invalid(data)
    
    #get only the valid values
    x1 = x_share[~data.mask]
    y1 = y_share[~data.mask]
    new_data = data[~data.mask]
    
    new_data = interpolate.griddata((x1, y1), new_data.ravel(), (x_share, y_share),method='linear')
    data = np.ma.MaskedArray(new_data, mask2)
    
    #reshape to a given output grid
    if (x_final is not None) and (y_final is not None):
        print("\tRegridding")
        print("\t - origin boundaries : x = [%d, %d], y = [%d, %d]" % (x_share.min(), x_share.max(), y_share.min(), y_share.max()))
        print("\t - target boundaries : x = [%d, %d], y = [%d, %d]" % (x_final.min(), x_final.max(), y_final.min(), y_final.max()))
        
        y_offset = 5e4 #I have to apply offset but I don't know why
        
        interpolator2 = RegularGridInterpolator((x_share[0,:], y_share[::-1,0] - y_offset), data.T, bounds_error=False)
        xx,yy = np.meshgrid(x_final,y_final)
        points_final = list(zip(xx.flatten(), yy.flatten()))
        
        data = interpolator2(points_final)
        data = data.reshape((len(x_final),len(y_final)))
        x, y = x_final, y_final
        
    else: 
        x, y = x_share, y_share
    print("\tDone.")
    return x, y, data

def import_mask(path_file, xres,yres):
    xres=xres*1000
    yres = yres*1000
    os.system('gdalwarp -tr '+str(xres)+' '+str(yres)+' -r average '+path_file+' '+path_file.split('.tif')[0]+'_reproj.tif') 
    dataset = gdal.Open(path_file.split('.tif')[0]+'_reproj.tif')
    #dataset = gdal.Open(path_file)

    if dataset is None:
        raise FileNotFoundError(f"File {path_file} cannot be found")
    
    # Récupère la première bande (en supposant que le masque soit en une seule bande)
    bande = dataset.GetRasterBand(1)
    
    # Lit les données de la bande sous forme de tableau NumPy
    mask_data = bande.ReadAsArray()
    
    # Récupère des informations supplémentaires (dimensions, géoréférencement, projection)
    largeur = dataset.RasterXSize
    hauteur = dataset.RasterYSize
    geotransform = dataset.GetGeoTransform()
    projection = dataset.GetProjection()
    
    # Ferme le dataset (libération des ressources)
    dataset = None
    
    return {
        'data': mask_data,
        'width': largeur,
        'height': hauteur,
        'geotransform': geotransform,
        'projection': projection
    }

def read_bedmachine(resolution=5):
    print('Reading bedmachine...')
    print('\t -using a %0.1f km resolution' % (resolution*0.5))
    bedmachine_filename = '../../../DATA/BEDROCK/BedMachineAntarctica-v03.nc'
    
    x, y, bed = netcdf.readgrid(bedmachine_filename, 'bed')
    bed = np.ma.masked_array(bed,mask2)
    surf = np.ma.masked_array(surf,mask2)
    thickness = np.ma.masked_array(thickness,mask2)
    
    version = bedmachine_filename.split('/')[-1][-6:-3]
    
    return x,y,bed,surf,thickness, version


def grounding(thickness, bed, ri=0.917, rw = 1.028, z_sl=0):
    cond = thickness*ri/rw + bed
    gmask = np.ones_like(cond)
    gmask[cond<0] = -1  #floating
    gmask[cond>0] = 1 #grounded
    return gmask


def new_floating_surf(thickness, bed, gmask, ri=0.917, rw = 1.028, zs_sl=0):
    surface_ground = bed + thickness
    surface_float = thickness*(1-ri/rw)
    grounded_mask = 1+((gmask-1)/2) 
    floating_mask = (-((gmask-1)/2))
    surface = surface_ground*grounded_mask + surface_float*floating_mask
    return surface


def get_GL(x, y, gmask):
    print('Extracting grounding line...')
    x_gl,y_gl = [],[]
    n, m = len(x), len(y)
    for i in range(n):
        for j in range(m):
            if gmask[i,j] == 1 and np.mean(gmask[i-1:i+2,j-1:j+2])!=1: #on regarde ici si le point i,j égal à 1 (grounded) cotoit un point flottant
                x_gl.append(x[j])
                y_gl.append(y[i])
        old_i = 0
        # if int((i/n+1e-2)*100) != old_i:
        #     sys.stdout.write('\r')
        #     sys.stdout.flush()
        #     # the exact output you're looking for:
        #     #sys.stdout.write("[%-20s] %d%%" % ('='*int(i/n*20+1), (i/n+1e-2)*100))
        #     sys.stdout.write(" [%d%%]" % ((i/n+1e-2)*100))
        #     old_i = int((i/n+1e-2)*100)
        #     #sys.stdout.write("[{:{}}] {:.1f}%".format("="*i, n-1, (100/(n-1)*i)))
        #     sys.stdout.flush()
    print('\n')
    return x_gl, y_gl


def create_directory(directory_path):
    """
    Create a directory if it does not exist.

    Parameters:
    - directory_path (str): The path of the directory to be created.
    """
    try:
        # Check if the directory already exists
        if not os.path.exists(directory_path):
            # Create the directory
            os.makedirs(directory_path)
            print(f"Directory '{directory_path}' created successfully.")
        else:
            print(f"Directory '{directory_path}' already exists.")
    except Exception as e:
        print(f"Error creating directory '{directory_path}': {str(e)}")

        
def compute_rate(year0, year1, power_rate):
    annual_rate = []  
    years = np.arange(year0, year1,1)
    for i,year in enumerate(years[:-1]):
        annual_rate.append(((years[i]-years[0])/(years[-1]-years[0]))*power_rate)
    rate = sum(annual_rate)/(years[-1]-years[0])
    return rate

#%% Loading files
x, y, bed, surf, thickness, version = read_bedmachine(resolution=int(resolution*2))
x_IS, y_IS, dhdt_IS = Build_IceSat2_ContinuousMap(x_final=x_fin, y_final=y_fin)



for targetyear in years : 
    print(targetyear)
    #make directory

    create_directory(f'{targetyear}_State')


    #%% build dhdt trend and corresponding dh
    #the trend is for the period 2003-2019, we assume it is similar over our period of reconstruction
    #extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)

    extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)]

    dhdt_IS = dhdt_IS[::-1,:]

    if dhdt_IS.data.shape != surf.mask.shape:
        #surf.mask = np.resize(surf.mask, dhdt_IS.data.shape)

        scale_factors = (surf.shape[0] / dhdt_IS.shape[0], surf.shape[1]/ dhdt_IS.shape[1])

        dhdt_IS = zoom(dhdt_IS, scale_factors,order = 1)
        print(f'dhdt_IS shape is : {dhdt_IS.shape}' )
        print(f'surf mask shape is : {surf.mask.shape}' )


        
def compute_rate(year0, year1, power_rate):
    annual_rate = []  
    years = np.arange(year0, year1,1)
    for i,year in enumerate(years[:-1]):
        annual_rate.append(((years[i]-years[0])/(years[-1]-years[0]))*power_rate)
    rate = sum(annual_rate)/(years[-1]-years[0])
    return rate








#%% Loading files
x, y, bed, surf, thickness, version = read_bedmachine(resolution=int(resolution*2))
x_IS, y_IS, dhdt_IS = Build_IceSat2_ContinuousMap(x_final=x, y_final=y)

for targetyear in years : 
    print(targetyear)




    #make directory

    create_directory(f'{targetyear}_State')
    
    #import the mask
    #mask_path = f'../../../DATA/Mask/mask_PIG_{targetyear}_mask_full_mask.tif'
    mask_dic=import_mask(f'../../../DATA/Mask/mask_PIG_{targetyear}_mask_full_mask.tif',xres=0.5,yres=0.5)
    #mask_dic=import_mask(f'../../../DATA/Mask/mask_PIG_{targetyear}_mask_full_mask.tif')
    mask_year=mask_dic['data']
    print(mask_year)
    print("Valeurs masque romain : "+str(np.unique(mask_year)))

    #mask_year=np.nan_to_num(x=mask_year, nan=10.0)
    #print(mask_year.shape)
    
    #%% build dhdt trend and corresponding dh
    #the trend is for the period 2003-2019, we assume it is similar over our period of reconstruction
    # extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)] 

    extent_IS = [np.min(x_fin), np.max(x_fin), np.min(y_fin), np.max(y_fin)]

    dhdt_IS = dhdt_IS[::-1,:]
    dhdt_IS = np.ma.masked_array(dhdt_IS.data, surf.mask)


    rate = compute_rate(targetyear, 2022, power_rate)
    dh_IS_targetyear = dhdt_IS*(2022-targetyear)*rate 
    dh_IS_2020 = dhdt_IS*(2022-2020)


    base = surf - thickness

    gmask = np.ones_like(base)
    gmask[base-bed>0] = -1

    # target year
    thickness_targetyear = thickness - np.nan_to_num(dh_IS_targetyear)
    gmask_targetyear = grounding(thickness_targetyear, bed)
    surf_targetyear = new_floating_surf(thickness_targetyear, bed, gmask_targetyear)
    bed_targetyear = surf_targetyear - thickness_targetyear
    bedrock_targetyear = bed

    #print(gmask_targetyear.shape)
    #gmask_targetyear=np.delete(gmask_targetyear,-1,axis=1)
    #gmask_targetyear=np.delete(gmask_targetyear,-1,axis=0)
    #difference_mask=mask_year-gmask_targetyear
    
    #new_gmask = gmask_targetyear.copy()

    #new_gmask[difference_mask==1]=0#le front a reculé
    #new_gmask[difference_mask==3]=-1 #le front a avance

    if subplot_mask :
        # Pine Island
        x_min, x_max = -1694616.000, -1506616.000
        y_min, y_max = -371234.000, -204734.000
    
        # Créer un masque pour les coordonnées x et y
        x_mask = (x >= x_min) & (x <= x_max)
        y_mask = (y >= y_min) & (y <= y_max)


        # Appliquer le masque aux variables surface
        x_mask_subset = x[x_mask]
