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

import sys
import os

#chemin vers mes modules
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__),'../MyModule'))

if module_path not in sys.path:
    sys.path.append(module_path)

#from Antarctica_Background import Plot_Antarctica, basin, scale,plot_GL
from ordering import GL_ordering
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
from scipy.ndimage import gaussian_filter
import sys
from scipy import interpolate
import os
import math
from osgeo import gdal
import xarray as xr

#-------------------------------------------------------------------------
# USER KEYWORDS
#-------------------------------------------------------------------------
plot_mask = False
subplot_mask = False
plot = False
MakeMesh = True #g√©n√©re automatiquement les fichiers contours 1 et 2 necessaire a  la cre©ation du mesh
ri = 0.917
rw = 1.028
z_sl = 0
years = [1995, 2000, 2005, 2011,  2014,  2017]
#years=[2017]
n=len(years)
#domaine final
resolution=5 #km

x_fin=np.array(list(range(-1694617, -1506616,int(resolution*(10**3)) )))
y_fin=np.array(list(range(-371235, -204734, int(resolution*(10**3)))))
x_min, x_max = -1694616.000, -1506616.000
y_min, y_max = -371234.000, -204734.000


power_rate = 1 # for long period, we can use dhdt rates that decrease with t
#-------------------------------------------------------------------------

#%%Read dHdT from IceSat1-2 (over the shelves)
def Build_dhdt(x_final=None, y_final=None, masked=False, year=None, resolution=5, i=None):
    print("Building IceSat dhdt...")
    #Ice Sheet
    
    filename  = '../../../DATA/BEDROCK/ANT_G1920V01_GroundedIceHeight_mean.nc'
    ds_gd = xr.open_dataset(filename)
    ds_gd_2015 = ds_gd.sel(time=str(2015))
    ds_gd_year = ds_gd.sel(time=str(year))
    #x_gd, y_gd = np.meshgrid(ds_gd_year['x'].values, ds_gd_year['y'].values, indexing='ij')
    #print("x_gd shape:", x_gd.shape)

    x_gd= ds_gd_year['x'].values
    y_gd= ds_gd_year['y'].values
    dhdt_gd_year = np.squeeze(ds_gd_year['dh'].values)
    dhdt_gd_2015 = np.squeeze(ds_gd_2015['dh'].values)
    #accrochez vous, calculs savants
    dhdt_gd = dhdt_gd_year - (dhdt_gd_2015) # year - 2015 = year - 2013 - (2015-2013)

    #Ice Shelf
    filename = '../../../DATA/BEDROCK/ANT_G1920V01_IceShelfMelt_mean.nc'
    ds_is=xr.open_dataset(filename)
    ds_is_year = ds_is.sel(time=str(year))
    ds_is_2015 = ds_is.sel(time=str(2015))

    #x_is, y_is = np.meshgrid(ds_is_year['x'].values, ds_is_year['y'].values, indexing='ij')
    x_is= ds_is_year['x'].values
    y_is= ds_is_year['y'].values
    h_is_year=np.squeeze(ds_is_year['thickness'].values)
    h_is_2015=np.squeeze(ds_is_2015['thickness'].values) # importation des epaisseurs
    dhdt_is = h_is_year - h_is_2015 # variation d'√©paisseur par rapport a 2014
    dhdt_is = dhdt_is * ((rw-ri)/rw)# variation d'altitude
    #dhdt_is =np.squeeze(ds_is_year['height_change'].values)

    #mega plot, pt 1
    ax=axes[0, i]  
    mp = ax.imshow(dhdt_gd, cmap='bwr', vmin = -20, vmax= 20)
    fig.colorbar(mp, ax=ax,label=f'dhdt_gd, {year}')
    ax=axes[1, i]  
    mp = ax.imshow(dhdt_is, cmap='bwr', vmin = -20, vmax= 20)
    fig.colorbar(mp, ax=ax,label=f'dhdt_is, {year}')
    ax=axes[2, i]  
    mp = ax.imshow(dhdt_gd, cmap='bwr', vmin = -20, vmax= 20)
    mp = ax.imshow(dhdt_is, cmap='bwr', vmin = -20, vmax= 20)
    fig.colorbar(mp, ax=ax,label=f'dhdt_gd + dhdt_is, {year}')


    x_gd_grid, y_gd_grid = np.meshgrid(x_gd, y_gd, indexing='xy')
    x_is_grid, y_is_grid = np.meshgrid(x_is, y_is, indexing='xy')

    smoothing=True
    if smoothing :
        #Create a mask of distance to the grounding line (in pixels)
        from scipy.ndimage import distance_transform_edt
        nan_mask = np.isnan(dhdt_gd)
        distances = distance_transform_edt(~nan_mask)

        dhdt_gd[distances<10]=np.nan #masking all data that are within 5 pixels of the GL
        #Below I am merging floating and grounded ice
        dhdt_gd[np.isnan(dhdt_gd)] = dhdt_is[np.isnan(dhdt_gd)]

        #Then I linearly interpolate between where I have nans
        # Step 1: Interpolate NaN values
        # Create a grid of indices
        #x, y = np.indices(dh.shape)
        valid_mask = ~np.isnan(dhdt_gd)  # Mask for non-NaN values

        # Get the coordinates and values of valid (non-NaN) points
        x_valid = x_gd_grid[valid_mask]
        y_valid = y_gd_grid[valid_mask]
        values_valid = dhdt_gd[valid_mask]

        # Get the coordinates of NaN values (for interpolation)
        x_nan = x_gd_grid[np.isnan(dhdt_gd)]
        y_nan = y_gd_grid[np.isnan(dhdt_gd)]

        # Use griddata for interpolation (linear interpolation)
        interpolated_values = griddata(
            (x_valid, y_valid), values_valid, (x_nan, y_nan), method='linear'
        )

        # Assign the interpolated values back to the original array
        dhdt_gd[np.isnan(dhdt_gd)] = interpolated_values

    ax=axes[3, i]  
    mp = ax.imshow(dhdt_gd, cmap='bwr', vmin = -20, vmax= 20)
    mp = ax.imshow(dhdt_is, cmap='bwr', vmin = -20, vmax= 20)
    fig.colorbar(mp, ax=ax,label=f'dhdt_gd + dhdt_is, after interpolation {year}')



    x_gd = x_gd[::resolution]
    y_gd = y_gd[::resolution]
    dhdt_gd = dhdt_gd[::resolution, ::resolution]
    x_is = y_is[::resolution]
    y_is = y_is[::resolution]
    dhdt_is = dhdt_is[::resolution, ::resolution] 


    mask_x = (x_gd >= x_min) & (x_gd <= x_max)
    mask_y = (y_gd >= y_min) & (y_gd <= y_max)
    x_gd = x_gd[mask_x]
    x_is = x_is[mask_x]
    y_gd = y_gd[mask_y]
    y_is = y_is[mask_y]
    dhdt_gd = dhdt_gd[np.ix_(mask_y, mask_x)]
    dhdt_is = dhdt_is[np.ix_(mask_y, mask_x)]
    
    x_gd_grid, y_gd_grid = np.meshgrid(x_gd, y_gd, indexing='xy')
    x_is_grid, y_is_grid = np.meshgrid(x_is, y_is, indexing='xy')


    #print("x_gd_grid shape:", x_gd.shape)
    #print("y_gd_grid shape:", y_gd.shape)
    #print("dhdt_gd shape:", dhdt_gd.shape)
    #print("x_is_grid shape:", x_is.shape)
    #print("y_is_grid shape:", y_is.shape)
    #print("dhdt_is shape:", dhdt_is.shape)

    #Interpolation of the shelf data on the ground data
    print("\tInterpolating over same grid...")
    lim = 1
    #x_share, y_share = x_gd[lim:-lim], y_gd[lim:-lim] 
    x_share_grid, y_share_grid = x_gd_grid, y_gd_grid
    # print('xground = %d' % x_share.min())
    points_share = np.column_stack((x_share_grid.ravel(), y_share_grid.ravel()))
    
    #make the interpolator
    #interpolator = RegularGridInterpolator((y_is, x_is), dhdt_is, bounds_error=False)
    #interpolation
    #data_shelf_share = interpolator(points_share)
    #data_shelf = data_shelf_share.reshape(x_share_grid.shape)
    #data_shelf=data_shelf.T

    

    #mask 
    mask = np.logical_and(np.isnan(dhdt_is), np.isnan(dhdt_gd))
    n, m = mask.shape
    
    mask2 = copy.copy(mask)
    for i in np.arange(n):
        for j in np.arange(m):
            if mask[i,j]:
                if False in mask[i-1:i+2,j-1:j+2]:
                    mask2[i,j] = False
                    
    data = np.nan_to_num(dhdt_is) + np.nan_to_num(dhdt_gd)
    data[data==0.0]=np.nan
    
    #mask NaNs
    data = np.ma.masked_invalid(data)
    print(' data shape ', data.shape)
    print(' x_share_grid ', x_share_grid.shape) 
    print(' y_share_grid ', y_share_grid.shape) 
    
    #get only the valid values
    x1 = x_share_grid[~data.mask]
    y1 = y_share_grid[~data.mask]
    new_data = data[~data.mask]
    
    new_data = interpolate.griddata((x1, y1), new_data.ravel(), (x_share_grid, y_share_grid),method='linear')
    
    data = np.ma.MaskedArray(new_data, mask2)
    
    #plt.imshow(data, cmap='bwr', vmin=-5, vmax=5)
    #plt.colorbar(label='var altitude avant interp')
    #plt.show()


    #print(' data shape ', data.shape)
    #
    #print("y_share_grid shape ",y_share_grid[0,:].shape )
    #print("y_share_grid shape",x_share_grid[:,0].shape )

        
    #reshape to a given output grid
    resize=True
    if (x_final is not None) and (y_final is not None) and resize:
        print("\tRegridding")
        print("\t - origin boundaries : x = [%d, %d], y = [%d, %d]" % (x_share_grid.min(), x_share_grid.max(), y_share_grid.min(), y_share_grid.max()))
        print("\t - target boundaries : x = [%d, %d], y = [%d, %d]" % (x_final.min(), x_final.max(), y_final.min(), y_final.max()))
        
        y_offset = 5e4 #I have to apply offset but I don't know why
       
        #print("y_share_grid",y_share_grid[:,0] )
        #print("x_share_grid",x_share_grid[0,:] )
        interpolator2 = RegularGridInterpolator((y_share_grid[:,0] , x_share_grid[0,:] ), data, bounds_error=False)
        xx,yy = np.meshgrid(x_final,y_final)
        #print(' xx et yy  ', xx,'\n', yy[::-1])


        points_final = list(zip(yy[::-1].flatten(), xx.flatten()))
        #print(' points_final shape ', len(points_final))

        data = interpolator2(points_final)
        data = data.reshape((len(y_final),len(x_final)))
        x, y = y_final, x_final
        
    else: 
        x, y = x_share_grid, y_share_grid
    print("\tDone.")

    #plt.imshow(data, cmap='bwr', vmin=-3,vmax=3)
    #plt.colorbar(label="variation d'altitude apr√©s interp n2")
    #plt.show()
    return x, y, data

def import_mask(path_file,xres,yres):
    xres=xres*1000
    yres = yres*1000
    os.system('gdalwarp -tr '+str(xres)+' '+str(yres)+' -r average '+path_file+' '+path_file.split('.tif')[0]+'_reproj.tif')
    dataset = gdal.Open(path_file.split('.tif')[0]+'_reproj.tif')
    
    if dataset is None:
        raise FileNotFoundError(f"File {path_file} cannot be found")
   
    # Obtenir la transformation gÈographique
    geo_transform = dataset.GetGeoTransform()
    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_res = geo_transform[1]
    y_res = geo_transform[5]  # Note : nÈgatif pour les rasters orientÈs "top-down"

    #print(f"ymax {y_max}")

    # Obtenir les dimensions
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
  
    #print(rows)  
    #print("ymin"+str(y_max+ (rows * y_res)))

    # Calculer les coordonnÈes x et y
    x = np.arange(x_min, x_min + cols*x_res, x_res)
    y = np.arange( y_max, y_max + (rows * y_res),  y_res)
    
    #print("x")
    #print(x)
    #print("y")
    #print(y)


    # RÈcupËre la premiËre bande (en supposant que le masque soit en une seule bande)
    bande = dataset.GetRasterBand(1)
    mask_data = bande.ReadAsArray()
    dataset = None

    #gmask = np.full_like(mask_data, np.nan, dtype=np.float32)

    return x, y, mask_data

def read_bedmachine(resolution=5, x_min=None,x_max=None,y_min=None, y_max=None):
    print('Reading bedmachine...')
    print('\t -using a %0.1f km resolution' % (resolution*0.5))
    bedmachine_filename = '../../../DATA/BEDROCK/BedMachineAntarctica-v03.nc'
    
    x, y, bed = netcdf.readgrid(bedmachine_filename, 'bed')
    x, y, surf = netcdf.readgrid(bedmachine_filename, 'surface')
    x, y, thickness = netcdf.readgrid(bedmachine_filename, 'thickness')
    x, y, mask = netcdf.readgrid(bedmachine_filename, 'mask')
    
    

    # Filtrer selon les limites spÈcifiÈes
    
    x_mask = (x >= x_min-500) & (x <= x_max)
    y_mask = (y >= y_min-500) & (y <= y_max)
    x = x[x_mask]
    y = y[y_mask]
    bed = bed[np.ix_(y_mask, x_mask)]
    surf = surf[np.ix_(y_mask, x_mask)]
    thickness = thickness[np.ix_(y_mask, x_mask)]
    mask = mask[np.ix_(y_mask, x_mask)]


    #x_inf = np.argmin(np.abs(x-np.min(x_fin)))-1
    #x_sup = np.argmin(np.abs(x-np.max(x_fin)))
    #y_sup = np.argmin(np.abs(y-np.min(y_fin)))
    #y_inf = np.argmin(np.abs(y-np.max(y_fin)))-1
    #print('x_inf ', x_inf, ' x_sup ',x_sup, ' y_min ', y_inf, ' y_sup ', y_sup)
    #print()
    
    y_new = y[::resolution]
    x_new = x[::resolution]
    print('len x ', len(x_new), ' len y', len(y_new))


    bed = bed[::resolution, ::resolution]
    surf = surf[::resolution, ::resolution]
    thickness = thickness[::resolution, ::resolution]
    mask=mask[::resolution, ::resolution]

    print('bed shape ', bed.shape)
    #plt.imshow(thickness, cmap='Purples', vmin=0)
    #plt.show()

    mask2 = np.zeros_like(mask)
    mask2[mask==0] = True
    mask2[mask>0.5] = False
    
    bed = np.ma.masked_array(bed,mask2)
    surf = np.ma.masked_array(surf,mask2)
    thickness = np.ma.masked_array(thickness,mask2)
    
    version = bedmachine_filename.split('/')[-1][-6:-3]
    
    return x_new,y_new,bed,surf,thickness, version


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
    n, m = len(y), len(x)
    print('valeurs gmask ', np.unique(gmask))
    for i in range(n):
        for j in range(m):
            if gmask[i,j] == 1 and np.mean(gmask[i-1:i+2,j-1:j+2])!=1:
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
    return y_gl, x_gl

def get_FL(x, y, gmask):
    print('Extracting front line...')
    x_fl,y_fl = [],[]
    n, m = len(y), len(x)
    for i in range(n-1):
        for j in range(m-1):
            if i > 0 and i < n-2 and j > 0 and j < m-2:
                if (mask_year[i,j] == 0 | np.isnan(mask_year[i,j])) and np.max(mask_year[i-1:i+2,j-1:j+2]) > 2:
                    x_fl.append(x[j])
                    y_fl.append(y[i])
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
    return x_fl, y_fl



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
x, y, bed, surf_2015, thickness_2015, version = read_bedmachine(resolution=int(resolution*2), x_min = np.min(x_fin),x_max = np.max(x_fin), y_min = np.min(y_fin), y_max = np.max(y_fin)  )
#x_IS, y_IS, dhdt_IS = Build_IceSat2_ContinuousMap(x_final=x, y_final=y)
fig, axes=plt.subplots(5, n, figsize=(5*n, 25), sharey= True)
i=0
for i, targetyear in enumerate(years) : 
    print(targetyear)
    

    x_IS,y_IS, dhdt_IS = Build_dhdt(x_final=x_fin, y_final=y_fin, year=targetyear, resolution=int(resolution*2),i=i)

    #make directory

    create_directory(f'{targetyear}_State')
    
    #import the mask
    
    x_mask_year, y_mask_year, mask_year = import_mask(f'../../../DATA/Mask/mask_PIG_{targetyear}_mask_full_mask_crop.tif',xres=resolution,yres=resolution)
    


    #x_mask_year_max, y_mask_year_max = mask_year.shape
    #x_mask_year = range(1,x_mask_year_max)
    #y_mask_year = range(1,y_mask_year_max)

    #%% build dhdt trend and corresponding dh
    #the trend is for the period 2003-2019, we assume it is similar over our period of reconstruction
    # extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)] 

    extent_IS = [np.min(x_fin), np.max(x_fin), np.min(y_fin), np.max(y_fin)]

    #dhdt_IS = dhdt_IS[::-1,:]
    #dhdt_IS = np.ma.masked_array(dhdt_IS.data, surf.mask)
    
    
    base = surf_2015 - thickness_2015
    gmask = np.ones_like(base)
    gmask[base-bed>0] = -1

    # target year
    thickness_targetyear = thickness_2015 + np.nan_to_num(dhdt_IS)
    gmask_targetyear = grounding(thickness_targetyear, bed)
    surf_targetyear = new_floating_surf(thickness_targetyear, bed, gmask_targetyear)
    bed_targetyear = surf_targetyear - thickness_targetyear
    bedrock_targetyear = bed


    ax=axes[4,i]
    th = ax.imshow(thickness_targetyear, cmap='Purples')
    fig.colorbar(th, ax=ax, label='thickness m')
    


    #plt.imshow(gmask_targetyear, cmap='inferno')
    #plt.colorbar(label='gmask')
    #plt.show()
    
    # Afficher le masque sous forme d'image
    #plt.figure(figsize=(10, 8))
    #plt.imshow(gmask_targetyear, cmap="Pastel1", interpolation="nearest")
    #plt.colorbar(label="Valeurs du masque")
    #plt.title("Apercu du masque")
    #plt.xlabel("Longitude (pixels)")
    #plt.ylabel("Latitude (pixels)")
    #plt.show()




    if subplot_mask :
        # Pine Island
        x_min, x_max = -1694616.000, -1506616.000
        y_min, y_max = -371234.000, -204734.000
    
        # CrÈer un masque pour les coordonnÈes x et y
        x_mask = (x >= x_min) & (x <= x_max)
        y_mask = (y >= y_min) & (y <= y_max)


        # Appliquer le masque aux variables surface
        x_mask_subset = x[x_mask]
        y_mask_subset = y[y_mask]
        diffmask_subset = difference_mask[np.ix_(y_mask, x_mask)]
        gmask_subset = gmask_targetyear[np.ix_(y_mask, x_mask)]
        newmask_subset = new_gmask[np.ix_(y_mask, x_mask)]


        # CrÈation des sous-graphes pour afficher les trois masques
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        # Affichage de chaque masque dans un sous-graphe
        im=axes[0].imshow(gmask_subset, cmap='Pastel1', vmin=-2,vmax=4)
        axes[0].set_title('Masque modellise (gmask_targetyear)')
        axes[0].axis('off')

        im=axes[1].imshow(diffmask_subset, cmap='Pastel1' , vmin=-2,vmax=4)
        axes[1].set_title('Masque difference')
        axes[1].axis('off')

        im=axes[2].imshow(newmask_subset, cmap='Pastel1', vmin=-2,vmax=4)
        axes[2].set_title('Nouveau masque')
        axes[2].axis('off')

        # Affichage du graphique complet
        cbar=fig.colorbar(im, ax=axes,orientation='horizontal', fraction=0.05, pad=0.1)
        cbar.set_label('Valeurs')
        plt.tight_layout()
        #plt.show()
        plt.savefig(f"{targetyear}_State/Masques_Is")

    if plot_mask :




        # Pine Island
        x_min, x_max = -1694616.000, -1506616.000
        y_min, y_max = -371234.000, -204734.000
    
        #print(f"x shape {x.shape}, y shape : {y.shape}, mask shape {mask_year.shape}")

        # CrÈer un masque pour les coordonnÈes x et y
        x_mask = (x_mask_year >= x_min) & (x_mask_year <= x_max)
        y_mask = (y_mask_year >= y_min) & (y_mask_year <= y_max)


        # Appliquer le masque aux variables surface
        x_mask_subset = x_mask_year[x_mask]
        y_mask_subset = y_mask_year[y_mask]
        mask_subset = mask_year[np.ix_(y_mask, x_mask)]
        #mask_subset = gmask_targetyear[y_mask, x_mask]


        plt.imshow(mask_subset, cmap='Pastel1', extent=[x_min, x_max, y_min, y_max], vmin=-1, vmax=4)
        plt.colorbar(label='Values')
        plt.title(f'Grouded mask for {targetyear}')
        plt.xlabel('Longitude')
        plt.ylabel('Latitude')
        plt.savefig(f"{targetyear}_State/masque_romain_500")
        plt.show()



    #2020
    #thickness_2020 = thickness - np.nan_to_num(dh_IS_2020)
    #gmask_2020 = grounding(thickness_2020, bed)
    #surf_2020 = new_floating_surf(thickness_2020, bed, gmask_2020)
    #bed_2020 = surf_2020 - thickness_2020
    #bedrock_2020 = bed

    #%%
    thickness_targetyear[thickness_targetyear<10]=10.
    thickness_targetyear[thickness_targetyear<10]=10.

    #%% Get grounding lines for target year and 2020
    GL = True
    if GL:
        #x_gl_2020, y_gl_2020 = get_GL(x,y,gmask_2020) 
        y_gl_targetyear, x_gl_targetyear = get_GL(x,y,gmask_targetyear)
        print('extrems de x GL trouves bruts', np.min(x_gl_targetyear), np.max(x_gl_targetyear))
        print('extrems de y GL trouves bruts', np.min(y_gl_targetyear), np.max(y_gl_targetyear))

        #np.savetxt('2020_State/GL_2020.txt', [x_gl_2020, y_gl_2020])
        np.savetxt(f'{targetyear}_State/GL_{targetyear}.txt', [x_gl_targetyear, y_gl_targetyear])
        #x_fl_2020, y_fl_2020 = get_FL(x,y,gmask_2020) 
        x_fl_targetyear, y_fl_targetyear = get_FL(x_mask_year,y_mask_year,mask_year)

        #np.savetxt('2020_State/FL_2020.txt', [x_fl_2020, y_fl_2020])
        np.savetxt(f'{targetyear}_State/FL_{targetyear}.txt', [x_fl_targetyear, y_fl_targetyear])


        if MakeMesh :
            x_gl_targetyear_fin =[]
            y_gl_targetyear_fin =[]
            x_fl_targetyear_fin =[]
            y_fl_targetyear_fin =[]

            with open(f'{targetyear}_State/contours_GL_{targetyear}_paolo.txt','w') as GL:
                o=0
                #x_gl_targetyear,y_gl_targetyear = GL_ordering(x_gl_targetyear,y_gl_targetyear)
                for i in range(len(x_gl_targetyear)):
                    if x_gl_targetyear[i] >= x_fin[0] and x_gl_targetyear[i] <= x_fin[-1] and y_gl_targetyear[i] >= y_fin[0] and y_gl_targetyear[i] <= y_fin[-1]:
                        o=o+1    
                        GL.write( str(x_gl_targetyear[i])+'\t'+str(y_gl_targetyear[i])+'\n')
                        x_gl_targetyear_fin.append(x_gl_targetyear[i])
                        y_gl_targetyear_fin.append(y_gl_targetyear[i])
                print('Point des GL trouv√©s', o)
                #x_gl_fin, y_gl_GL_fin = GL_ordering(x_gl_targetyear_fin, y_gl_targetyear_fin)
                #for i in range(len(x_gl_fin)):
                    #GL.write( str(x_gl_targetyear[i])+'\t'+str(y_gl_targetyear[i])+'\n')
            GL.close()

            print('Contours GL done, nombre de points : '+str(o))         
            with open(f'{targetyear}_State/contours_FL_{targetyear}_paolo.txt', 'w') as FL:
                print(len(x_fl_targetyear))
                o=0
                for i in range(len(x_fl_targetyear)):
                    #if x_fl_targetyear[i] >= x_fin[0] and x_fl_targetyear[i] <= x_fin[-1] and y_fl_targetyear[i] >= y_fin[0] and y_fl_targetyear[i] <= y_fin[-1]:
                    o=o+1
                    FL.write( str(x_fl_targetyear[i])+'\t'+str(y_fl_targetyear[i])+'\n')
                    x_fl_targetyear_fin.append(x_fl_targetyear[i])
                    y_fl_targetyear_fin.append(y_fl_targetyear[i])
            FL.close()
            print('contours FL done, nombre de points : '+str(o))


    else:
        x_gl_2020, y_gl_2020 = np.loadtxt('2020_State/GL_2020.txt')
        x_gl_targetyear, y_gl_targetyear = np.loadtxt(f'{targetyear}_State/GL_{targetyear}.txt')
    #%%
    from Antarctica_Background import Plot_Antarctica, basin, scale, plot_front, basemap_LIMA
    
    extent=[x_fin[0],x_fin[-1], y_fin[0], y_fin[-1]]
    #print(extent)

    #fig, ax = plt.subplots(figsize =(20,20))
    #ax.imshow(dhdt_IS*rate, cmap = 'seismic',extent=extent, vmin = -10, vmax = 10, alpha = 0.9, zorder = 1e1)
    #ax[0].scatter(x_gl_2020, y_gl_2020, c='black', s=0.1, linewidths=0, zorder= 1e5, extent=extent)
    #ax.scatter(x_gl_targetyear_fin, y_gl_targetyear_fin, c='orange', s=0.1, linewidths=0, zorder = 1e5)
    #ax.scatter(x_fl_targetyear_fin, y_fl_targetyear_fin, c='blue', s=0.1, linewidths=0, zorder = 1e5)

    #plt.savefig(f'Figures/PIG_dHdt_GL{targetyear}.png', bbox_inches='tight', pad_inches=0.1, transparent = True)


    if plot:
        extent = np.asarray(basin.PanAntarctic())
        fig, ax = Plot_Antarctica(nrows=1, ncols=1, basemap = 'light', GL = False, icefront = False, continental_shelf=0.0, extent=extent, figsize = (30,25))
        #basemap_LIMA_AMU(ax[0])
        cb = ax[0].imshow(dhdt_IS*rate, cmap = 'seismic',extent=extent_IS, vmin = -10, vmax = 10, alpha = 0.9, zorder = 1e1)
        ax[0].scatter(x_gl_2020, y_gl_2020, c='black', s=0.1, linewidths=0, zorder= 1e5)
        ax[0].scatter(x_gl_targetyear, y_gl_targetyear, c='orange', s=0.1, linewidths=0, zorder = 1e5)
        ax[0].scatter(x_fl_targetyear, y_fl_targetyear, c='blue', s=0.1, linewidths=0, zorder = 1e5)
        cbar = ax[0].cax.colorbar(cb)
        cbar = ax.cbar_axes[0].colorbar(cb)
        cbaxes = fig.add_axes([0.33, 0.35, 0.12, 0.013]) 
        cbar = fig.colorbar(cb, cax=cbaxes, ticks=[-10,0,10],orientation='horizontal')
        cbar.set_label('ice thickness rate of change', fontsize=16)
        cbar.ax.set_xticklabels(['-10', '0', r'   10 m a$^{-1}$'], fontsize=13)   
        ax[0].axis('off')

        #plot_GL(ax[0], '1996', color = 'green', lw = 0.3, zorder = 1e5)
    
        plt.savefig(f'Figures/Antarctica_dHdt_GL{targetyear}.pdf', bbox_inches='tight', pad_inches=0.1, transparent = True)
        plt.close('all')
    #%%
    #targetyear State
    print(f'x is {x.shape}')
    #print(x)
    print(f'y is {y.shape}')
    #print(y)
    print(f'surf is {surf_targetyear.shape}')
    print(f'thick is {thickness_targetyear.shape}')
    print(f'bed is {bedrock_targetyear.shape}')

    i = i+1
    dic = {'x':x, 'y':y, 'surface':surf_targetyear, 'bed':bedrock_targetyear, 'thickness':thickness_targetyear}
    netcdf.write(f'{targetyear}_State/BedMachine_IS_{targetyear}_{version}_paolo.nc' , dic)


plt.tight_layout()
plt.savefig("epaisseur_paolo_variation.png")
plt.show()
    #2020 State
    #dic = {'x':x, 'y':y, 'surface':surf_2020, 'bed':bedrock_2020, 'thickness':thickness_2020}
    #netcdf.write('2020_State/BedMachine_IS_2020_{version}.nc', dic)
