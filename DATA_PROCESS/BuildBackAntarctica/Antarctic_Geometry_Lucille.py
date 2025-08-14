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
from netCDF4 import Dataset, num2date
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

#-------------------------------------------------------------------------
# USER KEYWORDS
#-------------------------------------------------------------------------
plot_mask = False
subplot_mask = False
plot = False
MakeMesh = True #gÃ©nÃ©re automatiquement les fichiers contours 1 et 2 necessaire a  la cre©ation du mesh
ri = 0.917
rw = 1.028
z_sl = 0
#years = [1991, 2000, 2002, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016, 2017]
#years=[2008,2011,2015]
years=[  2018,2020, 2022]
#domaine final
x_fin=np.array(list(range(-1694617, -1506616, 100)))
y_fin=np.array(list(range(-371235, -204734, 100)))
resolution = 0.1 #km

#-------------------------------------------------------------------------


def import_mask(path_file,xres,yres):
    xres=xres*1000 #(km en m)
    yres = yres*1000
    os.system('gdalwarp -tr '+str(xres)+' '+str(yres)+' -r average '+path_file+' '+path_file.split('.tif')[0]+'_reproj.tif')
    dataset = gdal.Open(path_file.split('.tif')[0]+'_reproj.tif')
    
    if dataset is None:
        raise FileNotFoundError(f"File {path_file} cannot be found")
   
    # Obtenir la transformation géographique
    geo_transform = dataset.GetGeoTransform()
    x_min = geo_transform[0]
    y_max = geo_transform[3]
    x_res = geo_transform[1]
    y_res = geo_transform[5]  # Note : négatif pour les rasters orientés "top-down"

    print(f"ymax {y_max}")

    # Obtenir les dimensions
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
  
    print(rows)  
    print("ymin"+str(y_max+ (rows * y_res)))

    # Calculer les coordonnées x et y
    x = np.arange(x_min, x_min + cols*x_res, x_res)
    y = np.arange( y_max, y_max + (rows * y_res),  y_res)
    
    #print("x")
    #print(x)
    #print("y")
    #print(y)


    # Récupère la première bande (en supposant que le masque soit en une seule bande)
    bande = dataset.GetRasterBand(1)
    mask_data = bande.ReadAsArray()
    dataset = None

    #gmask = np.full_like(mask_data, np.nan, dtype=np.float32)

    return x, y, mask_data

def read_bedmachine(resolution=5):
    print('Reading bedmachine...')
    print('\t -using a %0.1f km resolution' % (resolution*0.5))
    bedmachine_filename = '../../../DATA/BEDROCK/BedMachineAntarctica-v03.nc'
    
    x, y, bed = netcdf.readgrid(bedmachine_filename, 'bed')
    x, y, surf = netcdf.readgrid(bedmachine_filename, 'surface')
    x, y, thickness = netcdf.readgrid(bedmachine_filename, 'thickness')
    x, y, mask = netcdf.readgrid(bedmachine_filename, 'mask')
    dx_native = np.abs(x[1] - x[0])
    step = int(round(resolution / dx_native))
    print('\t -Step is %0.1f' % (step))
    if step < 1:
        step = 1

    x_inf, x_sup = +0, -1
    y_inf, y_sup = +0, -1
    x = x[x_inf:x_sup:step]
    y = y[y_inf:y_sup:step]
    bed = bed[x_inf:x_sup:step, y_inf:y_sup:step]
    surf = surf[x_inf:x_sup:step, y_inf:y_sup:step]
    thickness = thickness[x_inf:x_sup:step, y_inf:y_sup:step]
    mask=mask[x_inf:x_sup:step, y_inf:y_sup:step]
    
    mask2 = np.zeros_like(mask)
    mask2[mask==0] = True
    mask2[mask>0.5] = False
    
    #bed = np.ma.masked_array(bed,mask2)
    #surf = np.ma.masked_array(surf,mask2)
    #thickness = np.ma.masked_array(thickness,mask2)
    
    version = bedmachine_filename.split('/')[-1][-6:-3]
    
    return x,y,bed,thickness, version


def read_elevation(resolution,year):
    print('Reading elevation...')
    print('\t -using a %0.1f km resolution' % (resolution*0.5))
    bedmachine_filename = '../../../DATA/ELEVATION/elevation_lucille_'+str(year)+'.nc'
    
    x, y, surf = netcdf.readgrid(bedmachine_filename, 'surf')
        
    dx_native = np.abs(x[1] - x[0]) 
    step = int(round(resolution / dx_native))
    if step < 1:
        step = 1

    x_inf, x_sup = +0, -1
    y_inf, y_sup = +0, -1
    x = x[x_inf:x_sup:step]
    y = y[y_inf:y_sup:step]
    surf = surf[x_inf:x_sup:step, y_inf:y_sup:step]


    return x,y,surf



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


def new_floating_thickness(surf, bed, gmask, ri=0.917, rw = 1.028, zs_sl=0):
    thickness_ground = surf - bed
    thickness_float = surf/(1-ri/rw)
    grounded_mask = 1+((gmask-1)/2)
    floating_mask = (-((gmask-1)/2))
    thickness = thickness_ground*grounded_mask + thickness_float*floating_mask
    return thickness


def get_GL(x, y, gmask):
    print('Extracting grounding line...')
    x_gl,y_gl = [],[]
    n, m = len(x), len(y)
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
    return x_gl, y_gl

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
x, y, bed, thickness, version = read_bedmachine(resolution=resolution*2)
#x_IS, y_IS, dhdt_IS = Build_IceSat2_ContinuousMap(x_final=x, y_final=y)

for targetyear in years : 
    print(targetyear)
    x_s, y_s, surf = read_elevation(resolution=resolution*2, year=targetyear)
    
    X, Y = np.meshgrid(x, y)
    interp_func = RegularGridInterpolator((y_s, x_s), surf, bounds_error=False, fill_value=np.nan)
    points = np.array([Y.ravel(), X.ravel()]).T  # shape (N, 2)
    surf_targetyear = interp_func(points).reshape(Y.shape)

    #make directory

    create_directory(f'{targetyear}_State')
    
    #import the mask
    #mask_path = f'../../../DATA/Mask/mask_PIG_{targetyear}_mask_full_mask.tif'
    print('RES',resolution)
    x_mask_year, y_mask_year, mask_year = import_mask(f'../../../DATA/Mask/mask_PIG_{targetyear}_crop.tif',xres=resolution,yres=resolution)
    print("mask romain shape : "+str(mask_year.shape))
    

    #x_mask_year_max, y_mask_year_max = mask_year.shape
    #x_mask_year = range(1,x_mask_year_max)
    #y_mask_year = range(1,y_mask_year_max)

    #%% build dhdt trend and corresponding dh
    #the trend is for the period 2003-2019, we assume it is similar over our period of reconstruction
    # extent_IS = [np.min(x_IS), np.max(x_IS), np.min(y_IS), np.max(y_IS)] 

    extent_IS = [np.min(x_fin), np.max(x_fin), np.min(y_fin), np.max(y_fin)]

    
    
    base = surf_targetyear - thickness
    gmask = np.ones_like(base)
    gmask[base-bed>0] = -1

    # target year
    thickness_targetyear = new_floating_thickness(surf_targetyear, bed, gmask)
    gmask_targetyear = grounding(thickness_targetyear, bed)
    bed_targetyear = surf_targetyear - thickness_targetyear
    bedrock_targetyear = bed

    
    # Afficher le masque sous forme d'image
    plt.figure(figsize=(10, 8))
    plt.imshow(gmask_targetyear, cmap="Pastel1", interpolation="nearest")
    plt.colorbar(label="Valeurs du masque")
    plt.title("Apercu du masque")
    plt.xlabel("Longitude (pixels)")
    plt.ylabel("Latitude (pixels)")
    #plt.show()




    if subplot_mask :
        # Pine Island
        x_min, x_max = -1694616.000, -1506616.000
        y_min, y_max = -371234.000, -204734.000
    
        # Créer un masque pour les coordonnées x et y
        x_mask = (x >= x_min) & (x <= x_max)
        y_mask = (y >= y_min) & (y <= y_max)


        # Appliquer le masque aux variables surface
        x_mask_subset = x[x_mask]
        y_mask_subset = y[y_mask]
        diffmask_subset = difference_mask[np.ix_(y_mask, x_mask)]
        gmask_subset = gmask_targetyear[np.ix_(y_mask, x_mask)]
        newmask_subset = new_gmask[np.ix_(y_mask, x_mask)]


        # Création des sous-graphes pour afficher les trois masques
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
    
        print(f"x shape {x.shape}, y shape : {y.shape}, mask shape {mask_year.shape}")

        # Créer un masque pour les coordonnées x et y
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



    

    #%%
    thickness_targetyear[thickness_targetyear<10]=10.
    thickness_targetyear[thickness_targetyear<10]=10.

    #%% Get grounding lines for target year and 2020
    GL = False
    if GL:
        #x_gl_2020, y_gl_2020 = get_GL(x,y,gmask_2020) 
        x_gl_targetyear, y_gl_targetyear = get_GL(x,y,gmask_targetyear) 
        #np.savetxt('2020_State/GL_2020_Lucille.txt', [x_gl_2020, y_gl_2020])
        np.savetxt(f'{targetyear}_State/GL_{targetyear}_Lucille.txt', [x_gl_targetyear, y_gl_targetyear])
        #x_fl_2020, y_fl_2020 = get_FL(x,y,gmask_2020) 
        x_fl_targetyear, y_fl_targetyear = get_FL(x_mask_year,y_mask_year,mask_year)

        #np.savetxt('2020_State/FL_2020.txt', [x_fl_2020, y_fl_2020])
        np.savetxt(f'{targetyear}_State/FL_{targetyear}_Lucille.txt', [x_fl_targetyear, y_fl_targetyear])


        if MakeMesh :
            x_gl_targetyear_fin =[]
            y_gl_targetyear_fin =[]
            x_fl_targetyear_fin =[]
            y_fl_targetyear_fin =[]

            with open(f'{targetyear}_State/contours_GL_{targetyear}_Lucille.txt','w') as GL:
                o=0
                #x_gl_targetyear,y_gl_targetyear = GL_ordering(x_gl_targetyear,y_gl_targetyear)
                for i in range(len(x_gl_targetyear)):
                    if x_gl_targetyear[i] >= x_fin[0] and x_gl_targetyear[i] <= x_fin[-1] and y_gl_targetyear[i] >= y_fin[0] and y_gl_targetyear[i] <= y_fin[-1]:
                        o=o+1    
                        GL.write( str(x_gl_targetyear[i])+'\t'+str(y_gl_targetyear[i])+'\n')
                        x_gl_targetyear_fin.append(x_gl_targetyear[i])
                        y_gl_targetyear_fin.append(y_gl_targetyear[i])
                #x_gl_fin, y_gl_GL_fin = GL_ordering(x_gl_targetyear_fin, y_gl_targetyear_fin)
                #for i in range(len(x_gl_fin)):
                    #GL.write( str(x_gl_targetyear[i])+'\t'+str(y_gl_targetyear[i])+'\n')
            GL.close()

            print('Contours GL done, nombre de points : '+str(o))         
            with open(f'{targetyear}_State/contours_FL_{targetyear}_Lucille.txt', 'w') as FL:
                #print(len(x_fl_targetyear))
                o=0
                for i in range(len(x_fl_targetyear)):
                    #if x_fl_targetyear[i] >= x_fin[0] and x_fl_targetyear[i] <= x_fin[-1] and y_fl_targetyear[i] >= y_fin[0] and y_fl_targetyear[i] <= y_fin[-1]:
                    print('pt FL trouvÃ©')
                    o=o+1
                    FL.write( str(x_fl_targetyear[i])+'\t'+str(y_fl_targetyear[i])+'\n')
                    x_fl_targetyear_fin.append(x_fl_targetyear[i])
                    y_fl_targetyear_fin.append(y_fl_targetyear[i])
            FL.close()
            print('contours FL done, nombre de points : '+str(o))


    #else:
        #x_gl_2020, y_gl_2020 = np.loadtxt('2020_State/GL_2020.txt')
        #x_gl_targetyear, y_gl_targetyear = np.loadtxt(f'{targetyear}_State/GL_{targetyear}_Lucille.txt')
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
    print(f'y is {y.shape}')
    print(f'surf is {surf_targetyear.shape}')
    print(f'thick is {thickness_targetyear.shape}')
    print(f'bed is {bedrock_targetyear.shape}')



    dic = {'x':x, 'y':y, 'surface':surf_targetyear, 'bed':bedrock_targetyear, 'thickness':thickness_targetyear}
    netcdf.write('../../../DATA/ELEVATION/Elevation_{year}_Lucille_100m.nc' , dic)

