#from fpdf import FPDF
from matplotlib import pyplot as plt
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from math import *
import rasterio


#ds = xr.open_dataset("../DATA/VELOCITY/velocity_data_1991_2023_v2.nc")
#ds = xr.open_dataset("../DATA/VELOCITY/ASE_TimeSeries_1973-2018.nc")





years= [1991,1995,2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017]
years = [1995, 2005, 2017]
#years  = ds['time'].values

background_image = "../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse_crop.tif"  # Remplacez par le chemin de votre fichier

with rasterio.open(background_image) as src:
    background = src.read(1)
    extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)


#Vx = ds['VX'].isel(time=i)
#Vy = ds['VY'].isel(time=i)
#V = np.sqrt(Vx **2 + Vy **2)

xmin, xmax, ymin, ymax = -1694617,-1550000,-371235, -250000

# Configuration de la figure avec sous-graphiques (3 variables x nombre d'années)

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True, sharex=True, layout='constrained')
a=0

# Boucle pour chaque année
for i, year in enumerate(years):
    # Import du masque
    print(year)
    ds = xr.open_dataset(f"../DATA/VELOCITY/ASE_TimeSeries_{year}-new.nc")
    if year < 1992 :
        print('rentre boucle')
        with rasterio.open('../DATA/Mask/mask_PIG_jer_1992.tif') as src:
            mask = src.read(1)  # Lecture de la première bande du masque
            mask_transform = src.transform  # Obtient la transformation du masque
            mask_bounds = src.bounds

    if year > 2017 :
        with rasterio.open('../DATA/Mask/mask_full_PIG_'+str(year)+'.tif') as src:
            mask = src.read(1)  # Lecture de la première bande du masque
            mask_transform = src.transform  # Obtient la transformation du masque
            mask_bounds = src.bounds

    else :

        with rasterio.open('../DATA/Mask/mask_PIG_jer_'+str(year)+'.tif') as src:
            mask = src.read(1)  # Lecture de la première bande du masque
            mask_transform = src.transform  # Obtient la transformation du masque
            mask_bounds = src.bounds


    #b=a%4
    #c=floor(i/4)
    ax = axes[i]
    #a=a+1

    # Extraction de la couche temporelle (z) correspondant à l'année
    data = ds['v']
    #Vx = ds['VX'].isel(YEAR1=year)
    #Vy = ds['VY'].isel(YEAR1=year)
    #data = np.sqrt(Vx **2 + Vy **2)


    subset = data.where(
        (data.x >= xmin) & (data.x <= xmax) &
        (data.y >= ymin) & (data.y <= ymax),
        drop=True  # pour ne garder que les valeurs dans l'extent
    )

    alpha_values = np.where(subset < 10e7, 1, 0.0)
    ax.imshow(background, extent=extent,cmap="gray", aspect='auto')
    #ax.imshow(mask, cmap='gray', alpha=0, extent=[mask_bounds.left, mask_bounds.right, mask_bounds.bottom, mask_bounds.top], zorder=10)
    # Affichage de la carte pour la variable sélectio
    im = ax.imshow(subset.values, extent=extent, vmin=0, vmax=3000 ,  cmap='bwr',alpha=alpha_values)

    # Superposer le masque en utilisant imshow
    


    # Définir les titres
    ax.set_title(f'Velocity in {year}', fontsize=20)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

#fig.colorbar(im,ax=axes,orientation='vertical',fraction=0.02,pad=0.04)

# Ajustement de l'espace entre les subplots
fig.colorbar(im, ax=axes[ :], location='right',shrink =0.7, pad=0.02, label="m.y-1")
#plt.show()
plt.savefig("Timeseries_velocity_small")
