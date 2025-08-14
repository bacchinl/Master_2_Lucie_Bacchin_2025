import xarray as xr
import numpy as np
import os

# Charger le fichier NetCDF
file_path = "../DATA/VELOCITY/velocity_data_1991_2023_v2.nc"  # Chemin du fichier NetCDF source
ds = xr.open_dataset(file_path)

# Créer un répertoire pour stocker les fichiers par année
output_dir = "../DATA/VELOCITY"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#years = ds['YEAR1'].values
#years = [ 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016, 2017]
years = [2020, 2022]
print(years)

# Boucle sur chaque année disponible dans la dimension 'time'
for year in years:

    index_year = ds['time'].values == year

    # Sélectionner les variables 'vx' et 'vy' pour l'année spécifique
    vx = ds['vx'].isel(time=index_year).squeeze().data
    vy = ds['vy'].isel(time=index_year).squeeze().data
    vx = np.nan_to_num(vx, nan=-1e08)
    vy = np.nan_to_num(vy, nan=-1e08)
    # Calculer la vitesse totale v = sqrt(vx^2 + vy^2)
    #v = np.sqrt(vx**2 + vy**2)
    v= ds['v'].isel(time=index_year).squeeze().data
    v = np.nan_to_num(v, nan=-1e08)
    
    # Créer un nouveau Dataset pour cette année
    ds_year = xr.Dataset(
        {
            "vx": (("y", "x"), vx,  {"_FillValue": -1e08}),
            "vy": (("y", "x"), vy,  {"_FillValue": -1e08}),
            "v": (("y", "x"), v,  {"_FillValue": -1e08})  # Ajouter la vitesse totale
        },
        coords={
            "x": ds['x'],
            "y": ds['y'],
            "time": year
        }
    )
    # Nom du fichier de sortie pour l'année
    output_file = os.path.join(output_dir, f"velocity_lucille_{int(year)}.nc")
    
    # Sauvegarder dans un nouveau fichier NetCDF
    ds_year.to_netcdf(output_file)
    print(f"Fichier sauvegardé pour l'année {int(year)} : {output_file}")

# Fermer le Dataset original
ds.close()
