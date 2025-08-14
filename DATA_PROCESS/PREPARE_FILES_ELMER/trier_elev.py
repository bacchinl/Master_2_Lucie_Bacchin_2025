import xarray as xr
import numpy as np
import os

# Charger le fichier NetCDF
file_path = "../DATA/ELEVATION/elevation_data_2022_2023_50m.nc"  # Chemin du fichier NetCDF source
ds = xr.open_dataset(file_path)

# Créer un répertoire pour stocker les fichiers par année
output_dir = "../DATA/ELEVATION"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

#years = ds['YEAR1'].values
#years = [ 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2016, 2017]
years = [2017]
print(years)

# Boucle sur chaque année disponible dans la dimension 'time'
for year in years:

    index_year = ds['time'].values == year

    # Sélectionner les variables 'vx' et 'vy' pour l'année spécifique
    surf = ds['h'].isel(time=index_year).squeeze().data
    print("got h") 
    
    surf = np.nan_to_num(surf, nan=0)
    
    # Créer un nouveau Dataset pour cette année
    ds_year = xr.Dataset(
        {
            "surf": (("y", "x"), surf,  {"_FillValue": 0}),
    
        },
        coords={
            "x": ds['x'],
            "y": ds['y'],
            "time": year
        }
    )
    # Nom du fichier de sortie pour l'année
    output_file = os.path.join(output_dir, f"elevation_lucille_{int(year)}.nc")
    
    # Sauvegarder dans un nouveau fichier NetCDF
    ds_year.to_netcdf(output_file)
    print(f"Fichier sauvegardé pour l'année {int(year)} : {output_file}")

# Fermer le Dataset original
ds.close()
