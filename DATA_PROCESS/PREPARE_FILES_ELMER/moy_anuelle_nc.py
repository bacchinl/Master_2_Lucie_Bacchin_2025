import xarray as xr

# Charger le fichier NetCDF
file_path = "../DATA/BEDROCK/ANT_G1920V01_IceShelfMelt_mean"
ds = xr.open_dataset(file_path+".nc")

# Vérifier les dimensions du fichier
print(ds['time'])

ds=ds.assign_coords(year=ds['time'].dt.year)
print("done")
print(ds['year'])


# Assurez-vous que 'time' est interprété comme une dimension temporelle
# Si 'time' est un objet datetime, c'est suffisant.
# Sinon, vous pouvez convertir explicitement :
if not ds['time'].dtype == 'datetime64[ns]':
    ds['time'] = xr.decode_cf(ds['time'])
print("conversion done")
# Ajouter une coordonnée "année" à partir de 'time'
ds = ds.assign_coords(year=ds['time'].dt.year)
print("ajout fait")


# Calculer la moyenne annuelle pour toutes les variables
annual_mean = ds.groupby("year").mean(dim="time")

# Sauvegarder le résultat dans un nouveau fichier NetCDF (optionnel)
output_path = str(file_path)+"_moy.nc"
annual_mean.to_netcdf(output_path)

print(f"Moyenne annuelle calculée et sauvegardée dans : {output_path}")
