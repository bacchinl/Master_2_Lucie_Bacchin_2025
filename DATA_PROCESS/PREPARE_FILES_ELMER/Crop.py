import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import h5py
from scipy.interpolate import griddata


x_min, x_max = -2172277, -1233381
y_min, y_max = -763723, 18030



def filter_netcdf(input_nc, output_nc, x_min, x_max, y_min, y_max):
    ds = xr.open_dataset(input_nc)

    # Filtrer les coordonnées x et y
    filtered_ds = ds.where((ds['x'] >= x_min) & (ds['x'] <= x_max) &
                           (ds['y'] >= y_min) & (ds['y'] <= y_max), drop=True)

    # Sauvegarder le NetCDF modifié
    filtered_ds.to_netcdf(output_nc)
    print(f"NetCDF filtré enregistré sous {output_nc}")

# Fonction pour filtrer un HDF5
def filter_hdf5(input_h5, output_h5, x_min, x_max, y_min, y_max):
    with h5py.File(input_h5, 'r') as f_in:
        # Lire les coordonnées x et y
        x_data = f_in['x'][:]
        y_data = f_in['y'][:]

        # Trouver les indices correspondant aux limites de x et y
        x_indices = np.where((x_data >= x_min) & (x_data <= x_max))[0]
        y_indices = np.where((y_data >= y_min) & (y_data <= y_max))[0]

        # Filtrer les données pour chaque dataset dans le fichier HDF5
        with h5py.File(output_h5, 'w') as f_out:
            for name in f_in:
                if name in ['x', 'y']:  # Filtrer les coordonnées
                    if name == 'x':
                        filtered_data = x_data[x_indices]
                    else:
                        filtered_data = y_data[y_indices]
                    f_out.create_dataset(name, data=filtered_data)
                else:  # Filtrer les autres données selon les indices
                    data = f_in[name][:]
                    if data.shape == (len(x_data), len(y_data)):
                        filtered_data = data[np.ix_(x_indices, y_indices)]
                        f_out.create_dataset(name, data=filtered_data)
                    else:
                        f_out.create_dataset(name, data=data)
            print(f"HDF5 filtré enregistré sous {output_h5}")

# Filtrer le NetCDF
filter_netcdf('../DATA/BedMachineAntarctica-v03.nc', 'BedMachineAntartica_crop.nc', x_min, x_max, y_min, y_max)

# Filtrer le premier fichier HDF5
filter_hdf5('../DATA/AIS_mass_change.h5', 'AIS_mass_change_crop.h5', x_min, x_max, y_min, y_max)

# Filtrer le second fichier HDF5
filter_hdf5('../DATA/ICE1_ICE2_AnIS_dHdt_2003_2018_R209_05KM_FLOAT_MASS_F2.h5', 'ICE1_ICE2_crop.h5', x_min, x_max, y_min, y_max)

