import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import rasterio
from pyproj import Transformer
from matplotlib.patches import Rectangle
from rasterio import transform

years_nc = [1995, 2000, 2005, 2011, 2017] #icesat
years_tif = range(1991,2017,1 )#paolo
years_nc_L = [2002, 2006, 2008, 2009, 2011, 2012, 2015, 2016, 2018]#lucille

nb_zones=3

x_pt_1 = -1593500-2500
y_pt_1 = -274500 -8000
x_pt_2 = -1591000+2000 
y_pt_2 = -255000+1000
x_pt_3 =(-1658500-1638700)/2
y_pt_3 = -328500 + 20000


dist = 8000 



zone_1 = {
    "x_min": x_pt_1-dist,
    "x_max": x_pt_1+dist,
    "y_max": y_pt_1-dist,
    "y_min": y_pt_1+dist
}

zone_2 = {
    "x_min": x_pt_2-dist,
    "x_max": x_pt_2+dist,
    "y_min": y_pt_2-dist,
    "y_max": y_pt_2+dist}

zone_3 = {
    "x_min": x_pt_3-dist,
    "x_max": x_pt_3+dist,
    "y_min": y_pt_3-dist,
    "y_max": y_pt_3+dist}



zones=[zone_1, zone_3]


h_nc_l_1=[0]*len(years_nc_L)
h_nc_1=[0]*len(years_nc)
h_tif_1=[]

h_nc_l_2=[0]*len(years_nc_L)
h_nc_2=[0]*len(years_nc)
h_tif_2=[]

h_nc_l_3=[0]*len(years_nc_L)
h_nc_3=[0]*len(years_nc)
h_tif_3=[]

with rasterio.open("../DATA/BACKGROUND_IMAGES/rema_mosaic_crop_georef.tif") as src:
    image = src.read(1)
    #print("CRS",src.crs, "EPSG", src.crs.to_epsg())
    #print("Transform",src.transform)

    transform = src.transform

    # Coordonnées raster pour zone 1
    row_min_1, col_min_1 = ~transform * (zone_1["x_min"], zone_1["y_min"])
    row_max_1, col_max_1 = ~transform * (zone_1["x_max"], zone_1["y_max"])

    # Coordonnées raster pour zone 2
    row_min_2, col_min_2 = ~transform * (zone_2["x_min"], zone_2["y_min"])
    row_max_2, col_max_2 = ~transform * (zone_2["x_max"], zone_2["y_max"])
    #print(row_tif,col_tif)
    row_min_3, col_min_3 = ~transform * (zone_3["x_min"], zone_3["y_min"])
    row_max_3, col_max_3 = ~transform * (zone_3["x_max"], zone_3["y_max"])


for i, year in enumerate(years_tif) :

    paolo_file_path = f"../DATA/BEDROCK/thickness_{year}_warp_ps_crop.tif"
    

    with rasterio.open(paolo_file_path) as src:
    # Lire les données en tant qu'array
        paolo_year = src.read(1)  # Lire la première bande
        #Convertir les coordonnées géographiques en coordonnées raster
        #print(paolo_year)

        for j, zone in enumerate(zones) :
            row_start, col_start = src.index(zone["x_min"], zone["y_max"])
            row_stop, col_stop   = src.index(zone["x_max"], zone["y_min"])
            row_start, row_stop = sorted([row_start, row_stop])
            col_start, col_stop = sorted([col_start, col_stop])
            #row_min, col_min = ~transform * (zone["x_min"], zone["y_min"])
            #row_max, col_max = ~transform * (zone["x_max"], zone["y_max"]) 

            window = src.read(1, window=((row_start, row_stop), (col_start, col_stop)))

            valid_values = window[~np.isnan(window)]
            #print(f"[TIF] Zone {j+1}, year {year} — valid pixels: {np.count_nonzero(~np.isnan(window))}/{window.size}")
            mean_val = np.nanmean(valid_values) #if valid_values.size > 0 else np.nan
            #print(j)
            if j == 0 :
                h_tif_1.append(mean_val)
            if j == 1 :
                h_tif_2.append(mean_val)
            if j == 2 :
                h_tif_3.append(mean_val)


for i, year in enumerate(years_nc) :

    # Charger le fichier NetCDF
    file_path = f"../DATA/BEDROCK/BedMachine_IS_{year}_v03.nc"
    
    ds = xr.open_dataset(file_path)
    
    for j, zone in enumerate(zones) :
        
        x_min, x_max = sorted([zone["x_min"], zone["x_max"]])
        y_min, y_max = sorted([zone["y_min"], zone["y_max"]])
        
        # Corriger l'ordre des slices en fonction du sens de la grille
        x_increasing = (ds.x[1] - ds.x[0]) > 0
        y_increasing = (ds.y[1] - ds.y[0]) > 0

        x_slice = slice(x_min, x_max) if x_increasing else slice(x_max, x_min)
        y_slice = slice(y_min, y_max) if y_increasing else slice(y_max, y_min)

        zone_data = ds['thickness'].sel(x=x_slice, y=y_slice)
        values = zone_data.values
        valid = values[~np.isnan(values) & (values > 10)]

        if valid.size > 0:
            mean_val = np.nanmean(valid).item()  # assure conversion propre
        else:
            mean_val = np.nan


        
        #mean_val = float(zone_data.mean().values)
        #mean_val = float(np.nanmean(zone_data.values))

        #print(f"x: {x_min} to {x_max}")

        #print(f"y: {float(ds.y.min())} to {float(ds.y.max())}")

        #print(year,zone,mean_val)
        if j == 0 :
            h_nc_1[i] = mean_val
        if j == 1 :
            h_nc_2[i] = mean_val
        if j == 2 :
            h_nc_2[i] = mean_val



    #print(pt_nc)
    ds.close()


for i, year in enumerate(years_nc_L) :

    # Charger le fichier NetCDF
    lucille_file_path = f"../DATA/ELEVATION/Bedmachine_elevation_Lucille_{year}.nc"
    #lucille_file_path = f"../DATA/ELEVATION/elevation_lucille_{year}.nc"
    ds=xr.open_dataset(lucille_file_path)

    

    for j, zone in enumerate(zones) :
        x_min, x_max = sorted([zone["x_min"], zone["x_max"]])
        y_min, y_max = sorted([zone["y_min"], zone["y_max"]])

        # Corriger l'ordre des slices en fonction du sens de la grille
        x_increasing = (ds.x[1] - ds.x[0]) > 0
        y_increasing = (ds.y[1] - ds.y[0]) > 0

        x_slice = slice(x_min, x_max) if x_increasing else slice(x_max, x_min)
        y_slice = slice(y_min, y_max) if y_increasing else slice(y_max, y_min)

        zone_data = ds['thickness'].sel(x=x_slice, y=y_slice)

        values = zone_data.values
        valid = values[~np.isnan(values)]

        if valid.size > 0:
            mean_val = np.nanmean(valid).item()  # assure conversion propre
        else:
            mean_val = np.nan
    
        
        #print(f"Zone {j+1}, year {year} — shape: {zone_data.shape}")
        #print(f"  min: {np.nanmin(zone_data.values)}, max: {np.nanmax(zone_data.values)}")
        #print(f"  NaN count: {np.isnan(zone_data.values).sum()}, total points: {zone_data.size}")


        if j == 0 :
            
            h_nc_l_1[i] = mean_val
        if j == 1 :
            
            h_nc_l_2[i] = mean_val
        if j == 2 :
            
            h_nc_l_3[i] = mean_val
    


    #print(pt_nc)
    ds.close()




print("Icesat zone 1: ", h_tif_1 )
print("PAOLO zone 1: ", h_nc_1 )
print("LUCILLE zone 1: ", h_nc_l_1 )

print("Icesat zone 2: ", h_tif_2 )
print("PAOLO zone 2: ", h_nc_2 )
print("LUCILLE zone 2: ", h_nc_l_2 )

fig, axs = plt.subplots (1, 3, figsize=(30,10))
# --- Subplot 1: Time series plot ---
axs[0].plot(years_nc, h_nc_1, marker="o", label=f"LASER_ALT ", color="green")
axs[0].plot(years_tif, h_tif_1, marker="o", label=f"RADAR_ALT ", color="red")
axs[0].plot(years_nc_L, h_nc_l_1, marker="o", label=f"STEREO_ALT ", color="magenta")
axs[0].set_xlabel("Time")
axs[0].set_ylabel("Thickness")
axs[0].set_title(f"Temporal evolution of mean thickess in Zone 1 (Center : {x_pt_1},{y_pt_1}, Side : {dist*2}) ")
axs[0].legend()
axs[0].grid(True)


# --- Subplot 1: Time series plot ---
axs[1].plot(years_nc, h_nc_2, marker="o", label=f"LASER_ALT ", color="green")
axs[1].plot(years_tif, h_tif_2, marker="o", label=f"RADAR_ALT", color="red")
axs[1].plot(years_nc_L, h_nc_l_2, marker="o", label=f"STEREO_ALT", color="magenta")
axs[1].set_xlabel("Time")
axs[1].set_ylabel("Thickness")
axs[1].set_title(f"Temporal evolution of mean thickess in Zone 2 (Center : {x_pt_2},{y_pt_2}, Side : {dist*2})")
axs[1].legend()
axs[1].grid(True)





# --- Subplot 2: Image with point location ---
axs[2].imshow(image, cmap='gray')
rect1 = Rectangle((row_min_1, col_min_1), row_max_1 - row_min_1, col_max_1 - col_min_1,
                  linewidth=2, edgecolor='blue', facecolor='none', label="Zone 1")
#rect2 = Rectangle((row_min_2, col_min_2), row_max_2 - row_min_2, col_max_2 - col_min_2,
#                  linewidth=2, edgecolor='black', facecolor='none', label="Zone 2")
rect3 = Rectangle((row_min_3, col_min_3), row_max_3 - row_min_3, col_max_3 - col_min_3,
                  linewidth=2, edgecolor='orange', facecolor='none', label="Zone 3")



axs[2].add_patch(rect1)
#axs[2].add_patch(rect2)
axs[2].add_patch(rect3)
axs[2].axis('off')
axs[2].set_title(f"Location of the zone")
axs[2].legend()

# --- Save the combined figure ---
plt.tight_layout()
plt.savefig(f"subplot_ev_tempo_3.png")

