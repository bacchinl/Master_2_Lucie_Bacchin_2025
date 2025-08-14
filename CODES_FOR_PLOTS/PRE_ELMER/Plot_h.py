import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import rasterio

year_init = 1995
years = [1995, 2000, 2005, 2011, 2014, 2017]

fig, axes = plt.subplots(2, len(years), figsize=(len(years)*5,10), sharey=True)


for i, year in enumerate(years) :

    # Charger le fichier NetCDF
    file_path_init = "./Python-Scripts-main/BuildBackAntarctica/"+str(year_init)+"_State/BedMachine_IS_"+str(year_init)+"_v03_paolo.nc"  # Remplacez par le chemin de votre fichier NetCDF
    file_path_paolo = "./Python-Scripts-main/BuildBackAntarctica/"+str(year)+"_State/BedMachine_IS_"+str(year)+"_v03_paolo.nc"
    file_path_BM = "./Python-Scripts-main/BuildBackAntarctica/"+str(year)+"_State/BedMachine_IS_"+str(year)+"_v03.nc"

    paolo_file_path = f"../DATA/BEDROCK/thickness_{year}_warp_ps_crop.tif"
    paolo_file_1995 = "../DATA/BEDROCK/thickness_1995_warp_ps_crop.tif"
    paolo_file_gd = f"../DATA/BEDROCK/ANT_G1920V01_GroundedIceHeight_mean.nc"


    ds_init = xr.open_dataset(file_path_init)
    ds = xr.open_dataset(file_path_paolo)
    ds_BM = xr.open_dataset(file_path_BM)
    ds_gd = xr.open_dataset(paolo_file_gd )
    ds_gd_year = ds_gd.sel(time=str(year))
    ds_gd_1995 = ds_gd.sel(time=str(year_init))

    with rasterio.open(paolo_file_path) as src1, rasterio.open(paolo_file_1995) as src2:
    # Lire les données en tant qu'array
        paolo_year = src1.read(1)  # Lire la première bande
        paolo_1995 = src2.read(1)


    #print(paolo_ds)

    GL_file = "Python-Scripts-main/BuildBackAntarctica/"+str(year)+"_State/contours_GL_"+str(year)+".txt" 
    FL_file ="Python-Scripts-main/BuildBackAntarctica/"+str(year)+"_State/contours_FL_"+str(year)+".txt"

    GL= np.loadtxt(GL_file, delimiter='\t')
    FL= np.loadtxt(FL_file, delimiter='\t')

    # Séparer les colonnes en coordonnées x et y
    x_GL = GL[:, 0]
    y_GL = GL[:, 1]

    x_FL = FL[:, 0]
    y_FL = FL[:, 1]

    #print(ds.data_vars)

    # Définir les limites du domaine
    xmin = -1694616.000
    xmax = -1506616.000
    ymin = -371234.000
    ymax =-204734.000

    # Sélectionner l'année 1995 dans la variable YEAR1
    #year_index = ds['YEAR1'].values == year

    # Vérifier si l'année 1995 est présente dans YEAR1
    Cond = True
    if Cond:
        # Sélectionner les coordonnées dans le domaine défini
        x_mask = (ds['x'] >= xmin) & (ds['x'] <= xmax)
        y_mask = (ds['y'] >= ymin) & (ds['y'] <= ymax)
    
        # Extraire VX et VY pour l'année 199 le domaine spatial défini
        h_p= ds['thickness'].where(x_mask &  y_mask, drop = True)
        h_BM= ds_BM['thickness'].where(x_mask &  y_mask, drop = True)

        h_init= ds_init['thickness'].sel(x=x_mask, y=y_mask)
        dhdt_gd_year = np.squeeze(ds_gd_year['dh']).values #.sel(x=x_mask, y=y_mask) 
        dhdt_gd_year_init = np.squeeze(ds_gd_1995['dh']).values #.sel(x=x_mask, y=y_mask)   
        



        dhdt_gd = dhdt_gd_year - dhdt_gd_year_init
        #dif_h = h - h_init
        

        #paolo

        dif_h_paolo = paolo_year - paolo_1995


        ax = axes[0,i]
        ctr=ax.pcolormesh(ds_BM['x'], ds_BM['y'], h_BM, vmin=0, vmax=1000,cmap="Purples")
        fig.colorbar(ctr, ax=ax, label="Epaisseur (m)")
        ax.set_title(f"Epaisseur en {year}, utilisée dans elme, avec BedMachine")

        # Ajouter les points du fichier texte au plot
        #ax.scatter(x_GL, y_GL, color='black', marker='o', label='GL', s=1)
        ax.scatter(x_FL, y_FL, color='black', marker='o', label='FL', s=1)
        ax.legend()
    
        ax = axes[1,i]
        #ctr2=ax.imshow( dif_h_paolo,cmap="bwr")
        ctr = ax.pcolormesh( x_mask, y_mask, h_paolo,vmin=0, vmax=1000, cmap='Purples')
        #ctr = ax.imshow( dhdt_gd,vmin=-100, vmax=100, cmap='bwr')

        fig.colorbar(ctr, ax=ax, label="Epaisseur (m)")
        ax.set_title(f"Epaisseur en {year}, utilisée dans elme, avec paolo datase ")
        #ax.set_title(f"variation d'epaisseur en {year}, pure")

        # Ajouter les points du fichier texte au plot
        #ax.scatter(x_GL, y_GL, color='black', marker='o', label='GL', s=1)
        ax.scatter(x_FL, y_FL, color='black', marker='o', label='FL', s=1)
        ax.legend()


    # Fermer le Dataset
    ds_init.close()
    ds.close()

#cbar =  fig.colorbar(ctr, ax=axes[0,i], orientation='vertical', fraction=0.15, pad=0.04)
#cbar =  fig.colorbar(ctr2, ax=axes[1,i], orientation='vertical', fraction=0.15, pad=0.04)

#cbar.set_label("Variation d'epaisseur (m)")


plt.tight_layout()
plt.savefig("comparaison-epaisseur-paolo.png")
plt.show()
