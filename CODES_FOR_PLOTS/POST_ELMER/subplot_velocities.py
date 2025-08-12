# faire un subplot avec les viscosité pour les 3 années, la différence de vitesse en pourcentage

import pyvista as pv
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from scipy.interpolate import griddata
import rasterio

#year2 - year1
#years = [2007, 2009, 2012, 2016, 2018, 2020, 2022]
#years = [1995, 2000, 2005, 2011, 2014, 2017]
#years = [2009,2012, 2016,2018]
years =[1995, 2005, 2017]
#years = [2007, 2016, 2022]
n = len(years)
#A = 2.4*(10**(-24)) # 0°C
A = 3.5*(10.**(-25.)) #-10°C
A = A*(10.**(18.))*(31536.*(10.**3.)) # convertir dans les unitées de sortie elmer
i_lambda=1

# charger le netcdf (thickness)
#netcdf_file =f"../../../../DATA/BEDROCK/BedMachine_IS_{years[0]}_v03.nc" 
#ds = xr.open_dataset(netcdf_file)

# charger le fond
background_image = "../../../../../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse_crop.tif"  # Remplacez par le chemin de votre fichier

with rasterio.open(background_image) as src:
    background = src.read(1)
    extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)


# Charger les fichiers VTU
meshes = []
for year in years:
    try:
        #mesh = pv.read(f"../Simu_{year}/OPTIM_OPT_{year}-{i_lambda}_L_t0002.pvtu")
        mesh = pv.read(f"../Simu_{year}/Mesh_{year}/OPTIM_OPT_{year}-{i_lambda}_t0002.pvtu")

        meshes.append(mesh)
    except Exception as e:
        print(f"Erreur lors du chargement du fichier pour {year}: {e}")
        meshes.append(None)

# Ouvre un fichier texte pour sauvegarder les statistiques globales

f = open("errors_Lucille_ds.txt", "w")
f.write("Year\t & Mean\t & Median\t & Standard deviation\t & RMSE\\ \n")  # En-tête


all_erreur = []



n= len(years)
fig, axes=plt.subplots(3, n, figsize=(5*n, 15), sharey=True, sharex=True, layout='constrained')
# Assurez-vous que les maillages sont de même taille et structure avant de comparer
for i, (mesh,year) in enumerate(zip(meshes,years)):

    for j in range(3):  # Deux rangées : viscosité et erreur
        ax = axes[j, i]
        ax.imshow(background, extent=extent,cmap="gray", aspect='auto')  # Ajustez l'extent selon vos coordonnées

    # Affichage de la différence "alpha_difference"
    x = mesh.points[:, 0]  # Coordonnées X
    y = mesh.points[:, 1]  # Coordonnées Y

    if "uobs" in mesh.point_data and "ssavelocity" in mesh.point_data:
        # Pourcentage d'erreur sur les vitesses
        uobs = mesh.point_data["uobs"]
        uobs_x = mesh.point_data["uobs 1"]
        uobs_y = mesh.point_data["uobs 2"]
        ssavelocity = mesh.point_data["ssavelocity"]
        umod_x =ssavelocity[:, 0]
        umod_y = ssavelocity[:, 1]
        umod = np.sqrt(umod_x**2+umod_y**2)
        mask = uobs > 0 # éviter de diviser par 0
    

        dif_v = np.zeros_like(uobs)
        dif_v[mask]= uobs[mask]-umod[mask]

        ######### STATISTICS STUFF #########

        
        erreur = dif_v[mask]
        q05 = np.percentile(erreur, 10)
        q95 = np.percentile(erreur, 90)
        erreur = erreur[(erreur >= q05) & (erreur <= q95)]
        erreur = erreur[np.abs(erreur) <= 10000]

        moyenne = np.mean(erreur)
        mediane = np.median(erreur)
        std = np.std(erreur)
        rmse = np.sqrt(np.mean(erreur**2))
        all_erreur.extend(erreur.tolist())
        
        f.write(f"{year} & {moyenne:.2f} & {mediane:.2f}&{std:.2f} & {rmse:.2f} \\  \n")


        #### PLOT STUFF  ######

        base_cmap=plt.cm.get_cmap("bwr")
        colors = [
                *base_cmap(np.linspace(0,1, 256)),
                (1,1,0),]
        custom_cmap_bwr = LinearSegmentedColormap.from_list("bwr_custom", colors)

        bounds = [-np.inf, -400000,-200000,0,200000, 400000 , np.inf]
        norm=BoundaryNorm(bounds, custom_cmap_bwr.N)



        ax=axes[0, i]  # Sélection de l'axe correspondant
        scatter = ax.scatter(x,y,c=uobs, cmap =custom_cmap_bwr,vmin=0,vmax=5000, marker='.', s=1, alpha=0.8)
        ax.set_title(f'Observed velocities for {year}')
        ax.set_xlabel('X')
        if i == len(years)-1:  # Ajouter un label Y uniquement pour le premier subplot
            ax.set_ylabel('Y')
            #fig.colorbar(scatter, ax=ax, orientation="vertical", label="m.y-1")
        #fig.colorbar(scatter, ax=axes[0, :], location='right', shrink=0.9, label="m.y-1")
        #fig.colorbar(scatter, ax=axes[0,:],location="right", orientation="vertical", label="m.y-1 " , shrink=0.8)



        ax=axes[1, i]  # Sélection de l'axe correspondant
        scatter = ax.scatter(x,y,c=umod, cmap ='bwr',vmin=0,vmax=5000, marker='.', s=1, alpha=0.8)
        ax.set_title(f'Modeled velocities for {year}')
        ax.set_xlabel('X')
        if i == len(years)-1:  # Ajouter un label Y uniquement pour le premier subplot
            ax.set_ylabel('Y')
            #fig.colorbar(scatter, ax=ax, orientation="vertical", label="m.y-1")
        #fig.colorbar(scatter, ax=axes[1, :], location='right', shrink=0.9, label="m.y-1")

        #ax=axes[2, i]  # Sélection de l'axe correspondant
        #scatter = ax.scatter(x,y,c=dif_v, cmap ='bwr',vmin=-2000,vmax=2000, marker='.', s=1, alpha=0.8)
        #x.set_title(f'Modeled velocities for {year}, with paolo dhdt')
        #ax.set_xlabel('X')
        #if i == 0:  # Ajouter un label Y uniquement pour le premier subplot
        #    ax.set_ylabel('Y')
        #fig.colorbar(scatter, ax=ax, orientation="vertical", label="m.y-1")


        
        pourc  = np.zeros_like(uobs)  # Tableau de différence initialisé à zéro
        pourc[mask]= 100 * ((umod[mask] - uobs[mask]) / uobs[mask])
    

        # Ajouter les données de différence au maillage
        mesh.point_data["por_erreur"] = pourc
        print("Le pourcentage d'erreur sur la vitesse est compris entre "+str(np.min(pourc))+" et "+str(np.max(pourc)))



        base_cmap=plt.cm.get_cmap("BrBG_r")
        colors = [
                (1,1,0),
                *base_cmap(np.linspace(0,1, 256)),
                (1,1,0)]
        custom_cmap_err = LinearSegmentedColormap.from_list("BrBG_r_custom", colors)

        bounds = [-np.inf, -10000,-50,0,50, 10000 , np.inf]
        norm=BoundaryNorm(bounds, custom_cmap_err.N)


        ax=axes[2, i]  # Sélection de l'axe correspondant
        scatter_err = ax.scatter(x,y,c=-dif_v, cmap =custom_cmap_err,vmin=-500, vmax=500, marker='.', s=1)
        label = 'X \n'
        label +=  f'Abs mean : {np.mean(abs(pourc)):.1f} %' 
        ax.set_title(f'Error for {year}')
        ax.set_xlabel(label)
        if i == len(years)-1:  # Ajouter un label Y uniquement pour le premier subplot
            ax.set_ylabel('Y')
            #fig.colorbar(scatter, ax=ax, orientation="vertical", label="Error (%)")
        #fig.colorbar(scatter, ax=axes[2, :], location='right', shrink=0.1, label="Error (%)")

all_erreur = np.array(all_erreur)
moyenne = np.mean(all_erreur)
mediane = np.median(all_erreur)
std = np.std(all_erreur)
rmse = np.sqrt(np.mean(all_erreur**2))


f.write(f"All year & {moyenne:.2f} & {mediane:.2f}&{std:.2f} & {rmse:.2f} \\ \n")
f.close()
fig.colorbar(scatter, ax=axes[0, :], location='right', pad=0.02, label="m.y-1")
fig.colorbar(scatter, ax=axes[1, :], location='right', pad=0.02, label="m.y-1")
fig.colorbar(scatter_err, ax=axes[2, :], location='right', pad=0.02, label="Error m.y-1")

#plt.tight_layout()
#plt.subplots_adjust(hspace=0.2,wspace=0.2)
plt.savefig(f"./Velocity_plots_dpi_300/Subplot_velocities_ICESat_red_light", dpi = 300)
#plt.savefig(f"./Subplot_velocities-paolo-lambda1e04.png", dpi = 600)

#plt.show()
