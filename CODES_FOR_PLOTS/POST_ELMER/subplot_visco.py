# faire un subplot avec les viscosité pour les 3 années, la différence de vitesse en pourcentage

import pyvista as pv
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
from scipy.interpolate import griddata
import rasterio

#year2 - year1
years = [1995,2000,2005,2011, 2014,2017]
#years=[2007, 2009,2012, 2016,2018, 2020, 2022]
#years=[1995,2005, 2017]
#years=[2007, 2016, 2022]
n = len(years)
#A = 2.4*(10**(-24)) # 0°C
A = 3.5*(10.**(-25.)) #-10°C
A = A*(10.**(18.))*(31536.*(10.**3.)) # convertir dans les unitées de sortie elmer
i_lambda=1

# charger le netcdf (thickness)
#netcdf_file =f"../../../../../DATA/BEDROCK/BedMachine_IS_{years[0]}_v03_paolo.nc" 
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
        mesh = pv.read(f"../Simu_{year}/Mesh_{year}_paolo_2/OPTIM_OPT_{year}-{i_lambda}_t0002.pvtu")
        meshes.append(mesh)
    except Exception as e:
        print(f"Erreur lors du chargement du fichier pour {year}: {e}")
        meshes.append(None)



n= len(years)
YlOrBr=mpl.colormaps['YlOrBr'].resampled(n)
colors_GL=YlOrBr(range(n))
print(colors_GL)
fig, axes=plt.subplots(1, n, figsize=(5*n, 5), sharey=True, layout='constrained')
# Assurez-vous que les maillages sont de même taille et structure avant de comparer
f= open("visco_stats_Paolo.txt","w")
f.write("Year \t & Mean \t & \t Median \t & Standard deviation\t & Min \t & Max\\ \n")
all_visco=[]

for i, (mesh,year) in enumerate(zip(meshes,years)):

    print(year)
    #GL_year='../../../../../DATA_PROCESS/Python-Scripts-main/BuildBackAntarctica/'+str(year)+'_State/contours_GL_'+str(year)+'.txt'
    #GL = np.loadtxt(GL_year, delimiter='\t')
    #x_GL=GL[:,0]
    #y_GL=GL[:,1]

    for j in range(1):  # Deux rangées : viscosité et erreur
        ax = axes[i]
        ax.imshow(background, extent=extent,cmap="gray", aspect='auto')  # Ajustez l'extent selon vos coordonnées
        
    
    if "alpha" in mesh.point_data:
        # Extraire les valeurs "alpha" des deux maillages après interpolation
        alpha = mesh.point_data["alpha"]
        
        eta = np.zeros_like(alpha)  # Tableau de différence initialisé à zéro
        mask = alpha!=0
        eta[mask] = alpha[mask]**2 # Différence aux endroits communs

        E = np.zeros_like(alpha) 
        D = np.zeros_like(alpha)  # Tableau de différence initialisé à zéro
        E[mask] = 2 * eta[mask] * A**(1/n)
        D[mask] = 1- (2*eta[mask]*(A **(1/n))) # Différence aux endroits communs
    
        # Ajouter les données de différence au maillage
        mesh.point_data["D"] = D
        
        # Affichage de la différence "alpha_difference"
        x = mesh.points[:, 0]  # Coordonnées X
        y = mesh.points[:, 1]  # Coordonnées Y

        grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        grid_z = griddata((x,y), eta, (grid_x, grid_y), method="linear")


        eta = eta[mask]
        q_min = np.percentile(eta,5)
        q_max = np.percentile(eta,95)
        eta_val = eta[(eta >= q_min) & (eta <= q_max) ]

        moyenne = np.mean(eta_val)
        medianne = np.median(eta_val)
        std = np.std(eta_val)
        Min = np.min(eta_val)
        Max= np.max(eta_val)
        all_visco.extend(eta_val.tolist())

        f.write(f"{year} \t & {moyenne:.2f} \t &  {medianne:.2f} \t & {std:.2f} \t & {Min:.2f} \t & {Max:.2f}\\ \n")

        # Ajouter les données de différence au maillage
        mesh.point_data["eta"] = eta

        # Affichage de la différence "alpha_difference"
        x = mesh.points[:, 0]  # Coordonnées X
        y = mesh.points[:, 1]  # Coordonnées Y

        grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        grid_z = griddata((x,y), eta, (grid_x, grid_y), method="linear")
     
        bu_gn_r = plt.cm.get_cmap("BuGn_r", 256)
        hot_r=plt.cm.get_cmap("hot_r", 256)
        
        colors_bu_gn = bu_gn_r(np.linspace(0, 1, 128))  # BuGn_r pour 0 à 0.2
        colors_hot = hot_r(np.linspace(0, 1, 1152))       # 1152 pour vmax 1, 3304 pour vmax=3 hot_r pour 0.2 à 1

        combined_colors = np.vstack((colors_bu_gn, colors_hot))
        custom_cmap = LinearSegmentedColormap.from_list("BuGn_Hot", combined_colors)

        # Créer un scatter plot des valeurs (ou utiliser un griddata si une interpolation est souhaitée)
        ax=axes[i]
        scatter = ax.scatter(x, y, c=eta, cmap=custom_cmap,vmin=0, vmax=3.5, marker='.', s=0.5)
        #scatter = ax.scatter(x, y, c=eta, cmap=custom_cmap,vmin=0, vmax=1, marker='.', s=0.5)
        ax.set_title(f'Viscosity in {year}')
        ax.set_xlabel('X')
        if i == 0:  # Ajouter un label Y uniquement pour le premier subplot
            ax.set_ylabel('Y')
        #fig.colorbar(scatter, ax=ax, orientation="vertical", label="Viscosity (eta)")

    #for j in range(2):  # Deux rangées : viscosité et erreur
        #for h in range(i, n):
            #ax = axes[j, h]
            #ax.scatter(x_GL, y_GL, color=colors_GL[i], marker='o', label='GL', s=0.5)  # Tracé des contours


all_visco=np.array(all_visco)
moy=np.mean(all_visco)
med=np.median(all_visco)
std=np.std(all_visco)
Min=np.min(all_visco)
Max=np.max(all_visco)
f.write(f"All years \t & {moy:.2f} \t &  {med:.2f} \t & {std:.2f} \t & {Min:.2f} \t & {Max:.2f}\\ \n")
f.close

fig.colorbar(scatter, ax=axes[:], location='right',pad=0.02, label="Viscosity (eta)")
plt.savefig(f"./Subplot_visco-.png", dpi = 100)
#plt.show()
