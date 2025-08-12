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
#years = [1995,2000,2005,2011, 2014,2017]
years=[2007, 2009,2012, 2016,2018, 2020, 2022]
#years=[1995,2005, 2017]
#years=[2007, 2016, 2022]
n = len(years)
#A = 2.4*(10**(-24)) # 0°C
A = 3.5*(10.**(-25.)) #-10°C
A = A*(10.**(18.))*(31536.*(10.**3.)) # convertir dans les unitées de sortie elmer
i_lambda=1


### Zone de définition des différentes viscosité
# Visco_GD, viscosité a la grounding line
y_lim_gd = -265800

# visco_SM, visco des deux shear margins
#pt_sud_east_SM 
x1_e, y1_e = -1597950, -331888
#pt_nord_east_SM 
x2_e, y2_e = -1578000, -271455

#pt_sud_west_SM 
x1_w, y1_w= (-1623740, -308246)

#pt_nord_west_SM 
x2_w, y2_w = (-1609760, -268751)


largeur_SM = 3000

##### center
#pt_sud_center
x1_c, y1_c= (x1_e+x1_w)/2, (y1_e+y1_w)/2 

#pt_nord_west_SM
x2_c, y2_c =(x2_e+x2_w)/2, (y2_e+y2_w)/2

largeur_center = abs(x2_c-x1_c)+ (largeur_SM/2) ## sécu


def distance_point_to_line(px, py, x1, y1, x2, y2):
    # Calcul de la distance d'un point (px, py) à la ligne (x1,y1)-(x2,y2)
    num = abs((y2 - y1)*px - (x2 - x1)*py + x2*y1 - y2*x1)
    den = np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    return num / den


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
        mesh = pv.read(f"../Simu_{year}/OPTIM_OPT_{year}-{i_lambda}_L_t0002.pvtu")
        #mesh = pv.read(f"../Simu_{year}/Mesh_{year}_paolo_2/OPTIM_OPT_{year}-{i_lambda}_t0002.pvtu")
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
f= open("visco_stats_Lucille.txt","w")
f.write("Year \t & Area & Mean \t & \t Median \t & Standard deviation\t & Min \t & Max\\ \n")
all_visco=[]

for i, (mesh,year) in enumerate(zip(meshes,years)):

    print(year)
    #GL_year='../../../../../DATA_PROCESS/Python-Scripts-main/BuildBackAntarctica/'+str(year)+'_State/contours_GL_'+str(year)+'.txt'
    #GL = np.loadtxt(GL_year, delimiter='\t')
    #x_GL=GL[:,0]
    #y_GL=GL[:,1]

    #for j in range(1):  # Deux rangées : viscosité et erreur
    #    ax = axes[i]
    #    ax.imshow(background, extent=extent,cmap="gray", aspect='auto')  # Ajustez l'extent selon vos coordonnées
        
    
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
        x_points = mesh.points[:, 0]  # Coordonnées X
        y_points = mesh.points[:, 1]  # Coordonnées Y

        #grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        #grid_z = griddata((x,y), eta, (grid_x, grid_y), method="linear")


        #### CALCUL DES VISCO ZONE #####

        dist_to_SM_est = distance_point_to_line(x_points, y_points, x1_e, y1_e, x2_e, y2_e)
        dist_to_SM_ouest = distance_point_to_line(x_points, y_points, x1_w, y1_w, x2_w, y2_w)

        # grounding line (GD)
        mask_GD = y_points > y_lim_gd

        # shear margins (SM), bande autour des deux lignes, distance < largeur_SM
        mask_SM_est = dist_to_SM_est < largeur_SM
        mask_SM_ouest = dist_to_SM_ouest < largeur_SM
        mask_SM = mask_SM_est | mask_SM_ouest

        # centre (entre les deux SM), sous la grounding line (y < y_lim_gd), hors SM
        dist_to_center = distance_point_to_line(x_points, y_points, x1_c, y1_c, x2_c, y2_c)
        mask_center = (~mask_SM) & (y_points < y_lim_gd) & (dist_to_center < largeur_center)

        # Est de la SM est
        mask_east_IS = (x_points > x1_e) & (~mask_SM_est) & (~mask_center) & (~mask_GD)

        # Ouest de la SM ouest
        mask_west_IS = (x_points < x1_w) & (~mask_SM_ouest) & (~mask_center)       

        # Extraire les viscosités dans chaque zone
        eta=eta[mask]
        visco_GD = eta[mask_GD]
        visco_SM = eta[mask_SM]
        visco_center = eta[mask_center]
        visco_east_IS = eta[mask_east_IS]
        visco_west_IS = eta[mask_west_IS]



        viscos = [eta, visco_GD, visco_SM, visco_center, visco_east_IS, visco_west_IS]
        zones = ["all", "GD", "SM", "Center", "East_IS", "West_IS"]
          
        for j, visco_zone in enumerate(viscos):
            if len(visco_zone) == 0:
                # Eviter erreur si zone vide
                print(visco_zone, " is empty for", year)
                continue

        
            q_min = np.percentile(visco_zone,5)
            q_max = np.percentile(visco_zone,95)
            eta_val = visco_zone[(visco_zone >= q_min) & (visco_zone <= q_max) ]

            moyenne = np.mean(eta_val)
            medianne = np.median(eta_val)
            std = np.std(eta_val)
            Min = np.min(eta_val)
            Max= np.max(eta_val)
            #print(j)
            if j == 0 :
                all_visco.extend(eta_val.tolist())

            f.write(f"{year} \t & {zones[j]} & {moyenne:.2f} \t &  {medianne:.2f} \t & {std:.2f} \t & {Min:.2f} \t & {Max:.2f}\\ \n")


        # Ajouter les données de différence au maillage
        #mesh.point_data["eta"] = eta

        # Affichage de la différence "alpha_difference"
        x = mesh.points[:, 0]  # Coordonnées X
        y = mesh.points[:, 1]  # Coordonnées Y

        


        #grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        #grid_z = griddata((x,y), eta, (grid_x, grid_y), method="linear")
     
        bu_gn_r = plt.cm.get_cmap("BuGn_r", 256)
        hot_r=plt.cm.get_cmap("hot_r", 256)
        
        colors_bu_gn = bu_gn_r(np.linspace(0, 1, 128))  # BuGn_r pour 0 à 0.2
        colors_hot = hot_r(np.linspace(0, 1, 1152))       # 1152 pour vmax 1, 3304 pour vmax=3 hot_r pour 0.2 à 1

        combined_colors = np.vstack((colors_bu_gn, colors_hot))
        custom_cmap = LinearSegmentedColormap.from_list("BuGn_Hot", combined_colors)

        # Créer un scatter plot des valeurs (ou utiliser un griddata si une interpolation est souhaitée)
        #ax=axes[i]
        #scatter = ax.scatter(x, y, c=eta, cmap=custom_cmap,vmin=0, vmax=3.5, marker='.', s=0.5)
        #scatter = ax.scatter(x, y, c=eta, cmap=custom_cmap,vmin=0, vmax=1, marker='.', s=0.5)
        #ax.set_title(f'Viscosity in {year}')
        #ax.set_xlabel('X')
        #if i == 0:  # Ajouter un label Y uniquement pour le premier subplot
        #    ax.set_ylabel('Y')
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


plt.figure(figsize=(10,8))

# Tous les points en gris clair en fond (optionnel)
plt.scatter(x, y, c='lightgrey', s=1, label='Other points')

# Points de chaque zone avec couleur différente
plt.scatter(x[mask_GD], y[mask_GD], c='red', s=5, label='Grounding Line')
plt.scatter(x[mask_SM], y[mask_SM], c='blue', s=5, label='Shear Margins')
plt.scatter(x[mask_center], y[mask_center], c='green', s=5, label='Center')
plt.scatter(x[mask_east_IS], y[mask_east_IS], c='orange', s=5, label='East IS')
plt.scatter(x[mask_west_IS], y[mask_west_IS], c='purple', s=5, label='West IS')

plt.xlabel('X')
plt.ylabel('Y')
plt.title('Visualization of Viscosity Zones')
plt.legend()
plt.axis('equal')
#plt.show()




#fig.colorbar(scatter, ax=axes[:], location='right',pad=0.02, label="Viscosity (eta)")
#plt.savefig(f"./Subplot_visco-.png", dpi = 100)
#plt.show()
