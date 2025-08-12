import pyvista as pv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.interpolate import griddata
import rasterio
print(matplotlib.__version__)


#year2 - year1
#years = [ 2007, 2009, 2012,  2016, 2018, 2020, 2022]
years = [1995, 2000, 2005, 2011, 2014, 2017]
#years = [1995, 2005, 2017]
#years=[2007, 2016, 2022]
n = 3
A = 2.4*(10**(-24))
Amax=2.4*(10**(-24)) # 0°Cunités Pa-3.s-1
Amin=3.5*(10**(-25)) # -10°

# convertir unitees elmer
Amax=Amax*10**(18)*31536*10**3
Amin=Amin*10**(18)*31536*10**3

i_lambda=1

# charger le fond

background_image_1990 ="../../../../../DATA/BACKGROUND_IMAGES/background_1990_crop.tif"
background_image_2000 ="../../../../../DATA/BACKGROUND_IMAGES/background_2000_crop.tif"
background_image_2010 ="../../../../../DATA/BACKGROUND_IMAGES/background_2010_crop.tif"
background_image_2015 = "../../../../../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse_crop.tif"  # Remplacez par le chemin de votre fichier





meshes =[]
for year in years :
    try:
        # Charger les fichiers VTU
        #mesh = pv.read(f"../Simu_{year}/OPTIM_OPT_{year}-{i_lambda}_L_t0002.pvtu") #, file_format="pvtu")
        mesh = pv.read(f"../Simu_{year}/Mesh_{year}_paolo_2/OPTIM_OPT_{year}-{i_lambda}_t0002.pvtu")
        meshes.append(mesh)
    except Exception as e:
        print(f"Erreur annee {year}")



fig,axes=plt.subplots(2,len(years),figsize=(5*(len(years)),10), sharey=True, layout='constrained')
# Assurez-vous que les maillages sont de même taille et structure avant de comparer
for i, (mesh,year) in enumerate(zip(meshes,years)):
    
    if year < 2000 : 
        with rasterio.open(background_image_1990) as src:
            background = src.read(1)
            extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)
    elif 2000 <= year <= 2015 :
        with rasterio.open(background_image_2000) as src:
            background = src.read(1)
            extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)
    elif 3005 < year <= 3015 :
        with rasterio.open(background_image_2010) as src:
            background = src.read(1)
            extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)
    elif 2015 < year  :
        with rasterio.open(background_image_2015) as src:
            background = src.read(1)
            extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)



    for j in range(2):  # Deux rangées : viscosité et erreur
        ax = axes[j, i]
        ax.imshow(background, extent=extent,cmap="gray", aspect='auto')  # Ajustez l'extent selon vos coordonnées


    if "alpha" in mesh.point_data and "mu" in mesh.point_data:
        # Extraire les valeurs "alpha" des deux maillages après interpolation
        alpha = mesh.point_data["alpha"]
        mu = mesh.point_data["mu"]
    
        # Calculer la différence uniquement là où les deux valeurs sont définies

        eta = np.zeros_like(alpha)  # Tableau de différence initialisé à zéro
        mask= alpha != 0
        eta[mask] = alpha[mask]**2 # Différence aux endroits communs
        print(type(alpha[mask]))    
            # Ajouter les données de différence au maillage
        mesh.point_data["eta"] = eta

        # Affichage de la différence "alpha_difference"
        x = mesh.points[:, 0]  # Coordonnées X
        y = mesh.points[:, 1]  # Coordonnées Y

        grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        grid_z = griddata((x,y), eta, (grid_x, grid_y), method="linear")




        # Calcul du enhancement factor

        Emin=np.zeros_like(alpha)
        Emax=np.zeros_like(alpha)
        Emin[mask]=1/((2*eta[mask])**n * Amax)
        Emax[mask]=1/((2*eta[mask])**n * Amin)
        print("Emin mean is "+str(np.mean(Emin)))
        print("Emawx mean is "+str(np.mean(Emax)))





        # Calcul du damage min, damage max

        Dmin=np.zeros_like(alpha)
        Dmax=np.zeros_like(alpha)
    
        Dmin[mask]=1-((1/Emin[mask])**(1/n))
        Dmax[mask]=1-((1/Emax[mask])**(1/n))
    
        
        print("Dmin mean is "+str(np.mean(Dmin[mask])))
        print("Dmawx mean is "+str(np.mean(Dmax[mask])))

        Dmin_selected_idx = np.where((Dmin <1) & (Dmin>0))
        Dmax_selected_idx = np.where((Dmax <1) & (Dmax>0))

        Dmin_selected = Dmin[Dmin_selected_idx]
        Dmax_selected = Dmax[Dmax_selected_idx]

        print("Dmin selected mean is "+str(np.mean(Dmin_selected)))
        print("Dmawx selected mean is "+str(np.mean(Dmax_selected)))


        alpha_values_Dmin = np.where(Dmin > 0.3, 1, 0.0)# alpha=1 si Dsuperieur à 0.3, invisible sinon
        
        ax=axes[0,i]
        scatter = ax.scatter(x, y, c=Dmin, cmap="Reds",vmin=0,vmax=1, s=0.15, alpha=alpha_values_Dmin)
        ax.set_title(f"Dmin - {year}")
        label ='X\n'
        label += f'Dmin Mean = {np.mean(Dmin_selected):.2f}'
        ax.set_xlabel(label)
        ax.set_ylabel("Y")
        #fig.colorbar(scatter, ax=ax, orientation="vertical", label="Dmin (units)")

        alpha_values_Dmax = np.where(Dmax > 0.3, 1, 0.0)


        ax=axes[1,i]
        scatter = ax.scatter(x, y, c=Dmax, cmap="Reds",vmin=0,vmax=1,s=0.15, alpha=alpha_values_Dmax)
        label ='X\n'
        label += f'Dmax Mean = {np.mean(Dmax_selected):.2f}'
        ax.set_title(f"Dmax - {year}")
        ax.set_xlabel(label)
        #fig.colorbar(scatter, ax=ax, orientation="vertical", label="Dmax (units)")


fig.colorbar(scatter, ax=axes[0,:], location='right',pad=0.02, label="Dmin (units)")
fig.colorbar(scatter, ax=axes[1,:], location='right',pad=0.02, label="Dmax (units)")

plt.savefig(f"Dmin-Dmax_Paolo_light.png", dpi=100)
#plt.show()





    # Calcul du damage

    #D = np.zeros_like(alpha)  # Tableau de différence initialisé à zéro
    #D[mask]= 1 - (2*(eta[mask]))*(A**(1/n))
    # Ajouter les données de différence au maillage
    #mesh.point_data["D"] = D
    #print("Le damage est compris entre "+str(np.min(D))+" et "+str(np.max(D)))

    # Enregistrer dans un nouveau fichier PVTU
    #mesh.save(f"Mesh_{year}/alpha_div_mu_{year}.vtu")
    #print(f"Le fichier '/alpha_div_mu_{year}.vtu' a été créé avec su dans Mesh_{year}.")    


    # Créer une figure avec Matplotlib
    #plt.figure(figsize=(8,6))
    
    # Créer un scatter plot des valeurs (ou utiliser un griddata si une interpolation est souhaitée)
    #sc = plt.scatter(x, y, c=D, cmap="Reds", marker='.', s=10)
    
    # Ajouter une barre de couleur
    #plt.colorbar(sc, label="D")
    
    # Ajouter des labels et un titre
    #plt.xlabel("longitude")
    #plt.ylabel("latitude")
    #plt.title("Damage (D)")
    
    # Afficher la figure
    #plt.show()
    #plt.savefig(f"Mesh_{year}/Damage_{year}.png")


#else:
    #print("Les données des points ne sont pas disponibles ou ne portent pas le même nom.")

