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
years=[2007]
#years = [2007,2009, 2012, 2016,2018, 2020, 2022]

n = len(years)
#A = 2.4*(10**(-24)) # 0°C
A = 3.5*(10.**(-25.)) #-10°C
A = A*(10.**(18.))*(31536.*(10.**3.)) # convertir dans les unitées de sortie elmer
i_lambda=1
rho_w=1028 #kg.m-3
rho_i=917 #kg.m-3
g =9.81 #m.s-2
year_to_sec = 3.15576e+7

# charger le netcdf (thickness)
#netcdf_file =f"../../../../../DATA/BEDROCK/BedMachine_IS_{years[0]}_v03_paolo.nc" 
#ds = xr.open_dataset(netcdf_file)

# charger le fond
background_image = "../../../../../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse_crop.tif"  # Remplacez par le chemin de votre fichier

with rasterio.open(background_image) as src:
    background = src.read(1)
    extent = (src.bounds.left, src.bounds.right,src.bounds.bottom,src.bounds.top)

#define color map as the one if Furst et al 
colors = [
    "#39A0CD","#54B6E0","#6BC9F2","#85D5F0", "#9BD5C5", "#B1D295", "#CACF67","#E1CE34","#F7CD19", "#E9BD1F", "#D0AA2C","#B79530"
]
col_butt = LinearSegmentedColormap.from_list("col_butt", colors, N=12)


# Charger les fichiers VTU
meshes = []
for year in years:
    try:
        mesh = pv.read(f"../Simu_{year}/STRESS_OPT_{year}_L_t0001.pvtu")
        meshes.append(mesh)
    except Exception as e:
        print(f"Erreur lors du chargement du fichier pour {year}: {e}")
        meshes.append(None)



n= len(years)


fig, ax=plt.subplots(1, n, figsize=(5*n, 4), sharey=True)
# Assurez-vous que les maillages sont de même taille et structure avant de comparer
for i, (mesh,year) in enumerate(zip(meshes,years)):

    print(year)
    #GL_year='../../../../Mesh_Generation/Contours/'+str(year)+'_State/contours_1_'+str(year)+'_paolo_Clean.txt'
    #GL = np.loadtxt(GL_year, delimiter=' ')
    #x_GL=GL[:,0]
    #y_GL=GL[:,1]
    x_1, y_1=-1582400,-277000
    x_2,y_2=-1608500,-271000
    n_x, n_y =-(y_2-y_1), x_2-x_1
    len_n = np.sqrt(n_x**2+n_y**2)
    n_x, n_y = n_x/len_n, n_y/len_n
    vec_n=np.array([n_x,n_y])
    print("vec_n", vec_n)
    
    x=mesh.points[:,0]
    y=mesh.points[:,1]

    
    #ax = axes[i]
    ax.imshow(background, extent=extent,cmap="gray", aspect='auto')  # Ajustez l'extent selon vos coordonnées
        
    
    if ("stress 1" in mesh.point_data) & ("stress 2" in mesh.point_data) & ("stress 4" in mesh.point_data) & ("h" in mesh.point_data) & ("uobs 1" in mesh.point_data)& ("uobs 2" in mesh.point_data):
        # Extraire les valeurs "alpha" des deux maillages après interpolation
        sigma_xx = mesh.point_data["stress 1"]
        sigma_yy = mesh.point_data["stress 2"]
        sigma_xy = mesh.point_data["stress 4"]
        h = mesh.point_data["h"]
        u1 = mesh.point_data["uobs 1"]
        u2 = mesh.point_data["uobs 2"]
        gdmask = mesh.point_data["groundedmask"]



        print("shape sigma_xx", np.shape(sigma_xx))    

        mask = sigma_xx != 0

        if np.sum(mask) == 0:
            raise ValueError("Mask selects zero points — check your condition!")


        #initialisation of new var
        sigma_n = np.zeros_like(sigma_xx)
        sigma_eig = np.zeros_like(sigma_xx)
        sigma_flow = np.zeros_like(sigma_xx)
        No = np.zeros_like(h)
        Kn = np.zeros_like(sigma_xx)
        lim_grounded = np.zeros_like(gdmask)
        u1_norm, u2_norm = np.zeros_like(u1), np.zeros_like(u2)
        vec_flow = np.zeros((len(mask), 2))



        #p[mask] = rho_w * g * ( (rho_i/rho_w)* h[mask])
        lim_grounded_indices = np.where(gdmask == 0)[0]

        indices = np.where(mask)[0]

        

        for idx in np.where(mask)[0]:
            sigma = np.array([
                [sigma_xx[idx], sigma_xy[idx]],
                [sigma_xy[idx], sigma_yy[idx]]
            ])
           

            norm = np.sqrt(u1[idx]**2 + u2[idx]**2)

            if norm > 0:
                u1_norm[idx] = u1[idx] / norm
                u2_norm[idx] = u2[idx] / norm
            else:
                u1_norm[idx], u2_norm[idx] = 0, 0  # ou np.nan, selon ce que tu veux faire
            
            vec = np.array([u1_norm[idx], u2_norm[idx]])
            vec_flow[idx] = vec


            eigvals, eigvecs = np.linalg.eigh(sigma)
            vec_eig=eigvecs[:,1]
            len_vec_eig = np.sqrt(vec_eig[0]**2+vec_eig[1]**2)
            vec_eig = vec_eig/len_vec_eig
            

            sigma_n[idx] = (vec_n @ sigma @ vec_n)  *10**6#elmer MPa
            sigma_eig[idx]= (vec_eig @ sigma @ vec_eig) *10**6 
            sigma_flow[idx]= (vec_flow[idx] @ sigma @ vec_flow[idx]) *10**6

            #if abs(sigma_eig[idx]) < p[idx] :  
            #    sigma_eig[idx] = p[idx]
        
            #if abs(sigma_eig[idx]) > 0.5e11 :
            #    sigma_eig[idx] = 1e9

        #print("max sigma_n", np.max(sigma_n))
        #print("min sigma_n", np.min(sigma_n))
        #print("max sigma_xx", np.max(sigma_xx))
        #print("min sigma_xx", np.min(sigma_xx))
        #print("vec_n.T",np.shape(vec_n.T),"sigma_n",np.shape(sigma_n),"vec_n",np.shape(vec_n))

        



        #Nb[mask] = 1 - (abs(p[mask])/ (abs(sigma_eig[mask])-abs(p[mask])))
        #Nb[mask] = 1 - (abs(p[mask])/ abs(sigma_eig[mask]))

        # definition of Nb according to Furst 2016
        No[mask] = (1/2) * rho_i * (1- (rho_i/rho_w)) * g * h[mask]
        #Kn[mask] = 1- (sigma_eig[mask]/No[mask])
        Kn[mask] = 1-((No[mask]/sigma_flow[mask]))
        #Kn[mask] = 1-np.exp(-np.abs(np.log((No[mask]/sigma_flow[mask]))))

        #std = np.std(Kn[mask])
        #mean=np.mean(Kn[mask])
        #Kn[mask] = (Kn[mask]-mean )/std # standardisation de Kn

        print("max Kn", np.nanmax(Kn))
        print("min Kn", np.nanmin(Kn))

        print("max No", np.max(No))
        print("min No", np.min(No))
        print("max sigma_flow", np.max(sigma_flow))
        print("min sigma_flow", np.min(sigma_flow))
        



        # Ajouter les données de différence au maillage
        mesh.point_data["Kn"] = Kn
        mesh.point_data["No"] = No

        # Affichage de la différence "alpha_difference"
        x = mesh.points[:, 0]  # Coordonnées X
        y = mesh.points[:, 1]  # Coordonnées Y
        x_lim = x[lim_grounded_indices]
        y_lim = y[lim_grounded_indices]


        grid_x, grid_y = np.mgrid[min(x):max(x):500j, min(y):max(y):500j]
        grid_z = griddata((x,y), Kn, (grid_x, grid_y), method="linear")
         
    
        # Créer un scatter plot des valeurs (ou utiliser un griddata si une interpolation est souhaitée)
        #ax=axes[0,i]
        #scatter = ax.scatter(x, y, c=p, cmap="Purples", marker='.', s=1)
        scatter = ax.scatter(x, y, c=Kn, cmap=col_butt, vmin = 0, vmax=1, marker='.', s=1, alpha=1)
        #scatter = ax.scatter(x, y, c=sigma_flow,cmap='inferno_r' ,vmin=1e6 ,vmax=15e6, marker='.', s=1)
        ax.set_title(f'Kn in {year}')
        ax.set_xlabel('X')
        if i == 0:  # Ajouter un label Y uniquement pour le premier subplot
            ax.set_ylabel('Y')
        fig.colorbar(scatter, ax=ax, orientation="vertical", label="buttressing number")
        #scatter = ax.scatter(x_lim, y_lim, cmap="cyan", marker='.', s=1)


        
        
    
    
    #for j in range(2):  # Deux rangées : viscosité et erreur
        #for h in range(i, n):
            #ax = axes[j, h]
            #ax.scatter(x_GL, y_GL, color=colors_GL[i], marker='o', label='GL', s=0.5)  # Tracé des contours
   

plt.tight_layout()
plt.subplots_adjust(hspace=0.2,wspace=0.2)
plt.savefig(f"./Subplot_buttressing_eigen1.png", dpi = 600)
plt.show()





U1,V1=eigvecs[:,0]
U2,V2=eigvecs[:,1]
print('vec_flow ',np.shape(vec_flow))
u1_norm, u2_norm = vec_flow[:,0], vec_flow[:,1]

step=500
X_sub = x[::step]
Y_sub = y[::step]
#U1_sub = U1[:step]  # U = x-component of eigenvector
#V1_sub = V1[:step]
#U2_sub = U2[:step]  # U = x-component of eigenvector
#V2_sub = V2[:step]


x_sub=x[::step]
y_sub=y[::step]
#u1_sub = u1_norm[::step, :: step]
#u2_sub = u2_norm[::step, ::step]



plt.figure()
plt.imshow(background, extent=extent,cmap="gray", aspect='auto')
#plt.streamplot(x, y, u1_norm, u2_norm, density=0.5)
plt.quiver(x[::step], y[::step], u1_norm[::step], u2_norm[::step], color='red',scale = 50, width=0.002, label='flow direction')
#plt.quiver(X_sub, Y_sub, U1, V1, color='red', label='σ₁ direction')
#plt.quiver(X_sub, Y_sub, U2, V2, color='blue', label='σ₂ direction')
plt.legend()
plt.title("Principal Stress Directions")
plt.axis('equal')
plt.savefig(f"Plot_eigenvecs.png")
#plt.show()


plt.figure(figsize=(10, 6))


Kn_valid = Kn[(Kn >= -1) & (Kn <= 2)]
total_points = len(Kn)
valid_points = len(Kn_valid)
percent_valid = (valid_points / total_points) * 100

sigma_valid = sigma_flow[(sigma_flow >= 1e6) & (sigma_flow <= 15e6)]

plt.hist(Kn_valid, bins=100, alpha=0.6, label='Kn', color='blue', density=True)
#plt.hist(sigma_valid/(10**6), bins=100, alpha=0.6, label='Sigma, Mpa', color='red', density=True)
#plt.hist(No[mask], bins=100, alpha=0.6, label='No', color='green', density=True)


plt.axvline(np.mean(Kn), color='blue', linestyle='--', label=f'mean ⟨σₙ⟩ {np.mean(Kn)}')
#plt.axvline(np.mean(sigma_flow/(10**6)), color='red', linestyle='--', label=f'mean sigma_flow (MPa) {np.mean(sigma_flow/(10**6))}')
#plt.axvline(np.mean(No[mask]), color='Green', linestyle='--', label=f'mean sigma_flow (MPa) {np.mean(No[mask])}')
plt.text(0.02, 0.95, f'Points entre 0 et 1 : {percent_valid:.1f}%', transform=plt.gca().transAxes,
         fontsize=12, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, edgecolor='gray'))



plt.xlabel('Valeur Kn')
plt.ylabel('Densité normale')
plt.title("Histogramme de la répartition statistique de Nb")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"Plot_distrib_Nb.png")
plt.show()





