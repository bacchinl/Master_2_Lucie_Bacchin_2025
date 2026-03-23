""" This code works with the plot_buttressing functions """


from turtle import speed

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import numpy as np
import pyvista as pv
import matplotlib.tri as mtri
from plot_buttressing import (
    plot_buttressing,
    plot_corrected,
    plot_stress,
    plot_strain,
    get_furst_colormap,
    test_var
)

small = 1e-10
year = 1995

# print ressource 
import resource
print("Max memory:", resource.getrlimit(resource.RLIMIT_AS))


calcul_eig_furst_way = False

beta = '1e-06'
## ------------------------------------------------- IMPORT DATA  ------------------------------------------

filename = f"../Pine_Island2/Mesh_{year}/OPTIM_{year}_OPT_{year}-1_t0002.pvtu"


## ------------------------------------------------- BUTTRESSING CALCULATION ------------------------------------------


def compute_buttressing_unstructured(filename, nn=3, rhoice=917, rhow=1028, grav=9.81, thiampl=1):
    """
    Calcule KT, KN1, KN2, KN3 directement sur un maillage VTU/PVTU triangulaire.
    AUCUNE reconstruction de grille → très faible RAM.
    """

    year_in_sec = 365.25 * 24 * 3600

    # --- lecture du maillage ---
    mesh = pv.read(filename)

    #print variables disponibles dans le maillage
    # print(mesh.point_data.keys())
    # noms plausibles — à ajuster selon le VTU
    vx = mesh.point_data["ssavelocity"][:,0]    # vx m/a
    vy = mesh.point_data["ssavelocity"][:,1]    # vy m/a
    eta = mesh.point_data["alpha"]**2  # viscosité optim, in MPa.a^1/3
    thi = mesh.point_data["h"]         # thickness in m
 
    pts = mesh.points
    x = pts[:,0]
    y = pts[:,1]
    

    # ---------------------------------------------------------
    # 1) CALCUL DES GRADIENTS PAR VTK → hyper rapide et léger
    # ---------------------------------------------------------
    
    # gradient of a velocity component is a 3D array with shape (n_points, 9) where the 9 columns correspond to the components of the gradient tensor:
    grad_v = mesh.compute_derivative(scalars="ssavelocity", gradient=True).point_data["gradient"]
    

    # grad_vx[:,0] = ∂vx/∂x, grad_vx[:,1] = ∂vx/∂y
    dudx = grad_v[:,0]
    dudy = grad_v[:,1]
    
    dvdx = grad_v[:,3]
    dvdy = grad_v[:,4]

    print(np.shape(dudx), np.shape(dudy), np.shape(dvdx), np.shape(dvdy))
    # ---------------------------------------------------------
    
    # 2) TENSEUR DE DÉFORMATION (epsilon)
    # ---------------------------------------------------------
    exx = dudx
    eyy = dvdy
    exy = 0.5 * (dudy + dvdx)
    ezz = -exx - eyy

    small = 1e-10
    eff = 0.5*(exx**2 + eyy**2 + ezz**2) + exy**2 + small

    # ---------------------------------------------------------
    # 3) DEVIATORIC STRESS (glen law)
    # ---------------------------------------------------------

    # factor = rf**(-1/nn) * eff**((1-nn)/(2*nn))
    eta_effective = eta * eff**((1-nn)/(2*nn))

    #ttxx = exx * factor
    #ttyy = eyy * factor
    #ttxy = exy * factor
    ttxx = 2 * eta_effective * exx
    ttyy = 2 * eta_effective * eyy
    ttxy = 2 * eta_effective * exy
    
    # ---------------------------------------------------------
    # 4) TENSEUR T MODIFIÉ (stress tensor in SSA) 
    # ---------------------------------------------------------
    t11 = 2 * ttxx + ttyy
    t12 = ttxy
    t21 = ttxy
    t22 = ttxx + 2* ttyy


    # ---------------------------------------------------------
    # 5) EIGENVALUES
    # ---------------------------------------------------------
    
    bb = -(t11 + t22)
    cc = t11*t22 - t21*t12

    disc = np.sqrt(bb**2 - 4*cc)

    eig1 = (-bb + disc) / 2 #in MPa
    eig2 = (-bb - disc) / 2 #in MPa

    
    # ---------------------------------------------------------
    # 5bis) EIGENVECTOR associated with eig2
    # ---------------------------------------------------------
    eig2_vx = t12
    eig2_vy = eig2 - t11

    # normalisation
    norm = np.sqrt(eig2_vx**2 + eig2_vy**2) + 1e-12

    if calcul_eig_furst_way:
        norm = np.sqrt(1**2 + eig2_vy**2) + 1e-12
        eig2 = np.sqrt(eig2_vx**2 + eig2_vy**2)
    
    eig2_vx /= norm
    eig2_vy /= norm

    # ---------------------------------------------------------
    # 6) DIRECTIONS / PROJECTIONS (KN1)
    # ---------------------------------------------------------
    speed = np.sqrt(vx**2 + vy**2)
    vnorm = np.maximum(speed, 1e-12)

    n1 = vx / vnorm
    n2 = vy / vnorm

    m1 = -n2
    m2 =  n1

    tau = n1*(t11*n1 + t12*n2) + n2*(t12*n1 + t22*n2)
    psi = m1*(t11*m1 + t12*m2) + m2*(t12*m1 + t22*m2)

    # ---------------------------------------------------------
    # 7) BUTTRESSING RATIOS
    # ---------------------------------------------------------
    tau0 = 0.5*rhoice*grav*thi*(1 - rhoice/rhow) 
    conversion = 1e-6 #convert to MPa
    tau0 = tau0 * conversion
    
    KT  = psi / tau0
    KN1 = 1 - (tau  / tau0)
    KN2 = 1 - (eig1 / tau0)
    KN3 = 1 - (eig2 / tau0)

    return x, y, vx, vy, t11, t12, t21, t22, exx, eyy, exy, KT, KN1, KN2, KN3, eig1, eig2, eig2_vx, eig2_vy


x, y, vx, vy, t11, t12, t21, t22, exx, eyy, exy, KT, KN1, KN2, KN3, eig1, eig2, eig2_vx, eig2_vy = compute_buttressing_unstructured(filename)

## tests 
def test_var(var, name):
    """
    Affiche les statistiques de base d'une variable :
      - min
      - max
      - mean
      - median
      - nombre de NaN
    """
    print(f"\n===== {name} =====")
    print(f"  min     : {np.nanmin(var)}")
    print(f"  max     : {np.nanmax(var)}")
    print(f"  mean    : {np.nanmean(var)}")
    print(f"  median  : {np.nanmedian(var)}")
    print(f"  NaNs    : {np.isnan(var).sum()}")
    print("====================")

test_var(KT,  "KT  (transverse)")
test_var(KN1, "KN1 (along-flow)")
test_var(KN2, "KN2 (eigenvalue 1)")
test_var(KN3, "KN3 (eigenvalue 2)")


#------------------------------
#   PLOT
#------------------------------
# palette comme dans furst et al 


cmap = get_furst_colormap()

#plot_buttressing(x, y, KT, KN1, KN2, KN3, color_fuerst, year)

fig, ax = plot_buttressing(
    x, y,
    KT, KN1, KN2, KN3,
    cmap,     
    year
)
plt.savefig(f"./Figures/buttressing_corrected_{year}_FRT_{beta}", dpi = 300)
plt.close()


fig, ax = plot_strain(x, y, exx, exy, exy, eyy, year, beta)
plt.savefig(f"./Figures/strain_{year}_FRT_{beta}", dpi = 300)
plt.close()

fig, ax = plot_stress(x, y, eig2_vx, eig2_vy, year, beta)
plt.title(f"Stress for year {year}, friction param of {beta}")
plt.savefig(f"./Figures/stress_{year}_FRT_{beta}", dpi = 300)
plt.close()

# fig, ax = plot_eigenvec(
#     x, y, eig2_vx, eig2_vy, year, beta
# )
# plt.savefig(f"./Figures/eigenvec_{year}_FRT_{beta}", dpi = 300)
# plt.show()

