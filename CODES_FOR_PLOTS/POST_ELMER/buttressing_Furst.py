""" This code works with the plot_buttressing functions """


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
    plot_stress_simple,
    plot_eigenvec,
    get_furst_colormap,
    test_var
)

small = 1e-10
year = 2016
year_to_sec = 3.15576e+7

# print ressource 
import resource
print("Max memory:", resource.getrlimit(resource.RLIMIT_AS))


calcul_eig_furst_way = False

beta = '1e-06'
## ------------------------------------------------- IMPORT DATA  ------------------------------------------

filename = f"../Simu_{year}/OPTIM_{year}_L_FRT_{beta}_t0002.pvtu"


## ------------------------------------------------- BUTTRESSING CALCULATION ------------------------------------------


def compute_buttressing_unstructured(filename, nn=3, rhoice=917, rhow=1028, grav=9.81, thiampl=1):
    """
    Calcule KT, KN1, KN2, KN3 directement sur un maillage VTU/PVTU triangulaire.
    AUCUNE reconstruction de grille → très faible RAM.
    """

    # --- lecture du maillage ---
    mesh = pv.read(filename)

    # noms plausibles — à ajuster selon le VTU
    vx = mesh.point_data["uobs 1"]    # vx m/a
    vy = mesh.point_data["uobs 2"]    # vy m/a
    eta = mesh.point_data["alpha"]**2  # viscosité optim, stockée en racine ds pvtu
    thi = mesh.point_data["h"]         # thicknessi m
 
    pts = mesh.points
    x = pts[:,0]
    y = pts[:,1]
    


    # ---------------------------------------------------------
    # 1) CALCUL DES GRADIENTS PAR VTK → hyper rapide et léger
    # ---------------------------------------------------------
    grad_vx = mesh.compute_derivative(scalars="uobs 1", gradient=True).point_data["gradient"]
    grad_vy = mesh.compute_derivative(scalars="uobs 2", gradient=True).point_data["gradient"]

    # grad_vx[:,0] = ∂vx/∂x, grad_vx[:,1] = ∂vx/∂y
    dudx = grad_vx[:,0]
    dudy = grad_vx[:,1]

    dvdx = grad_vy[:,0]
    dvdy = grad_vy[:,1]

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

    rf = eta  # viscosité, stockée en racine dans pvtu
    #factor = rf**(-1/nn) * eff**((1-nn)/(2*nn))
    #factor = rf * eff**((1-nn)/(2*nn))

    #ttxx = exx * factor
    #ttyy = eyy * factor
    #ttxy = exy * factor
    ttxx = 2 * eta * exx
    ttyy = 2 * eta * eyy
    ttxy = 2 * eta * exy
    # ---------------------------------------------------------
    # 4) TENSEUR T MODIFIÉ (stress tensor in SSA) 
    # ---------------------------------------------------------
    t11 = 2*ttxx + ttyy
    t12 = ttxy
    t21 = ttxy
    t22 = ttxx + 2*ttyy

    # ---------------------------------------------------------
    # 5) EIGENVALUES
    # ---------------------------------------------------------
    
    bb = -(t11 + t22)
    cc = t11*t22 - t21*t12

    disc = np.sqrt(bb**2 - 4*cc)

    eig1 = (-bb + disc) / 2 #* 1e6 to Pa ?
    eig2 = (-bb - disc) / 2 #* 1e6

    
    # ---------------------------------------------------------
    # 5bis) EIGENVECTOR associated with eig2
    # ---------------------------------------------------------
    eig2_vx = t12
    eig2_vy = eig2 - t11

    # normalisation
    norm = np.sqrt(eig2_vx**2 + eig2_vy**2) + 1e-12
    #norm = np.sqrt(1**2 + eig2_vy**2) + 1e-12
    eig2_vx /= norm
    eig2_vy /= norm

    if calcul_eig_furst_way :
        eig2 = np.sqrt(eig2_vx**2 + eig2_vy**2)

    # ---------------------------------------------------------
    # 6) DIRECTIONS / PROJECTIONS (KN1)
    # ---------------------------------------------------------
    vnorm = np.sqrt(vx**2 + vy**2) + 1e-12
    n1 = vx / vnorm
    n2 = vy / vnorm

    tau = ((t11*n1 + t12*n2)*n1 +
           (t21*n1 + t22*n2)*n2)

    psi = ((t11*n2 - t12*n1)*n2 -
           (t21*n2 - t22*n1)*n1)

    # ---------------------------------------------------------
    # 7) BUTTRESSING RATIOS
    # ---------------------------------------------------------
    tau0 = 0.5*rhoice*grav*thi*(1 - rhoice/rhow) 

    KT  = psi / tau0
    KN1 = 1 - tau  / tau0
    KN2 = 1 - eig1 / tau0
    KN3 = 1 - (eig2 / tau0) *10e6

    return x, y, vx, vy, t11, t12, t21, t22, KT, KN1, KN2, KN3, eig1, eig2, eig2_vx, eig2_vy




x, y, vx, vy, t11, t12, t21, t22, KT, KN1, KN2, KN3, eig1, eig2, eig2_vx, eig2_vy = compute_buttressing_unstructured(filename)

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

fig, ax = plot_corrected(
    x, y,
    KN3,
    vx,vy,
    cmap,        # ta palette blue–brown ou custom
    title="KN3 – Corrected Version - Friction at {beta}"
)
plt.savefig(f"../BUTTRESSING/buttressing_corrected_{year}_FRT_{beta}", dpi = 300)
plt.show()


fig, ax = plot_stress(
    x, y,t11, t12, t21, t22, year, beta
)
#plt.title(f"Stress for year {year}, friction param of {beta}")
plt.savefig(f"../STRESS/stress_{year}_FRT_{beta}", dpi = 300)
plt.show()


fig, ax = plot_stress_simple(
    x, y,t11, t12, t21, t22, year, beta
)
plt.savefig(f"../STRESS/stress_simple_{year}_FRT_{beta}", dpi = 300)
plt.show()

fig, ax = plot_eigenvec(
    x, y, eig2_vx, eig2_vy, year, beta
)
plt.savefig(f"../STRESS/eigenvec_{year}_FRT_{beta}", dpi = 300)
plt.show()


