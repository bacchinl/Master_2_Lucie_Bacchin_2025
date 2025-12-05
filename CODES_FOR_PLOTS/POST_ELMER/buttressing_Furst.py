import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm
import numpy as np
import pyvista as pv
import matplotlib.tri as mtri

small = 1e-10
year = 2016
year_to_sec = 3.15576e+7

# print ressource 
import resource
print("Max memory:", resource.getrlimit(resource.RLIMIT_AS))

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
    vx = mesh.point_data["uobs 1"]     # vx m/a
    vy = mesh.point_data["uobs 2"]     # vy m/a
    eta = mesh.point_data["alpha"]**2 # viscosité optim, stockée en racine ds pvtu
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
    # 2) TENSEUR DE DÉFORMATION
    # ---------------------------------------------------------
    exx = dudx
    eyy = dvdy
    exy = 0.5 * (dudy + dvdx)
    ezz = -exx - eyy

    small = 1e-10
    eff = 0.5*(exx**2 + eyy**2 + ezz**2) + exy**2 + small

    # ---------------------------------------------------------
    # 3) DEVIATORIC STRESS
    # ---------------------------------------------------------
    rf = eta  # viscosité, stockée en racine dans pvtu
    factor = rf**(-1/nn) * eff**((1-nn)/(2*nn))

    ttxx = exx * factor
    ttyy = eyy * factor
    ttxy = exy * factor

    # ---------------------------------------------------------
    # 4) TENSEUR T MODIFIÉ
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

    eig1 = (-bb + disc) / 2
    eig2 = (-bb - disc) / 2

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
    tau0 = 0.5*rhoice*grav*thiampl*thi*(1 - rhoice/rhow)

    KT  = psi / tau0
    KN1 = 1 - tau  / tau0
    KN2 = 1 - eig1 / tau0
    KN3 = 1 - eig2 / tau0

    return x, y, vx, vy, KT, KN1, KN2, KN3
#------------------------------
#   PLOT
#------------------------------
# palette comme dans furst et al 

colors = [
    "#39A0CD","#54B6E0","#6BC9F2","#85D5F0", "#9BD5C5", "#B1D295", "#CACF67","#E1CE34","#F7CD19", "#E9BD1F", "#D0AA2C","#B79530"
]
color_fuerst = LinearSegmentedColormap.from_list("color_fuerst", colors, N=12)

def plot_buttressing(x, y, KT, KN1, KN2, KN3, color_fuerst, year):

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    maps = [
        (KT,  "KT – transverse"),
        (KN1, "KN1 – along-flow"),
        (KN2, "KN2 – eigenvalue 1"),
        (KN3, "KN3 – eigenvalue 2"),
    ]

    for ax, (K, title) in zip(axes.ravel(), maps):

        scatter = ax.scatter(
            x, y,
            c=K,
            cmap=color_fuerst,
            vmin=0, vmax=1,
            marker='.',
            s=1,
            alpha=1
        )

        ax.set_title(f"{title} ({year})")
        ax.set_xlabel("X")

        # label Y uniquement sur la colonne de gauche
        if ax in axes[:,0]:
            ax.set_ylabel("Y")

        fig.colorbar(scatter, ax=ax, label="buttressing number")


    plt.tight_layout()
    plt.show()
    plt.savefig(f"./Subplot_buttressing_Furst_{year}", dpi = 300)


def plot_corrected(x, y, KN3, vx, vy, color_f, title="Corrected KN3"):
    """
    Make it look like Furst et al 2016
    - KN3 * 1 / sign(vm)
    - niveaux [-10 -2 -0.5 0:0.1:1 2 4 100]
    - colormap blue–brown (déjà définie dans color_f)
    - contourf non linéaire, style MATLAB
    """

     # ---------------------------
    # Calcul de vmELM = magnitude du vecteur vitesse
    # ---------------------------
    vm = np.sqrt(vx**2 + vy**2)
    vm[vm == 0] = 1e-12  # éviter division par zéro

    # ---------------------------
    # TRICK MATLAB
    # ---------------------------
    Kn_plot = (KN3 / vm)*100    # pas de shelfmask car tu es sur un seul ice shelf
    test_var(Kn_plot, "Kn_plot (KN corrected)")
    # ---------------------------
    # same contour levels
    # ---------------------------
    levels = [-10, -2, -0.5] + list(np.arange(0, 1.1, 0.1)) + [2, 4, 100]

    # ---------------------------
    # figure
    # ---------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    scatter = ax.scatter(
        x, y,
        c=Kn_plot,
        cmap=color_f,
        marker='.',
        s=1,
        alpha=1,
        vmin=0, vmax=1  # tu peux adapter la plage
    )
    # ---------------------------
    # colorbar
    # ---------------------------
    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Corrected Kn")

    # ---------------------------
    # labels, title
    # ---------------------------
    ax.set_aspect("equal")
    ax.set_title(title)
    ax.set_xlabel("X (m)")
    ax.set_ylabel("Y (m)")

    return fig, ax



##### LAUNCH #####

x, y, vx, vy, KT, KN1, KN2, KN3 = compute_buttressing_unstructured(filename)

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




#plot_buttressing(x, y, KT, KN1, KN2, KN3, color_fuerst, year)

fig, ax = plot_corrected(
    x, y,
    KN3,
    vx,vy,
    color_fuerst,        # ta palette blue–brown ou custom
    title="KNflow – Corrected Version - Friction at {beta}"
)
#plt.savefig(f"./buttressing/buttressing_corrected_{year}_FRT_{beta}", dpi = 300)
plt.show()
