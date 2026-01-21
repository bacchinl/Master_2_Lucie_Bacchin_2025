# plot_stress.py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap




# ===============================
# Colormap Furst et al.
# ===============================
def get_furst_colormap():
    colors = [
        "#39A0CD", "#54B6E0", "#6BC9F2", "#85D5F0",
        "#9BD5C5", "#B1D295", "#CACF67", "#E1CE34",
        "#F7CD19", "#E9BD1F", "#D0AA2C", "#B79530"
    ]
    return LinearSegmentedColormap.from_list(
        "color_fuerst", colors, N=12
    )

# ===============================
# Plot buttressing classique
# ===============================

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
    return fig, axes


def plot_stress(x, y, t11, t12, t21, t22, year, beta):

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    maps = [
        (t11,  "Extension compression(x axis)"),
        (t12, "Shearing stress txy=tyx"),
        (t21, "Shearing stress tyx=txy"),
        (t22, "Extension compression(y axis)"),
    ]

    for ax, (K, title) in zip(axes.ravel(), maps):

        scatter = ax.scatter(
            x, y,
            c=K,
            cmap="seismic",
            vmin=-1, vmax=1,
            marker='.',
            s=1,
            alpha=1
        )

        ax.set_title(f"{title} ({year})")
        ax.set_xlabel("X")

        # label Y uniquement sur la colonne de gauche
        if ax in axes[:,0]:
            ax.set_ylabel("Y")

        fig.colorbar(scatter, ax=ax, label="Pa")

    plt.suptitle(f"Stress for year {year}, friction param of {beta}") 
    plt.tight_layout()
    return fig, axes

# ===============================
# Plot corrigé (Furst-like)
# ===============================



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


