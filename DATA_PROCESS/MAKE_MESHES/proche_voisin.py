#afinn de faire le mesh, il est necessaire de réarranger les points constituants le front et ceux constituants la GL dans un ordre défini. Réglez les "seuil_distance pour supprimer pinning points et artefacts.

import numpy as np
import sys
year=sys.argv[1]

#try : 
#    sys.argv[1]
#except : 
#    print('error : try python proche_voisin.py year ')
#    exit()

#year = sys.argv[1]

def read_points(file_path):
    """
    Lire les points depuis un fichier texte avec deux colonnes (x, y).
    """
    points = []
    with open(file_path, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split('\t'))  # Supposons que les coordonnées sont séparées par des tabulations
            points.append((x, y))
    return points

def write_points(file_path, points):
    """
    Écrire les points triés dans un fichier texte.
    """
    with open(file_path, 'w') as f:
        for x, y in points:
            f.write(f"{x}\t{y}\n")

def find_nearest_neighbor(current_point, points):
    """
    Trouver le point le plus proche parmi les points restants.
    """
    min_dist = float('inf')
    nearest_point = None
    nearest_index = -1
    for i, point in enumerate(points):
        dist = np.sqrt((current_point[0] - point[0])**2 + (current_point[1] - point[1])**2)
        if dist < min_dist:
            min_dist = dist
            nearest_point = point
            nearest_index = i
    return nearest_point, nearest_index

def sort_points(points, seuil):
    """
    Trier les points en commençant par le point le plus bas à gauche et en utilisant l'algorithme du plus proche voisin.
    """
    # Trouver le point de départ (celui en bas à gauche)
    start_point = min(points, key=lambda p: (p[1], p[0]))  # Priorité à y, puis à x
    sorted_points = [start_point]
    points.remove(start_point)

    current_point = start_point
    while points:
        nearest_point, nearest_index = find_nearest_neighbor(current_point, points)
        dist = np.sqrt((current_point[0] - nearest_point[0])**2 + (current_point[1] - nearest_point[1])**2)
        
        # Ne prendre en compte que les points dans une distance raisonnable
        if dist <= seuil_distance:
            sorted_points.append(nearest_point)
            current_point = nearest_point
            points.pop(nearest_index)  # Retirer le point visité
        else:
            points.pop(nearest_index)  # Si trop loin, on l'ignore quand même 

    return sorted_points

# Exemple d'utilisation
input_file = f'./{year}_State/contours_GL_{year}_Lucille.txt'  # Remplace par le chemin de ton fichier de coordonnées
output_file = f'../../../ELMERICE/Mesh_Generation/Contours/{year}_State/contours_1_{year}_Lucille.txt'
#output_file = f'./{year}_State/contours_1_{year}_Lucille.txt'

#input_file = './../../../ELMERICE/Mesh_Generation/Contours/Contour_1_inv.txt'  # Remplace par le chemin de ton fichier de coordonnées
#output_file = './../../../ELMERICE/Mesh_Generation/Contours/Contour_1_inv.txt'


# Lire les points depuis le fichier texte
points = read_points(input_file)

# Définir une distance seuil au-dessus de laquelle les points sont ignorés
seuil_distance = 80000  # 20000 en m

# Trier les points
sorted_points = sort_points(points, seuil_distance)

# Écrire les points triés dans un nouveau fichier
write_points(output_file, sorted_points)

# Exemple d'utilisation
input_file = f'./{year}_State/contours_FL_{year}_Lucille.txt'  # Remplace par le chemin de ton fichier de coordonnées
output_file = f'../../../ELMERICE/Mesh_Generation/Contours/{year}_State/contours_2_{year}_Lucille.txt'

#output_file = f'./{year}_State/contours_2_{year}_Lucille.txt'

# Lire les points depuis le fichier texte
points = read_points(input_file)

# Définir une distance seuil au-dessus de laquelle les points sont ignorés
seuil_distance = 80000  # en m

# Trier les points
sorted_points = sort_points(points, seuil_distance)

# Écrire les points triés dans un nouveau fichier
write_points(output_file, sorted_points)


