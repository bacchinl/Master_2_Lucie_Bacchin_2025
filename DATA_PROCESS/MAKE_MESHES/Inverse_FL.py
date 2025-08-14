import sys

year=2000
def inverser_lignes(input_file, output_file):
    """
    Lit un fichier texte, inverse les lignes et écrit le résultat dans un nouveau fichier.
    
    :param input_file: Chemin du fichier texte à inverser.
    :param output_file: Chemin du fichier texte où enregistrer les lignes inversées.
    """
    with open(input_file, 'r') as f:
        lignes = f.readlines()  # Lire toutes les lignes du fichier
        
    lignes_inversées = lignes[::-1]  # Inverser l'ordre des lignes
    
    with open(output_file, 'w') as f:
        f.writelines(lignes_inversées)  # Écrire les lignes inversées dans le fichier de sortie

# Exemple d'utilisation
#input_file = f'{year}_State/contours_2_{year}_paolo_Clean.txt'  # Remplace par ton fichier d'entrée
#output_file = f'{year}_State/contours_2_{year}_paolo_Clean.txt'  # Fichier de sortie avec les lignes inversées
input_file = f'2.txt'
output_file = f'2.txt'


inverser_lignes(input_file, output_file)

print(f"Les lignes du fichier '{input_file}' ont été inversées et enregistrées dans '{output_file}'.")



