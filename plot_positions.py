"""
plot_positions - Tracé de la carte de positionnement des nuages d'hydrogène neutre dans la Galaxie. Un modèle récent
de bras spiraux a été utilisé dans la carte de la Galaxie. Il s'agit du modèle publié dans l'article
"Trigonometric Parallaxes Of High-Mass Star Forming Regions: Our View Of The Milky Way" (M. J. Reid et al.).
Voir https://arxiv.org/abs/1910.03357.
"""

# Copyright (C) 2026  Philippe DEVERCHERE
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D

# Soleil à (0, 8.5) en coordonnées catésiennes (kpc). Le centre galactique est à (0, 0).
# La vitesse radiale du Soleil est de 220 km/s
V0 = 220.0
R0 = 8.5

#-----------------------------------------------------------------------------------------------------------------------
def get_spiral_arms():
    """
    Retourne un modèle des coordonnées (X, Y) des portions des bras principaux dans le quadrant I selon un modèle
    dérivé de Reid et al. (2014). L'affichage est strictement découpé et limité au Quadrant I (X >= 0 et Y <= 8.5).
    """
    arms = {}

    # Paramètres cinématiques stricts de Reid et al. (2014). La signification des paramètres est la suivante:
    # - R_ref (Rayon de référence, en kpc) : Il s'agit du rayon galactocentrique du bras à une position angulaire
    #   donnée définie par beta_ref. C'est le point d'ancrage radial du bras dans le disque galactique.
    # - beta_ref (Azimut de référence, en degrés) : C'est l'angle galactocentrique (mesuré depuis l'axe du Soleil)
    #   auquel le bras se trouve à la distance R_ref. Dans le modèle, ce paramètre positionne le début du tracé
    #   du bras.
    # - psi (angle de tangage, en degrés) : Cet angle définit l'ouverture de la spirale. Un angle faible
    #   (ex: 5° à 8°) signifie que le bras est très serré et s'enroule lentement autour du centre.
    # - Un angle élevé (ex: 15° à 20°) signifie que le bras est très ouvert et s'éloigne rapidement du centre
    #   galactique.
    # - ranges (Plages angulaires, en degrés) : Ce paramètre définit les bornes de l'azimut lesquelles le bras
    #   doit être tracé. Il permet de limiter physiquement la longueur du segment tracé sur la carte, évitant
    #   ainsi que le bras ne s'enroule indéfiniment autour de la Galaxie ou ne traverse des zones non désirées.
    models = {
        # Le bras Scutum-Centaurus-OSC possède deux segments distincts dans le Quadrant I :
        # la partie interne et la partie externe au tour suivant.
        '3kpc': {'R_ref': 3.6, 'beta_ref': 25, 'psi': 5, 'ranges': [[-60, 180]]},
        'Scutum-Centaurus-OSC-1': {'R_ref': 5.0, 'beta_ref': 15.6, 'psi': 11, 'ranges': [[-30, 80]]},
        'Scutum-Centaurus-OSC-2': {'R_ref': 5.3, 'beta_ref': 6, 'psi': 11.1, 'ranges': [[-370, -170]]},
        'Sagittarius-Carina': {'R_ref': 6.2, 'beta_ref': 15.6, 'psi': 6.9, 'ranges': [[-30, 180]]},
        'Local': {'R_ref': 8.4, 'beta_ref': 8.9, 'psi': 12.8, 'ranges': [[-20, 35]]},
        'Perseus': {'R_ref': 9, 'beta_ref': 30.2, 'psi': 12.4, 'ranges': [[-30, 180]]},
        'Norma-Outer': {'R_ref': 12, 'beta_ref': 30.2, 'psi': 12.4, 'ranges': [[-30, 180]]}
    }
    
    for name, m in models.items():
        x_list, y_list = [], []
        
        for b_min, b_max in m['ranges']:
            # Génération d'une séquence de points très dense pour une courbe lisse
            beta_deg = np.linspace(b_min, b_max, 500)
            beta = np.radians(beta_deg)
            beta_ref_rad = np.radians(m['beta_ref'])
            psi_rad = np.radians(m['psi'])
            
            # Équation de la spirale logarithmique
            r = m['R_ref'] * np.exp(-(beta - beta_ref_rad) * np.tan(psi_rad))
            
            # Conversion cartésienne
            x = r * np.sin(beta)
            y = r * np.cos(beta)
            
            # Masque strict limitant géométriquement au Quadrant I
            mask = (x >= 0) & (y <= 8.5)
            
            # Le bras de 3 kpc est limité à la région interne
            if name == '3kpc':
                mask = mask & (r >= 2.0) & (r <= 4.0)
                
            # Extraction exclusive des points valides formant une ligne continue
            x_seg = x[mask]
            y_seg = y[mask]
            
            if len(x_seg) > 0:
                # Injection d'un NaN pour séparer proprement les segments multiples 
                # d'un même bras (le cas de Scutum-Centaurus-OSC) et éviter les lignes droites
                if len(x_list) > 0:
                    x_list.append(np.nan)
                    y_list.append(np.nan)
                x_list.extend(x_seg)
                y_list.extend(y_seg)
                
        if len(x_list) > 0:
            arms[name] = (np.array(x_list), np.array(y_list))
            
    return arms
    
#-----------------------------------------------------------------------------------------------------------------------
def get_min_distance_to_arms(x, y, arms):
    """
    Calcule la distance géométrique minimale entre un point (x,y) et les bras spiraux,
    en ignorant les valeurs NaN utilisées pour le formatage du tracé.
    """
    min_distance = float('inf')
    for name, (arm_x, arm_y) in arms.items():
        # Calcul vectoriel
        distances = np.sqrt((arm_x - x)**2 + (arm_y - y)**2)
        
        # Filtrage des NaN pour ne conserver que les distances réelles
        valid_distances = distances[~np.isnan(distances)]
        
        if len(valid_distances) > 0:
            min_dist_to_this_arm = np.min(valid_distances)
            if min_dist_to_this_arm < min_distance:
                min_distance = min_dist_to_this_arm
                
    return min_distance

#-----------------------------------------------------------------------------------------------------------------------
def select_distance_by_arms(d1, d2, l_rad, arms, threshold_kpc=1.0):
    """
    Sélectionne la distance (Near ou Far) en fonction de la proximité aux bras spiraux.
    Le seuil par défaut est fixé à 1.0 kpc.
    """
    # Identification stricte de Near (plus petite distance) et Far (plus grande distance)
    d_near = min(d1, d2)
    d_far = max(d1, d2)
    
    # Projection cartésienne du point Near
    x_near = d_near * np.sin(l_rad)
    y_near = 8.5 - d_near * np.cos(l_rad)
    
    # Projection cartésienne du point Far
    x_far = d_far * np.sin(l_rad)
    y_far = 8.5 - d_far * np.cos(l_rad)
    
    # Distances respectives aux bras spiraux modélisés
    dist_arms_near = get_min_distance_to_arms(x_near, y_near, arms)
    dist_arms_far = get_min_distance_to_arms(x_far, y_far, arms)
    
    # Application de la règle probabiliste
    if dist_arms_far <= threshold_kpc and dist_arms_far < dist_arms_near:
        return d_far
    else:
        return d_near
        
#-----------------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Tracé de la carte galactique à partir d'un fichier de synthèse CSV.")
    parser.add_argument("csv_file", help="Fichier CSV de synthèse")
    parser.add_argument("-l", "--labels", action="store_true", help="Agffichage des étiquettes (Longitude, V_lsr) à côté des points")
    parser.add_argument("-s", "--sight", action="store_true", help="Tracé des lignes de visée")
    parser.add_argument("-n", "--near_far", action="store", help="Choix pour les positions Near-Far: a (auto), n (near) ou f (far)")
    parser.add_argument("-m", "--median_arm", action="store_true", help="Tracé des lignes médianes des bras spiraux")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_file, sep=';')
    
    l_rad = np.radians(df['Longitude'])
    
    # Génération du modèle des bras spiraux pour la comparaison
    arms_model = get_spiral_arms()
    
    # Calcul de la distance la plus probable dans le cadre d'un couple Near Far
    df['Selected_Dist'] = np.nan
    for i, row in df.iterrows():
        if not np.isnan(row['Dist_2_kpc']):
            l2 = np.radians(row['Longitude'])
            d1 = row['Dist_1_kpc']
            d2 = row['Dist_2_kpc']
            
            # Sélection géométrique (seuil ajustable, ex: 1.0 kpc)
            df.at[i, 'Selected_Dist'] = select_distance_by_arms(d1, d2, l2, arms_model, threshold_kpc=1)
        else:
            df.at[i, 'Selected_Dist'] = row['Dist_1_kpc']
        
    # On utilise le modèle de distance Near / Far demandé par l'utilisateur
    near_far = 0
    if args.near_far:
        if args.near_far == 'n':
            near_far = 1
        elif args.near_far == 'f':
            near_far = 2
    if near_far == 0:
        d = df['Selected_Dist']
        d = d.fillna(df['Dist_1_kpc'])
    elif near_far == 1:
        d = df['Dist_1_kpc']
        d = d.fillna(df['Dist_2_kpc'])
    else:
        d = df['Dist_2_kpc']
        d = d.fillna(df['Dist_1_kpc'])
    
    # Conversion coordonnées polaires (l, d) en cartésiennes (X, Y)
    # Pour que l=0 pointe vers le centre galactique (0,0) : X = d * sin(l), Y = 8.5 - d * cos(l)
    x = d * np.sin(l_rad)
    y = 8.5 - d * np.cos(l_rad)

    # --- Création de la figure ---
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # 1. Tracé des bras
    if args.median_arm:
        for name, (a_x, a_y) in arms_model.items():
            plt.plot(a_x, a_y, label=f"Bras {name}", color='black', linewidth=1, alpha=0.5)

    # --- Chargement de l'image de fond ---
    try:
        # 'extent' définit les limites du graphique (gauche, droite, bas, haut) en kpc
        # Ces valeurs sont ajustées pour que l'image colle parfaitement aux axes du graphique
        img = mpimg.imread('bras_piraux.jpg')
        offset_x = -1.19
        offset_y = -1.06
        scale = 17.4
        ax.imshow(img, extent=[-scale+offset_x, scale+offset_x, -scale+offset_y, scale+offset_y], aspect='equal', alpha=0.8, zorder=0)
    except FileNotFoundError:
        print("Image 'bras_piraux.jpg' introuvable.")
        
    # Tracé du Soleil et du Centre Galactique
    ax.scatter(0, 8.5, color='black', marker='*', s=80, label='Soleil', zorder=10)
    ax.scatter(0, 0, color='black', marker='x', s=100, label='Centre Galactique', zorder=10)

    # --- Tracé des lignes de visée ---
    if args.sight:
        for i in df.index:
            # Filtre sur les quadrants I, II et III et vérification de validité des coordonnées
            if df.at[i, 'Quadrant'] in [1, 2, 3] and pd.notna(x[i]) and pd.notna(y[i]):
                # Tracé d'une ligne entre le Soleil (0, 8.5) et le point mesuré (x[i], y[i])
                ax.plot([0, x[i]], [8.5, y[i]], color='gray', linestyle='--', linewidth=0.3, alpha=0.5, zorder=1)
                
    # On crée un point noir fictif pour la légende
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Mesures HI', markerfacecolor='black', markersize=6),
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='black', markersize=14, label='Soleil'),
        plt.Line2D([0], [0], marker='x', color='black', markersize=8, label='Centre Galactique', linestyle='None')
    ]
    
    # Tracé des points mesurés
    color_quadrant = {1: '#d62728', 2: '#ff7f0e', 3: '#2ca02c', 4: '#1f77b4'}
    ax.scatter(x, y, c=df['Quadrant'].map(color_quadrant), s=25, edgecolors='white', linewidth=0.5, zorder=5)

    # --- Tracé des étiquettes (Longitude, V_lsr) ---
    if args.labels:
        for i in df.index:
            # Vérification que le point possède des coordonnées valides
            if pd.notna(x[i]) and pd.notna(y[i]):
                lon = int(round(df.at[i, 'Longitude']))
                v_lsr = df.at[i, 'V_lsr']
                # Ajout d'un léger décalage (+0.2) pour ne pas masquer le point
                ax.text(x[i] + 0.2, y[i] + 0.2, f"{lon}°, {v_lsr}", 
                        fontsize=7, color='black', alpha=0.8, zorder=8)
                        
    # --- Tracé des limites de quadrants ---
    # Lignes croisées au niveau du Soleil (0, 8.5)
    ax.axhline(y=8.5, color='gray', linestyle='-', linewidth=1.5, zorder=2, alpha=0.8)
    ax.axvline(x=0, color='gray', linestyle='-', linewidth=1.5, zorder=2, alpha=0.8)

    # Ajout des labels (I, II, III, IV) avec un fond blanc semi-transparent pour la lisibilité
    bbox_props = dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='none')
    ax.text(13, -12, 'I', fontsize=16, color=color_quadrant[1], fontweight='bold', ha='center', va='center', zorder=3, bbox=bbox_props)
    ax.text(13, 13, 'II', fontsize=16, color=color_quadrant[2], fontweight='bold', ha='center', va='center', zorder=3, bbox=bbox_props)
    ax.text(-13, 13, 'III', fontsize=16, color=color_quadrant[3], fontweight='bold', ha='center', va='center', zorder=3, bbox=bbox_props)
    ax.text(-13, -12, 'IV', fontsize=16, color=color_quadrant[4], fontweight='bold', ha='center', va='center', zorder=3, bbox=bbox_props)

    # --- Tracé des labels des bras spiraux ---
    # Les coordonnées et rotations sont ajustées itérativement
    bras_labels = [
        {'nom': 'Local arm', 'x': -5.5, 'y': 7.9, 'rot': 17},
        {'nom': 'Perseus arm', 'x': -6.5, 'y': 9.2, 'rot': 17},
        {'nom': 'Sagittarius–Carina arm', 'x': 0.0, 'y': -13.2, 'rot': -21},
        {'nom': 'Norma–Outer arm', 'x': -8.5, 'y': 10.7, 'rot': 17},
        {'nom': 'Scutum–Centaurus–OSC arm', 'x': 3.0, 'y': 17.5, 'rot': -35},
        {'nom': '3kpc arm', 'x': -2.5, 'y': 1.6, 'rot': 62}
    ]

    for bras in bras_labels:
        ax.text(bras['x'], bras['y'], bras['nom'], 
                fontsize=11, fontweight='normal', color='dimgray', 
                rotation=bras['rot'], ha='center', va='center', zorder=4,
                bbox=dict(boxstyle='round,pad=0.1', facecolor='white', alpha=0.5, edgecolor='none'))
                
    ax.set_xlim(-15, 15)
    ax.set_ylim(-15, 15)
    ax.set_xlabel('X (kpc)')
    ax.set_ylabel('Y (kpc)')
    plt.title('Positions des nuages HI dans la Voie lactée', fontsize=16, fontweight='normal')
    plt.axis('equal')
    plt.grid(True, linestyle='--', alpha=0.5)
    ax.legend(handles=legend_elements, loc='best')
    plt.tight_layout()
    
    print("Affichage de la carte galactique...")

    plt.show()

    # Recherche de points aberrants pour les signaler
    mask_aberrant = (df['Quadrant'] == 2) & (y < 0)
    lignes_aberrantes = df[mask_aberrant]
    
    if not lignes_aberrantes.empty:
        print("Données du ou des points aberrants :")
        print(lignes_aberrantes[['Longitude', 'Quadrant', 'V_lsr', 'Dist_1_kpc', 'Dist_2_kpc', 'Selected_Dist']])

if __name__ == "__main__":
    main()