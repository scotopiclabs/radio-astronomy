"""
plot_speed_curve - Production de la courbe de rotation de la Galaxie. On utilise les données des points tangents
disponibles: le rayon tangeant R_tangent_kpc et la vitesse orbitale tangeante V_orb_tangent.
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
from matplotlib.ticker import MultipleLocator

#-----------------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Tracé de la courbe de rotation de la Galaxie à partir des points tangents.")
    parser.add_argument("csv_file", help="Fichier CSV de synthèse (ex: synthese_observations.csv)")
    parser.add_argument("-d", "--degree", type=int, default=3, help="Degré de la régression polynomiale (défaut: 3)")
    parser.add_argument("-l", "--labels", action="store_true", help="Affichage des étiquettes de longitude à côté des points")
    parser.add_argument("-e", "--error", action="store_true", help="Affichage des barres d'erreur")
    args = parser.parse_args()

    print(f"Chargement du fichier : {args.csv_file}")
    try:
        df = pd.read_csv(args.csv_file, sep=';', skipinitialspace=True)
    except FileNotFoundError:
        raise SystemExit(f"Erreur : Le fichier '{args.csv_file}' est introuvable.")

    df.columns = df.columns.str.strip()
    
    # Conversion des virgules en points le cas échéant et forçage numérique
    cols = ['R_tangent_kpc', 'V_orb_tangent']
    for col in cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col].astype(str).str.replace(',', '.'), errors='coerce')

    # Filtrage : on ne garde que la ligne correspondant au pic tangent (la seule renseignée pour chaque longitude)
    df_tangent = df.dropna(subset=['R_tangent_kpc', 'V_orb_tangent']).copy()
            
    if df_tangent.empty:
        raise SystemExit("Erreur : Aucune donnée de point tangent trouvée dans le CSV (vérifiez que les quadrants I/IV sont présents).")

    print(f"Extraction de {len(df_tangent)} points tangents.")

    R = df_tangent['R_tangent_kpc'].values
    V = df_tangent['V_orb_tangent'].values

    # Tri des données pour que la ligne de régression polynomiale soit propre
    sort_idx = np.argsort(R)
    R_sorted = R[sort_idx]
    V_sorted = V[sort_idx]

    degree = min(args.degree, len(np.unique(R)) - 1)
    coeffs = np.polyfit(R_sorted, V_sorted, degree)
    poly = np.poly1d(coeffs)

    # --- Tracé de la figure ---
    
    plt.figure(figsize=(10, 6))
        
    # Tracé de la ligne de régression
    x_plot = np.linspace(R_sorted.min(), R_sorted.max(), 300)
    plt.plot(x_plot, poly(x_plot), color='gray', linestyle='-', linewidth=1, 
             label=f'Régression polynomiale (degré {degree})', zorder=4)

    # Calcul de la vitesse sur la courbe de régression pour le rayon maximum
    r_max = R_sorted.max()
    v_max_calc = poly(r_max)

    # Affichage d'un marqueur et de l'étiquette de valeur
    plt.scatter(r_max, v_max_calc, color='red', marker='s', s=20, zorder=6, label='Vitesse max (régression)')
    plt.text(r_max - 1.5, 10, f"V({r_max:.1f} kpc) = {v_max_calc:.1f} km/s", 
             color='red', fontsize=10, fontweight='bold', zorder=7,
             bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
                 
    # Tracé des barres d'erreur si demandé
    if args.error:
        # Erreur asymptotique (dispersion intrinsèque du gaz HI en km/s)
        err_base = 10.0
        
        # Terme d'erreur supplémentaire modélisant les mouvements non circulaires de la barre
        # L'amplitude (60.0) et la constante de décroissance (3 kpc) sont ajustables
        erreur_V = err_base + 60.0 * np.exp(-R / 3)
    
        # Tracé des points de données avec barres d'erreur
        plt.errorbar(R, V, yerr=erreur_V, fmt='o', color='lightblue', markeredgecolor='black', 
                     markersize=6, ecolor='gray', elinewidth=1, capsize=3,
                     label='Points tangents mesurés', zorder=3)
    else:
        # Tracé des points de données sans barre d'erreur
        plt.scatter(R, V, color='lightblue', edgecolor='black', s=40, 
                    label='Points tangents mesurés', zorder=3)

    # Tracé des étiquettes si demandé
    if args.labels:
        for i, row in df_tangent.iterrows():
            r_val = row['R_tangent_kpc']
            v_val = row['V_orb_tangent']
            
            try:
                lon = int(round(float(row['Longitude'])))
            except ValueError:
                lon = row['Longitude']
                
            v_lsr = row['V_lsr']
            
            # Décalage de +0.1 en X et +1.5 en Y pour la lisibilité
            #plt.text(r_val + 0.08, v_val + 1.5, f"{lon}°, {v_lsr}", 
                     #fontsize=8, color='black', alpha=0.8, zorder=5)
            plt.text(r_val + 0.06, v_val - 6, f"{lon}°", 
                     fontsize=9, color='black', alpha=0.8, zorder=5)
                     
    plt.xlabel('Rayon galactique $R$ (kpc)', fontsize=13)
    plt.ylabel('Vitesse orbitale $V(R)$ (km/s)', fontsize=13)
    plt.title('Courbe de rotation galactique (points tangents)', fontsize=16, fontweight='normal')
    
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', alpha=0.3)
    plt.legend(loc='best', fontsize=11)
    plt.tight_layout()
    
    # Affichage
    print("Affichage de la courbe...")
    plt.show()

if __name__ == "__main__":
    main()