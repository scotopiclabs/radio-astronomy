"""
plot_spectrum - Tracé d'un spectre vitesse radiale / température à partir d'un fichier de données
issu d'un radio-télescope EU-HOU.
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

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, NullFormatter
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

V0 = 220.0
R0 = 8.5

#-----------------------------------------------------------------------------------------------------------------------
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

#-----------------------------------------------------------------------------------------------------------------------
def multi_gaussian(x, *params):
    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        y += gaussian(x, params[i], params[i+1], params[i+2])
    return y

#-----------------------------------------------------------------------------------------------------------------------
def estimate_signal_window(v, t):
    """Estime automatiquement la fenêtre de signal en excluant le bruit de fond."""
    idx = np.argsort(v)
    v_sorted = v.values[idx] if hasattr(v, 'values') else v[idx]
    t_sorted = t.values[idx] if hasattr(t, 'values') else t[idx]

    # On utilise les 15% aux extrémités pour définir le bruit de référence
    n = len(v_sorted)
    margin = max(int(n * 0.15), 10) 
    
    noise_v = np.concatenate([v_sorted[:margin], v_sorted[-margin:]])
    noise_t = np.concatenate([t_sorted[:margin], t_sorted[-margin:]])
    
    # Soustraction basique pour aplanir lors de l'estimation
    p = np.polyfit(noise_v, noise_t, 1)
    t_flat = t_sorted - np.polyval(p, v_sorted)
    
    # Calcul du seuil de détection (RMS)
    noise_std = np.std(np.concatenate([t_flat[:margin], t_flat[-margin:]]))
    
    # Détection du signal (seuil à 3*sigma)
    signal_mask = t_flat > (3 * noise_std)
    if not np.any(signal_mask):
        return -30.0, 30.0
        
    v_signal = v_sorted[signal_mask]
    return np.min(v_signal) - 15.0, np.max(v_signal) + 15.0

#-----------------------------------------------------------------------------------------------------------------------
def get_galactic_distances(v_lsr, l_rad):
    """Calcule la distance cinématique pour la galaxie interne (2 solutions) ou externe (1 solution).
    Args:
        v_lsr   Vitesse locale standard de repos --> vitesse radiale résultant de la différence de rotation
                entre le nuage et le Soleil.
        l_rad   Longitude galactique de la ligne de visée en radians
    Returns:
        Liste constituée de deux distances (quadrants I et IV) ou d'une seule distance (quadrants II et III). On
        retourne None si l'on est très proche des longitudes 0 et 180 pour lesquelles la vitesse radiale est nulle
        (i.e. on ne peut pas déterminer une distance)
    """
    # Éviter la division par zéro aux longitudes 0 et 180
    if abs(np.sin(l_rad)) < 0.05:
        return None
            
    # R est la distance au centre galactique
    # Formule : V_lsr = V0 * sin(l) * (R0/R - 1)
    # => R = R0 / (V_lsr / (V0 * sin(l)) + 1)
    denom = (v_lsr / (V0 * np.sin(l_rad))) + 1
    if denom <= 0: return None
    R = R0 / denom
    
    # --- Cas 1 : Galaxie Externe (R >= R0) pour les quadrants II et III ---
    # Il n'y a qu'une seule distance possible
    if R >= R0:
        d = R0 * np.cos(l_rad) + np.sqrt(max(0, R**2 - R0**2 * np.sin(l_rad)**2))
        return [round(d, 2)]
    
    # --- Cas 2 : Galaxie Interne (R < R0) pour les quadrants I et IV ---
    # Ambiguïté Near/Far : deux distances possibles
    else:
        det = R**2 - (R0**2 * (np.sin(l_rad)**2))
        if det < 0: return None
        d1 = R0 * np.cos(l_rad) - np.sqrt(det)
        d2 = R0 * np.cos(l_rad) + np.sqrt(det)
        return sorted([round(d1, 2), round(d2, 2)])

#-----------------------------------------------------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Analyse d'un spectre radio HI (21 cm) à partir de données brutes")
    parser.add_argument("measure_file", help="Fichier contenant la mesure au format EU-HOU")
    parser.add_argument("mask_min", type=float, help="Vitesse du début de la zone de signal (km/s)")
    parser.add_argument("mask_max", type=float, help="Vitesse de la fin de la zone de signal (km/s)")
    parser.add_argument("offset", type=float, help="Offset en vitesse à appliquer au spectre (km/s)")
    parser.add_argument("display", type=int, help="Affichage de la courbe produite (1=oui, 0=non)")
    parser.add_argument("-a", "--add", type=str, help="Pics à ajouter, séparés par des virgules (ex: 2.6,21.4)")
    parser.add_argument("-r", "--remove", type=str, help="Index des pics à supprimer, séparés par des virgules (ex: 1,3)")
    parser.add_argument("-l", "--legend", type=str, help="Position de la légende (ur, ul, lr, ll, cr, cl)")
    parser.add_argument("-p", "--peak_height", type=float, help="Ratio de hauteur pour détecter un pic (défaut 0.4)")
    parser.add_argument("-b", "--blank", action="store_true", help="Ne pas enregistrer le résulat dans le fichier CSV de sortie")
    args = parser.parse_args()

    basename = os.path.basename(args.measure_file)
    
    # Expression régulière pour cibler et extraire les composants du nom de fichier
    pattern = r"^(.*?)-(\d{4}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2}_\d{2})-(-?\d+)-(-?\d+)-SRT-(.*?)-3m\.fits\.csv$"
    match = re.match(pattern, basename)
    
    if match:
        # Récupération des différents blocs
        obs = match.group(1)
        raw_date = match.group(2)
        lon_str = match.group(3)
        lat_str = match.group(4)
        location = match.group(5)
    
        # Conversion de la longitude (-180/180 vers 0/360)
        # L'opérateur modulo % en Python gère automatiquement les nombres négatifs (ex: -10 % 360 = 350)
        galac_lon = float(lon_str) % 360.0
        l_rad = np.radians(galac_lon)
        
        # Assignation de la latitude
        galac_lat = float(lat_str)
        
        # Formatage de la date "yyyy-mm-dd hh:mm:ss"
        # On découpe la chaîne brute et on ignore le dernier élément (les centièmes/millisecondes)
        y, mo, d, h, mi, s, ms = raw_date.split('_')
        utc_datetime = f"{y}-{mo}-{d} {h}:{mi}:{s}"
    else:
        raise SystemExit(f"Le fichier '{basename}' ne respecte pas le format attendu.")
        
    # --- Lecture des données de mesure ---
    print(f"\nFichier de mesures: {basename}")
    df_mesure = pd.read_csv(args.measure_file, sep=';', skipinitialspace=True)
    df_mesure.columns = df_mesure.columns.str.strip()
    if 'Freq' in df_mesure.columns:
        df_mesure = df_mesure.rename(columns={'Freq': 'Frequency'})
    exposure = int(df_mesure['ObsTime'].iloc[0])
    
    C = 299792.458
    F0 = 1420.405751
    vitesse_brute = C * (F0 - df_mesure['Frequency']) / F0
    vlsr_correction = df_mesure['Vlsr'].iloc[0]
    vitesse_mesure = vitesse_brute + vlsr_correction + args.offset

    # --- Détermination du quadrant de la longitude galactique ---
    if galac_lon >= 0 and galac_lon < 90:
        quadrant = 'I'
    elif galac_lon >= 90 and galac_lon < 180:
        quadrant = 'II'
    elif galac_lon >= 180 and galac_lon < 270:
        quadrant = 'III'
    elif galac_lon >= 270 and galac_lon < 360:
        quadrant = 'IV'
    else:
        raise SystemExit('Longitude galactique incorrecte: {galac_lon}')
    print(f"Longitude galactique: {galac_lon}")
    print(f"Quadrant: {quadrant}")

    # --- Estimation automatique de la fenêtre de signal ---
    est_min, est_max = estimate_signal_window(vitesse_mesure, df_mesure['Temperature'])
    
    if args.mask_min == 0.0 and args.mask_max == 0.0:
        mask_min, mask_max = est_min, est_max
        print(f"Utilisation du masque automatique : [{mask_min:.1f}, {mask_max:.1f}] km/s")
    else:
        mask_min, mask_max = args.mask_min, args.mask_max

    # --- Soustraction de la ligne de base ---
    temperature_brute = df_mesure['Temperature']
    masque_bruit = (vitesse_mesure < mask_min) | (vitesse_mesure > mask_max)
    coefficients = np.polyfit(vitesse_mesure[masque_bruit], temperature_brute[masque_bruit], 2)
    temperature_corrigee = temperature_brute - np.polyval(coefficients, vitesse_mesure)

    # --- Identification des pics ---
    height = 0.4
    if args.peak_height:
        height = float(args.peak_height)
    peaks, props = find_peaks(temperature_corrigee, height=np.max(temperature_corrigee)*height, distance=5)
    p0 = []
    
    # Sécurisation de l'accès par position via Numpy
    v_arr = vitesse_mesure.values if hasattr(vitesse_mesure, 'values') else vitesse_mesure
    t_arr = temperature_corrigee.values if hasattr(temperature_corrigee, 'values') else temperature_corrigee

    for peak_idx in peaks:
        v_peak = v_arr[peak_idx]
        # Filtrage 1 : on ne garde que les pics situés dans la fenêtre estimée/déclarée
        if mask_min <= v_peak <= mask_max:
            # amplitude, centre, largeur (sigma) initiale de 5 km/s
            p0.extend([t_arr[peak_idx], v_peak, 5.0]) 
    
    y_fit = np.zeros_like(vitesse_mesure)
    main_peak_center = 0
    popt = []
    
    if len(p0) > 0:
        try:
            # Augmentation du maxfev pour laisser plus de temps à l'algo de converger
            popt_brut, _ = curve_fit(multi_gaussian, vitesse_mesure, temperature_corrigee, p0=p0, maxfev=10000)
            
            # Filtrage 2 : on supprime les gaussiennes qui auraient "dérivé" hors de la fenêtre pendant le calcul
            for i in range(len(popt_brut)//3):
                centre_fit = popt_brut[i*3 + 1]
                if mask_min <= centre_fit <= mask_max:
                    popt.extend(popt_brut[i*3 : i*3+3])
                    
            if len(popt) > 0:
                y_fit = multi_gaussian(vitesse_mesure, *popt)
                print(f"Ajustement réussi : {len(popt)//3} gaussienne(s) valide(s) détectée(s) dans la fenêtre.")
            else:
                print("Aucune gaussienne valide n'a été conservée dans la fenêtre après ajustement.")
        except Exception as e:
            print(f"Erreur lors de l'ajustement gaussien : {e}")
    else:
        print("Aucun pic significatif détecté dans la fenêtre pour l'ajustement.")

    # --- Extraction des sommets réels de la courbe ajustée (y_fit) ---
    final_peaks_v = []
    if len(popt) > 0:
        # On cherche les maxima locaux directement sur la courbe mathématique générée
        peaks_fit, _ = find_peaks(y_fit)
        v_arr = vitesse_mesure.values if hasattr(vitesse_mesure, 'values') else vitesse_mesure
        final_peaks_v = [v_arr[idx] for idx in peaks_fit]

    # Supression éventuelle de certains pics
    if args.remove:
        indices_to_remove = sorted([int(x) - 1 for x in args.remove.split(',')], reverse=True)
        for idx in indices_to_remove:
            if 0 <= idx < len(final_peaks_v):
                del final_peaks_v[idx : idx + 1]
                print(f"Pic n°{idx+1} supprimé de la liste finale.")

    # ajout éventuel de nouveaux pics
    if args.add:
        peaks_to_add = sorted([float(x) for x in args.add.split(',')], reverse=False)
        for peak in peaks_to_add:
            final_peaks_v.append(peak)

    # --- 2. Traitement de la référence LAB ---
    ref_file = f"reference_{galac_lon:03.0f}_{galac_lat:.0f}.txt"
    if not os.path.exists(ref_file):
        raise SystemExit(f"Spectre de référence {ref_file} non disponible")
    vitesse_ref, temperature_ref = [], []
    with open(ref_file, 'r') as f:
        in_lab = False
        for line in f:
            if line.startswith('%%LAB'): 
                in_lab = True
            elif line.startswith('%%'): 
                in_lab = False
            elif in_lab and not line.startswith('%'):
                parts = line.split()
                if len(parts) == 2:
                    vitesse_ref.append(float(parts[0]))
                    temperature_ref.append(float(parts[1]))
    
    # --- Calculs pour la console ---
    ref_peak = vitesse_ref[np.argmax(temperature_ref)]
    main_peak_center = vitesse_mesure[np.argmax(temperature_corrigee)]
    
    # Calcul du facteur d'étalonnage par intégration (aire sous la courbe)
    v_ref_arr = np.array(vitesse_ref)
    t_ref_arr = np.array(temperature_ref)
    
    # On restreint l'intégration à la fenêtre contenant le signal (pour ignorer le bruit lointain)
    mask_int_mes = (vitesse_mesure >= mask_min) & (vitesse_mesure <= mask_max)
    mask_int_ref = (v_ref_arr >= mask_min) & (v_ref_arr <= mask_max)
    
    # np.trapezoid calcule l'intégrale
    # On utilise abs() car selon l'ordre de l'axe X (croissant ou décroissant), l'intégrale peut être négative.
    integrale_mes = abs(np.trapezoid(temperature_corrigee[mask_int_mes], vitesse_mesure[mask_int_mes]))
    integrale_ref = abs(np.trapezoid(t_ref_arr[mask_int_ref], v_ref_arr[mask_int_ref]))
    
    facteur_etalonnage = integrale_ref / integrale_mes
    
    print(f"Fenêtre de signal estimée : [{est_min:.1f}, {est_max:.1f}] km/s")
    print(f"Sommet Référence: {ref_peak:.2f} km/s")
    print(f"Sommet Mesure: {main_peak_center:.2f} km/s")
    print(f"Offset suggéré pour alignement: {(ref_peak - main_peak_center):.2f} km/s")
    print(f"Facteur d'étalonnage global : x{facteur_etalonnage:.2f}")
    
    for i, v_c in enumerate(final_peaks_v):
        dist = get_galactic_distances(v_c, l_rad)

        V_reelle = None
        if abs(np.sin(l_rad)) >= 0.05:
            denom = (v_c / (V0 * np.sin(l_rad))) + 1
            if denom > 0:
                R = R0 / denom
                V_reelle = R * ((v_c / (R0 * np.sin(l_rad))) + (V0 / R0))
                        
        # Gestion dynamique de l'affichage des distances
        if dist:
            if len(dist) == 1:
                d_str = f"{dist[0]:.2f} kpc"
            else:
                d_str = f"{dist[0]:.2f} & {dist[1]:.2f} kpc"
        else:
            d_str = "Hors modèle"

        if V_reelle is not None:
            print(f"Pic {i+1} : V_lsr={v_c:.1f} km/s -> Distances : {d_str} | V_orb={V_reelle:.1f} km/s")
        else:
            print(f"Pic {i+1} : V_lsr={v_c:.1f} km/s -> Distances : {d_str}")            

    # --- Calcul du point tangent (Quadrants I et IV) ---
    V_tangent = None
    R_tangent = None
    v_ext = None
    if len(final_peaks_v) > 0 and quadrant in ['I', 'IV']:
        v_centers = final_peaks_v
        
        if quadrant == 'I':
            v_ext = max(v_centers)
            V_tangent = v_ext + V0 * np.sin(l_rad)
            R_tangent = R0 * np.sin(l_rad)
            print(f"Point tangent: quadrant I (Vmax = {v_ext:.1f} km/s) -> R = {R_tangent:.2f} kpc | V(R) = {V_tangent:.1f} km/s")
            
        elif quadrant == 'IV':
            v_ext = min(v_centers)
            V_tangent = -(v_ext + V0 * np.sin(l_rad))
            R_tangent = R0 * abs(np.sin(l_rad))
            print(f"Point tangent: quadrant IV (Vmin = {v_ext:.1f} km/s) -> R = {R_tangent:.2f} kpc | V(R) = {V_tangent:.1f} km/s")

    # --- Sauvegarde des résultats pour traitement ultérieur ---
    quadrant_dict = {'I': 1, 'II': 2, 'III': 3, 'IV': 4}
    quad_int = quadrant_dict.get(quadrant, 0)
    
    rows_to_save = []
    if len(final_peaks_v) == 0:
        # Aucun pic trouvé : on sauvegarde l'observation à vide
        rows_to_save.append({
            'Longitude': round(galac_lon, 2),
            'Quadrant': quad_int,
            'R_tangent_kpc': '',
            'V_orb_tangent': '',
            'R_kpc': '',
            'V_lsr': '',
            'Dist_1_kpc': '',
            'Dist_2_kpc': ''
        })
    else:
        # Un ou plusieurs pics : on crée une ligne par pic
        for v_c in final_peaks_v:
            
            # On vérifie si le pic en cours de traitement est le pic "extrême" (pour la courbe de rotation)
            is_tangent = (v_ext is not None) and (v_c == v_ext)
            
            # Calcul classique avec équation cinématique (maintient les nuages dans leurs bras)
            dist = get_galactic_distances(v_c, l_rad)
            
            if dist is not None:
                # Le nuage a une (ou deux) solution(s) géométrique(s) valide(s)
                d1 = round(dist[0], 2) if len(dist) > 0 else ''
                d2 = round(dist[1], 2) if len(dist) > 1 else ''
            else:
                # Hors modèle (vitesse dépassant la limite théorique due à la turbulence)
                if is_tangent:
                    # On "sauve" ce pic en le plaçant au point tangent géométrique
                    d_tangent = R0 * np.cos(l_rad)
                    d1 = round(d_tangent, 2)
                    d2 = ''
                else:
                    d1 = ''
                    d2 = ''

            # 2. On affecte les métadonnées tangentes UNIQUEMENT à la ligne de ce pic 
            # (Ces colonnes ne servent qu'à la courbe de rotation, pas à la carte)
            if is_tangent:
                val_R_tangent = round(R_tangent, 2) if R_tangent is not None else ''
                val_V_tangent = round(V_tangent, 2) if V_tangent is not None else ''
            else:
                val_R_tangent = ''
                val_V_tangent = ''

            # Calcul du Rayon R classique pour le nuage
            R_calc = None
            if abs(np.sin(l_rad)) >= 0.05:
                denom = (v_c / (V0 * np.sin(l_rad))) + 1
                if denom > 0:
                    R_calc = R0 / denom            
                    
            rows_to_save.append({
                'Longitude': round(galac_lon, 2),
                'Quadrant': quad_int,
                'R_tangent_kpc': val_R_tangent,
                'V_orb_tangent': val_V_tangent,
                'R_kpc': round(R_calc, 2) if R_calc is not None else '',
                'V_lsr': round(v_c, 2),
                'Dist_1_kpc': d1,
                'Dist_2_kpc': d2
            })
                        
    # Écriture en mode "append" (ajout à la fin) avec point-virgule comme séparateur
    df_save = pd.DataFrame(rows_to_save)
    if not args.blank:
        csv_filename = "synthese_observations.csv"
        file_exists = os.path.isfile(csv_filename)
        df_save.to_csv(csv_filename, mode='a', index=False, header=not file_exists, sep=';')
        print(f"[+] Données sauvegardées dans le fichier : {csv_filename}")
    print(f"Données extraites:")
    print(df_save.to_string(index=False))
        
    # --- Affichage du spectre ---
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    # Axe Gauche : Référence
    ax1.set_xlabel('Vitesse $V_{LSR}$ (km/s)', fontsize=13, fontweight='normal')
    ax1.set_ylabel('Température (K) - Référence LAB', color='tab:orange', fontsize=13)
    ligne_ref, = ax1.plot(vitesse_ref, temperature_ref, '-', color='tab:orange', linewidth=1.3, label='Référence (LAB)')
    ax1.tick_params(axis='y', labelcolor='tab:orange')
    ax1.set_xlim(vitesse_mesure.min(), vitesse_mesure.max())
    
    # Grille et ticks secondaires
    ax1.xaxis.set_minor_locator(MultipleLocator(10))
    ax1.xaxis.set_minor_formatter(NullFormatter())
    ax1.minorticks_on()
    ax1.grid(which='minor', linestyle=':', alpha=0.6)
    plt.grid(True, linestyle='--', alpha=0.5)

    # Axe Droit : Mesure
    ax2 = ax1.twinx()
    ligne_mesure, = ax2.plot(vitesse_mesure, temperature_corrigee, '-', color='tab:blue', linewidth=1.3, label='Mesure (Corrigée)')
    ligne_fit, = ax2.plot(vitesse_mesure, y_fit, '--', color='green', linewidth=1.0, label='Gaussiennes')
    ax2.set_ylabel('Intensité (K) - Mesure', color='tab:blue', fontsize=13)
    ax2.tick_params(axis='y', labelcolor='tab:blue')

    # Alignement proportionnel des zéros
    y1_min, y1_max = ax1.get_ylim()
    y2_min, y2_max = ax2.get_ylim()
    z1 = -y1_min / (y1_max - y1_min) if (y1_max - y1_min) != 0 else 0 
    z2 = -y2_min / (y2_max - y2_min) if (y2_max - y2_min) != 0 else 0
    z_max = max(z1, z2)
    
    if z_max < 1:
        ax1.set_ylim(-y1_max * (z_max / (1 - z_max)), y1_max)
        ax2.set_ylim(-y2_max * (z_max / (1 - z_max)), y2_max)

    ax2.axhline(0, color='black', linewidth=0.8, alpha=0.5)
    
    # Annotation des pics
    delta_temp = max(temperature_corrigee) - min(temperature_corrigee)
    if len(final_peaks_v) > 0:
        for v_c in final_peaks_v:
            # Ligne verticale courte pour marquer le pic
            ax2.vlines(v_c, 0, delta_temp/35, color='green', ls='-')
            # Texte du pic
            ax2.text(v_c, delta_temp/27, f'{v_c:.1f}', 
                     rotation=90, horizontalalignment='center', verticalalignment='bottom')
            
    # Légende et finitions
    legend_pos = 'upper right'
    if args.legend:
        if args.legend == 'ul':
            legend_pos = 'upper left'
        elif args.legend == 'lr':
            legend_pos = 'lower right'
        elif args.legend == 'll':
            legend_pos = 'lower left'
        elif args.legend == 'cr':
            legend_pos = 'center right'
        elif args.legend == 'cl':
            legend_pos = 'center left'
    ax1.legend([ligne_ref, ligne_mesure, ligne_fit], ['Référence (LAB)', 'Mesure (corrigée)', 'Gaussiennes'], loc=legend_pos)
    plt.title(f'Spectre Radio - Hydrogène Neutre (21 cm) - l={galac_lon:.1f}° b={galac_lat:.1f}°', fontsize=14, fontweight='bold')
    
    # Ajout des métadonnées et facteur d'étalonnage
    meta_str = f"{location} / {utc_datetime} / {exposure}s"
    ax1.text(0.02, 0.95, meta_str, transform=ax1.transAxes, fontsize=11, fontstyle='italic', va='bottom')
    bbox_props = dict(boxstyle="round,pad=0.4", fc="white", ec="gray", lw=1, alpha=0.9)
    ax1.text(0.02, 0.92, f"Facteur d'étalonnage global : x{facteur_etalonnage:.2f}", 
             transform=ax1.transAxes, fontsize=10, va='top', bbox=bbox_props)
    
    plt.tight_layout()
    
    # Sauvegarde du spectre
    plt.savefig(f'./{basename}.png')
    
    if args.display:
        plt.show()

if __name__ == "__main__":
    main()