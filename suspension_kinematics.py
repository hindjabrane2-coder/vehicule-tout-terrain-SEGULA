"""
Cinematique de la suspension a double triangulation -- vehicule tout-terrain.
Calcul du debattement, carrossage, variation de voie, centre instantane de rotation,
centre de roulis, anti-plongee, anti-cabrage, rapport de mouvement et bump steer.
Projet SEGULA Technologies.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import csv
import os

# Geometrie de la suspension (mm) -- plan transversal (x = lateral, y = vertical)
# Les cotes correspondent a un demi-train (cote droit vu de l'arriere).
CHASSIS_SUP_INT = np.array([150.0, 350.0])   # pivot interieur triangle superieur
CHASSIS_SUP_EXT = np.array([350.0, 360.0])   # rotule superieure (porte-fusee)
CHASSIS_INF_INT = np.array([120.0, 150.0])   # pivot interieur triangle inferieur
CHASSIS_INF_EXT = np.array([380.0, 160.0])   # rotule inferieure (porte-fusee)

L_SUP = np.linalg.norm(CHASSIS_SUP_EXT - CHASSIS_SUP_INT)
L_INF = np.linalg.norm(CHASSIS_INF_EXT - CHASSIS_INF_INT)
DIST_ROTULES_NOM = np.linalg.norm(CHASSIS_SUP_EXT - CHASSIS_INF_EXT)

# Decalage du moyeu par rapport au milieu des deux rotules
MOYEU_OFFSET = np.array([30.0, -80.0])

# Point d'attache ressort sur le triangle inferieur (pour calcul du motion ratio)
RESSORT_ATTACHE = np.array([250.0, 155.0])
RESSORT_ANCRAGE_CHASSIS = np.array([250.0, 400.0])

# Point de biellette de direction sur le porte-fusee (pour calcul bump steer)
BIELLETTE_ROTULE = np.array([370.0, 200.0])
BIELLETTE_CREMAILLERE = np.array([200.0, 210.0])

# Demi-voie nominale
DEMI_VOIE = 700.0  # mm

# Hauteur du centre de gravite
H_CDG = 500.0  # mm

# Empattement (pour anti-dive/anti-squat)
EMPATTEMENT = 2200.0  # mm
DIST_AXE_AVANT = 1100.0  # distance du CdG a l'axe avant

# Geometrie vue de cote pour l'anti-plongee (plan longitudinal)
# Points projetes sur le plan lateral (x = longitudinal, y = vertical)
PIVOT_INF_AVANT = np.array([0.0, 150.0])
PIVOT_INF_ARRIERE = np.array([350.0, 155.0])
PIVOT_SUP_AVANT = np.array([50.0, 350.0])
PIVOT_SUP_ARRIERE = np.array([300.0, 345.0])
# Fonctions geometriques de base

def position_rotule(pivot, longueur, angle):
    """Position de la rotule exterieure pour un angle de bras donne."""
    return pivot + longueur * np.array([np.cos(angle), np.sin(angle)])
def angle_nominal(pivot_int, pivot_ext):
    """Angle nominal d'un bras a partir de ses deux extremites."""
    delta = pivot_ext - pivot_int
    return np.arctan2(delta[1], delta[0])
def intersection_droites(p1, d1, p2, d2):
    """
    Intersection de deux droites definies par un point et un vecteur directeur.
    Retourne None si les droites sont paralleles.
    """
    det = d1[0] * d2[1] - d1[1] * d2[0]
    if abs(det) < 1e-10:
        return None
    dp = p2 - p1
    t = (dp[0] * d2[1] - dp[1] * d2[0]) / det
    return p1 + t * d1
# Resolution cinematique principale

def resoudre_suspension(dz, sup_int=None, sup_ext=None, inf_int=None, inf_ext=None):
    """
    Resout la cinematique pour un deplacement vertical dz (mm).
    dz > 0 = compression (bosse), dz < 0 = detente (creux).
    On peut passer des pivots personnalises pour la comparaison de geometries.
    Retourne un dict avec toutes les grandeurs calculees.
    """
    si = sup_int if sup_int is not None else CHASSIS_SUP_INT
    se = sup_ext if sup_ext is not None else CHASSIS_SUP_EXT
    ii = inf_int if inf_int is not None else CHASSIS_INF_INT
    ie = inf_ext if inf_ext is not None else CHASSIS_INF_EXT

    l_sup = np.linalg.norm(se - si)
    l_inf = np.linalg.norm(ie - ii)
    dist_rot = np.linalg.norm(se - ie)

    a_sup_0 = angle_nominal(si, se)
    a_inf_0 = angle_nominal(ii, ie)

    milieu_nom = (se + ie) / 2
    moyeu_nom = milieu_nom + MOYEU_OFFSET

    def erreur(d_angle_inf):
        rot_inf = position_rotule(ii, l_inf, a_inf_0 + d_angle_inf)

        def err_sup(d_angle_sup):
            rot_sup = position_rotule(si, l_sup, a_sup_0 + d_angle_sup)
            return np.linalg.norm(rot_sup - rot_inf) - dist_rot

        try:
            das = brentq(err_sup, -0.5, 0.5)
        except ValueError:
            return 1e6

        rot_sup = position_rotule(si, l_sup, a_sup_0 + das)
        milieu = (rot_sup + rot_inf) / 2
        moyeu = milieu + MOYEU_OFFSET
        return moyeu[1] - moyeu_nom[1] - dz

    try:
        d_angle = brentq(erreur, -0.4, 0.4)
    except ValueError:
        return None

    # recalcul final avec l'angle trouve
    rot_inf = position_rotule(ii, l_inf, a_inf_0 + d_angle)

    def err_sup(d_angle_sup):
        rot_sup = position_rotule(si, l_sup, a_sup_0 + d_angle_sup)
        return np.linalg.norm(rot_sup - rot_inf) - dist_rot

    d_angle_sup = brentq(err_sup, -0.5, 0.5)
    rot_sup = position_rotule(si, l_sup, a_sup_0 + d_angle_sup)

    milieu = (rot_sup + rot_inf) / 2
    moyeu = milieu + MOYEU_OFFSET

    # carrossage : angle de la fusee (droite rotule inf -> rotule sup) par rapport a la verticale
    fusee = rot_sup - rot_inf
    carrossage = np.degrees(np.arctan2(fusee[0], fusee[1]))

    # variation de voie laterale
    delta_voie = moyeu[0] - moyeu_nom[0]

    # centre instantane de rotation (intersection des prolongements des deux bras)
    dir_sup = rot_sup - si
    dir_inf = rot_inf - ii
    ci = intersection_droites(si, dir_sup, ii, dir_inf)

    # centre de roulis = intersection entre la droite CI-contact_pneu et l'axe du vehicule
    # le contact pneu est sous le moyeu, au niveau du sol (y = 0)
    contact_pneu = np.array([moyeu[0] + MOYEU_OFFSET[0], 0.0])
    centre_roulis_h = None
    if ci is not None:
        dir_cr = contact_pneu - ci
        if abs(dir_cr[0]) > 1e-6:
            t_cr = (0.0 - ci[0]) / dir_cr[0]
            centre_roulis_h = ci[1] + t_cr * dir_cr[1]

    return {
        'dz': dz,
        'moyeu': moyeu,
        'rot_sup': rot_sup,
        'rot_inf': rot_inf,
        'carrossage': carrossage,
        'delta_voie': delta_voie,
        'centre_instantane': ci,
        'centre_roulis_h': centre_roulis_h,
        'd_angle_inf': d_angle,
        'd_angle_sup': d_angle_sup,
    }
def calculer_motion_ratio(dz, pas=0.5):
    """
    Rapport de mouvement : course du ressort / course de la roue.
    On evalue le deplacement du point d'attache ressort sur le triangle inferieur.
    """
    a_inf_0 = angle_nominal(CHASSIS_INF_INT, CHASSIS_INF_EXT)
    bras_attache = np.linalg.norm(RESSORT_ATTACHE - CHASSIS_INF_INT)
    angle_attache = angle_nominal(CHASSIS_INF_INT, RESSORT_ATTACHE)
    offset_angle = angle_attache - a_inf_0

    res_ref = resoudre_suspension(dz)
    res_voisin = resoudre_suspension(dz + pas)
    if res_ref is None or res_voisin is None:
        return None

    # position de l'attache ressort pour les deux configurations
    pos_ref = position_rotule(CHASSIS_INF_INT, bras_attache, a_inf_0 + res_ref['d_angle_inf'] + offset_angle)
    pos_voisin = position_rotule(CHASSIS_INF_INT, bras_attache, a_inf_0 + res_voisin['d_angle_inf'] + offset_angle)

    # course ressort = variation de distance entre attache et ancrage chassis
    d_ref = np.linalg.norm(RESSORT_ANCRAGE_CHASSIS - pos_ref)
    d_voisin = np.linalg.norm(RESSORT_ANCRAGE_CHASSIS - pos_voisin)

    delta_ressort = d_voisin - d_ref
    delta_roue = pas
    if abs(delta_roue) < 1e-10:
        return 0.0
    return delta_ressort / delta_roue
def calculer_bump_steer(dz):
    """
    Variation de l'angle de braquage (toe) due au debattement.
    On suit le deplacement de la rotule de biellette de direction
    qui est solidaire du porte-fusee.
    """
    res = resoudre_suspension(dz)
    res_0 = resoudre_suspension(0.0)
    if res is None or res_0 is None:
        return None

    # la biellette de direction tourne avec la fusee
    # on calcule l'angle de rotation du porte-fusee par rapport a la position nominale
    fusee_0 = res_0['rot_sup'] - res_0['rot_inf']
    fusee = res['rot_sup'] - res['rot_inf']

    angle_0 = np.arctan2(fusee_0[1], fusee_0[0])
    angle = np.arctan2(fusee[1], fusee[0])
    rotation_fusee = angle - angle_0

    # position de la rotule de biellette (solidaire du porte-fusee)
    centre_fusee_0 = (res_0['rot_sup'] + res_0['rot_inf']) / 2
    centre_fusee = (res['rot_sup'] + res['rot_inf']) / 2
    offset_biel = BIELLETTE_ROTULE - centre_fusee_0

    # rotation de l'offset
    c, s = np.cos(rotation_fusee), np.sin(rotation_fusee)
    offset_rot = np.array([c * offset_biel[0] - s * offset_biel[1],
                           s * offset_biel[0] + c * offset_biel[1]])
    biel_pos = centre_fusee + offset_rot

    # angle de la biellette par rapport a la cremaillere
    dir_biel = biel_pos - BIELLETTE_CREMAILLERE
    dir_biel_0 = BIELLETTE_ROTULE - BIELLETTE_CREMAILLERE
    toe = np.degrees(np.arctan2(dir_biel[1], dir_biel[0]) -
                     np.arctan2(dir_biel_0[1], dir_biel_0[0]))
    return toe
def anti_plongee_pourcentage():
    """
    Pourcentage d'anti-plongee (anti-dive) pour le train avant.
    On prolonge les bras vus de cote pour trouver le point de convergence,
    puis on calcule la tangente de l'angle par rapport a l'horizontale.
    """
    dir_inf = PIVOT_INF_ARRIERE - PIVOT_INF_AVANT
    dir_sup = PIVOT_SUP_ARRIERE - PIVOT_SUP_AVANT

    ci_lateral = intersection_droites(PIVOT_INF_AVANT, dir_inf, PIVOT_SUP_AVANT, dir_sup)
    if ci_lateral is None:
        return 0.0

    # angle entre la droite CI-contact_sol et l'horizontale
    contact_sol = np.array([0.0, 0.0])
    dx = ci_lateral[0] - contact_sol[0]
    dy = ci_lateral[1] - contact_sol[1]
    if abs(dx) < 1e-6:
        return 100.0

    tan_theta = dy / dx
    # pourcentage anti-dive = tan(theta) * (empattement / h_cdg) * repartition_freinage
    repartition_av = DIST_AXE_AVANT / EMPATTEMENT  # fraction du freinage sur l'avant
    anti_dive = tan_theta * (EMPATTEMENT * repartition_av / H_CDG) * 100.0
    return anti_dive
def anti_cabrage_pourcentage():
    """
    Pourcentage d'anti-cabrage (anti-squat) pour le train arriere.
    Meme methode que l'anti-plongee mais avec la geometrie arriere inversee.
    """
    # pour le train arriere, on inverse la direction longitudinale
    dir_inf = PIVOT_INF_AVANT - PIVOT_INF_ARRIERE
    dir_sup = PIVOT_SUP_AVANT - PIVOT_SUP_ARRIERE

    ci_lat = intersection_droites(PIVOT_INF_ARRIERE, dir_inf, PIVOT_SUP_ARRIERE, dir_sup)
    if ci_lat is None:
        return 0.0

    dx = ci_lat[0]
    dy = ci_lat[1]
    if abs(dx) < 1e-6:
        return 100.0

    tan_theta = dy / abs(dx)
    repartition_ar = 1.0 - DIST_AXE_AVANT / EMPATTEMENT
    anti_squat = tan_theta * (EMPATTEMENT * repartition_ar / H_CDG) * 100.0
    return anti_squat
def derivees_cinematiques(debattements, carrossages, variations_voie):
    """
    Taux de variation du carrossage et de la voie en fonction du debattement.
    On utilise les differences finies centrees.
    """
    d_camber = np.gradient(carrossages, debattements)
    d_voie = np.gradient(variations_voie, debattements)
    return d_camber, d_voie
def balayage_cinematique(dz_min=-80, dz_max=80, n_points=100,
                         sup_int=None, sup_ext=None,
                         inf_int=None, inf_ext=None):
    """
    Balayage complet de la cinematique sur la plage de debattement.
    Retourne des tableaux numpy pour chaque grandeur.
    """
    debattements = np.linspace(dz_min, dz_max, n_points)
    carrossages = np.full(n_points, np.nan)
    variations_voie = np.full(n_points, np.nan)
    centres_roulis = np.full(n_points, np.nan)
    motion_ratios = np.full(n_points, np.nan)
    bump_steer_arr = np.full(n_points, np.nan)

    for i, dz in enumerate(debattements):
        res = resoudre_suspension(dz, sup_int, sup_ext, inf_int, inf_ext)
        if res is not None:
            carrossages[i] = res['carrossage']
            variations_voie[i] = res['delta_voie']
            if res['centre_roulis_h'] is not None:
                centres_roulis[i] = res['centre_roulis_h']

        mr = calculer_motion_ratio(dz)
        if mr is not None:
            motion_ratios[i] = mr

        bs = calculer_bump_steer(dz)
        if bs is not None:
            bump_steer_arr[i] = bs

    return {
        'debattements': debattements,
        'carrossages': carrossages,
        'variations_voie': variations_voie,
        'centres_roulis': centres_roulis,
        'motion_ratios': motion_ratios,
        'bump_steer': bump_steer_arr,
    }
def exporter_csv(donnees, fichier='cinematique_suspension.csv'):
    """Exporte les courbes cinematiques dans un fichier CSV."""
    with open(fichier, 'w', newline='') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(['Debattement_mm', 'Carrossage_deg', 'Var_voie_mm',
                         'Centre_roulis_mm', 'Motion_ratio', 'Bump_steer_deg'])
        for i in range(len(donnees['debattements'])):
            writer.writerow([
                f"{donnees['debattements'][i]:.2f}",
                f"{donnees['carrossages'][i]:.4f}",
                f"{donnees['variations_voie'][i]:.4f}",
                f"{donnees['centres_roulis'][i]:.2f}",
                f"{donnees['motion_ratios'][i]:.4f}",
                f"{donnees['bump_steer'][i]:.4f}",
            ])
    print(f"Export CSV : {fichier}")
def sensibilite_pivot(pivot_nom, index_coord, delta_mm=5.0):
    """
    Analyse de sensibilite : on deplace un pivot de +/- delta_mm
    et on mesure l'effet sur la courbe de carrossage.
    pivot_nom : 'sup_int', 'sup_ext', 'inf_int', 'inf_ext'
    index_coord : 0 pour x (lateral), 1 pour y (vertical)
    """
    pivots = {
        'sup_int': CHASSIS_SUP_INT.copy(),
        'sup_ext': CHASSIS_SUP_EXT.copy(),
        'inf_int': CHASSIS_INF_INT.copy(),
        'inf_ext': CHASSIS_INF_EXT.copy(),
    }

    resultats = {}
    for signe, label in [(-1, f'-{delta_mm}mm'), (0, 'nominal'), (1, f'+{delta_mm}mm')]:
        p = pivots.copy()
        p[pivot_nom] = pivots[pivot_nom].copy()
        p[pivot_nom][index_coord] += signe * delta_mm
        data = balayage_cinematique(
            sup_int=p['sup_int'], sup_ext=p['sup_ext'],
            inf_int=p['inf_int'], inf_ext=p['inf_ext']
        )
        resultats[label] = data
    return resultats
# Visualisation 3D simplifiee (projection sur 3 positions)

def dessiner_suspension(ax, dz, couleur='b', label_prefix=''):
    """Dessine la suspension dans une position donnee sur un axe matplotlib."""
    res = resoudre_suspension(dz)
    if res is None:
        return

    si = CHASSIS_SUP_INT
    ri_sup = res['rot_sup']
    ii = CHASSIS_INF_INT
    ri_inf = res['rot_inf']
    moyeu = res['moyeu']

    lbl = f"{label_prefix}dz={dz:+.0f}mm" if label_prefix else f"dz={dz:+.0f}mm"
    ax.plot([si[0], ri_sup[0]], [si[1], ri_sup[1]], color=couleur, linewidth=2, label=lbl)
    ax.plot([ii[0], ri_inf[0]], [ii[1], ri_inf[1]], color=couleur, linewidth=2)
    ax.plot([ri_sup[0], ri_inf[0]], [ri_sup[1], ri_inf[1]], color=couleur, linewidth=1.5, linestyle='--')
    ax.plot(*moyeu, 'o', color=couleur, markersize=6)

    # marquer le centre instantane
    if res['centre_instantane'] is not None:
        ci = res['centre_instantane']
        if -500 < ci[0] < 1000 and -200 < ci[1] < 800:
            ax.plot(*ci, 'x', color=couleur, markersize=8)
# Programme principal

if __name__ == '__main__':
    # calcul anti-plongee et anti-cabrage
    ad = anti_plongee_pourcentage()
    asq = anti_cabrage_pourcentage()
    print(f"Anti-plongee (avant) : {ad:.1f} %")
    print(f"Anti-cabrage (arriere) : {asq:.1f} %")
    print()

    # balayage cinematique complet
    donnees = balayage_cinematique()
    deb = donnees['debattements']
    cam = donnees['carrossages']
    dv = donnees['variations_voie']
    cr = donnees['centres_roulis']
    mr = donnees['motion_ratios']
    bs = donnees['bump_steer']

    # derivees
    d_cam, d_voie = derivees_cinematiques(deb, cam, dv)

    # export CSV
    exporter_csv(donnees)

    # valeurs cles
    idx_0 = np.argmin(np.abs(deb))
    idx_p80 = np.argmin(np.abs(deb - 80))
    idx_m80 = np.argmin(np.abs(deb + 80))
    print(f"Carrossage nominal         : {cam[idx_0]:.2f} deg")
    print(f"Carrossage a +80mm (bosse) : {cam[idx_p80]:.2f} deg")
    print(f"Carrossage a -80mm (creux) : {cam[idx_m80]:.2f} deg")
    print(f"Centre de roulis nominal   : {cr[idx_0]:.1f} mm")
    print(f"Motion ratio nominal       : {mr[idx_0]:.3f}")
    print(f"Bump steer a +50mm         : {bs[np.argmin(np.abs(deb - 50))]:.3f} deg")
    print()

    # Comparaison geometrie stock vs modifiee (bras superieur 10mm plus long)
    sup_ext_mod = CHASSIS_SUP_EXT + np.array([10.0, 0.0])
    donnees_mod = balayage_cinematique(sup_ext=sup_ext_mod)

    # Figures
    fig = plt.figure(figsize=(16, 14))

    # 1 -- carrossage
    ax1 = fig.add_subplot(3, 3, 1)
    ax1.plot(deb, cam, 'b-', linewidth=2, label='Stock')
    ax1.plot(deb, donnees_mod['carrossages'], 'g--', linewidth=1.5, label='Modifie')
    ax1.set_xlabel('Debattement (mm)')
    ax1.set_ylabel('Carrossage (deg)')
    ax1.set_title('Carrossage')
    ax1.axhline(0, color='k', linestyle=':', alpha=0.3)
    ax1.axvline(0, color='k', linestyle=':', alpha=0.3)
    ax1.legend(fontsize=8)
    ax1.grid(True)

    # 2 -- variation de voie
    ax2 = fig.add_subplot(3, 3, 2)
    ax2.plot(deb, dv, 'r-', linewidth=2, label='Stock')
    ax2.plot(deb, donnees_mod['variations_voie'], 'g--', linewidth=1.5, label='Modifie')
    ax2.set_xlabel('Debattement (mm)')
    ax2.set_ylabel('Variation voie (mm)')
    ax2.set_title('Variation de voie')
    ax2.axhline(0, color='k', linestyle=':', alpha=0.3)
    ax2.legend(fontsize=8)
    ax2.grid(True)

    # 3 -- centre de roulis
    ax3 = fig.add_subplot(3, 3, 3)
    ax3.plot(deb, cr, 'm-', linewidth=2)
    ax3.set_xlabel('Debattement (mm)')
    ax3.set_ylabel('Hauteur centre roulis (mm)')
    ax3.set_title('Centre de roulis')
    ax3.grid(True)

    # 4 -- motion ratio
    ax4 = fig.add_subplot(3, 3, 4)
    ax4.plot(deb, mr, 'c-', linewidth=2)
    ax4.set_xlabel('Debattement (mm)')
    ax4.set_ylabel('Motion ratio')
    ax4.set_title('Rapport de mouvement')
    ax4.axhline(np.nanmean(mr), color='k', linestyle='--', alpha=0.4, label=f'Moy={np.nanmean(mr):.3f}')
    ax4.legend(fontsize=8)
    ax4.grid(True)

    # 5 -- bump steer
    ax5 = fig.add_subplot(3, 3, 5)
    ax5.plot(deb, bs, color='orange', linewidth=2)
    ax5.set_xlabel('Debattement (mm)')
    ax5.set_ylabel('Angle de pincement (deg)')
    ax5.set_title('Bump steer')
    ax5.axhline(0, color='k', linestyle=':', alpha=0.3)
    ax5.grid(True)

    # 6 -- derivees
    ax6 = fig.add_subplot(3, 3, 6)
    ax6.plot(deb, d_cam, 'b-', linewidth=1.5, label='dCamber/dz (deg/mm)')
    ax6_b = ax6.twinx()
    ax6_b.plot(deb, d_voie, 'r--', linewidth=1.5, label='dVoie/dz (mm/mm)')
    ax6.set_xlabel('Debattement (mm)')
    ax6.set_ylabel('dCarrossage/dz', color='b')
    ax6_b.set_ylabel('dVoie/dz', color='r')
    ax6.set_title('Taux de variation')
    ax6.grid(True)

    # 7 -- dessin de la suspension a 3 positions
    ax7 = fig.add_subplot(3, 3, 7)
    dessiner_suspension(ax7, 0, couleur='k', label_prefix='')
    dessiner_suspension(ax7, 60, couleur='b', label_prefix='')
    dessiner_suspension(ax7, -60, couleur='r', label_prefix='')
    ax7.set_xlabel('x lateral (mm)')
    ax7.set_ylabel('y vertical (mm)')
    ax7.set_title('Geometrie (3 positions)')
    ax7.set_aspect('equal')
    ax7.legend(fontsize=7)
    ax7.grid(True)

    # 8 -- sensibilite du pivot sup_int en y
    ax8 = fig.add_subplot(3, 3, 8)
    sens = sensibilite_pivot('sup_int', 1, delta_mm=8.0)
    for label, data in sens.items():
        ax8.plot(data['debattements'], data['carrossages'], linewidth=1.5, label=label)
    ax8.set_xlabel('Debattement (mm)')
    ax8.set_ylabel('Carrossage (deg)')
    ax8.set_title('Sensibilite pivot sup. int. (y)')
    ax8.legend(fontsize=7)
    ax8.grid(True)

    # 9 -- synthese carrossage + voie
    ax9 = fig.add_subplot(3, 3, 9)
    ax9.plot(deb, cam, 'b-', linewidth=2, label='Carrossage')
    ax9_b = ax9.twinx()
    ax9_b.plot(deb, dv, 'r--', linewidth=1.5, label='Var. voie')
    ax9.set_xlabel('Debattement (mm)')
    ax9.set_ylabel('Carrossage (deg)', color='b')
    ax9_b.set_ylabel('Variation voie (mm)', color='r')
    ax9.set_title('Synthese')
    ax9.annotate(f'AD={ad:.0f}% / AS={asq:.0f}%', xy=(0.02, 0.95),
                 xycoords='axes fraction', fontsize=7, va='top',
                 bbox=dict(boxstyle='round', fc='wheat', alpha=0.5))
    ax9.grid(True)

    fig.suptitle('Cinematique suspension double triangulation -- SEGULA', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('suspension_kinematics.png', dpi=150)
    print("Figure sauvegardee : suspension_kinematics.png")
