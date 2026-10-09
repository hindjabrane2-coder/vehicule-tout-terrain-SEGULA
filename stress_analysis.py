"""
Analyse des contraintes sur chassis tubulaire -- vehicule tout-terrain.
Verification Von Mises, fatigue, flambement, concentration de contraintes,
optimisation de masse et assemblage par elements finis simplifie.
Projet SEGULA Technologies.
"""
import numpy as np
import matplotlib.pyplot as plt

# Materiaux disponibles
MATERIAUX = {
    'S355': {
        'E': 210e3,           # module d'Young (MPa)
        'nu': 0.3,
        'sigma_y': 355,       # limite elastique (MPa)
        'sigma_u': 510,       # resistance a la traction (MPa)
        'sigma_f': 240,       # limite d'endurance en flexion alternee (MPa)
        'rho': 7850,          # masse volumique (kg/m^3)
        'alpha': 12e-6,       # coefficient de dilatation (1/K)
    },
    'E355': {
        'E': 210e3,
        'nu': 0.3,
        'sigma_y': 355,
        'sigma_u': 580,
        'sigma_f': 260,
        'rho': 7850,
        'alpha': 12e-6,
    },
    '4130_chromoly': {
        'E': 205e3,
        'nu': 0.29,
        'sigma_y': 460,
        'sigma_u': 670,
        'sigma_f': 310,
        'rho': 7850,
        'alpha': 11.2e-6,
    },
}

# Materiau par defaut
MAT = MATERIAUX['S355']

# Sections de tubes candidates
TUBES = [
    {'label': '42.4x2.0', 'd_ext': 42.4, 'ep': 2.0},
    {'label': '42.4x2.6', 'd_ext': 42.4, 'ep': 2.6},
    {'label': '48.3x2.6', 'd_ext': 48.3, 'ep': 2.6},
    {'label': '48.3x3.2', 'd_ext': 48.3, 'ep': 3.2},
    {'label': '50.8x2.6', 'd_ext': 50.8, 'ep': 2.6},
    {'label': '50.8x3.2', 'd_ext': 50.8, 'ep': 3.2},
    {'label': '60.3x3.2', 'd_ext': 60.3, 'ep': 3.2},
]

# Parametres du vehicule
MASSE = 450           # kg
EMPATTEMENT = 2200    # mm
VOIE = 1400           # mm
H_CDG = 500           # mm

# Facteurs de concentration de contraintes pour les soudures
KT_SOUDURE_BOUT = 1.5    # soudure bout a bout correcte
KT_SOUDURE_ANGLE = 2.0   # soudure d'angle non meulée
KT_SOUDURE_MEULE = 1.3   # soudure meulée et controlee
# Proprietes geometriques

def section_tube(d_ext, ep):
    """Proprietes geometriques d'un tube circulaire."""
    d_int = d_ext - 2 * ep
    aire = np.pi / 4 * (d_ext**2 - d_int**2)
    inertie = np.pi / 64 * (d_ext**4 - d_int**4)
    module_flexion = inertie / (d_ext / 2)
    j_polaire = 2 * inertie  # moment d'inertie polaire (tube mince)
    rayon_giration = np.sqrt(inertie / aire)
    masse_lin = aire * 7850 / 1e6  # kg/mm
    return {
        'aire': aire,
        'inertie': inertie,
        'module': module_flexion,
        'j_polaire': j_polaire,
        'r_giration': rayon_giration,
        'd_ext': d_ext,
        'd_int': d_int,
        'ep': ep,
        'masse_lin': masse_lin,
    }
def section_plaque(largeur, epaisseur):
    """Proprietes d'une plaque rectangulaire (gousset, renfort)."""
    aire = largeur * epaisseur
    inertie = largeur * epaisseur**3 / 12
    module_flexion = inertie / (epaisseur / 2)
    return {
        'aire': aire,
        'inertie': inertie,
        'module': module_flexion,
        'largeur': largeur,
        'epaisseur': epaisseur,
    }
# Contraintes elementaires

def contrainte_traction(force, section):
    """Contrainte normale en traction/compression (MPa)."""
    return force / section['aire']
def contrainte_flexion(moment, section):
    """Contrainte maximale en flexion (MPa). Moment en N.mm."""
    return abs(moment) / section['module']
def contrainte_torsion(couple, section):
    """Contrainte de cisaillement maximale en torsion (MPa)."""
    return abs(couple) * (section['d_ext'] / 2) / section['j_polaire']
def von_mises(sigma_x, sigma_y, tau_xy):
    """Contrainte equivalente de Von Mises (etat plan de contrainte)."""
    return np.sqrt(sigma_x**2 - sigma_x * sigma_y + sigma_y**2 + 3 * tau_xy**2)
def von_mises_simple(sigma, tau):
    """Von Mises simplifie : traction + cisaillement."""
    return np.sqrt(sigma**2 + 3 * tau**2)
def coefficient_securite(sigma_eq, sigma_lim):
    """Coefficient de securite statique."""
    if sigma_eq < 1e-10:
        return float('inf')
    return sigma_lim / sigma_eq
# Cas de charge

def cas_statique():
    """Poids reparti sur le chassis -- poutre sur deux appuis."""
    poids = MASSE * 9.81
    q = poids / EMPATTEMENT
    m_flex = q * EMPATTEMENT**2 / 8
    return {'nom': 'Statique', 'force_axiale': 0, 'moment': m_flex, 'couple': 0}
def cas_virage(g_lateral=1.5):
    """Virage a acceleration laterale donnee."""
    poids = MASSE * 9.81
    f_lat = MASSE * g_lateral * 9.81
    m_roulis = f_lat * H_CDG
    m_flex = poids * EMPATTEMENT / 8
    return {'nom': f'Virage {g_lateral}g', 'force_axiale': 0, 'moment': m_flex, 'couple': m_roulis}
def cas_freinage(g_long=1.0):
    """Freinage -- transfert de charge longitudinal."""
    poids = MASSE * 9.81
    f_frein = MASSE * g_long * 9.81
    m_flex = poids * EMPATTEMENT / 8 + f_frein * H_CDG / 2
    return {'nom': f'Freinage {g_long}g', 'force_axiale': f_frein / 2, 'moment': m_flex, 'couple': 0}
def cas_torsion():
    """Torsion diagonale -- une roue en l'air."""
    poids = MASSE * 9.81
    couple_torsion = poids / 4 * VOIE / 2
    return {'nom': 'Torsion diag.', 'force_axiale': 0, 'moment': 0, 'couple': couple_torsion}
def cas_choc(facteur=3.0):
    """Choc vertical avec facteur multiplicatif."""
    stat = cas_statique()
    return {
        'nom': f'Choc x{facteur:.0f}',
        'force_axiale': 0,
        'moment': stat['moment'] * facteur,
        'couple': 0,
    }
def cas_combine():
    """Combinaison freinage + virage (racine de la somme des carres)."""
    frein = cas_freinage(0.8)
    virage = cas_virage(1.0)
    m_comb = np.sqrt(frein['moment']**2 + virage['moment']**2)
    c_comb = virage['couple']
    f_comb = frein['force_axiale']
    return {'nom': 'Combine fr+vir', 'force_axiale': f_comb, 'moment': m_comb, 'couple': c_comb}
TOUS_CAS = [cas_statique(), cas_virage(), cas_freinage(), cas_torsion(),
            cas_choc(), cas_combine()]
# Analyse d'un tube sous un cas de charge

def analyser_tube(sec, cas_charge, mat, kt=1.0):
    """
    Calcul complet des contraintes pour un tube donne sous un cas de charge.
    kt : facteur de concentration (soudure).
    """
    sigma_flex = contrainte_flexion(cas_charge['moment'], sec)
    sigma_axe = contrainte_traction(cas_charge['force_axiale'], sec)
    tau = contrainte_torsion(cas_charge['couple'], sec)

    sigma_total = (sigma_flex + abs(sigma_axe)) * kt
    tau_total = tau * kt
    vm = von_mises_simple(sigma_total, tau_total)
    cs = coefficient_securite(vm, mat['sigma_y'])

    return {
        'sigma_flex': sigma_flex,
        'sigma_axe': sigma_axe,
        'tau': tau,
        'sigma_total': sigma_total,
        'tau_total': tau_total,
        'von_mises': vm,
        'cs': cs,
    }
# Flambement (charge critique d'Euler)

def charge_critique_euler(sec, longueur, k_flamb=1.0, E=210e3):
    """
    Charge critique de flambement pour un poteau comprime.
    k_flamb : coefficient de longueur effective (0.5 a 2.0 selon les conditions d'appui).
    """
    l_eff = k_flamb * longueur
    p_crit = np.pi**2 * E * sec['inertie'] / l_eff**2
    return p_crit
def verifier_flambement(sec, force_compression, longueur, k_flamb=1.0, E=210e3):
    """Verifie si le tube flambe sous la charge donnee."""
    p_crit = charge_critique_euler(sec, longueur, k_flamb, E)
    cs = p_crit / abs(force_compression) if abs(force_compression) > 0 else float('inf')
    return {'p_crit': p_crit, 'cs_flamb': cs, 'longueur': longueur}
# Fatigue (diagramme de Goodman)

def goodman_sigma_admissible(sigma_m, sigma_u, sigma_f):
    """
    Contrainte alternee admissible selon le critere de Goodman.
    sigma_m : contrainte moyenne (MPa)
    sigma_u : resistance a la traction (MPa)
    sigma_f : limite d'endurance (MPa)
    """
    if sigma_m >= sigma_u:
        return 0.0
    return sigma_f * (1 - sigma_m / sigma_u)
def duree_de_vie_sn(sigma_a, sigma_f, n_coude=1e6, exposant=5):
    """
    Estimation de la duree de vie (nombre de cycles) a partir d'une loi de Basquin simplifiee.
    sigma_a : contrainte alternee (MPa)
    sigma_f : limite d'endurance (MPa)
    n_coude : nombre de cycles au coude de la courbe S-N
    """
    if sigma_a <= sigma_f:
        return float('inf')
    return n_coude * (sigma_f / sigma_a) ** exposant
def analyse_fatigue(sigma_max, sigma_min, mat):
    """Analyse fatigue complete pour un cycle de contrainte."""
    sigma_m = (sigma_max + sigma_min) / 2
    sigma_a = (sigma_max - sigma_min) / 2
    sigma_adm = goodman_sigma_admissible(sigma_m, mat['sigma_u'], mat['sigma_f'])
    n_cycles = duree_de_vie_sn(sigma_a, mat['sigma_f'])
    cs_fatigue = sigma_adm / sigma_a if sigma_a > 0 else float('inf')
    return {
        'sigma_m': sigma_m,
        'sigma_a': sigma_a,
        'sigma_adm_goodman': sigma_adm,
        'n_cycles': n_cycles,
        'cs_fatigue': cs_fatigue,
    }
# Contrainte thermique

def contrainte_thermique(delta_T, mat):
    """Contrainte thermique dans un element bride (dilatation empechee)."""
    return mat['E'] * mat['alpha'] * delta_T
# Gousset de renfort

def analyse_gousset(force, largeur_gousset, ep_gousset, longueur_soudure, mat):
    """
    Verification d'un gousset (plaque de renfort) soude au chassis.
    On verifie la contrainte dans la plaque et dans le cordon de soudure.
    """
    sec_plaque = section_plaque(largeur_gousset, ep_gousset)
    sigma_plaque = abs(force) / sec_plaque['aire']

    # cordon de soudure : gorge = 0.7 * epaisseur
    gorge = 0.7 * ep_gousset
    aire_soudure = gorge * longueur_soudure
    tau_soudure = abs(force) / aire_soudure if aire_soudure > 0 else float('inf')

    cs_plaque = mat['sigma_y'] / sigma_plaque if sigma_plaque > 0 else float('inf')
    # contrainte admissible soudure = 0.6 * sigma_y en cisaillement
    cs_soudure = (0.6 * mat['sigma_y']) / tau_soudure if tau_soudure > 0 else float('inf')

    return {
        'sigma_plaque': sigma_plaque,
        'tau_soudure': tau_soudure,
        'cs_plaque': cs_plaque,
        'cs_soudure': cs_soudure,
    }
# Matrice de rigidite poutre 2D (element fini simplifie)

def matrice_rigidite_poutre(E, I, A, L):
    """
    Matrice de rigidite 6x6 d'un element poutre 2D (Euler-Bernoulli).
    DDL par noeud : (u, v, theta) -- axial, transversal, rotation.
    """
    k_axial = E * A / L
    k_flex = E * I / L**3

    ke = np.zeros((6, 6))
    # rigidite axiale
    ke[0, 0] = k_axial
    ke[0, 3] = -k_axial
    ke[3, 0] = -k_axial
    ke[3, 3] = k_axial

    # rigidite en flexion
    ke[1, 1] = 12 * k_flex
    ke[1, 2] = 6 * k_flex * L
    ke[1, 4] = -12 * k_flex
    ke[1, 5] = 6 * k_flex * L

    ke[2, 1] = 6 * k_flex * L
    ke[2, 2] = 4 * k_flex * L**2
    ke[2, 4] = -6 * k_flex * L
    ke[2, 5] = 2 * k_flex * L**2

    ke[4, 1] = -12 * k_flex
    ke[4, 2] = -6 * k_flex * L
    ke[4, 4] = 12 * k_flex
    ke[4, 5] = -6 * k_flex * L

    ke[5, 1] = 6 * k_flex * L
    ke[5, 2] = 2 * k_flex * L**2
    ke[5, 4] = -6 * k_flex * L
    ke[5, 5] = 4 * k_flex * L**2

    return ke
def assembler_chassis_simplifie(sec, mat):
    """
    Assemblage d'un chassis simplifie a 3 poutres (5 noeuds) :
    - poutre avant (noeud 0 -> 1)
    - longeron (noeud 1 -> 2 -> 3)
    - poutre arriere (noeud 3 -> 4)
    Retourne la matrice globale et les correspondances.
    """
    E = mat['E']
    I = sec['inertie']
    A = sec['aire']

    # longueurs des elements
    L_traverse = VOIE / 2         # demi-traverse
    L_longeron = EMPATTEMENT / 2  # demi-longeron

    elements = [
        (0, 1, L_traverse),  # traverse avant
        (1, 2, L_longeron),  # longeron avant
        (2, 3, L_longeron),  # longeron arriere
        (3, 4, L_traverse),  # traverse arriere
    ]

    n_noeuds = 5
    n_ddl = 3 * n_noeuds
    K_glob = np.zeros((n_ddl, n_ddl))

    for n1, n2, L in elements:
        ke = matrice_rigidite_poutre(E, I, A, L)
        ddl = [3*n1, 3*n1+1, 3*n1+2, 3*n2, 3*n2+1, 3*n2+2]
        for i in range(6):
            for j in range(6):
                K_glob[ddl[i], ddl[j]] += ke[i, j]

    return K_glob, elements
def resoudre_chassis(sec, mat, forces_noeuds):
    """
    Resolution du systeme simplifie.
    forces_noeuds : dict {noeud: (Fx, Fy, Mz)}
    Conditions limites : noeuds 0 et 4 encastres.
    """
    K, elements = assembler_chassis_simplifie(sec, mat)
    n_ddl = K.shape[0]

    F = np.zeros(n_ddl)
    for noeud, (fx, fy, mz) in forces_noeuds.items():
        F[3*noeud] = fx
        F[3*noeud + 1] = fy
        F[3*noeud + 2] = mz

    # conditions limites : noeuds 0 et 4 bloques
    ddl_bloques = [0, 1, 2, 12, 13, 14]
    ddl_libres = [i for i in range(n_ddl) if i not in ddl_bloques]

    K_red = K[np.ix_(ddl_libres, ddl_libres)]
    F_red = F[ddl_libres]

    U_red = np.linalg.solve(K_red, F_red)
    U = np.zeros(n_ddl)
    for i, ddl in enumerate(ddl_libres):
        U[ddl] = U_red[i]

    # reactions aux appuis
    R = K @ U - F
    return U, R
# Optimisation de masse

def optimiser_tube(cas_charges, mat, cs_min=2.0):
    """
    Trouve le tube le plus leger qui respecte le coefficient de securite minimal
    sous tous les cas de charge.
    """
    meilleur = None
    for tube_def in TUBES:
        sec = section_tube(tube_def['d_ext'], tube_def['ep'])
        cs_min_trouve = float('inf')
        for cas in cas_charges:
            res = analyser_tube(sec, cas, mat, kt=KT_SOUDURE_BOUT)
            cs_min_trouve = min(cs_min_trouve, res['cs'])
        if cs_min_trouve >= cs_min:
            if meilleur is None or sec['masse_lin'] < meilleur['masse_lin']:
                meilleur = {
                    'label': tube_def['label'],
                    'section': sec,
                    'cs_min': cs_min_trouve,
                    'masse_lin': sec['masse_lin'],
                }
    return meilleur
# Programme principal

if __name__ == '__main__':
    sec_ref = section_tube(48.3, 2.6)
    print(f"Tube de reference : 48.3x2.6 mm")
    print(f"  Aire          : {sec_ref['aire']:.1f} mm2")
    print(f"  Inertie       : {sec_ref['inertie']:.0f} mm4")
    print(f"  Module flexion: {sec_ref['module']:.0f} mm3")
    print(f"  Rayon giration: {sec_ref['r_giration']:.1f} mm")
    print(f"  Masse lineique: {sec_ref['masse_lin']*1000:.2f} g/mm")
    print()

    # Analyse sous tous les cas de charge
    print("=" * 70)
    print(f"{'Cas':<16} {'sigma_f':>8} {'sigma_a':>8} {'tau':>8} {'VM':>8} {'CS':>6}")
    print("=" * 70)
    resultats_tous = []
    for cas in TOUS_CAS:
        r = analyser_tube(sec_ref, cas, MAT, kt=KT_SOUDURE_BOUT)
        r['nom'] = cas['nom']
        resultats_tous.append(r)
        flag = " (!)" if r['cs'] < 1.5 else ""
        print(f"{cas['nom']:<16} {r['sigma_flex']:8.1f} {r['sigma_axe']:8.1f} "
              f"{r['tau']:8.1f} {r['von_mises']:8.1f} {r['cs']:6.2f}{flag}")
    print()

    # Comparaison des materiaux
    print("--- Comparaison materiaux (cas choc x3) ---")
    cas_ref = cas_choc()
    for nom_mat, props in MATERIAUX.items():
        r = analyser_tube(sec_ref, cas_ref, props, kt=KT_SOUDURE_BOUT)
        print(f"  {nom_mat:<16} VM={r['von_mises']:.1f} MPa  CS={r['cs']:.2f}")
    print()

    # Flambement
    print("--- Flambement (entretoise 600 mm, articule-articule) ---")
    longueur_entretoise = 600  # mm
    for tube_def in TUBES:
        sec_t = section_tube(tube_def['d_ext'], tube_def['ep'])
        flamb = verifier_flambement(sec_t, MASSE * 9.81, longueur_entretoise, k_flamb=1.0)
        print(f"  {tube_def['label']:<12}  Pcrit={flamb['p_crit']:.0f} N  CS={flamb['cs_flamb']:.1f}")
    print()

    # Fatigue
    print("--- Analyse fatigue (cycle statique / choc x2) ---")
    sigma_stat = analyser_tube(sec_ref, cas_statique(), MAT, kt=KT_SOUDURE_ANGLE)['von_mises']
    sigma_choc = analyser_tube(sec_ref, cas_choc(2.0), MAT, kt=KT_SOUDURE_ANGLE)['von_mises']
    fat = analyse_fatigue(sigma_choc, sigma_stat, MAT)
    print(f"  Sigma moyen    : {fat['sigma_m']:.1f} MPa")
    print(f"  Sigma alterne  : {fat['sigma_a']:.1f} MPa")
    print(f"  Goodman adm.   : {fat['sigma_adm_goodman']:.1f} MPa")
    print(f"  CS fatigue     : {fat['cs_fatigue']:.2f}")
    if fat['n_cycles'] == float('inf'):
        print(f"  Duree de vie   : endurance illimitee")
    else:
        print(f"  Duree de vie   : {fat['n_cycles']:.0f} cycles")
    print()

    # Contrainte thermique
    delta_T = 60  # variation de temperature en service (degC)
    sig_th = contrainte_thermique(delta_T, MAT)
    print(f"--- Contrainte thermique (deltaT={delta_T} degC) ---")
    print(f"  Sigma thermique : {sig_th:.1f} MPa")
    print()

    # Gousset de renfort
    force_gousset = MASSE * 9.81 * 1.5 / 4  # charge sur un noeud en virage
    gousset = analyse_gousset(force_gousset, 60, 4, 80, MAT)
    print("--- Gousset 60x4 mm, cordon 80 mm ---")
    print(f"  Sigma plaque  : {gousset['sigma_plaque']:.1f} MPa  CS={gousset['cs_plaque']:.2f}")
    print(f"  Tau soudure   : {gousset['tau_soudure']:.1f} MPa  CS={gousset['cs_soudure']:.2f}")
    print()

    # Resolution chassis simplifie
    print("--- Modele poutre simplifie (5 noeuds) ---")
    # charge au centre du chassis (noeud 2)
    poids_total = MASSE * 9.81
    U, R = resoudre_chassis(sec_ref, MAT, {2: (0, -poids_total / 2, 0)})
    print(f"  Fleche au centre : {U[3*2+1]:.3f} mm")
    print(f"  Rotation noeud 1 : {np.degrees(U[3*1+2]):.4f} deg")
    print()

    # Optimisation de masse
    opt = optimiser_tube(TOUS_CAS, MAT, cs_min=2.0)
    if opt:
        print(f"--- Tube optimal (CS >= 2.0) ---")
        print(f"  {opt['label']}  masse={opt['masse_lin']*1000:.2f} g/mm  CS_min={opt['cs_min']:.2f}")
    else:
        print("  Aucun tube ne satisfait CS >= 2.0 sous tous les cas")
        opt_15 = optimiser_tube(TOUS_CAS, MAT, cs_min=1.5)
        if opt_15:
            print(f"  Avec CS >= 1.5 : {opt_15['label']}  CS_min={opt_15['cs_min']:.2f}")
    print()

    # Figures
    fig = plt.figure(figsize=(16, 12))

    # 1 -- Von Mises par cas de charge
    ax1 = fig.add_subplot(2, 3, 1)
    noms = [r['nom'] for r in resultats_tous]
    vms = [r['von_mises'] for r in resultats_tous]
    couleurs_vm = ['#2196F3' if r['cs'] >= 2.0 else '#FF9800' if r['cs'] >= 1.5 else '#F44336'
                   for r in resultats_tous]
    ax1.barh(noms, vms, color=couleurs_vm)
    ax1.axvline(x=MAT['sigma_y'], color='r', linestyle='--', label=f"Re = {MAT['sigma_y']} MPa")
    ax1.set_xlabel('Von Mises (MPa)')
    ax1.set_title('Contraintes par cas de charge')
    ax1.legend(fontsize=8)
    ax1.grid(True, axis='x')

    # 2 -- Coefficients de securite (carte de couleur)
    ax2 = fig.add_subplot(2, 3, 2)
    tubes_labels = [t['label'] for t in TUBES]
    cas_noms = [c['nom'] for c in TOUS_CAS]
    cs_matrice = np.zeros((len(TUBES), len(TOUS_CAS)))
    for i, tube_def in enumerate(TUBES):
        sec_t = section_tube(tube_def['d_ext'], tube_def['ep'])
        for j, cas in enumerate(TOUS_CAS):
            r = analyser_tube(sec_t, cas, MAT, kt=KT_SOUDURE_BOUT)
            cs_matrice[i, j] = min(r['cs'], 5.0)  # plafonner l'affichage
    im = ax2.imshow(cs_matrice, aspect='auto', cmap='RdYlGn', vmin=0.5, vmax=4.0)
    ax2.set_xticks(range(len(cas_noms)))
    ax2.set_xticklabels(cas_noms, rotation=45, ha='right', fontsize=7)
    ax2.set_yticks(range(len(tubes_labels)))
    ax2.set_yticklabels(tubes_labels, fontsize=8)
    ax2.set_title('CS par tube et cas de charge')
    fig.colorbar(im, ax=ax2, label='Coeff. securite')

    # 3 -- Decomposition contraintes (tube de reference)
    ax3 = fig.add_subplot(2, 3, 3)
    x_pos = np.arange(len(TOUS_CAS))
    sig_f = [r['sigma_flex'] for r in resultats_tous]
    sig_a = [abs(r['sigma_axe']) for r in resultats_tous]
    tau_arr = [r['tau'] for r in resultats_tous]
    w = 0.25
    ax3.bar(x_pos - w, sig_f, w, label='Flexion', color='#42A5F5')
    ax3.bar(x_pos, sig_a, w, label='Axiale', color='#66BB6A')
    ax3.bar(x_pos + w, tau_arr, w, label='Torsion', color='#FFA726')
    ax3.set_xticks(x_pos)
    ax3.set_xticklabels(cas_noms, rotation=45, ha='right', fontsize=7)
    ax3.set_ylabel('Contrainte (MPa)')
    ax3.set_title('Decomposition des contraintes')
    ax3.legend(fontsize=8)
    ax3.grid(True, axis='y')

    # 4 -- Courbe S-N
    ax4 = fig.add_subplot(2, 3, 4)
    sigma_range = np.linspace(MAT['sigma_f'] * 0.8, MAT['sigma_u'] * 0.8, 200)
    n_vie = [duree_de_vie_sn(s, MAT['sigma_f']) for s in sigma_range]
    n_vie = np.array([n if n != float('inf') else 1e8 for n in n_vie])
    ax4.semilogy(sigma_range, n_vie, 'b-', linewidth=2)
    ax4.axhline(1e6, color='k', linestyle=':', alpha=0.3)
    ax4.axvline(MAT['sigma_f'], color='g', linestyle='--', label=f"Endurance {MAT['sigma_f']} MPa")
    ax4.set_xlabel('Contrainte alternee (MPa)')
    ax4.set_ylabel('Nombre de cycles')
    ax4.set_title('Courbe S-N (Basquin)')
    ax4.legend(fontsize=8)
    ax4.grid(True, which='both')

    # 5 -- Diagramme de Goodman
    ax5 = fig.add_subplot(2, 3, 5)
    sm_range = np.linspace(0, MAT['sigma_u'], 100)
    sa_adm = [goodman_sigma_admissible(sm, MAT['sigma_u'], MAT['sigma_f']) for sm in sm_range]
    ax5.plot(sm_range, sa_adm, 'b-', linewidth=2, label='Limite Goodman')
    ax5.fill_between(sm_range, 0, sa_adm, alpha=0.15, color='blue')
    ax5.plot(fat['sigma_m'], fat['sigma_a'], 'ro', markersize=8, label='Point de fonctionnement')
    ax5.set_xlabel('Contrainte moyenne (MPa)')
    ax5.set_ylabel('Contrainte alternee (MPa)')
    ax5.set_title('Diagramme de Goodman')
    ax5.legend(fontsize=8)
    ax5.grid(True)

    # 6 -- Masse lineique vs CS minimal
    ax6 = fig.add_subplot(2, 3, 6)
    for tube_def in TUBES:
        sec_t = section_tube(tube_def['d_ext'], tube_def['ep'])
        cs_min_t = float('inf')
        for cas in TOUS_CAS:
            r = analyser_tube(sec_t, cas, MAT, kt=KT_SOUDURE_BOUT)
            cs_min_t = min(cs_min_t, r['cs'])
        masse_g = sec_t['masse_lin'] * 1000
        couleur = 'green' if cs_min_t >= 2.0 else 'orange' if cs_min_t >= 1.5 else 'red'
        ax6.plot(masse_g, cs_min_t, 'o', color=couleur, markersize=10)
        ax6.annotate(tube_def['label'], (masse_g, cs_min_t),
                     textcoords='offset points', xytext=(5, 5), fontsize=7)
    ax6.axhline(2.0, color='g', linestyle='--', alpha=0.5, label='CS = 2.0')
    ax6.axhline(1.5, color='orange', linestyle='--', alpha=0.5, label='CS = 1.5')
    ax6.set_xlabel('Masse lineique (g/mm)')
    ax6.set_ylabel('CS minimal')
    ax6.set_title('Optimisation masse / securite')
    ax6.legend(fontsize=8)
    ax6.grid(True)

    fig.suptitle(f'Analyse structurelle chassis tubulaire -- SEGULA', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('stress_analysis.png', dpi=150)
    print("Figure sauvegardee : stress_analysis.png")
