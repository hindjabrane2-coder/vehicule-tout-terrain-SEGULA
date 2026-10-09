"""
Dimensionnement ressort, amortisseur et barre anti-roulis -- vehicule tout-terrain.
Calcul de la frequence propre, du taux d'amortissement, de la transmissibilite,
de l'acceleration RMS (confort ISO 2631) et du comportement dynamique.
Projet SEGULA Technologies.
"""
import numpy as np
import matplotlib.pyplot as plt

# Parametres du vehicule
MASSE_TOTALE = 450          # kg (vehicule + 2 passagers)
MASSE_SUSPENDUE = 350       # kg
MASSE_NON_SUSPENDUE = 25    # kg par roue (moyeu, triangle, frein)
MASSE_PAR_ROUE = MASSE_SUSPENDUE / 4
EMPATTEMENT = 2200          # mm
REPARTITION_AVANT = 0.52   # fraction de la masse sur l'essieu avant
VOIE_AVANT = 1400           # mm
VOIE_ARRIERE = 1380         # mm
H_CDG = 500                 # mm
DEBATTEMENT_COMP = 80       # mm (course en compression)
DEBATTEMENT_DET = 80        # mm (course en detente)
RATIO_LEVIER = 0.70         # rapport de mouvement (motion ratio)

# Frequences cibles
F_CIBLE_AV = 1.3            # Hz -- avant (confort)
F_CIBLE_AR = 1.5            # Hz -- arriere (un peu plus raide pour l'aplomb)

# Materiau du ressort : 50CrV4 (acier a ressort classique)
G_RESSORT = 78500           # module de cisaillement (MPa)
SIGMA_ADM = 800             # contrainte admissible statique (MPa)
SIGMA_END = 550             # limite d'endurance en cisaillement corrigee (MPa)
RHO_ACIER = 7850            # kg/m^3
# Raideur et frequence propre

def raideur_roue(masse, freq):
    """Raideur au niveau de la roue pour obtenir la frequence propre souhaitee (N/mm)."""
    return masse * (2 * np.pi * freq)**2 / 1000
def raideur_ressort(k_roue, ratio):
    """Raideur du ressort physique, en tenant compte du rapport de levier."""
    return k_roue / ratio**2
def frequence_propre(k_roue, masse):
    """Frequence propre non amortie (Hz)."""
    return np.sqrt(k_roue * 1000 / masse) / (2 * np.pi)
def frequences_tangage_roulis(k_av, k_ar, m_av, m_ar):
    """
    Frequences naturelles de tangage et de roulis (modele 2 DDL decouples).
    k en N/mm par roue, m en kg par essieu.
    """
    # roulis : les deux roues d'un meme essieu travaillent en opposition
    # on prend la frequence de rebond vertical comme approximation
    f_rebond_av = frequence_propre(k_av, m_av / 2)
    f_rebond_ar = frequence_propre(k_ar, m_ar / 2)

    # tangage : couplage avant/arriere
    k_tot = (2 * k_av + 2 * k_ar) * 1000  # N/m
    m_tot = m_av + m_ar
    inertie_tangage = m_av * (EMPATTEMENT * REPARTITION_AVANT / 1000)**2 + \
                      m_ar * (EMPATTEMENT * (1 - REPARTITION_AVANT) / 1000)**2

    # raideur de tangage (difference de raideur avant/arriere)
    l_av = EMPATTEMENT * REPARTITION_AVANT  # mm
    l_ar = EMPATTEMENT * (1 - REPARTITION_AVANT)
    k_tangage = 2 * k_av * (l_av / 1000)**2 + 2 * k_ar * (l_ar / 1000)**2  # N.m/rad (approx)
    f_tangage = np.sqrt(k_tangage * 1000 / inertie_tangage) / (2 * np.pi)

    return {
        'f_rebond_av': f_rebond_av,
        'f_rebond_ar': f_rebond_ar,
        'f_tangage': f_tangage,
    }
# Dimensionnement du ressort helicoldal

def facteur_wahl(C):
    """Facteur de correction de Wahl pour la contrainte dans le fil."""
    return (4 * C - 1) / (4 * C - 4) + 0.615 / C
def dimensionner_ressort(k_ressort, charge_stat, d_fil_init=10.0, C_cible=6.0):
    """
    Dimensionnement geometrique du ressort helicoldal.
    k_ressort : raideur du ressort (N/mm)
    charge_stat : charge statique sur le ressort (N)
    C_cible : rapport d'enroulement vise (D_moyen / d_fil)
    """
    d_fil = d_fil_init
    D_moyen = C_cible * d_fil
    C = C_cible
    Kw = facteur_wahl(C)

    # charge maximale = charge statique + compression complete
    course_ressort = DEBATTEMENT_COMP / RATIO_LEVIER
    F_max = charge_stat + k_ressort * course_ressort

    # verification de la contrainte sous F_max
    tau_max = Kw * 8 * F_max * D_moyen / (np.pi * d_fil**3)

    if tau_max > SIGMA_ADM:
        # recalculer le diametre de fil necessaire
        d_fil = (Kw * 8 * F_max * C_cible * d_fil / (np.pi * SIGMA_ADM))**(1/3)
        d_fil = np.ceil(d_fil * 2) / 2  # arrondi au demi-mm superieur
        D_moyen = C_cible * d_fil
        C = D_moyen / d_fil
        Kw = facteur_wahl(C)
        tau_max = Kw * 8 * F_max * D_moyen / (np.pi * d_fil**3)

    # nombre de spires actives
    n_actives = G_RESSORT * d_fil**4 / (8 * D_moyen**3 * k_ressort)

    # hauteurs
    n_total = n_actives + 2  # 2 spires mortes (extremites meulees)
    hauteur_bloc = n_total * d_fil  # hauteur a spires jointives
    fleche_stat = charge_stat / k_ressort
    fleche_max = F_max / k_ressort
    L0 = hauteur_bloc + fleche_max + 5.0  # 5 mm de marge avant clash

    # jeu au bloc
    jeu_bloc = L0 - fleche_max - hauteur_bloc

    # masse du ressort
    longueur_fil = np.pi * D_moyen * n_total
    volume_fil = np.pi / 4 * d_fil**2 * longueur_fil  # mm^3
    masse_ressort = volume_fil * RHO_ACIER / 1e9  # kg

    return {
        'd_fil': d_fil,
        'D_moyen': D_moyen,
        'D_ext': D_moyen + d_fil,
        'C': C,
        'Kw': Kw,
        'n_actives': n_actives,
        'n_total': n_total,
        'L0': L0,
        'hauteur_bloc': hauteur_bloc,
        'k': k_ressort,
        'tau_max': tau_max,
        'fleche_stat': fleche_stat,
        'fleche_max': fleche_max,
        'jeu_bloc': jeu_bloc,
        'masse': masse_ressort,
    }
def verifier_clash(ressort):
    """Verifie que le ressort ne va pas a bloc en compression maximale."""
    if ressort['jeu_bloc'] <= 0:
        return False, f"CLASH : jeu negatif ({ressort['jeu_bloc']:.1f} mm)"
    if ressort['jeu_bloc'] < 3.0:
        return True, f"Jeu faible ({ressort['jeu_bloc']:.1f} mm), prevoir butee"
    return True, f"OK (jeu = {ressort['jeu_bloc']:.1f} mm)"
def duree_vie_ressort(tau_max, tau_min=0):
    """Estimation du nombre de cycles avant rupture en fatigue."""
    tau_a = (tau_max - tau_min) / 2
    if tau_a <= SIGMA_END * 0.5:
        return float('inf')
    # loi de Basquin simplifiee
    return 1e7 * (SIGMA_END * 0.5 / tau_a) ** 8
# Ressort a raideur progressive (bi-taux avec helper spring)

def force_ressort_progressif(deplacement, k1, k2, seuil):
    """
    Force d'un systeme bi-taux.
    k1 : raideur basse (debut de course), N/mm
    k2 : raideur haute (apres le seuil), N/mm
    seuil : deplacement ou le helper spring se ferme (mm)
    """
    if deplacement <= seuil:
        return k1 * deplacement
    return k1 * seuil + k2 * (deplacement - seuil)
def courbe_progressive(k1, k2, seuil, course_max):
    """Genere la courbe force-deplacement du systeme bi-taux."""
    x = np.linspace(0, course_max, 200)
    f = np.array([force_ressort_progressif(xi, k1, k2, seuil) for xi in x])
    return x, f
# Butee de compression (bump stop)

def force_butee(deplacement, debut_contact, k_butee, exposant=2.5):
    """
    Force de la butee elastomere (loi puissance).
    debut_contact : deplacement a partir duquel la butee touche (mm)
    k_butee : raideur initiale (N/mm)
    exposant : non-linearite de l'elastomere (typiquement 2 a 3)
    """
    if deplacement <= debut_contact:
        return 0.0
    penetration = deplacement - debut_contact
    return k_butee * penetration**exposant / 10**(exposant - 1)
# Amortisseur

def amortissement_critique(k, m):
    """Coefficient d'amortissement critique (N.s/m)."""
    return 2 * np.sqrt(k * 1000 * m)
def force_amortisseur(vitesse, c_comp, c_det):
    """
    Force d'amortissement asymetrique (compression / detente).
    vitesse en m/s, force en N.
    """
    if vitesse >= 0:
        return -c_comp * vitesse
    return -c_det * vitesse
def dimensionner_amortisseur(k_roue, masse, zeta_cible=0.35):
    """
    Dimensionnement de l'amortisseur.
    zeta_cible : taux d'amortissement vise (0.25 a 0.4 pour tout-terrain).
    """
    c_crit = amortissement_critique(k_roue, masse)
    c_total = zeta_cible * c_crit

    # repartition compression / detente (typiquement 30/70)
    rapport_cd = 0.30
    c_comp = c_total * rapport_cd / (rapport_cd + (1 - rapport_cd))
    c_det = c_total * (1 - rapport_cd) / (rapport_cd + (1 - rapport_cd))

    # diametre de tige (estimation pour la pression de service)
    pression_service = 25  # bar (gaz)
    force_max = c_det * 1.5  # vitesse de piston max ~ 1.5 m/s (estimation)
    d_piston = 36  # mm (standard)
    d_tige = 16    # mm (standard)

    return {
        'c_total': c_total,
        'c_comp': c_comp,
        'c_det': c_det,
        'c_crit': c_crit,
        'zeta': zeta_cible,
        'd_piston': d_piston,
        'd_tige': d_tige,
    }
# Transmissibilite et confort

def transmissibilite(f_excitation, f_propre, zeta):
    """
    Rapport d'amplitude (transmissibilite) en fonction de la frequence d'excitation.
    """
    r = f_excitation / f_propre
    num = 1 + (2 * zeta * r)**2
    den = (1 - r**2)**2 + (2 * zeta * r)**2
    return np.sqrt(num / den)
def psd_route(f, v_kmh=40, classe='C'):
    """
    Densite spectrale de puissance (PSD) du profil de route (ISO 8608).
    f : frequence spatiale -> on convertit en temporelle via la vitesse.
    v_kmh : vitesse du vehicule (km/h)
    classe : 'A' (autoroute) a 'E' (piste)
    """
    # valeurs de reference a 0.1 cycle/m (Gd(n0))
    gd_ref = {'A': 2, 'B': 8, 'C': 32, 'D': 128, 'E': 512}  # mm^2 / (cycle/m) * 1e-6
    gd0 = gd_ref.get(classe, 32) * 1e-6  # m^2 / (cycle/m)
    n0 = 0.1  # cycle/m (frequence de reference)

    v_ms = v_kmh / 3.6
    # frequence spatiale correspondant a la frequence temporelle f
    n = f / v_ms  # cycle/m
    # PSD du deplacement de la route
    psd = gd0 * (n / n0)**(-2) * v_ms  # m^2/Hz (PSD temporelle)
    return psd
def acceleration_rms(k_roue, c, masse, v_kmh=40, classe='C', f_max=30):
    """
    Acceleration RMS de la masse suspendue (m/s^2).
    Critere de confort ISO 2631 : < 0.315 m/s^2 pour "pas inconfortable".
    """
    f_propre = frequence_propre(k_roue, masse)
    zeta = c / amortissement_critique(k_roue, masse)

    freqs = np.linspace(0.5, f_max, 500)
    df = freqs[1] - freqs[0]

    acc_psd = np.zeros_like(freqs)
    for i, f in enumerate(freqs):
        T = transmissibilite(f, f_propre, zeta)
        S_route = psd_route(f, v_kmh, classe)
        # PSD de l'acceleration = (2*pi*f)^2 * T^2 * PSD_deplacement * (2*pi*f)^2
        # on simplifie : acc = -omega^2 * x
        omega = 2 * np.pi * f
        acc_psd[i] = omega**4 * T**2 * S_route

    rms = np.sqrt(np.sum(acc_psd) * df)
    return rms, freqs, acc_psd
# Barre anti-roulis

def raideur_barre_anti_roulis(d_barre, longueur_bras, longueur_torsion, G=78500):
    """
    Raideur de la barre anti-roulis (N/mm au niveau de la roue).
    d_barre : diametre de la barre pleine (mm)
    longueur_bras : longueur du levier (mm)
    longueur_torsion : longueur de la partie torsion (mm)
    """
    J = np.pi * d_barre**4 / 32
    k_torsion = G * J / longueur_torsion  # N.mm/rad
    k_roue = k_torsion / longueur_bras**2  # N/mm
    return k_roue
def dimensionner_bar(objectif_k_roll, longueur_bras, longueur_torsion):
    """Trouve le diametre de barre pour la raideur de roulis souhaitee."""
    # k_roue = G * pi * d^4 / (32 * L_torsion * L_bras^2)
    # d^4 = k_roue * 32 * L_torsion * L_bras^2 / (G * pi)
    d4 = objectif_k_roll * 32 * longueur_torsion * longueur_bras**2 / (78500 * np.pi)
    d = d4**0.25
    return d
# Comparaison de types de ressorts

def raideur_lame(E, b, h, L, n_lames):
    """
    Raideur d'un ressort a lames (faisceau parabolique simplifie).
    b : largeur (mm), h : epaisseur d'une lame (mm), L : demi-longueur (mm)
    """
    I_lame = b * h**3 / 12
    k = 3 * E * I_lame * n_lames / (2 * L**3)
    return k
def raideur_pneumatique(pression_bar, aire_piston_cm2, volume_cm3, gamma=1.4):
    """
    Raideur d'un ressort pneumatique (processus adiabatique).
    pression en bar, aire en cm^2, volume en cm^3.
    """
    p = pression_bar * 1e5  # Pa
    A = aire_piston_cm2 * 1e-4  # m^2
    V = volume_cm3 * 1e-6  # m^3
    k = gamma * p * A**2 / V  # N/m
    return k / 1000  # N/mm
# Simulation temporelle (Euler ameliore)

def reponse_echelon(k, c, m, x0=30.0, t_max=2.0, dt=0.0005):
    """
    Reponse a un echelon (obstacle de x0 mm).
    Integration par la methode de Heun (Euler ameliore).
    """
    omega_n = np.sqrt(k * 1000 / m)
    zeta = c / (2 * m * omega_n)

    n_pts = int(t_max / dt)
    t = np.linspace(0, t_max, n_pts)
    x = np.zeros(n_pts)
    v = np.zeros(n_pts)
    x[0] = x0

    for i in range(n_pts - 1):
        a1 = (-k * 1000 * x[i] - c * v[i]) / m  # m/s^2, mais x en mm -> corriger
        # on travaille en mm et en s
        a1 = (-k * x[i] - c * v[i] / 1000) / (m / 1000)  # mm/s^2
        v_pred = v[i] + a1 * dt
        x_pred = x[i] + v[i] * dt

        a2 = (-k * x_pred - c * v_pred / 1000) / (m / 1000)
        v[i+1] = v[i] + 0.5 * (a1 + a2) * dt
        x[i+1] = x[i] + 0.5 * (v[i] + v_pred) * dt

    return t, x, omega_n, zeta
# Programme principal

if __name__ == '__main__':
    # Masses par essieu
    m_av = MASSE_SUSPENDUE * REPARTITION_AVANT
    m_ar = MASSE_SUSPENDUE * (1 - REPARTITION_AVANT)
    m_roue_av = m_av / 2
    m_roue_ar = m_ar / 2

    charge_stat_av = m_roue_av * 9.81  # N
    charge_stat_ar = m_roue_ar * 9.81

    # raideur au niveau de la roue
    k_roue_av = raideur_roue(m_roue_av, F_CIBLE_AV)
    k_roue_ar = raideur_roue(m_roue_ar, F_CIBLE_AR)

    # raideur du ressort physique
    k_ressort_av = raideur_ressort(k_roue_av, RATIO_LEVIER)
    k_ressort_ar = raideur_ressort(k_roue_ar, RATIO_LEVIER)

    print("=" * 60)
    print("DIMENSIONNEMENT SUSPENSION -- VEHICULE TOUT-TERRAIN SEGULA")
    print("=" * 60)
    print()
    print(f"{'Parametre':<30} {'Avant':>12} {'Arriere':>12}")
    print("-" * 56)
    print(f"{'Masse par roue (kg)':<30} {m_roue_av:12.1f} {m_roue_ar:12.1f}")
    print(f"{'Charge statique (N)':<30} {charge_stat_av:12.1f} {charge_stat_ar:12.1f}")
    print(f"{'Frequence cible (Hz)':<30} {F_CIBLE_AV:12.2f} {F_CIBLE_AR:12.2f}")
    print(f"{'K roue (N/mm)':<30} {k_roue_av:12.2f} {k_roue_ar:12.2f}")
    print(f"{'K ressort (N/mm)':<30} {k_ressort_av:12.2f} {k_ressort_ar:12.2f}")
    print()

    # dimensionnement des ressorts
    ressort_av = dimensionner_ressort(k_ressort_av, charge_stat_av / RATIO_LEVIER)
    ressort_ar = dimensionner_ressort(k_ressort_ar, charge_stat_ar / RATIO_LEVIER)

    for nom, res in [('AVANT', ressort_av), ('ARRIERE', ressort_ar)]:
        ok, msg = verifier_clash(res)
        n_vie = duree_vie_ressort(res['tau_max'])
        print(f"--- Ressort {nom} ---")
        print(f"  Fil            : {res['d_fil']:.1f} mm")
        print(f"  D moyen        : {res['D_moyen']:.0f} mm")
        print(f"  D exterieur    : {res['D_ext']:.0f} mm")
        print(f"  Spires actives : {res['n_actives']:.1f}")
        print(f"  Longueur libre : {res['L0']:.0f} mm")
        print(f"  Hauteur bloc   : {res['hauteur_bloc']:.0f} mm")
        print(f"  Tau max        : {res['tau_max']:.0f} MPa (adm: {SIGMA_ADM} MPa)")
        print(f"  Masse          : {res['masse']:.2f} kg")
        print(f"  Clash          : {msg}")
        if n_vie == float('inf'):
            print(f"  Duree de vie   : illimitee")
        else:
            print(f"  Duree de vie   : {n_vie:.0f} cycles")
        print()

    # frequences de tangage et roulis
    freq_modes = frequences_tangage_roulis(k_roue_av, k_roue_ar, m_av, m_ar)
    print(f"Frequence rebond avant  : {freq_modes['f_rebond_av']:.2f} Hz")
    print(f"Frequence rebond arriere: {freq_modes['f_rebond_ar']:.2f} Hz")
    print(f"Frequence tangage       : {freq_modes['f_tangage']:.2f} Hz")
    print()

    # amortisseur
    amort_av = dimensionner_amortisseur(k_roue_av, m_roue_av, zeta_cible=0.35)
    amort_ar = dimensionner_amortisseur(k_roue_ar, m_roue_ar, zeta_cible=0.30)
    print(f"--- Amortisseur avant (zeta={amort_av['zeta']}) ---")
    print(f"  C critique : {amort_av['c_crit']:.1f} N.s/m")
    print(f"  C total    : {amort_av['c_total']:.1f} N.s/m")
    print(f"  C comp     : {amort_av['c_comp']:.1f} N.s/m")
    print(f"  C detente  : {amort_av['c_det']:.1f} N.s/m")
    print()

    # barre anti-roulis
    k_bar = 8.0  # N/mm objectif au niveau de la roue
    d_bar = dimensionner_bar(k_bar, 200, 800)
    k_bar_verif = raideur_barre_anti_roulis(d_bar, 200, 800)
    print(f"--- Barre anti-roulis ---")
    print(f"  Objectif     : {k_bar:.1f} N/mm")
    print(f"  Diametre     : {d_bar:.1f} mm")
    print(f"  K verifie    : {k_bar_verif:.2f} N/mm")
    print()

    # comparaison de types de ressorts
    k_lame = raideur_lame(210e3, 60, 6, 500, 4)
    k_pneum = raideur_pneumatique(8, 30, 500)
    print(f"--- Comparaison types de ressorts ---")
    print(f"  Helicoldal    : {k_ressort_av:.2f} N/mm")
    print(f"  Lames (4x60x6, L=500) : {k_lame:.2f} N/mm")
    print(f"  Pneumatique (8bar, 30cm2, 500cm3) : {k_pneum:.2f} N/mm")
    print()

    # acceleration RMS (confort)
    rms_val, _, _ = acceleration_rms(k_roue_av, amort_av['c_total'], m_roue_av,
                                     v_kmh=40, classe='C')
    print(f"Acceleration RMS (40 km/h, piste C) : {rms_val:.3f} m/s^2")
    if rms_val < 0.315:
        print(f"  -> Confort correct (ISO 2631 : < 0.315)")
    elif rms_val < 0.63:
        print(f"  -> Legerement inconfortable (ISO 2631)")
    else:
        print(f"  -> Inconfortable (ISO 2631)")
    print()

    # Figures
    fig = plt.figure(figsize=(16, 14))

    # 1 -- Reponse a un obstacle pour differents amortissements
    ax1 = fig.add_subplot(3, 3, 1)
    for z in [0.15, 0.25, 0.35, 0.50, 0.70]:
        c = z * amort_av['c_crit']
        t, x, wn, zeta_calc = reponse_echelon(k_roue_av, c, m_roue_av)
        ax1.plot(t, x, label=f'zeta = {z:.2f}')
    ax1.set_xlabel('Temps (s)')
    ax1.set_ylabel('Deplacement (mm)')
    ax1.set_title('Reponse obstacle (30 mm)')
    ax1.legend(fontsize=7)
    ax1.grid(True)
    ax1.axhline(0, color='k', linewidth=0.5)

    # 2 -- Raideur vs frequence
    ax2 = fig.add_subplot(3, 3, 2)
    freqs_plot = np.linspace(0.5, 3.0, 100)
    k_plot_av = [raideur_roue(m_roue_av, f) for f in freqs_plot]
    k_plot_ar = [raideur_roue(m_roue_ar, f) for f in freqs_plot]
    ax2.plot(freqs_plot, k_plot_av, 'b-', linewidth=2, label='Avant')
    ax2.plot(freqs_plot, k_plot_ar, 'r--', linewidth=2, label='Arriere')
    ax2.axvspan(1.0, 2.0, alpha=0.1, color='green', label='Plage cible')
    ax2.plot(F_CIBLE_AV, k_roue_av, 'bo', markersize=8)
    ax2.plot(F_CIBLE_AR, k_roue_ar, 'rs', markersize=8)
    ax2.set_xlabel('Frequence propre (Hz)')
    ax2.set_ylabel('Raideur roue (N/mm)')
    ax2.set_title('Raideur necessaire')
    ax2.legend(fontsize=8)
    ax2.grid(True)

    # 3 -- Force-deplacement (lineaire + progressif + butee)
    ax3 = fig.add_subplot(3, 3, 3)
    dep = np.linspace(0, DEBATTEMENT_COMP + 20, 200)
    f_lin = k_roue_av * dep
    # systeme progressif : k1 = 0.7*k, k2 = 1.8*k, seuil a 50 mm
    k1_prog = 0.7 * k_roue_av
    k2_prog = 1.8 * k_roue_av
    f_prog = np.array([force_ressort_progressif(d, k1_prog, k2_prog, 50) for d in dep])
    f_but = np.array([force_butee(d, 70, 5.0) for d in dep])
    ax3.plot(dep, f_lin, 'b-', linewidth=2, label='Lineaire')
    ax3.plot(dep, f_prog, 'g--', linewidth=2, label='Progressif')
    ax3.plot(dep, f_but, 'r:', linewidth=2, label='Butee')
    ax3.plot(dep, f_prog + f_but, 'k-', linewidth=1.5, alpha=0.6, label='Prog. + butee')
    ax3.set_xlabel('Deplacement (mm)')
    ax3.set_ylabel('Force (N)')
    ax3.set_title('Caracteristiques force-deplacement')
    ax3.legend(fontsize=7)
    ax3.grid(True)

    # 4 -- Transmissibilite
    ax4 = fig.add_subplot(3, 3, 4)
    f_exc = np.linspace(0.1, 20, 500)
    f_n = frequence_propre(k_roue_av, m_roue_av)
    for z in [0.15, 0.25, 0.35, 0.50]:
        T = [transmissibilite(fe, f_n, z) for fe in f_exc]
        ax4.semilogy(f_exc, T, label=f'zeta={z:.2f}')
    ax4.axhline(1, color='k', linestyle=':', alpha=0.3)
    ax4.axvline(f_n, color='gray', linestyle='--', alpha=0.3, label=f'fn={f_n:.1f} Hz')
    ax4.set_xlabel('Frequence excitation (Hz)')
    ax4.set_ylabel('Transmissibilite')
    ax4.set_title('Transmissibilite')
    ax4.legend(fontsize=7)
    ax4.grid(True, which='both')

    # 5 -- PSD route et acceleration
    ax5 = fig.add_subplot(3, 3, 5)
    rms_40, freqs_psd, acc_psd_40 = acceleration_rms(k_roue_av, amort_av['c_total'],
                                                      m_roue_av, 40, 'C')
    rms_80, freqs_psd, acc_psd_80 = acceleration_rms(k_roue_av, amort_av['c_total'],
                                                      m_roue_av, 80, 'C')
    ax5.semilogy(freqs_psd, acc_psd_40, 'b-', label=f'40 km/h (RMS={rms_40:.2f})')
    ax5.semilogy(freqs_psd, acc_psd_80, 'r--', label=f'80 km/h (RMS={rms_80:.2f})')
    ax5.set_xlabel('Frequence (Hz)')
    ax5.set_ylabel('PSD acceleration (m/s^2)^2/Hz')
    ax5.set_title('Densite spectrale acceleration')
    ax5.legend(fontsize=7)
    ax5.grid(True, which='both')

    # 6 -- Force amortisseur (courbe vitesse)
    ax6 = fig.add_subplot(3, 3, 6)
    vitesses = np.linspace(-2, 2, 200)
    forces_amort = [force_amortisseur(v, amort_av['c_comp'], amort_av['c_det']) for v in vitesses]
    ax6.plot(vitesses, forces_amort, 'b-', linewidth=2)
    ax6.set_xlabel('Vitesse piston (m/s)')
    ax6.set_ylabel('Force (N)')
    ax6.set_title('Caracteristique amortisseur')
    ax6.axhline(0, color='k', linewidth=0.5)
    ax6.axvline(0, color='k', linewidth=0.5)
    ax6.annotate(f'Comp: {amort_av["c_comp"]:.0f} N.s/m', xy=(-0.5, 200), fontsize=8)
    ax6.annotate(f'Det: {amort_av["c_det"]:.0f} N.s/m', xy=(0.5, -400), fontsize=8)
    ax6.grid(True)

    # 7 -- Comparaison de types de ressorts
    ax7 = fig.add_subplot(3, 3, 7)
    types_noms = ['Helicoldal', 'Lames', 'Pneumatique']
    types_k = [k_ressort_av, k_lame, k_pneum]
    types_couleurs = ['#1976D2', '#388E3C', '#F57C00']
    ax7.bar(types_noms, types_k, color=types_couleurs)
    ax7.set_ylabel('Raideur (N/mm)')
    ax7.set_title('Comparaison types de ressorts')
    ax7.grid(True, axis='y')

    # 8 -- Raideur barre anti-roulis en fonction du diametre
    ax8 = fig.add_subplot(3, 3, 8)
    d_range = np.linspace(15, 30, 50)
    k_bar_range = [raideur_barre_anti_roulis(d, 200, 800) for d in d_range]
    ax8.plot(d_range, k_bar_range, 'b-', linewidth=2)
    ax8.plot(d_bar, k_bar_verif, 'ro', markersize=8, label=f'd={d_bar:.1f} mm')
    ax8.set_xlabel('Diametre barre (mm)')
    ax8.set_ylabel('K roulis (N/mm)')
    ax8.set_title('Barre anti-roulis')
    ax8.legend(fontsize=8)
    ax8.grid(True)

    # 9 -- Tableau recapitulatif (sous forme de texte dans un subplot)
    ax9 = fig.add_subplot(3, 3, 9)
    ax9.axis('off')
    tableau = [
        ['Parametre', 'Avant', 'Arriere'],
        ['Freq. propre (Hz)', f'{F_CIBLE_AV}', f'{F_CIBLE_AR}'],
        ['K roue (N/mm)', f'{k_roue_av:.1f}', f'{k_roue_ar:.1f}'],
        ['K ressort (N/mm)', f'{k_ressort_av:.1f}', f'{k_ressort_ar:.1f}'],
        ['Fil (mm)', f'{ressort_av["d_fil"]:.1f}', f'{ressort_ar["d_fil"]:.1f}'],
        ['D moyen (mm)', f'{ressort_av["D_moyen"]:.0f}', f'{ressort_ar["D_moyen"]:.0f}'],
        ['L0 (mm)', f'{ressort_av["L0"]:.0f}', f'{ressort_ar["L0"]:.0f}'],
        ['Zeta amort.', f'{amort_av["zeta"]}', f'{amort_ar["zeta"]}'],
        ['Masse ress. (kg)', f'{ressort_av["masse"]:.2f}', f'{ressort_ar["masse"]:.2f}'],
    ]
    table = ax9.table(cellText=tableau[1:], colLabels=tableau[0],
                      cellLoc='center', loc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)
    ax9.set_title('Recapitulatif', pad=20)

    fig.suptitle('Dimensionnement suspension -- SEGULA Technologies', fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('spring_sizing.png', dpi=150)
    print("Figure sauvegardee : spring_sizing.png")
