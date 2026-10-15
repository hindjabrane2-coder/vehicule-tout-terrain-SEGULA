# Vehicule Tout-Terrain SEGULA

Outils de calcul pour la conception d'un vehicule tout-terrain dans le cadre d'un projet chez SEGULA Technologies. Le code couvre trois domaines : cinematique de suspension, analyse structurelle du chassis et dimensionnement des organes de suspension.

Projet ecrit en C# pur, sans dependances externes. Les solveurs numeriques (bisection, elimination de Gauss, integration de Heun) sont implementes directement.


## Structure

```
src/
    CinematiqueSuspension.cs    Cinematique double triangulation, centre de roulis, bump steer
    AnalyseContraintes.cs       Contraintes chassis tubulaire, fatigue, flambement, FEA poutre
    DimensionnementRessort.cs   Ressort, amortisseur, barre anti-roulis, confort ISO 2631
```


## Cinematique de suspension

Modelisation d'une suspension a double triangulation (double wishbone) dans le plan transversal. La resolution utilise une methode de bisection imbriquee pour trouver la position d'equilibre de la fusee a chaque hauteur de debattement.

Grandeurs calculees sur la plage de debattement :
- Carrossage (angle de la fusee par rapport a la verticale)
- Variation de voie laterale
- Centre instantane de rotation (intersection des prolongements des bras)
- Hauteur du centre de roulis
- Rapport de mouvement (course ressort / course roue)
- Bump steer (variation de pincement due au debattement)
- Anti-plongee et anti-cabrage (geometrie vue de cote)
- Analyse de sensibilite des pivots

Geometrie nominale :
- Pivot superieur interieur : (150, 350) mm
- Pivot inferieur interieur : (120, 150) mm
- Demi-voie : 700 mm
- Empattement : 2200 mm


## Analyse structurelle

Verification d'un chassis tubulaire sous plusieurs cas de charge : statique, virage (1.5g), freinage (1.0g), torsion diagonale, choc vertical (x3) et combinaison freinage + virage.

Fonctionnalites :
- Proprietes de section (tube circulaire et plaque rectangulaire)
- Contrainte de Von Mises (etat plan de contrainte)
- Charge critique de flambement (Euler)
- Fatigue par diagramme de Goodman et courbe S-N (Basquin)
- Contrainte thermique
- Verification de goussets et cordons de soudure
- Modele elements finis simplifie a 5 noeuds (poutre Euler-Bernoulli 2D, 3 DDL par noeud)
- Resolution par elimination de Gauss avec pivot partiel
- Optimisation de masse (tube le plus leger satisfaisant le coefficient de securite)

Materiaux disponibles : S355, E355, 4130 chromoly. Catalogue de tubes de 42.4x2.0 a 60.3x3.2 mm.


## Dimensionnement des organes de suspension

Calcul complet des ressorts, amortisseurs et barres anti-roulis pour un vehicule de 450 kg.

Ressort helicoldal :
- Raideur au niveau de la roue a partir de la frequence propre cible
- Conversion raideur roue / raideur ressort par le rapport de levier
- Dimensionnement geometrique (diametre de fil, D moyen, nombre de spires, longueur libre)
- Facteur de correction de Wahl
- Verification du clash (spires jointives)
- Duree de vie en fatigue

Amortisseur :
- Amortissement critique et taux d'amortissement
- Repartition compression / detente (30/70)
- Caracteristique force-vitesse asymetrique

Confort et dynamique :
- Transmissibilite en fonction de la frequence d'excitation
- Densite spectrale de puissance du profil de route (ISO 8608, classes A a E)
- Acceleration RMS de la masse suspendue (critere ISO 2631)
- Simulation temporelle par methode de Heun (reponse a un obstacle)

Autres types de ressorts :
- Ressort a raideur progressive (bi-taux avec helper spring)
- Butee de compression (loi puissance)
- Ressort a lames (faisceau parabolique)
- Ressort pneumatique (processus adiabatique)
- Barre anti-roulis (barre de torsion pleine)


## Parametres du vehicule

| Parametre | Valeur |
|-----------|--------|
| Masse totale | 450 kg |
| Masse suspendue | 350 kg |
| Empattement | 2200 mm |
| Voie avant | 1400 mm |
| Hauteur CdG | 500 mm |
| Frequence cible avant | 1.3 Hz |
| Frequence cible arriere | 1.5 Hz |
| Debattement | +/- 80 mm |
| Motion ratio | 0.70 |


## Compilation

```bash
csc -out:VehiculeToutTerrain.exe src/*.cs
```

Necessite le SDK .NET ou Mono. Aucune dependance NuGet.
