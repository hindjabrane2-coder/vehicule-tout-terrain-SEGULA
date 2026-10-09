# Vehicule tout-terrain -- SEGULA Technologies

Stage ingenieur chez SEGULA Technologies (avril - aout 2024). Conception mecanique d'un vehicule tout-terrain avec focus sur le systeme de suspension et l'analyse structurelle du chassis.

## Contenu

- `suspension_kinematics.py` : modelisation cinematique de la suspension a double triangulation, calcul du debattement, du carrossage et de la voie en fonction de la course verticale
- `stress_analysis.py` : analyse des contraintes sur le chassis tubulaire, verification au critere de Von Mises sous charges statiques et dynamiques (choc, virage)
- `spring_sizing.py` : dimensionnement du ressort helicoldal et de l'amortisseur, calcul de la frequence propre et du taux d'amortissement

## Contexte

Le projet portait sur la conception d'un vehicule tout-terrain leger pour usage recreatif. Mon travail a couvert la modelisation cinematique du train avant (suspension a double triangulation), le dimensionnement des elements de suspension, et la verification structurelle du chassis par elements finis sous ANSYS. Les modeles CAO ont ete realises sous SolidWorks.

## Dependances

- Python 3.8+
- numpy, matplotlib, scipy
