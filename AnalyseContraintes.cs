using System;
using System.Collections.Generic;

namespace VehiculeToutTerrain
{
    // Analyse des contraintes sur chassis tubulaire
    // Von Mises, fatigue (Goodman), flambement (Euler), concentration de contraintes,
    // assemblage par elements finis simplifie (poutre Euler-Bernoulli)

    public static class AnalyseContraintes
    {
        // Materiaux

        public class Materiau
        {
            public double E;           // module d'Young (MPa)
            public double Nu;
            public double SigmaY;      // limite elastique (MPa)
            public double SigmaU;      // resistance a la traction (MPa)
            public double SigmaF;      // limite d'endurance flexion alternee (MPa)
            public double Rho;         // masse volumique (kg/m^3)
            public double Alpha;       // coefficient de dilatation (1/K)
        }

        public static readonly Dictionary<string, Materiau> Materiaux = new Dictionary<string, Materiau>
        {
            ["S355"] = new Materiau
            {
                E = 210e3, Nu = 0.3, SigmaY = 355, SigmaU = 510,
                SigmaF = 240, Rho = 7850, Alpha = 12e-6
            },
            ["E355"] = new Materiau
            {
                E = 210e3, Nu = 0.3, SigmaY = 355, SigmaU = 580,
                SigmaF = 260, Rho = 7850, Alpha = 12e-6
            },
            ["4130_chromoly"] = new Materiau
            {
                E = 205e3, Nu = 0.29, SigmaY = 460, SigmaU = 670,
                SigmaF = 310, Rho = 7850, Alpha = 11.2e-6
            },
        };

        public static Materiau Mat = Materiaux["S355"];

        // Sections de tubes

        public class TubeDefinition
        {
            public string Label;
            public double DExt;
            public double Ep;
        }

        public static readonly TubeDefinition[] Tubes =
        {
            new TubeDefinition { Label = "42.4x2.0", DExt = 42.4, Ep = 2.0 },
            new TubeDefinition { Label = "42.4x2.6", DExt = 42.4, Ep = 2.6 },
            new TubeDefinition { Label = "48.3x2.6", DExt = 48.3, Ep = 2.6 },
            new TubeDefinition { Label = "48.3x3.2", DExt = 48.3, Ep = 3.2 },
            new TubeDefinition { Label = "50.8x2.6", DExt = 50.8, Ep = 2.6 },
            new TubeDefinition { Label = "50.8x3.2", DExt = 50.8, Ep = 3.2 },
            new TubeDefinition { Label = "60.3x3.2", DExt = 60.3, Ep = 3.2 },
        };

        // Parametres du vehicule
        public const double Masse = 450;
        public const double EmpattementChassis = 2200;
        public const double Voie = 1400;
        public const double HCdg = 500;

        // Facteurs de concentration de contraintes
        public const double KtSoudureBout = 1.5;
        public const double KtSoudureAngle = 2.0;
        public const double KtSoudureMeule = 1.3;

        // Proprietes geometriques d'un tube

        public class SectionTube
        {
            public double Aire;
            public double Inertie;
            public double Module;       // module de flexion
            public double JPolaire;
            public double RGiration;
            public double DExt;
            public double DInt;
            public double Ep;
            public double MasseLin;     // kg/mm
        }

        public static SectionTube CalculerSectionTube(double dExt, double ep)
        {
            double dInt = dExt - 2 * ep;
            double aire = Math.PI / 4 * (dExt * dExt - dInt * dInt);
            double inertie = Math.PI / 64 * (Math.Pow(dExt, 4) - Math.Pow(dInt, 4));
            double module = inertie / (dExt / 2);
            double jPolaire = 2 * inertie;
            double rGiration = Math.Sqrt(inertie / aire);
            double masseLin = aire * 7850 / 1e6;

            return new SectionTube
            {
                Aire = aire, Inertie = inertie, Module = module,
                JPolaire = jPolaire, RGiration = rGiration,
                DExt = dExt, DInt = dInt, Ep = ep, MasseLin = masseLin,
            };
        }

        // Section plaque (gousset, renfort)

        public class SectionPlaque
        {
            public double Aire;
            public double Inertie;
            public double Module;
            public double Largeur;
            public double Epaisseur;
        }

        public static SectionPlaque CalculerSectionPlaque(double largeur, double epaisseur)
        {
            double aire = largeur * epaisseur;
            double inertie = largeur * Math.Pow(epaisseur, 3) / 12;
            return new SectionPlaque
            {
                Aire = aire, Inertie = inertie, Module = inertie / (epaisseur / 2),
                Largeur = largeur, Epaisseur = epaisseur,
            };
        }

        // Contraintes elementaires

        public static double ContrainteTraction(double force, SectionTube sec)
        {
            return force / sec.Aire;
        }

        public static double ContrainteFlexion(double moment, SectionTube sec)
        {
            return Math.Abs(moment) / sec.Module;
        }

        public static double ContrainteTorsion(double couple, SectionTube sec)
        {
            return Math.Abs(couple) * (sec.DExt / 2) / sec.JPolaire;
        }

        public static double VonMises(double sigmaX, double sigmaY, double tauXy)
        {
            return Math.Sqrt(sigmaX * sigmaX - sigmaX * sigmaY + sigmaY * sigmaY
                             + 3 * tauXy * tauXy);
        }

        public static double VonMisesSimple(double sigma, double tau)
        {
            return Math.Sqrt(sigma * sigma + 3 * tau * tau);
        }

        public static double CoefficientSecurite(double sigmaEq, double sigmaLim)
        {
            if (sigmaEq < 1e-10) return double.PositiveInfinity;
            return sigmaLim / sigmaEq;
        }

        // Cas de charge

        public class CasCharge
        {
            public string Nom;
            public double ForceAxiale;
            public double Moment;
            public double Couple;
        }

        public static CasCharge CasStatique()
        {
            double poids = Masse * 9.81;
            double q = poids / EmpattementChassis;
            double mFlex = q * EmpattementChassis * EmpattementChassis / 8;
            return new CasCharge { Nom = "Statique", ForceAxiale = 0, Moment = mFlex, Couple = 0 };
        }

        public static CasCharge CasVirage(double gLateral = 1.5)
        {
            double poids = Masse * 9.81;
            double fLat = Masse * gLateral * 9.81;
            double mRoulis = fLat * HCdg;
            double mFlex = poids * EmpattementChassis / 8;
            return new CasCharge
            {
                Nom = $"Virage {gLateral}g",
                ForceAxiale = 0, Moment = mFlex, Couple = mRoulis
            };
        }

        public static CasCharge CasFreinage(double gLong = 1.0)
        {
            double poids = Masse * 9.81;
            double fFrein = Masse * gLong * 9.81;
            double mFlex = poids * EmpattementChassis / 8 + fFrein * HCdg / 2;
            return new CasCharge
            {
                Nom = $"Freinage {gLong}g",
                ForceAxiale = fFrein / 2, Moment = mFlex, Couple = 0
            };
        }

        public static CasCharge CasTorsion()
        {
            double poids = Masse * 9.81;
            double coupleTorsion = poids / 4 * Voie / 2;
            return new CasCharge
            {
                Nom = "Torsion diag.",
                ForceAxiale = 0, Moment = 0, Couple = coupleTorsion
            };
        }

        public static CasCharge CasChoc(double facteur = 3.0)
        {
            var stat = CasStatique();
            return new CasCharge
            {
                Nom = $"Choc x{facteur:F0}",
                ForceAxiale = 0, Moment = stat.Moment * facteur, Couple = 0
            };
        }

        public static CasCharge CasCombine()
        {
            var frein = CasFreinage(0.8);
            var virage = CasVirage(1.0);
            double mComb = Math.Sqrt(frein.Moment * frein.Moment + virage.Moment * virage.Moment);
            return new CasCharge
            {
                Nom = "Combine fr+vir",
                ForceAxiale = frein.ForceAxiale, Moment = mComb, Couple = virage.Couple
            };
        }

        public static CasCharge[] TousLesCas()
        {
            return new[] { CasStatique(), CasVirage(), CasFreinage(),
                           CasTorsion(), CasChoc(), CasCombine() };
        }

        // Analyse d'un tube sous un cas de charge

        public class ResultatAnalyse
        {
            public string Nom;
            public double SigmaFlex;
            public double SigmaAxe;
            public double Tau;
            public double SigmaTotal;
            public double TauTotal;
            public double VonMisesVal;
            public double Cs;
        }

        public static ResultatAnalyse AnalyserTube(SectionTube sec, CasCharge cas,
            Materiau mat, double kt = 1.0)
        {
            double sigmaFlex = ContrainteFlexion(cas.Moment, sec);
            double sigmaAxe = ContrainteTraction(cas.ForceAxiale, sec);
            double tau = ContrainteTorsion(cas.Couple, sec);

            double sigmaTotal = (sigmaFlex + Math.Abs(sigmaAxe)) * kt;
            double tauTotal = tau * kt;
            double vm = VonMisesSimple(sigmaTotal, tauTotal);
            double cs = CoefficientSecurite(vm, mat.SigmaY);

            return new ResultatAnalyse
            {
                Nom = cas.Nom,
                SigmaFlex = sigmaFlex, SigmaAxe = sigmaAxe, Tau = tau,
                SigmaTotal = sigmaTotal, TauTotal = tauTotal,
                VonMisesVal = vm, Cs = cs,
            };
        }

        // Flambement (Euler)

        public class ResultatFlambement
        {
            public double PCrit;
            public double CsFlamb;
            public double Longueur;
        }

        public static double ChargeCritiqueEuler(SectionTube sec, double longueur,
            double kFlamb = 1.0, double E = 210e3)
        {
            double lEff = kFlamb * longueur;
            return Math.PI * Math.PI * E * sec.Inertie / (lEff * lEff);
        }

        public static ResultatFlambement VerifierFlambement(SectionTube sec,
            double forceCompression, double longueur, double kFlamb = 1.0, double E = 210e3)
        {
            double pCrit = ChargeCritiqueEuler(sec, longueur, kFlamb, E);
            double cs = Math.Abs(forceCompression) > 0 ?
                pCrit / Math.Abs(forceCompression) : double.PositiveInfinity;
            return new ResultatFlambement { PCrit = pCrit, CsFlamb = cs, Longueur = longueur };
        }

        // Fatigue (Goodman)

        public class ResultatFatigue
        {
            public double SigmaM;
            public double SigmaA;
            public double SigmaAdmGoodman;
            public double NCycles;
            public double CsFatigue;
        }

        public static double GoodmanSigmaAdmissible(double sigmaM, double sigmaU, double sigmaF)
        {
            if (sigmaM >= sigmaU) return 0.0;
            return sigmaF * (1 - sigmaM / sigmaU);
        }

        public static double DureeDeVieSn(double sigmaA, double sigmaF,
            double nCoude = 1e6, double exposant = 5)
        {
            if (sigmaA <= sigmaF) return double.PositiveInfinity;
            return nCoude * Math.Pow(sigmaF / sigmaA, exposant);
        }

        public static ResultatFatigue AnalyseFatigue(double sigmaMax, double sigmaMin, Materiau mat)
        {
            double sigmaM = (sigmaMax + sigmaMin) / 2;
            double sigmaA = (sigmaMax - sigmaMin) / 2;
            double sigmaAdm = GoodmanSigmaAdmissible(sigmaM, mat.SigmaU, mat.SigmaF);
            double nCycles = DureeDeVieSn(sigmaA, mat.SigmaF);
            double csFatigue = sigmaA > 0 ? sigmaAdm / sigmaA : double.PositiveInfinity;

            return new ResultatFatigue
            {
                SigmaM = sigmaM, SigmaA = sigmaA,
                SigmaAdmGoodman = sigmaAdm, NCycles = nCycles, CsFatigue = csFatigue,
            };
        }

        // Contrainte thermique

        public static double ContrainteThermique(double deltaT, Materiau mat)
        {
            return mat.E * mat.Alpha * deltaT;
        }

        // Gousset de renfort

        public class ResultatGousset
        {
            public double SigmaPlaque;
            public double TauSoudure;
            public double CsPlaque;
            public double CsSoudure;
        }

        public static ResultatGousset AnalyseGousset(double force, double largeurGousset,
            double epGousset, double longueurSoudure, Materiau mat)
        {
            var secPlaque = CalculerSectionPlaque(largeurGousset, epGousset);
            double sigmaPlaque = Math.Abs(force) / secPlaque.Aire;

            double gorge = 0.7 * epGousset;
            double aireSoudure = gorge * longueurSoudure;
            double tauSoudure = aireSoudure > 0 ?
                Math.Abs(force) / aireSoudure : double.PositiveInfinity;

            double csPlaque = sigmaPlaque > 0 ?
                mat.SigmaY / sigmaPlaque : double.PositiveInfinity;
            double csSoudure = tauSoudure > 0 && tauSoudure < double.PositiveInfinity ?
                (0.6 * mat.SigmaY) / tauSoudure : double.PositiveInfinity;

            return new ResultatGousset
            {
                SigmaPlaque = sigmaPlaque, TauSoudure = tauSoudure,
                CsPlaque = csPlaque, CsSoudure = csSoudure,
            };
        }

        // Matrice de rigidite poutre 2D (Euler-Bernoulli, 6 DDL)

        public static double[,] MatriceRigiditePoutre(double E, double I, double A, double L)
        {
            double kAxial = E * A / L;
            double kFlex = E * I / (L * L * L);

            var ke = new double[6, 6];

            // rigidite axiale
            ke[0, 0] = kAxial;
            ke[0, 3] = -kAxial;
            ke[3, 0] = -kAxial;
            ke[3, 3] = kAxial;

            // rigidite en flexion
            ke[1, 1] = 12 * kFlex;
            ke[1, 2] = 6 * kFlex * L;
            ke[1, 4] = -12 * kFlex;
            ke[1, 5] = 6 * kFlex * L;

            ke[2, 1] = 6 * kFlex * L;
            ke[2, 2] = 4 * kFlex * L * L;
            ke[2, 4] = -6 * kFlex * L;
            ke[2, 5] = 2 * kFlex * L * L;

            ke[4, 1] = -12 * kFlex;
            ke[4, 2] = -6 * kFlex * L;
            ke[4, 4] = 12 * kFlex;
            ke[4, 5] = -6 * kFlex * L;

            ke[5, 1] = 6 * kFlex * L;
            ke[5, 2] = 2 * kFlex * L * L;
            ke[5, 4] = -6 * kFlex * L;
            ke[5, 5] = 4 * kFlex * L * L;

            return ke;
        }

        // Assemblage chassis simplifie (5 noeuds, 4 elements poutre)

        public static double[,] AssemblerChassisSimplifie(SectionTube sec, Materiau mat)
        {
            double E = mat.E;
            double I = sec.Inertie;
            double A = sec.Aire;

            double lTraverse = Voie / 2;
            double lLongeron = EmpattementChassis / 2;

            int[][] elements = {
                new[] { 0, 1 }, new[] { 1, 2 }, new[] { 2, 3 }, new[] { 3, 4 }
            };
            double[] longueurs = { lTraverse, lLongeron, lLongeron, lTraverse };

            int nNoeuds = 5;
            int nDdl = 3 * nNoeuds;
            var kGlob = new double[nDdl, nDdl];

            for (int e = 0; e < 4; e++)
            {
                int n1 = elements[e][0], n2 = elements[e][1];
                double L = longueurs[e];
                var ke = MatriceRigiditePoutre(E, I, A, L);

                int[] ddl = { 3 * n1, 3 * n1 + 1, 3 * n1 + 2,
                              3 * n2, 3 * n2 + 1, 3 * n2 + 2 };
                for (int i = 0; i < 6; i++)
                    for (int j = 0; j < 6; j++)
                        kGlob[ddl[i], ddl[j]] += ke[i, j];
            }

            return kGlob;
        }

        // Resolution du systeme (elimination de Gauss avec pivot partiel)

        static double[] ResoudreSystemeLineaire(double[,] A, double[] b)
        {
            int n = b.Length;
            var aug = new double[n, n + 1];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++) aug[i, j] = A[i, j];
                aug[i, n] = b[i];
            }

            // elimination avec pivot partiel
            for (int k = 0; k < n; k++)
            {
                int maxRow = k;
                double maxVal = Math.Abs(aug[k, k]);
                for (int i = k + 1; i < n; i++)
                {
                    if (Math.Abs(aug[i, k]) > maxVal)
                    {
                        maxVal = Math.Abs(aug[i, k]);
                        maxRow = i;
                    }
                }

                if (maxRow != k)
                {
                    for (int j = k; j <= n; j++)
                    {
                        double tmp = aug[k, j];
                        aug[k, j] = aug[maxRow, j];
                        aug[maxRow, j] = tmp;
                    }
                }

                for (int i = k + 1; i < n; i++)
                {
                    double factor = aug[i, k] / aug[k, k];
                    for (int j = k; j <= n; j++)
                        aug[i, j] -= factor * aug[k, j];
                }
            }

            // substitution arriere
            var x = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                x[i] = aug[i, n];
                for (int j = i + 1; j < n; j++)
                    x[i] -= aug[i, j] * x[j];
                x[i] /= aug[i, i];
            }

            return x;
        }

        // Resolution chassis

        public class ResultatChassis
        {
            public double[] Deplacements;
            public double[] Reactions;
        }

        public static ResultatChassis ResoudreChassis(SectionTube sec, Materiau mat,
            Dictionary<int, double[]> forcesNoeuds)
        {
            var K = AssemblerChassisSimplifie(sec, mat);
            int nDdl = 15; // 5 noeuds * 3 DDL

            var F = new double[nDdl];
            foreach (var kv in forcesNoeuds)
            {
                int noeud = kv.Key;
                F[3 * noeud] = kv.Value[0];
                F[3 * noeud + 1] = kv.Value[1];
                F[3 * noeud + 2] = kv.Value[2];
            }

            // conditions limites : noeuds 0 et 4 bloques (DDL 0,1,2,12,13,14)
            int[] ddlBloques = { 0, 1, 2, 12, 13, 14 };
            var bloqueSet = new HashSet<int>(ddlBloques);
            var ddlLibres = new List<int>();
            for (int i = 0; i < nDdl; i++)
                if (!bloqueSet.Contains(i)) ddlLibres.Add(i);

            int nLibre = ddlLibres.Count;
            var kRed = new double[nLibre, nLibre];
            var fRed = new double[nLibre];

            for (int i = 0; i < nLibre; i++)
            {
                fRed[i] = F[ddlLibres[i]];
                for (int j = 0; j < nLibre; j++)
                    kRed[i, j] = K[ddlLibres[i], ddlLibres[j]];
            }

            double[] uRed = ResoudreSystemeLineaire(kRed, fRed);

            var U = new double[nDdl];
            for (int i = 0; i < nLibre; i++)
                U[ddlLibres[i]] = uRed[i];

            // reactions : R = K * U - F
            var R = new double[nDdl];
            for (int i = 0; i < nDdl; i++)
            {
                double somme = 0;
                for (int j = 0; j < nDdl; j++)
                    somme += K[i, j] * U[j];
                R[i] = somme - F[i];
            }

            return new ResultatChassis { Deplacements = U, Reactions = R };
        }

        // Optimisation de masse

        public class ResultatOptimisation
        {
            public string Label;
            public SectionTube Section;
            public double CsMin;
            public double MasseLin;
        }

        public static ResultatOptimisation OptimiserTube(CasCharge[] casCharges,
            Materiau mat, double csMin = 2.0)
        {
            ResultatOptimisation meilleur = null;

            foreach (var tubeDef in Tubes)
            {
                var sec = CalculerSectionTube(tubeDef.DExt, tubeDef.Ep);
                double csMinTrouve = double.PositiveInfinity;

                foreach (var cas in casCharges)
                {
                    var res = AnalyserTube(sec, cas, mat, KtSoudureBout);
                    csMinTrouve = Math.Min(csMinTrouve, res.Cs);
                }

                if (csMinTrouve >= csMin)
                {
                    if (meilleur == null || sec.MasseLin < meilleur.MasseLin)
                    {
                        meilleur = new ResultatOptimisation
                        {
                            Label = tubeDef.Label,
                            Section = sec,
                            CsMin = csMinTrouve,
                            MasseLin = sec.MasseLin,
                        };
                    }
                }
            }

            return meilleur;
        }
    }
}
