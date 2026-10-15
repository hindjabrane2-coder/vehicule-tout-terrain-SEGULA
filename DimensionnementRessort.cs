using System;

namespace VehiculeToutTerrain
{
    // Dimensionnement ressort, amortisseur et barre anti-roulis
    // Frequence propre, taux d'amortissement, transmissibilite,
    // acceleration RMS (confort ISO 2631), simulation temporelle (Heun)

    public static class DimensionnementRessort
    {
        // Parametres du vehicule

        public const double MasseTotale = 450;
        public const double MasseSuspendue = 350;
        public const double MasseNonSuspendue = 25;
        public const double MasseParRoue = MasseSuspendue / 4;
        public const double Empattement = 2200;
        public const double RepartitionAvant = 0.52;
        public const double VoieAvant = 1400;
        public const double VoieArriere = 1380;
        public const double HCdg = 500;
        public const double DebattementComp = 80;
        public const double DebattementDet = 80;
        public const double RatioLevier = 0.70;

        // Frequences cibles
        public const double FCibleAv = 1.3;
        public const double FCibleAr = 1.5;

        // Materiau du ressort (50CrV4)
        public const double GRessort = 78500;
        public const double SigmaAdm = 800;
        public const double SigmaEnd = 550;
        public const double RhoAcier = 7850;

        // Raideur et frequence propre

        public static double RaideurRoue(double masse, double freq)
        {
            double omega = 2 * Math.PI * freq;
            return masse * omega * omega / 1000;
        }

        public static double RaideurRessort(double kRoue, double ratio)
        {
            return kRoue / (ratio * ratio);
        }

        public static double FrequencePropre(double kRoue, double masse)
        {
            return Math.Sqrt(kRoue * 1000 / masse) / (2 * Math.PI);
        }

        // Frequences de tangage et roulis

        public class FrequencesModes
        {
            public double FRebondAv;
            public double FRebondAr;
            public double FTangage;
        }

        public static FrequencesModes FrequencesTangageRoulis(
            double kAv, double kAr, double mAv, double mAr)
        {
            double fRebondAv = FrequencePropre(kAv, mAv / 2);
            double fRebondAr = FrequencePropre(kAr, mAr / 2);

            double lAv = Empattement * RepartitionAvant;
            double lAr = Empattement * (1 - RepartitionAvant);

            double inertieTangage = mAv * Math.Pow(lAv / 1000, 2)
                                  + mAr * Math.Pow(lAr / 1000, 2);
            double kTangage = 2 * kAv * Math.Pow(lAv / 1000, 2)
                            + 2 * kAr * Math.Pow(lAr / 1000, 2);
            double fTangage = Math.Sqrt(kTangage * 1000 / inertieTangage) / (2 * Math.PI);

            return new FrequencesModes
            {
                FRebondAv = fRebondAv, FRebondAr = fRebondAr, FTangage = fTangage
            };
        }

        // Dimensionnement du ressort helicoldal

        public static double FacteurWahl(double C)
        {
            return (4 * C - 1) / (4 * C - 4) + 0.615 / C;
        }

        public class Ressort
        {
            public double DFil;
            public double DMoyen;
            public double DExt;
            public double C;
            public double Kw;
            public double NActives;
            public double NTotal;
            public double L0;
            public double HauteurBloc;
            public double K;
            public double TauMax;
            public double FlecheStat;
            public double FlecheMax;
            public double JeuBloc;
            public double Masse;
        }

        public static Ressort DimensionnerRessort(double kRessort, double chargeStat,
            double dFilInit = 10.0, double cCible = 6.0)
        {
            double dFil = dFilInit;
            double dMoyen = cCible * dFil;
            double Cval = cCible;
            double Kw = FacteurWahl(Cval);

            double courseRessort = DebattementComp / RatioLevier;
            double fMax = chargeStat + kRessort * courseRessort;

            double tauMax = Kw * 8 * fMax * dMoyen / (Math.PI * dFil * dFil * dFil);

            if (tauMax > SigmaAdm)
            {
                dFil = Math.Pow(Kw * 8 * fMax * cCible * dFil / (Math.PI * SigmaAdm), 1.0 / 3);
                dFil = Math.Ceiling(dFil * 2) / 2;
                dMoyen = cCible * dFil;
                Cval = dMoyen / dFil;
                Kw = FacteurWahl(Cval);
                tauMax = Kw * 8 * fMax * dMoyen / (Math.PI * dFil * dFil * dFil);
            }

            double nActives = GRessort * Math.Pow(dFil, 4) /
                              (8 * Math.Pow(dMoyen, 3) * kRessort);
            double nTotal = nActives + 2;
            double hauteurBloc = nTotal * dFil;
            double flecheStat = chargeStat / kRessort;
            double flecheMax = fMax / kRessort;
            double L0 = hauteurBloc + flecheMax + 5.0;
            double jeuBloc = L0 - flecheMax - hauteurBloc;

            double longueurFil = Math.PI * dMoyen * nTotal;
            double volumeFil = Math.PI / 4 * dFil * dFil * longueurFil;
            double masse = volumeFil * RhoAcier / 1e9;

            return new Ressort
            {
                DFil = dFil, DMoyen = dMoyen, DExt = dMoyen + dFil,
                C = Cval, Kw = Kw, NActives = nActives, NTotal = nTotal,
                L0 = L0, HauteurBloc = hauteurBloc, K = kRessort,
                TauMax = tauMax, FlecheStat = flecheStat, FlecheMax = flecheMax,
                JeuBloc = jeuBloc, Masse = masse,
            };
        }

        public static (bool ok, string message) VerifierClash(Ressort ressort)
        {
            if (ressort.JeuBloc <= 0)
                return (false, $"CLASH : jeu negatif ({ressort.JeuBloc:F1} mm)");
            if (ressort.JeuBloc < 3.0)
                return (true, $"Jeu faible ({ressort.JeuBloc:F1} mm), prevoir butee");
            return (true, $"OK (jeu = {ressort.JeuBloc:F1} mm)");
        }

        public static double DureeVieRessort(double tauMax, double tauMin = 0)
        {
            double tauA = (tauMax - tauMin) / 2;
            if (tauA <= SigmaEnd * 0.5) return double.PositiveInfinity;
            return 1e7 * Math.Pow(SigmaEnd * 0.5 / tauA, 8);
        }

        // Ressort a raideur progressive (bi-taux)

        public static double ForceRessortProgressif(double deplacement,
            double k1, double k2, double seuil)
        {
            if (deplacement <= seuil) return k1 * deplacement;
            return k1 * seuil + k2 * (deplacement - seuil);
        }

        public static (double[] x, double[] f) CourbeProgressive(
            double k1, double k2, double seuil, double courseMax)
        {
            int n = 200;
            var x = new double[n];
            var f = new double[n];
            for (int i = 0; i < n; i++)
            {
                x[i] = i * courseMax / (n - 1);
                f[i] = ForceRessortProgressif(x[i], k1, k2, seuil);
            }
            return (x, f);
        }

        // Butee de compression

        public static double ForceButee(double deplacement, double debutContact,
            double kButee, double exposant = 2.5)
        {
            if (deplacement <= debutContact) return 0.0;
            double penetration = deplacement - debutContact;
            return kButee * Math.Pow(penetration, exposant) / Math.Pow(10, exposant - 1);
        }

        // Amortisseur

        public static double AmortissementCritique(double k, double m)
        {
            return 2 * Math.Sqrt(k * 1000 * m);
        }

        public static double ForceAmortisseur(double vitesse, double cComp, double cDet)
        {
            if (vitesse >= 0) return -cComp * vitesse;
            return -cDet * vitesse;
        }

        public class Amortisseur
        {
            public double CTotal;
            public double CComp;
            public double CDet;
            public double CCrit;
            public double Zeta;
            public int DPiston;
            public int DTige;
        }

        public static Amortisseur DimensionnerAmortisseur(double kRoue, double masse,
            double zetaCible = 0.35)
        {
            double cCrit = AmortissementCritique(kRoue, masse);
            double cTotal = zetaCible * cCrit;

            double rapportCd = 0.30;
            double cComp = cTotal * rapportCd;
            double cDet = cTotal * (1 - rapportCd);

            return new Amortisseur
            {
                CTotal = cTotal, CComp = cComp, CDet = cDet,
                CCrit = cCrit, Zeta = zetaCible,
                DPiston = 36, DTige = 16,
            };
        }

        // Transmissibilite

        public static double Transmissibilite(double fExcitation, double fPropre, double zeta)
        {
            double r = fExcitation / fPropre;
            double num = 1 + Math.Pow(2 * zeta * r, 2);
            double den = Math.Pow(1 - r * r, 2) + Math.Pow(2 * zeta * r, 2);
            return Math.Sqrt(num / den);
        }

        // PSD route (ISO 8608)

        public static double PsdRoute(double f, double vKmh = 40, char classe = 'C')
        {
            double gd0;
            switch (classe)
            {
                case 'A': gd0 = 2e-6; break;
                case 'B': gd0 = 8e-6; break;
                case 'D': gd0 = 128e-6; break;
                case 'E': gd0 = 512e-6; break;
                default: gd0 = 32e-6; break;
            }

            double n0 = 0.1;
            double vMs = vKmh / 3.6;
            double n = f / vMs;
            return gd0 * Math.Pow(n / n0, -2) * vMs;
        }

        // Acceleration RMS (confort ISO 2631)

        public class ResultatConfort
        {
            public double Rms;
            public double[] Frequences;
            public double[] AccPsd;
        }

        public static ResultatConfort AccelerationRms(double kRoue, double c, double masse,
            double vKmh = 40, char classe = 'C', double fMax = 30)
        {
            double fPropre = FrequencePropre(kRoue, masse);
            double zeta = c / AmortissementCritique(kRoue, masse);

            int nPts = 500;
            var freqs = new double[nPts];
            var accPsd = new double[nPts];
            double df = (fMax - 0.5) / (nPts - 1);

            for (int i = 0; i < nPts; i++)
            {
                freqs[i] = 0.5 + i * df;
                double T = Transmissibilite(freqs[i], fPropre, zeta);
                double sRoute = PsdRoute(freqs[i], vKmh, classe);
                double omega = 2 * Math.PI * freqs[i];
                accPsd[i] = Math.Pow(omega, 4) * T * T * sRoute;
            }

            double somme = 0;
            for (int i = 0; i < nPts; i++) somme += accPsd[i];
            double rms = Math.Sqrt(somme * df);

            return new ResultatConfort { Rms = rms, Frequences = freqs, AccPsd = accPsd };
        }

        // Barre anti-roulis

        public static double RaideurBarreAntiRoulis(double dBarre, double longueurBras,
            double longueurTorsion, double G = 78500)
        {
            double J = Math.PI * Math.Pow(dBarre, 4) / 32;
            double kTorsion = G * J / longueurTorsion;
            return kTorsion / (longueurBras * longueurBras);
        }

        public static double DimensionnerBarre(double objectifKRoll,
            double longueurBras, double longueurTorsion)
        {
            double d4 = objectifKRoll * 32 * longueurTorsion * longueurBras * longueurBras
                        / (78500 * Math.PI);
            return Math.Pow(d4, 0.25);
        }

        // Ressort a lames

        public static double RaideurLame(double E, double b, double h, double L, int nLames)
        {
            double iLame = b * h * h * h / 12;
            return 3 * E * iLame * nLames / (2 * L * L * L);
        }

        // Ressort pneumatique

        public static double RaideurPneumatique(double pressionBar, double airePistonCm2,
            double volumeCm3, double gamma = 1.4)
        {
            double p = pressionBar * 1e5;
            double A = airePistonCm2 * 1e-4;
            double V = volumeCm3 * 1e-6;
            double k = gamma * p * A * A / V;
            return k / 1000;
        }

        // Simulation temporelle (Heun)

        public class ReponseTemporelle
        {
            public double[] Temps;
            public double[] Deplacement;
            public double OmegaN;
            public double Zeta;
        }

        public static ReponseTemporelle ReponseEchelon(double k, double c, double m,
            double x0 = 30.0, double tMax = 2.0, double dt = 0.0005)
        {
            double omegaN = Math.Sqrt(k * 1000 / m);
            double zeta = c / (2 * m * omegaN);

            int nPts = (int)(tMax / dt);
            var t = new double[nPts];
            var x = new double[nPts];
            var v = new double[nPts];
            x[0] = x0;

            for (int i = 0; i < nPts; i++) t[i] = i * dt;

            for (int i = 0; i < nPts - 1; i++)
            {
                double a1 = (-k * x[i] - c * v[i] / 1000) / (m / 1000);
                double vPred = v[i] + a1 * dt;
                double xPred = x[i] + v[i] * dt;

                double a2 = (-k * xPred - c * vPred / 1000) / (m / 1000);
                v[i + 1] = v[i] + 0.5 * (a1 + a2) * dt;
                x[i + 1] = x[i] + 0.5 * (v[i] + vPred) * dt;
            }

            return new ReponseTemporelle
            {
                Temps = t, Deplacement = x, OmegaN = omegaN, Zeta = zeta
            };
        }
    }
}
