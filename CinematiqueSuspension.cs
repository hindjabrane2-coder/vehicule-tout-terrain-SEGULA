using System;
using System.Collections.Generic;

namespace VehiculeToutTerrain
{
    // Cinematique de la suspension a double triangulation
    // Debattement, carrossage, variation de voie, centre instantane, centre de roulis,
    // anti-plongee, anti-cabrage, rapport de mouvement et bump steer

    public static class CinematiqueSuspension
    {
        // Geometrie de la suspension (mm) -- plan transversal (x = lateral, y = vertical)

        public static readonly double[] ChassisSupInt = { 150.0, 350.0 };
        public static readonly double[] ChassisSupExt = { 350.0, 360.0 };
        public static readonly double[] ChassisInfInt = { 120.0, 150.0 };
        public static readonly double[] ChassisInfExt = { 380.0, 160.0 };

        public static readonly double[] MoyeuOffset = { 30.0, -80.0 };

        public static readonly double[] RessortAttache = { 250.0, 155.0 };
        public static readonly double[] RessortAncragesChassis = { 250.0, 400.0 };

        public static readonly double[] BieletteRotule = { 370.0, 200.0 };
        public static readonly double[] BieletteCremaillere = { 200.0, 210.0 };

        public const double DemiVoie = 700.0;
        public const double HCdg = 500.0;
        public const double Empattement = 2200.0;
        public const double DistAxeAvant = 1100.0;

        // Geometrie vue de cote (plan longitudinal)
        public static readonly double[] PivotInfAvant = { 0.0, 150.0 };
        public static readonly double[] PivotInfArriere = { 350.0, 155.0 };
        public static readonly double[] PivotSupAvant = { 50.0, 350.0 };
        public static readonly double[] PivotSupArriere = { 300.0, 345.0 };

        // Resultat cinematique

        public class ResultatCinematique
        {
            public double Dz;
            public double[] Moyeu;
            public double[] RotSup;
            public double[] RotInf;
            public double Carrossage;
            public double DeltaVoie;
            public double[] CentreInstantane;
            public double? CentreRoulisH;
            public double DAngleInf;
            public double DAngleSup;
        }

        // Resultats de balayage

        public class DonneesCinematiques
        {
            public double[] Debattements;
            public double[] Carrossages;
            public double[] VariationsVoie;
            public double[] CentresRoulis;
            public double[] MotionRatios;
            public double[] BumpSteer;
        }

        // Fonctions geometriques

        static double Norme(double[] a, double[] b)
        {
            double dx = b[0] - a[0], dy = b[1] - a[1];
            return Math.Sqrt(dx * dx + dy * dy);
        }

        static double[] PositionRotule(double[] pivot, double longueur, double angle)
        {
            return new double[] {
                pivot[0] + longueur * Math.Cos(angle),
                pivot[1] + longueur * Math.Sin(angle)
            };
        }

        static double AngleNominal(double[] pivotInt, double[] pivotExt)
        {
            return Math.Atan2(pivotExt[1] - pivotInt[1], pivotExt[0] - pivotInt[0]);
        }

        static double[] IntersectionDroites(double[] p1, double[] d1, double[] p2, double[] d2)
        {
            double det = d1[0] * d2[1] - d1[1] * d2[0];
            if (Math.Abs(det) < 1e-10) return null;

            double dpx = p2[0] - p1[0], dpy = p2[1] - p1[1];
            double t = (dpx * d2[1] - dpy * d2[0]) / det;
            return new double[] { p1[0] + t * d1[0], p1[1] + t * d1[1] };
        }

        // Recherche de zero par bisection (remplace scipy.optimize.brentq)

        static double Bisection(Func<double, double> f, double a, double b,
                                double tol = 1e-12, int maxIter = 200)
        {
            double fa = f(a), fb = f(b);
            if (fa * fb > 0) throw new InvalidOperationException("Pas de changement de signe");

            for (int i = 0; i < maxIter; i++)
            {
                double mid = (a + b) / 2;
                double fm = f(mid);
                if (Math.Abs(fm) < tol || (b - a) / 2 < tol) return mid;

                if (fa * fm < 0) { b = mid; fb = fm; }
                else { a = mid; fa = fm; }
            }
            return (a + b) / 2;
        }

        // Resolution cinematique principale

        public static ResultatCinematique ResoudreSuspension(double dz,
            double[] supInt = null, double[] supExt = null,
            double[] infInt = null, double[] infExt = null)
        {
            double[] si = supInt ?? ChassisSupInt;
            double[] se = supExt ?? ChassisSupExt;
            double[] ii = infInt ?? ChassisInfInt;
            double[] ie = infExt ?? ChassisInfExt;

            double lSup = Norme(si, se);
            double lInf = Norme(ii, ie);
            double distRot = Norme(se, ie);

            double aSupNom = AngleNominal(si, se);
            double aInfNom = AngleNominal(ii, ie);

            double[] milieuNom = { (se[0] + ie[0]) / 2, (se[1] + ie[1]) / 2 };
            double[] moyeuNom = { milieuNom[0] + MoyeuOffset[0], milieuNom[1] + MoyeuOffset[1] };

            double dAngle;
            try
            {
                dAngle = Bisection(dAngleInf =>
                {
                    double[] rotInfTest = PositionRotule(ii, lInf, aInfNom + dAngleInf);

                    double das;
                    try
                    {
                        das = Bisection(dAngleSup =>
                        {
                            double[] rotSupTest = PositionRotule(si, lSup, aSupNom + dAngleSup);
                            return Norme(rotSupTest, rotInfTest) - distRot;
                        }, -0.5, 0.5);
                    }
                    catch { return 1e6; }

                    double[] rotSupFinal = PositionRotule(si, lSup, aSupNom + das);
                    double[] milieu = { (rotSupFinal[0] + rotInfTest[0]) / 2,
                                        (rotSupFinal[1] + rotInfTest[1]) / 2 };
                    double[] moyeu = { milieu[0] + MoyeuOffset[0], milieu[1] + MoyeuOffset[1] };
                    return moyeu[1] - moyeuNom[1] - dz;
                }, -0.4, 0.4);
            }
            catch { return null; }

            // recalcul final
            double[] rotInf = PositionRotule(ii, lInf, aInfNom + dAngle);

            double dAngleSup2 = Bisection(da =>
            {
                double[] rs = PositionRotule(si, lSup, aSupNom + da);
                return Norme(rs, rotInf) - distRot;
            }, -0.5, 0.5);

            double[] rotSup = PositionRotule(si, lSup, aSupNom + dAngleSup2);
            double[] milieuFinal = { (rotSup[0] + rotInf[0]) / 2, (rotSup[1] + rotInf[1]) / 2 };
            double[] moyeuFinal = { milieuFinal[0] + MoyeuOffset[0], milieuFinal[1] + MoyeuOffset[1] };

            // carrossage
            double fuseX = rotSup[0] - rotInf[0];
            double fuseY = rotSup[1] - rotInf[1];
            double carrossage = Math.Atan2(fuseX, fuseY) * 180.0 / Math.PI;

            // variation de voie
            double deltaVoie = moyeuFinal[0] - moyeuNom[0];

            // centre instantane de rotation
            double[] dirSup = { rotSup[0] - si[0], rotSup[1] - si[1] };
            double[] dirInf = { rotInf[0] - ii[0], rotInf[1] - ii[1] };
            double[] ci = IntersectionDroites(si, dirSup, ii, dirInf);

            // centre de roulis
            double? centreRoulisH = null;
            if (ci != null)
            {
                double[] contactPneu = { moyeuFinal[0] + MoyeuOffset[0], 0.0 };
                double[] dirCr = { contactPneu[0] - ci[0], contactPneu[1] - ci[1] };
                if (Math.Abs(dirCr[0]) > 1e-6)
                {
                    double tCr = (0.0 - ci[0]) / dirCr[0];
                    centreRoulisH = ci[1] + tCr * dirCr[1];
                }
            }

            return new ResultatCinematique
            {
                Dz = dz,
                Moyeu = moyeuFinal,
                RotSup = rotSup,
                RotInf = rotInf,
                Carrossage = carrossage,
                DeltaVoie = deltaVoie,
                CentreInstantane = ci,
                CentreRoulisH = centreRoulisH,
                DAngleInf = dAngle,
                DAngleSup = dAngleSup2,
            };
        }

        // Rapport de mouvement (motion ratio)

        public static double? CalculerMotionRatio(double dz, double pas = 0.5)
        {
            double aInfNom = AngleNominal(ChassisInfInt, ChassisInfExt);
            double brasAttache = Norme(ChassisInfInt, RessortAttache);
            double angleAttache = AngleNominal(ChassisInfInt, RessortAttache);
            double offsetAngle = angleAttache - aInfNom;

            var resRef = ResoudreSuspension(dz);
            var resVoisin = ResoudreSuspension(dz + pas);
            if (resRef == null || resVoisin == null) return null;

            double[] posRef = PositionRotule(ChassisInfInt, brasAttache,
                aInfNom + resRef.DAngleInf + offsetAngle);
            double[] posVoisin = PositionRotule(ChassisInfInt, brasAttache,
                aInfNom + resVoisin.DAngleInf + offsetAngle);

            double dRef = Norme(RessortAncragesChassis, posRef);
            double dVoisin = Norme(RessortAncragesChassis, posVoisin);

            double deltaRessort = dVoisin - dRef;
            if (Math.Abs(pas) < 1e-10) return 0.0;
            return deltaRessort / pas;
        }

        // Bump steer

        public static double? CalculerBumpSteer(double dz)
        {
            var res = ResoudreSuspension(dz);
            var res0 = ResoudreSuspension(0.0);
            if (res == null || res0 == null) return null;

            double[] fusee0 = { res0.RotSup[0] - res0.RotInf[0], res0.RotSup[1] - res0.RotInf[1] };
            double[] fusee = { res.RotSup[0] - res.RotInf[0], res.RotSup[1] - res.RotInf[1] };

            double angle0 = Math.Atan2(fusee0[1], fusee0[0]);
            double angle = Math.Atan2(fusee[1], fusee[0]);
            double rotationFusee = angle - angle0;

            double[] centreFusee0 = { (res0.RotSup[0] + res0.RotInf[0]) / 2,
                                       (res0.RotSup[1] + res0.RotInf[1]) / 2 };
            double[] centreFusee = { (res.RotSup[0] + res.RotInf[0]) / 2,
                                      (res.RotSup[1] + res.RotInf[1]) / 2 };
            double[] offsetBiel = { BieletteRotule[0] - centreFusee0[0],
                                    BieletteRotule[1] - centreFusee0[1] };

            double c = Math.Cos(rotationFusee), s = Math.Sin(rotationFusee);
            double[] offsetRot = { c * offsetBiel[0] - s * offsetBiel[1],
                                   s * offsetBiel[0] + c * offsetBiel[1] };
            double[] bielPos = { centreFusee[0] + offsetRot[0], centreFusee[1] + offsetRot[1] };

            double[] dirBiel = { bielPos[0] - BieletteCremaillere[0],
                                 bielPos[1] - BieletteCremaillere[1] };
            double[] dirBiel0 = { BieletteRotule[0] - BieletteCremaillere[0],
                                  BieletteRotule[1] - BieletteCremaillere[1] };

            double toe = (Math.Atan2(dirBiel[1], dirBiel[0]) -
                          Math.Atan2(dirBiel0[1], dirBiel0[0])) * 180.0 / Math.PI;
            return toe;
        }

        // Anti-plongee (anti-dive)

        public static double AntiPlongeePourcentage()
        {
            double[] dirInf = { PivotInfArriere[0] - PivotInfAvant[0],
                                PivotInfArriere[1] - PivotInfAvant[1] };
            double[] dirSup = { PivotSupArriere[0] - PivotSupAvant[0],
                                PivotSupArriere[1] - PivotSupAvant[1] };

            double[] ci = IntersectionDroites(PivotInfAvant, dirInf, PivotSupAvant, dirSup);
            if (ci == null) return 0.0;

            double dx = ci[0], dy = ci[1];
            if (Math.Abs(dx) < 1e-6) return 100.0;

            double tanTheta = dy / dx;
            double repartitionAv = DistAxeAvant / Empattement;
            return tanTheta * (Empattement * repartitionAv / HCdg) * 100.0;
        }

        // Anti-cabrage (anti-squat)

        public static double AntiCabragePourcentage()
        {
            double[] dirInf = { PivotInfAvant[0] - PivotInfArriere[0],
                                PivotInfAvant[1] - PivotInfArriere[1] };
            double[] dirSup = { PivotSupAvant[0] - PivotSupArriere[0],
                                PivotSupAvant[1] - PivotSupArriere[1] };

            double[] ci = IntersectionDroites(PivotInfArriere, dirInf, PivotSupArriere, dirSup);
            if (ci == null) return 0.0;

            double dx = ci[0], dy = ci[1];
            if (Math.Abs(dx) < 1e-6) return 100.0;

            double tanTheta = dy / Math.Abs(dx);
            double repartitionAr = 1.0 - DistAxeAvant / Empattement;
            return tanTheta * (Empattement * repartitionAr / HCdg) * 100.0;
        }

        // Derivees cinematiques (differences finies centrees)

        public static void DerivesCinematiques(double[] debattements, double[] carrossages,
            double[] variationsVoie, out double[] dCamber, out double[] dVoie)
        {
            int n = debattements.Length;
            dCamber = new double[n];
            dVoie = new double[n];

            for (int i = 0; i < n; i++)
            {
                if (i == 0)
                {
                    double h = debattements[1] - debattements[0];
                    dCamber[i] = (carrossages[1] - carrossages[0]) / h;
                    dVoie[i] = (variationsVoie[1] - variationsVoie[0]) / h;
                }
                else if (i == n - 1)
                {
                    double h = debattements[n - 1] - debattements[n - 2];
                    dCamber[i] = (carrossages[n - 1] - carrossages[n - 2]) / h;
                    dVoie[i] = (variationsVoie[n - 1] - variationsVoie[n - 2]) / h;
                }
                else
                {
                    double h = debattements[i + 1] - debattements[i - 1];
                    dCamber[i] = (carrossages[i + 1] - carrossages[i - 1]) / h;
                    dVoie[i] = (variationsVoie[i + 1] - variationsVoie[i - 1]) / h;
                }
            }
        }

        // Balayage cinematique complet

        public static DonneesCinematiques BalayageCinematique(double dzMin = -80, double dzMax = 80,
            int nPoints = 100, double[] supInt = null, double[] supExt = null,
            double[] infInt = null, double[] infExt = null)
        {
            var donnees = new DonneesCinematiques
            {
                Debattements = new double[nPoints],
                Carrossages = new double[nPoints],
                VariationsVoie = new double[nPoints],
                CentresRoulis = new double[nPoints],
                MotionRatios = new double[nPoints],
                BumpSteer = new double[nPoints],
            };

            for (int i = 0; i < nPoints; i++)
            {
                double dz = dzMin + i * (dzMax - dzMin) / (nPoints - 1);
                donnees.Debattements[i] = dz;
                donnees.Carrossages[i] = double.NaN;
                donnees.VariationsVoie[i] = double.NaN;
                donnees.CentresRoulis[i] = double.NaN;
                donnees.MotionRatios[i] = double.NaN;
                donnees.BumpSteer[i] = double.NaN;

                var res = ResoudreSuspension(dz, supInt, supExt, infInt, infExt);
                if (res != null)
                {
                    donnees.Carrossages[i] = res.Carrossage;
                    donnees.VariationsVoie[i] = res.DeltaVoie;
                    if (res.CentreRoulisH.HasValue)
                        donnees.CentresRoulis[i] = res.CentreRoulisH.Value;
                }

                double? mr = CalculerMotionRatio(dz);
                if (mr.HasValue) donnees.MotionRatios[i] = mr.Value;

                double? bs = CalculerBumpSteer(dz);
                if (bs.HasValue) donnees.BumpSteer[i] = bs.Value;
            }

            return donnees;
        }

        // Export CSV

        public static void ExporterCsv(DonneesCinematiques donnees, string fichier)
        {
            using (var sw = new System.IO.StreamWriter(fichier))
            {
                sw.WriteLine("Debattement_mm;Carrossage_deg;Var_voie_mm;Centre_roulis_mm;Motion_ratio;Bump_steer_deg");
                for (int i = 0; i < donnees.Debattements.Length; i++)
                {
                    sw.WriteLine($"{donnees.Debattements[i]:F2};{donnees.Carrossages[i]:F4};" +
                        $"{donnees.VariationsVoie[i]:F4};{donnees.CentresRoulis[i]:F2};" +
                        $"{donnees.MotionRatios[i]:F4};{donnees.BumpSteer[i]:F4}");
                }
            }
        }

        // Analyse de sensibilite

        public static Dictionary<string, DonneesCinematiques> SensibilitePivot(
            string pivotNom, int indexCoord, double deltaMm = 5.0)
        {
            var pivots = new Dictionary<string, double[]>
            {
                ["sup_int"] = (double[])ChassisSupInt.Clone(),
                ["sup_ext"] = (double[])ChassisSupExt.Clone(),
                ["inf_int"] = (double[])ChassisInfInt.Clone(),
                ["inf_ext"] = (double[])ChassisInfExt.Clone(),
            };

            var resultats = new Dictionary<string, DonneesCinematiques>();

            foreach (var (signe, label) in new[] { (-1, $"-{deltaMm}mm"), (0, "nominal"), (1, $"+{deltaMm}mm") })
            {
                var p = new Dictionary<string, double[]>();
                foreach (var kv in pivots) p[kv.Key] = (double[])kv.Value.Clone();
                p[pivotNom][indexCoord] += signe * deltaMm;

                var data = BalayageCinematique(
                    supInt: p["sup_int"], supExt: p["sup_ext"],
                    infInt: p["inf_int"], infExt: p["inf_ext"]);
                resultats[label] = data;
            }

            return resultats;
        }
    }
}
