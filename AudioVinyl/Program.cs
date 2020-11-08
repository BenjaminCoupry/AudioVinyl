using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

namespace VinylAudio
{
    class Program
    {
        static void Main(string[] args)
        {
            Tuple<double[], int> W1 = FromWav("D:/lab/Audio/daft1.wav");
            Random r = new Random(45645);
            //double[] vinyled = EffetVinyle(W1.Item1, W1.Item2, 6, 0.025, ref r);
            //double[] bruited = InstabiliteAmplitude(vinyled, W1.Item2, 0.5, 0.6, ref r);
            //Console.WriteLine("...");
            //bruited = EffetCrepitement(bruited, W1.Item2, 5, 0.5, 0.001, 40, 900, ref r);
            //double[] filtred = FiltrerMinMax(bruited, 100,400, 15000, W1.Item2);
            //double[] bruited = EcrasementBoum(W1.Item1, W1.Item2, 0.1, 5);
            //Variation perlin de l'intensite sonore
            VersWav(W1.Item2, "D:/lab/Audio/test_Shifted.wav",ShiftFreqence(W1.Item1,W1.Item2,500,30,10000,100));
            //VersWav(W1.Item2, "D:/lab/Audio/test_f.wav", filtred);
        }
        static void VersWav(int fe, string Path, double[] inputs)
        {
            uint numsamples = (uint)inputs.GetLength(0);
            ushort numchannels = 1;
            ushort samplelength = 2; // in bytes
            uint samplerate = (uint)fe;

            FileStream f = new FileStream(Path, FileMode.Create);
            BinaryWriter wr = new BinaryWriter(f);
            wr.Write(System.Text.Encoding.ASCII.GetBytes("RIFF"));
            wr.Write((int)(36 + numsamples * numchannels * samplelength));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("WAVE"));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("fmt "));
            wr.Write(16);
            wr.Write((ushort)1);
            wr.Write((ushort)numchannels);
            wr.Write((int)samplerate);
            wr.Write((int)(samplerate * samplelength * numchannels));
            wr.Write((ushort)(samplelength * numchannels));
            wr.Write((ushort)(8 * samplelength));
            wr.Write(System.Text.Encoding.ASCII.GetBytes("data"));
            wr.Write((int)(numsamples * samplelength));
            for (int i = 0; i < numsamples; i++)
            {
                wr.Write(BitConverter.GetBytes((short)(inputs[i] * (double)short.MaxValue)));
            }
        }
        static Tuple<double[], int> FromWav(string Path)
        {
            FileStream f = new FileStream(Path, FileMode.Open);
            BinaryReader wr = new BinaryReader(f);
            //RIFF
            wr.ReadChars(4);
            //36 + numsamples * numchannels * samplelength
            int v1 = wr.ReadInt32();
            //WAVE
            wr.ReadChars(4);
            //"fmt "
            wr.ReadChars(4);
            //16
            wr.ReadInt32();
            //1
            wr.ReadUInt16();
            //num chanels
            ushort channels = wr.ReadUInt16(); ;
            //Sample rate
            int samplerate = wr.ReadInt32();
            //samplerate * samplelength * numchannels
            int v2 = wr.ReadInt32();
            //samplelength * numchannels
            ushort v3 = wr.ReadUInt16();
            //8 * samplelength
            ushort v4 = wr.ReadUInt16();
            //data
            wr.ReadChars(4);
            //numsamples * samplelength
            int v5 = wr.ReadInt32();
            //numSample
            int numsamples = (v5 / (v4/8));
            double[] retour = new double[numsamples/channels];
            for (int i = 0; i < numsamples/channels; i++)
            {
                double moy = 0;
                for (int j = 0; j < channels; j++)
                {
                    short val = wr.ReadInt16();
                    moy += (double)val / (double)short.MaxValue;
                }
                retour[i] = moy / channels;
            }
            return new Tuple<double[], int>(retour, samplerate);
        }

        static double[] Tune(double[] inputs, double fe, double f0, double fmin)
        {
            //TODO (penser a utiliser fmaxintensite.
            //Couper en morceaux, trouver fmaxintensite pour chaque
            
            return new double[0];
        }
        static double[] ShiftFreqence(double[] inputs,double fe, double df,double fmin, double fmax, int ordre)
        {
            //Eloigne le spectre de dF de l'origine
            int n = inputs.GetLength(0);
            double[] retour = new double[n];
            double f0 = 0;
            int n_ = 0;
            double[] travail = FiltrerMinMax(inputs, ordre, fmin, fmax, fe);
            while (f0 < fmax - df)
            {
                n_++;
                Console.WriteLine("f0:" + f0);
                //Selectionner un morceau de spectre de largeur 2df et le deplacer
                double fbas = Math.Max(0, f0 - df);
                double fhaut = Math.Max(0, f0 + df);
                double[] aShift = FiltrerMinMax(travail, ordre, fbas, fhaut, fe);
                double[] shifted = TransposerFreq(aShift, df, fe);
                fbas = Math.Max(0, f0);
                fhaut = Math.Max(0, f0 + 2*df);
                double[] cleaned = FiltrerMinMax(shifted, ordre, fbas, fhaut, fe);
                //Ajouter ce morceau au resultat
                VersWav((int)fe, "D:/lab/Audio/temp/"+n_+".wav", cleaned);
                Console.WriteLine(cleaned.Max());
                for(int i=0;i<cleaned.Length;i++)
                {
                    double added = cleaned[i];
                    retour[i] += added;
                }
                f0 += 2 * Math.Abs(df);
            }
            return FiltrerMinMax(retour,ordre,Math.Max(0,df),fmax,fe);
        }
        static double[] EffetCrepitement(double[] inputs, int fe, double F_Crepit, double Amp_Crepit,double T_crepit,double fmin_bruit,double fmax_bruit, ref Random r)
        {
            int nbCrepit = (int)(F_Crepit * inputs.Length / fe);
            double[] ret = (double[])inputs.Clone();

            for (int i=0;i<nbCrepit;i++)
            {
                double T = (inputs.Length / fe) * r.NextDouble();
                double duree = T_crepit * r.NextDouble();
                double freq = fmin_bruit + (fmax_bruit - fmin_bruit) * r.NextDouble();
                double amplitude = Amp_Crepit * r.NextDouble();
                int i0 = (int)(T * fe);
                int ni = (int)(duree * fe);
                int istart = Math.Max(0, Math.Min(inputs.Length - 1, i0 - ni));
                int ifin = Math.Max(0, Math.Min(inputs.Length - 1, i0 + ni));
                for (int u=istart;u<=ifin;u++)
                {
                    double t = u / (double)fe;
                    double val = Math.Cos(2.0 * Math.PI * freq * (t - T));
                    double fen = 0.5 * (1 + Math.Cos(Math.PI * Math.Abs((t - T) / duree)));
                    double bruit = fen * val * amplitude;
                    ret[u] += bruit;
                }

            }
            return ret;
        }
        static double[] EffetVinyle(double[] inputs, int fe, double F_distorsion, double Amp_distorsion, ref Random r)
        {
            double TDistorsion = 1 / F_distorsion;
            int n = inputs.Count();
            double dt = 1.0 / (double)fe;
            double TempsReel = 0;
            double TempsLu = 0;
            double distorsion = 0;
            List<double> Temps = new List<double>();
            double R0 = 0;
            double R1 = 2.0*(r.NextDouble()-0.5) * Amp_distorsion * dt;
            double Tchang = TDistorsion;
            int nprog=0;
            while(TempsLu*fe<n)
            {
                double k = (nprog*dt) / TDistorsion;
                double ksmooth = 0.5 * (1.0 - Math.Cos(Math.PI * k));
                distorsion = ksmooth * R1 + (1.0 - ksmooth) * R0;
                if(TempsReel>Tchang)
                {
                    nprog = 0;
                    Tchang += TDistorsion;
                    R0 = R1;
                    R1 = 2.0*(r.NextDouble()-0.5) * Amp_distorsion * dt;
                }
                nprog++;
                TempsReel += dt;
                TempsLu += dt + distorsion;
                Temps.Add(TempsLu);
            }
            Console.WriteLine(TempsLu);
            return ParcoursTemporel(inputs, fe, Temps.ToArray());
        }
        static double[] InstabiliteAmplitude(double[] inputs, int fe, double F_distorsion, double Amp_distorsion, ref Random r)
        {
            double A0 = 1.0 / (1.0 + Amp_distorsion);
            double TDistorsion = 1 / F_distorsion;
            int n = inputs.Count();
            double dt = 1.0 / (double)fe;
            double TempsReel = 0;
            double distorsion = 0;
            double R0 = 0;
            double R1 = 2.0 * (r.NextDouble() - 0.5) * Amp_distorsion;
            double Tchang = TDistorsion;
            int nprog = 0;
            double[] resultat = new double[n];
            int nlire = 0;
            while (nlire < n)
            {
                double k = (nprog * dt) / TDistorsion;
                double ksmooth = 0.5 * (1.0 - Math.Cos(Math.PI * k));
                distorsion = ksmooth * R1 + (1.0 - ksmooth) * R0;
                if (TempsReel > Tchang)
                {
                    nprog = 0;
                    Tchang += TDistorsion;
                    R0 = R1;
                    R1 = 2.0 * (r.NextDouble() - 0.5) * Amp_distorsion;
                }
                resultat[nlire] = inputs[nlire] * (1.0 + distorsion)*A0;
                nlire++;
                nprog++;
                TempsReel += dt;
            }
            return resultat;
        }
        static double[] ParcoursTemporel(double[] inputs, int fe,double[] temps)
        {
            double[] resultat = new double[temps.Length];
            int n = inputs.Length;
            for(int i =0;i<temps.Length;i++)
            {
                double t = temps[i];
                int echant = Math.Min(n-1,Math.Max(0,(int)(t*fe)));
                resultat[i] = inputs[echant];
            }
            return resultat;
        }
        static double[] Discretiser(double[] inputs, int nbNivEchant)
        {
            double[] result = new double[inputs.Length];
            for (int i = 0; i < inputs.Length; i++)
            {
                double v = inputs[i];
                result[i] = 1.0/(2.0*nbNivEchant)+Math.Floor(v * nbNivEchant ) / (double)nbNivEchant;
            }
            return result;
        }
        static double[] ReechantillonerDiscret(double[] inputs, int fe, int nouv_fe)
        {
            int ratioEch = fe / nouv_fe;
            double[] result = new double[inputs.Length];
            int u0 = 0;
            double val = inputs[0];
            for (int i = 0; i < inputs.Length; i++)
            {
                result[i] = val;
                u0++;
                if (u0 == ratioEch)
                {
                    u0 = 0;
                    val = inputs[i];
                }
            }
            return result;
        }
        static double[] Filtrer(double[] x, int Ordre, double f0, double deltaf, double fe)
        {
            //Filtre le signal x par bande centrée sur f0 et de largeur deltaf, avec un ordre donné, sachant la frequence d'echantillonage fe 
            double[] filtre = new double[Ordre];
            for (int i = 0; i < Ordre; i++)
            {
                double t = ((double)i - (Ordre / 2.0)) / fe;
                double cos = Math.Cos(2.0 * Math.PI * f0 * t);
                double sincard = 2.0*sinc(Math.PI * deltaf * t) * (deltaf / fe);
                filtre[i] = cos * sincard;
            }
            return CentrerConv(Conv(x, filtre), x.Length);
        }
        static double[] FiltrerMinMax(double[] x, int Ordre, double fmin, double fmax, double fe)
        {
            double m0 = Math.Max(0,Math.Min(fmin,fmax));
            double m1 = Math.Max(0,Math.Max(fmin, fmax));
            return Filtrer(x, Ordre,((m0+m1)/2.0),m1-m0, fe);
        }
        static double[] TransposerFreq(double[] x, double fp, double fe)
        {
            //Transpose le signal x sur fp et -fp, sachant fe la frequence d'echantillonage
            double[] ret = copiedb(x);
            for (int i = 0; i < ret.Length; i++)
            {
                double t = ((double)i - (ret.Length / 2.0)) / fe;
                double cos = Math.Cos(2.0 * Math.PI * fp * t); ;
                ret[i] = ret[i] * cos;
            }
            return ret
;
        }
        static double[] DSP(double[] x, int N)
        {
            //Calcul la dsp de x, en N points, entre 0 et fe
            int n = x.GetLength(0);
            double[] retour = new double[N/2];
            for (int i =0; i < N/2; i++)
            {
                if (i % 100 == 0)
                    Console.WriteLine(i + "/" + (N/2));
                double re = 0;
                double im = 0;
                for (int k = 0; k < n; k++)
                {
                    re += x[k] * Math.Cos(-2.0 * Math.PI * (double)k * (double)i / (double)N);
                    im += x[k] * Math.Sin(-2.0 * Math.PI * (double)k * (double)i / (double)N);
                }
                retour[i] = (re * re + im * im) / ((double)N);
            }

            return retour;
        }
        static double fMaxIntesnsite(double[] x, double fe, double fmin, int resolution)
        {
            //Renvoie la frequence au dela de fmin qui est le plus grand pic de dsp
            double[] DSp = DSP(x, nexpPow2(resolution));
            double dspmax = 0;
            double fmax = fmin;
            for(int i=0;i<DSp.Length;i++)
            {
                double f = (i / (double)DSp.Length) * fe;
                if(f>fmin)
                {
                    if(DSp[i]>dspmax)
                    {
                        dspmax = DSp[i];
                        fmax = f;
                    }
                }
            }
            return fmax;
        }
        static double[] multiplierFreq(double[]x, int fe, double kmult)
        {
            List<double> temps = new List<double>();
            double t = 0;
            double T = x.Length / (double)fe;
            while(t<T)
            {
                temps.Add(t);
                t += kmult * (1.0 / fe);
            }
            return ParcoursTemporel(x, fe, temps.ToArray());
        }
        
        //TODO : echo


        //Outils
        public static int nexpPow2(int N)
        {
            return (int)Math.Pow(2, Math.Floor(Math.Log(N, 2)) + 1);
        }
        public static double sinc(double x)
        {
            //Sinus cardinal de x
            if (x == 0)
            {
                return 1;
            }
            else
            {
                return Math.Sin(x) / x;
            }
        }
        public static double[] Conv(double[] u, double[] v)
        {
            //Convolue u et v
            int n = u.Length;
            int p = v.Length;
            double[] res = new double[n + p - 1];
            for (int i = 0; i < n + p + -1; i++)
            {
                double s = 0;
                for (int j = 0; j < p; j++)
                {
                    int z = i - p + 1 + j;
                    if (z < n && z >= 0)
                    {
                        s += v[j] * u[z];
                    }
                }
                res[i] = s;
            }
            return res;
        }
        public static double[] CentrerConv(double[] c, int n)
        {
            //Centre le retour d'un produit de convolution pour qu'il ait une taille n
            double[] retour = new double[n];
            int delta = (c.Length - n) / 2;
            for (int i = 0; i < n; i++)
            {
                retour[i] = c[i + delta];
            }
            return retour;
        }
        private static double[] copiedb(double[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }
    }
}
