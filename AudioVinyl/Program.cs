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
            double[] vinyled = EffetVinyle(W1.Item1, W1.Item2, 6, 0.025, ref r);
            double[] bruited = InstabiliteAmplitude(vinyled, W1.Item2, 0.5, 0.6, ref r);
            Console.WriteLine("...");
            bruited = EffetCrepitement(bruited, W1.Item2, 5, 0.5, 0.001, 40, 900, ref r);
            
            //double[] bruited = EcrasementBoum(W1.Item1, W1.Item2, 0.1, 5);
            //Variation perlin de l'intensite sonore
            VersWav(W1.Item2, "D:/lab/Audio/test_ar.wav", bruited);
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
        static double[] EcrasementBoum(double[] inputs, int fe, double T_fenetre, double lambda)
        {
            int nfenetre = (int)(T_fenetre * fe);
            int n = inputs.Length;
            double[] result = new double[n];
            int ifin = Math.Max(0, Math.Min(inputs.Length - 1, nfenetre));
            double vmax = 0;
            int indiceMax = 0;
            for (int k = 0; k <= ifin; k++)
            {
                if (Math.Abs(inputs[k]) > vmax)
                {
                    vmax = Math.Abs(inputs[k]);
                    indiceMax = k;
                }
            }
            for (int i=0;i<n;i++)
            {
                ifin = Math.Max(0, Math.Min(inputs.Length - 1, i+nfenetre));
                if(Math.Abs(inputs[ifin])>vmax)
                {
                    vmax = Math.Abs(inputs[ifin]);
                    indiceMax = ifin;
                }
                if(indiceMax<i-nfenetre)
                {
                    int istart = Math.Max(0, Math.Min(inputs.Length - 1, i-nfenetre));
                    for (int k = istart; k <= ifin; k++)
                    {
                        if (Math.Abs(inputs[k]) > vmax)
                        {
                            vmax = Math.Abs(inputs[k]);
                            indiceMax = k;
                        }
                    }
                }
                double val = inputs[i];
                double prop = val / vmax;
                double nval =Math.Sin(0.5 * Math.PI * prop);
                result[i] = Math.Sign(nval)*Math.Pow(Math.Abs(nval),lambda);
                if(i%1000==0)
                Console.WriteLine((double)i / n);
            }
            return result;
        }
        //A reparer
        static double[] IntensiteSonore(double[] inputs, int fe, double T_fenetre)
        {
            int nfenetre = (int)(T_fenetre * fe);
            double[] retour = new double[inputs.Length];

            int ifin = Math.Max(0, Math.Min(inputs.Length - 1, nfenetre));
            double moy = 0;
            for (int k = 0; k <= nfenetre-1; k++)
            {
                ifin = Math.Max(0, Math.Min(inputs.Length - 1, k + nfenetre));
                int istart = Math.Max(0, Math.Min(inputs.Length - 1, k - nfenetre));
                moy = 0;
                for(int u=istart;u<=ifin;u++)
                {
                    moy += Math.Abs(retour[u]);
                }
                moy = moy / (ifin - istart + 1);
                retour[k] = moy;
            }
            
            for(int i=nfenetre;i<inputs.Length-nfenetre;i++)
            {
                ifin = Math.Max(0, Math.Min(inputs.Length - 1, i + nfenetre));
                int istart = Math.Max(0, Math.Min(inputs.Length - 1, i - nfenetre));
                moy += (Math.Abs(inputs[ifin])-Math.Abs(inputs[istart])) / nfenetre;
                retour[i] = moy;
            }

            for (int k = inputs.Length - nfenetre; k <= inputs.Length-nfenetre - 1; k++)
            {
                ifin = Math.Max(0, Math.Min(inputs.Length - 1, k + nfenetre));
                int istart = Math.Max(0, Math.Min(inputs.Length - 1, k - nfenetre));
                moy = 0;
                for (int u = istart; u <= ifin; u++)
                {
                    moy += Math.Abs(retour[u]);
                }
                moy = moy / (ifin - istart + 1);
                retour[k] = moy;
            }
            return retour;
        }


    }
}
