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
            Tuple<double[], int> W1 = FromWav("D:/lab/Audio/test.wav");
            VersWav(W1.Item2, "D:/lab/Audio/test_ar.wav", W1.Item1);
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
            int numsamples = 8 * (v5 / v4);
            double[] retour = new double[numsamples];
            for (int i = 0; i < numsamples; i++)
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
    }
}
