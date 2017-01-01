/****************************************************************************
*
* NAME: PitchShift.cs
* VERSION: 1.2
* HOME URL: http://www.dspdimension.com
* KNOWN BUGS: none
*
* SYNOPSIS: Routine for doing pitch shifting while maintaining
* duration using the Short Time Fourier Transform.
*
* DESCRIPTION: The routine takes a pitchShift factor value which is between 0.5
* (one octave down) and 2. (one octave up). A value of exactly 1 does not change
* the pitch. numSampsToProcess tells the routine how many samples in indata[0...
* numSampsToProcess-1] should be pitch shifted and moved to outdata[0 ...
* numSampsToProcess-1]. The two buffers can be identical (ie. it can process the
* data in-place). fftFrameSize defines the FFT frame size used for the
* processing. Typical values are 1024, 2048 and 4096. It may be any value <=
* MAX_FRAME_LENGTH but it MUST be a power of 2. sampleOverlap is the STFT
* oversampling factor which also determines the overlap between adjacent STFT
* frames. It should at least be 4 for moderate scaling ratios. A value of 32 is
* recommended for best quality. sampleRate takes the sample rate for the signal 
* in unit Hz, ie. 44100 for 44.1 kHz audio. The data passed to the routine in 
* indata[] should be in the range [-1.0, 1.0), which is also the output range 
* for the data, make sure you scale the data accordingly (for 16bit signed integers
* you would have to divide (and multiply) by 32768). 
*
* COPYRIGHT 1999-2006 Stephan M. Bernsee <smb [AT] dspdimension [DOT] com>
*
* 						The Wide Open License (WOL)
*
* Permission to use, copy, modify, distribute and sell this software and its
* documentation for any purpose is hereby granted without fee, provided that
* the above copyright notice and this license appear in all source copies. 
* THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF
* ANY KIND. See http://www.dspguru.com/wol.htm for more information.
*
*****************************************************************************/

/****************************************************************************
*
* This code was converted to C# by Michael Knight
* madmik3 at gmail dot com. 
* http://sites.google.com/site/mikescoderama/ 
* 
*****************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;

namespace VisualGraph
{
    public class PitchShifter
    {
        public static void PitchShift(float pitchShift, int sampleRate, float[] indata)
        {
            if (!IsPowerOfTwo(indata.Length))
            {
                throw new ArgumentException(@"indata must be a power of two", nameof(indata));
            }

            var fftBuffer = indata.Zip(
                Enumerable.Repeat(0f, indata.Length), 
                (real, imaginary) => new Complex(real, imaginary)).ToArray();
            
            ShortTimeFourierTransform(fftBuffer, FftDirection.Forward);

            var bins = CalculateBins(sampleRate, fftBuffer);

            var dcOffset = bins[0].Magnitude / fftBuffer.Length;

            var shiftedBins = PitchShiftBins(pitchShift, bins);

            var newBuffer = SynthesizeFft(sampleRate, shiftedBins);
            
            ShortTimeFourierTransform(newBuffer, FftDirection.Inverse);

            var factor = newBuffer.Length / 2;

            for (var i = 0; i < fftBuffer.Length; i++)
            {
                indata[i] = newBuffer[i].Real / factor - dcOffset;
            }
        }

        private static Complex[] SynthesizeFft(float sampleRate, IList<Bin> bins)
        {
            const float expectedPhaseDifference = (float) (2*Math.PI);
            var fftBuffer = new Complex[bins.Count * 2];
            var frequencyPerBin = sampleRate / fftBuffer.Length;

            var phase = 0f;
          
            for (var i = 0; i < bins.Count; i++)
            {
                var tmp = bins[i].Frequency;

                tmp -= i*frequencyPerBin;

                tmp /= frequencyPerBin;

                tmp *= (float) (2*Math.PI);

                tmp += i*expectedPhaseDifference;

                phase += tmp;

                var real = (float)(bins[i].Magnitude * Math.Cos(phase));
                var imaginary = (float)(bins[i].Magnitude * Math.Sin(phase));

                fftBuffer[i] = new Complex(real, imaginary);
            }

            for (var i = bins.Count; i < fftBuffer.Length; i++)
            {
                fftBuffer[i] = new Complex(0f, 0f);
            }

            return fftBuffer;
        }

        private static Bin[] PitchShiftBins(float pitchShift, IList<Bin> bins)
        {
            const int binFactor = 2;

            var shiftedBins = new Bin[bins.Count * binFactor];
            
            for (var i = 0; i < bins.Count; i++)
            {
                var bin = bins[i];
                shiftedBins[i] = new Bin(bin.Frequency * pitchShift, bin.Magnitude * binFactor);
            }

            for (var i = bins.Count; i < shiftedBins.Length; i++)
            {
                shiftedBins[i] = new Bin(0f, 0f);
            }

            return shiftedBins;
        }

        private static Bin[] CalculateBins(int sampleRate, IList<Complex> fftBuffer)
        {
            const float expectedPhaseDifference = (float) (2*Math.PI);

            var frequencyPerBin = sampleRate / (float) fftBuffer.Count;
            var lastPhase = 0f;
            var bins = new Bin[fftBuffer.Count / 2];

            for (var i = 0; i < bins.Length; i++)
            {
                var magnitude = CalculateMagnitude(fftBuffer[i]);
                var phase = CalculatePhase(fftBuffer[i]);

                var tmp = phase - lastPhase;
                lastPhase = phase;

                tmp -= i*expectedPhaseDifference;

                var qpd = (long) (tmp/Math.PI);
                if (qpd >= 0) qpd += qpd & 1;
                else qpd -= qpd & 1;
                tmp -= (float) (Math.PI*qpd);

                tmp /= expectedPhaseDifference;

                tmp = i*frequencyPerBin + tmp*frequencyPerBin;
                
                bins[i] = new Bin(tmp, magnitude);
            }

            return bins;
        }

        private static float CalculatePhase(Complex c)
        {
            return (float)Math.Atan2(c.Imaginary, c.Real);
        }

        private static float CalculateMagnitude(Complex c)
        {
            return (float)Math.Sqrt(c.Real * c.Real + c.Imaginary * c.Imaginary);
        }

        private static void ShortTimeFourierTransform(IList<Complex> fftBuffer, FftDirection fftDirection)
        {
            var buffer = fftBuffer.Select(x => new[] {x.Real, x.Imaginary}).SelectMany(x => x).ToArray();

            ShortTimeFourierTransform(buffer, fftBuffer.Count, fftDirection == FftDirection.Forward ? -1 : 1);

            for (var i = 0; i < fftBuffer.Count; i++)
            {
                fftBuffer[i] = new Complex(buffer[i*2], buffer[i*2+1]);
            }
        }

        private static void ShortTimeFourierTransform(IList<float> fftBuffer, int fftFrameSize, int sign)
        {
            for (var i = 2; i < 2 * fftFrameSize - 2; i += 2)
            {
                var j = 0;
                for (var bitm = 2; bitm < 2 * fftFrameSize; bitm <<= 1)
                {
                    if ((i & bitm) != 0) j++;
                    j <<= 1;
                }

                if (i >= j)
                {
                    continue;
                }

                var temp = fftBuffer[i];
                fftBuffer[i] = fftBuffer[j];
                fftBuffer[j] = temp;
                temp = fftBuffer[i + 1];
                fftBuffer[i + 1] = fftBuffer[j + 1];
                fftBuffer[j + 1] = temp;
            }

            var max = (int)(Math.Log(fftFrameSize) / Math.Log(2.0d) + 0.5d);

            for (int k = 0, le = 2; k < max; k++)
            {
                le <<= 1;
                var le2 = le >> 1;
                var ur = 1.0F;
                var ui = 0.0F;
                var arg = (float)Math.PI / (le2 >> 1);
                var wr = (float)Math.Cos(arg);
                var wi = (float)(sign * Math.Sin(arg));
                for (var j = 0; j < le2; j += 2)
                {
                    float tr;
                    for (var i = j; i < 2 * fftFrameSize; i += le)
                    {
                        tr = fftBuffer[i + le2] * ur - fftBuffer[i + le2 + 1] * ui;
                        var ti = fftBuffer[i + le2] * ui + fftBuffer[i + le2 + 1] * ur;
                        fftBuffer[i + le2] = fftBuffer[i] - tr;
                        fftBuffer[i + le2 + 1] = fftBuffer[i + 1] - ti;
                        fftBuffer[i] += tr;
                        fftBuffer[i + 1] += ti;
                    }
                    tr = ur * wr - ui * wi;
                    ui = ur * wi + ui * wr;
                    ur = tr;
                }
            }
        }

        private static bool IsPowerOfTwo(int x)
        {
            return (x & (x - 1)) == 0;
        }
    }

    internal class Complex
    {
        public Complex(float real, float imaginary)
        {
            Real = real;
            Imaginary = imaginary;
        }

        public float Real { get; }
        public float Imaginary { get; }

        public override string ToString()
        {
            return $"({Real}, {Imaginary})";
        }
    }

    internal struct Bin
    {
        public Bin(float frequency, float magnitude)
        {
            Frequency = frequency;
            Magnitude = magnitude;
        }

        public float Frequency { get; }
        public float Magnitude { get; }

        public override string ToString()
        {
            return $"({Frequency}hz, {Magnitude})";
        }
    }


    internal enum FftDirection
    {
        Forward,
        Inverse
    }
}