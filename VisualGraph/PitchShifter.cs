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
        private const int SampleRate = 44100;

        private const int MaxFrameLength = 2048;
        private static readonly float[] InFifo = new float[MaxFrameLength];
        private static readonly float[] OutFifo = new float[MaxFrameLength];
        private static readonly float[] FftWorkspace = new float[2 * MaxFrameLength];
        private static readonly float[] LastPhase = new float[MaxFrameLength / 2 + 1];
        private static readonly float[] SumPhase = new float[MaxFrameLength / 2 + 1];
        private static readonly float[] OutputAccumilator = new float[2 * MaxFrameLength];
        private static readonly float[] AnalyzedFrequency = new float[MaxFrameLength];
        private static readonly float[] AnalyzedMagnitude = new float[MaxFrameLength];
        private static readonly float[] SynthesizedFrequency = new float[MaxFrameLength];
        private static readonly float[] SynthesizedMagnitude = new float[MaxFrameLength];

        public static void PitchShift(float pitchShift, int numSampsToProcess,
            int sampleRate, float[] indata)
        {
            PitchShift2(pitchShift, numSampsToProcess, 10, sampleRate, indata);
        }

        private static void PitchShift2(float pitchShift, int numSampsToProcess, int sampleOverlap, int sampleRate, IList<float> data)
        {
            var fftBuffer = data.Zip(
                Enumerable.Repeat(0f, data.Count), 
                (real, imaginary) => new Complex(real, imaginary)).ToArray();
            
            ShortTimeFourierTransform(fftBuffer, FftDirection.Forward);

            var bins = CalculateBins(SampleRate, fftBuffer);

            var dcOffset = bins[0].Magnitude / fftBuffer.Length;

            var shiftedBins = PitchShiftBins(pitchShift, bins);

            var newBuffer = SynthesizeFft(SampleRate, shiftedBins);
            
            ShortTimeFourierTransform(newBuffer, FftDirection.Inverse);

            var factor = (newBuffer.Length / 2f);

            for (var i = 0; i < fftBuffer.Length; i++)
            {
                data[i] = newBuffer[i].Real / factor - dcOffset;
            }
        }

        private static Complex[] SynthesizeFft(float sampleRate, Bin[] bins)
        {
            const float expectedPhaseDifference = (float) (2*Math.PI);
            var fftBuffer = new Complex[bins.Length * 2];
            var frequencyPerBin = (sampleRate / (float)fftBuffer.Length);

            var phase = 0f;
          
            for (var i = 0; i < bins.Length; i++)
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

            for (var i = bins.Length; i < fftBuffer.Length; i++)
            {
                fftBuffer[i] = new Complex(0f, 0f);
            }

            return fftBuffer;
        }

        private static Bin[] PitchShiftBins(float pitchShift, Bin[] bins)
        {
            return bins.Select(x => new Bin(x.Frequency*pitchShift, x.Magnitude)).ToArray();
            
            //var shiftedBins = new Bin[bins.Length];

            //for (var i = 0; i < shiftedBins.Length; i++)
            //{
            //    shiftedBins[i] = new Bin(0f, 0f);
            //}

            //for (var i = 0; i < bins.Length; i++)
            //{
            //    var index = (int)(i * pitchShift);

            //    if (index < bins.Length)
            //    {
            //        shiftedBins[i].Magnitude = bins[i].Magnitude;
            //        shiftedBins[i].Frequency = bins[i].Frequency * pitchShift;
            //    }
            //}

            //return shiftedBins;
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

        public static void PitchShift(float pitchShift, int numSampsToProcess, int fftFrameSize,
            int sampleOverlap, int sampleRate, float[] indata)
        {
            var outdata = indata;
            /* set up some handy variables */
            var fftFrameSize2 = fftFrameSize / 2;
            var stepSize = fftFrameSize / sampleOverlap;
            var freqPerBin = sampleRate / (float)fftFrameSize;
            var expct = 2.0 * Math.PI * stepSize / fftFrameSize;
            var inFifoLatency = fftFrameSize - stepSize;
            var rover = inFifoLatency;

            /* main processing loop */
            for (var i = 0; i < numSampsToProcess; i++)
            {
                /* As long as we have not yet collected enough data just read in */
                InFifo[rover] = indata[i];
                outdata[i] = OutFifo[rover - inFifoLatency];
                rover++;

                /* now we have enough data for processing */
                if (rover >= fftFrameSize)
                {
                    rover = inFifoLatency;

                    /* do windowing and re,im interleave */
                    for (var k = 0; k < fftFrameSize; k++)
                    {
                        FftWorkspace[2 * k] = InFifo[k] * CalculateWindow(k / (float)fftFrameSize);
                        FftWorkspace[2 * k + 1] = 0f;
                    }
                    
                    /* ***************** ANALYSIS ******************* */
                    /* do transform */
                    ShortTimeFourierTransform(FftWorkspace, fftFrameSize, -1);

                    /* this is the analysis step */
                    double magn;
                    double phase;
                    double tmp;
                    for (var k = 0; k <= fftFrameSize2; k++)
                    {
                        /* de-interlace FFT buffer */
                        double real = FftWorkspace[2 * k];
                        double imag = FftWorkspace[2 * k + 1];

                        /* compute magnitude and phase */
                        magn = 2.0 * Math.Sqrt(real * real + imag * imag);
                        phase = Math.Atan2(imag, real);

                        /* compute phase difference */
                        tmp = phase - LastPhase[k];
                        LastPhase[k] = (float)phase;

                        /* subtract expected phase difference */
                        tmp -= k * expct;

                        /* map delta phase into +/- Pi interval */
                        var qpd = (long)(tmp / Math.PI);
                        if (qpd >= 0) qpd += qpd & 1;
                        else qpd -= qpd & 1;
                        tmp -= Math.PI * qpd;
                        
                        /* get deviation from bin frequency from the +/- Pi interval */
                        tmp = sampleOverlap * tmp / (2.0 * Math.PI);

                        /* compute the k-th partials' true frequency */
                        tmp = k * freqPerBin + tmp * freqPerBin;

                        /* store magnitude and true frequency in analysis arrays */
                        AnalyzedMagnitude[k] = (float)magn;
                        AnalyzedFrequency[k] = (float)tmp;
                    }

                    /* ***************** PROCESSING ******************* */
                    /* this does the actual pitch shifting */
                    for (var k = 0; k < fftFrameSize; k++)
                    {
                        SynthesizedMagnitude[k] = 0f;
                        SynthesizedFrequency[k] = 0f;
                    }

                    for (var k = 0; k <= fftFrameSize2; k++)
                    {
                        var index = (int)(k * pitchShift);

                        if (index <= fftFrameSize2)
                        {
                            SynthesizedMagnitude[index] += AnalyzedMagnitude[k];
                            SynthesizedFrequency[index] = AnalyzedFrequency[k] * pitchShift;
                        }
                    }

                    /* ***************** SYNTHESIS ******************* */
                    /* this is the synthesis step */
                    for (var k = 0; k <= fftFrameSize2; k++)
                    {
                        /* get magnitude and true frequency from synthesis arrays */
                        magn = SynthesizedMagnitude[k];
                        tmp = SynthesizedFrequency[k];

                        /* subtract bin mid frequency */
                        tmp -= k * freqPerBin;

                        /* get bin deviation from freq deviation */
                        tmp /= freqPerBin;

                        /* take sampleOverlap into account */
                        tmp = 2.0 * Math.PI * tmp / sampleOverlap;

                        /* add the overlap phase advance back in */
                        tmp += k * expct;

                        /* accumulate delta phase to get bin phase */
                        SumPhase[k] += (float)tmp;
                        phase = SumPhase[k];

                        /* get real and imag part and re-interleave */
                        FftWorkspace[2 * k] = (float)(magn * Math.Cos(phase));
                        FftWorkspace[2 * k + 1] = (float)(magn * Math.Sin(phase));
                    }

                    /* zero negative frequencies */
                    for (var k = fftFrameSize + 2; k < 2 * fftFrameSize; k++) FftWorkspace[k] = 0.0F;

                    /* do inverse transform */
                    ShortTimeFourierTransform(FftWorkspace, fftFrameSize, 1);

                    /* do windowing and add to output accumulator */
                    for (var k = 0; k < fftFrameSize; k++)
                    {
                        var window = CalculateWindow(k / (float)fftFrameSize);
                        OutputAccumilator[k] += 2 * window * FftWorkspace[2 * k] / (fftFrameSize2 * sampleOverlap);
                    }
                    for (var k = 0; k < stepSize; k++) OutFifo[k] = OutputAccumilator[k];

                    /* shift accumulator */
                    for (var k = 0; k < fftFrameSize; k++)
                    {
                        OutputAccumilator[k] = OutputAccumilator[k + stepSize];
                    }

                    /* move input FIFO */
                    for (var k = 0; k < inFifoLatency; k++) InFifo[k] = InFifo[k + stepSize];
                }
            }
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

        private static float CalculateWindow(float interpolation)
        {
            return -0.5f * (float)Math.Cos(2.0 * Math.PI * interpolation) + 0.5f;
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

    internal class Bin
    {
        public Bin(float frequency, float magnitude)
        {
            Frequency = frequency;
            Magnitude = magnitude;
        }

        public float Frequency { get; set; }
        public float Magnitude { get; set; }

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