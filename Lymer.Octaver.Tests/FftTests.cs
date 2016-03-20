using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography;
using MathNet.Numerics.IntegralTransforms;
using NUnit.Framework;

namespace Lymer.Octaver.Tests
{
    [TestFixture]
    class FftTests
    {
        private const int NumberOfSamples = 1024;
        private const double SampleFrequency = 1024D;
        private static readonly double[] Scales = Fourier.FrequencyScale(NumberOfSamples, SampleFrequency);

        [Test]
        public void DoStuff()
        {
            var factor = Math.Sqrt(NumberOfSamples) / 2;

            var samples = GenerateSamples(new[] { 440d}, SampleFrequency, NumberOfSamples);

            var complexNumbers = samples.Select(x => new Complex(x, 0D)).ToArray();
            
            Fourier.Forward(complexNumbers);

            var fft = complexNumbers.ToArray();

            var pitchShiftedFft = PitchShift(fft);

            var magnitudes = fft.Select(x => x.Magnitude / factor).ToArray();
            var pitchShiftedMangitudes = pitchShiftedFft.Select(x => x.Magnitude / factor).ToArray();

            Fourier.Inverse(pitchShiftedFft);
            var inverted = complexNumbers.Select(x => x.Real).ToArray();

            var pairs = samples.Zip(inverted, (d, d1) => new {In = d, Out = d1});

            foreach (var pair in pairs)
            {
                Console.WriteLine(pair.Out);

                Assert.AreEqual(pair.In, pair.Out, 0.00000001);
            }
        }

        private static Complex[] PitchShift(Complex[] samples)
        {
            var output = new Complex[samples.Length];

            for (int i = 0, il = (samples.Length - 1) / 2; i < il; i++)
            {
                output[i] = (samples[i*2] + samples[i*2 + 1]) / 2;
            }

            return output;
        }

        private static IEnumerable<double> GenerateSamples(double[] signalFrequencies, double sampleFrequency, int numberOfSamples)
        {
            const double radiansPerPhase = Math.PI*2;
            var samplingPeriod = 1D / sampleFrequency;

            for (var i = 0; i < numberOfSamples; i++)
            {
                yield return signalFrequencies.Sum(x => Math.Sin(samplingPeriod*i*radiansPerPhase*x));
            }
        }
    }
}
