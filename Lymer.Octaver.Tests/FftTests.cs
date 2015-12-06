using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.IntegralTransforms;
using NUnit.Framework;

namespace Lymer.Octaver.Tests
{
    [TestFixture]
    class FftTests
    {
        [Test]
        public void DoStuff()
        {
            var signalFrequency = 440D;
            var sampleFrequency = 44100D;
            var numberOfSamples = 1024 * 8;

            var samples = GenerateSamples(signalFrequency, sampleFrequency, numberOfSamples);

            var complexNumbers = samples.Select(x => new Complex(x, 0D)).ToArray();

            Fourier.Forward(complexNumbers);

            var fft = complexNumbers.ToArray();

            var magnitudes = fft.Select(x => x.Magnitude);

            Fourier.Inverse(complexNumbers);
            var inverted = complexNumbers.Select(x => x.Real).ToArray();

            var pairs = samples.Zip(inverted, (d, d1) => new {In = d, Out = d1});

            foreach (var pair in pairs)
            {
                Assert.AreEqual(pair.In, pair.Out, 0.00000001);
            }
        }

        private static double[] CalculateAmplitude(Complex[] complexNumbers, double length)
        {
            return complexNumbers.Select(x => Math.Sqrt(Square(x.Real) + Square(x.Imaginary)) / length).ToArray();
        }

        private static double Square(double d)
        {
            return d*d;
        }

        private static IEnumerable<double> GenerateSamples(double signalFrequency, double sampleFrequency, int numberOfSamples)
        {
            const double radiansPerPhase = Math.PI*2;
            var samplingPeriod = 1D / sampleFrequency;

            for (var i = 0; i < numberOfSamples; i++)
            {
                yield return Math.Sin(samplingPeriod*i*radiansPerPhase*signalFrequency);
            }
        }
    }
}
