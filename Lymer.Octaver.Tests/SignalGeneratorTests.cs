using System;
using NUnit.Framework;

namespace Lymer.Octaver.Tests
{
    class SignalGeneratorTests
    {
        [TestFixture]
        public class WhenTellingSignalGeneratorToGetSample
        {
            [Test]
            [TestCase(1d, 0.0D)]
            [TestCase(1d, 0.5D)]
            [TestCase(1d, 1.0D)]
            [TestCase(1d, 1.5D)]
            [TestCase(2d, 0.0D)]
            [TestCase(2d, 0.25D)]
            [TestCase(2d, 0.5D)]
            [TestCase(2d, 0.75D)]
            [TestCase(2d, 1.0D)]
            [TestCase(2d, 1.25D)]
            public void ShouldReturnZeroSample(double frequency, double offset)
            {
                var sample = SignalGenerator.GetSample(frequency, offset);

                Assert.LessOrEqual(Math.Abs(sample), Math.Pow(0.1D, 14));
            }

            [Test]
            [TestCase(1d, 0.25D)]
            [TestCase(1d, 1.25D)]
            [TestCase(2d, 0.125D)]
            [TestCase(2d, 0.625D)]
            public void ShouldReturnOneSample(double frequency, double offset)
            {
                var sample = SignalGenerator.GetSample(frequency, offset);

                Assert.GreaterOrEqual(sample, 1D - Math.Pow(0.1D, 14));                
            }

            [Test]
            [TestCase(1d, 0.75D)]
            [TestCase(1d, 1.75D)]
            [TestCase(2d, 0.375D)]
            [TestCase(2d, 0.875d)]
            public void ShouldReturnNegativeOneSample(double frequency, double offset)
            {
                var sample = SignalGenerator.GetSample(frequency, offset);

                Assert.LessOrEqual(sample, -1D + Math.Pow(0.1D, 14));
            }
        }
    }

    internal static class SignalGenerator
    {
        private const double RadiansPerPhase = Math.PI*2D;
        
        public static double GetSample(double frequency, double timeOffsetInSeconds)
        {
            return Math.Sin(frequency * timeOffsetInSeconds * RadiansPerPhase);
        }
    }
}
