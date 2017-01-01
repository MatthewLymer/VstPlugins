using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace VisualGraph
{
    public partial class Form1 : Form
    {
        private static readonly double[] SignalFrequencies = { 256 };
        private const int SamplingFrequency = 44100;
        private const int NumberOfSamples = 512;
        private const float PitchShiftAmount = 0.5f;

        public Form1()
        {
            InitializeComponent();

            chart1.Series.Clear();

            var dry = chart1.Series.Add("dry");
            var wet = chart1.Series.Add("wet");

            dry.ChartType = SeriesChartType.Line;
            wet.ChartType = SeriesChartType.Line;

            var drySamples = GenerateSamples(SignalFrequencies, SamplingFrequency, NumberOfSamples).ToArray();

            for (var i = 0; i < drySamples.Length; i++)
            {
                dry.Points.AddXY(i, drySamples[i].YValues[0]);       
            }

            var firstWetSamples = drySamples.Select(x => (float)x.YValues[0]).Take(NumberOfSamples / 2).ToArray();
            var secondWetSamples = drySamples.Select(x => (float)x.YValues[0]).Skip(NumberOfSamples / 2).ToArray();

            PitchShifter.PitchShift(PitchShiftAmount, SamplingFrequency, firstWetSamples);
            PitchShifter.PitchShift(PitchShiftAmount, SamplingFrequency, secondWetSamples);

            var allWetSamples = firstWetSamples.Concat(secondWetSamples).ToArray();

            for (var i = 0; i < allWetSamples.Length; i++)
            {
                wet.Points.AddXY(i, allWetSamples[i]);
            }
        }

        private static IEnumerable<DataPoint> GenerateSamples(double[] signalFrequencies, double sampleFrequency, int numberOfSamples)
        {
            var samplingPeriod = 1D / sampleFrequency;

            for (var i = 0; i < numberOfSamples; i++)
            {
                var xValue = samplingPeriod * i;
                var yValue = signalFrequencies.Sum(f => Math.Sin(xValue * (Math.PI * 2) * f) * (1f / signalFrequencies.Length));
                yield return new DataPoint(xValue, yValue);
            }
        }
    }
}
