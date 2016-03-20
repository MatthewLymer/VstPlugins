using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;

namespace VisualGraph
{
    public partial class Form1 : Form
    {
        private const int SamplingFrequency = 44100;
        private const int NumberOfSamples = 128;
        private const double RadiansPerPhase = Math.PI * 2;

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            chart1.Series.Clear();

            var dry = chart1.Series.Add("dry");
            var wet = chart1.Series.Add("wet");

            dry.ChartType = SeriesChartType.Line;
            wet.ChartType = SeriesChartType.Line;

            var drySamples = GenerateSamples(new[] { 880d}, SamplingFrequency, NumberOfSamples).ToArray();

            for (var i = 0; i < drySamples.Length; i++)
            {
                dry.Points.AddXY(i, drySamples[i].YValues[0]);

                // dry.Points.Add(sample);                
            }

            var wetSamples = drySamples.Select(x => (float) x.YValues[0]).ToArray();

            PitchShifter.PitchShift(1, NumberOfSamples, SamplingFrequency, wetSamples);
            
            for (var i = 0; i < wetSamples.Length; i++)
            {
                wet.Points.AddXY(i, wetSamples[i]);
            }
        }

        private static IEnumerable<DataPoint> GenerateSamples(double[] signalFrequencies, double sampleFrequency, int numberOfSamples)
        {
            var samplingPeriod = 1D / sampleFrequency;

            for (var i = 0; i < numberOfSamples; i++)
            {
                var xValue = samplingPeriod * i;
                var yValue = signalFrequencies.Sum(x => Math.Sin(xValue * RadiansPerPhase * x) /*+ Math.Sin(xValue * RadiansPerPhase * x * 2)*/);
                yield return new DataPoint(xValue, yValue);
            }
        }
    }
}
