using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics.IntegralTransforms;

namespace VisualGraph
{
    public partial class Form1 : Form
    {
        private const int SamplingFrequency = 44100;
        private const int NumberOfSamples = SamplingFrequency / 220;
        private const double RadiansPerPhase = Math.PI * 2;

        public Form1()
        {
            InitializeComponent();
        }

        private void Form1_Load(object sender, EventArgs e)
        {
            chart1.Series.Clear();

            var dry = chart1.Series.Add("dry");
            var oc1 = chart1.Series.Add("oct-1");
            var oc2 = chart1.Series.Add("oct-2");

            dry.ChartType = SeriesChartType.Line;
            oc1.ChartType = SeriesChartType.Line;
            oc2.ChartType = SeriesChartType.Line;

            var drySamples = GenerateSamples(new[] {440d}, SamplingFrequency, NumberOfSamples).ToArray();

            foreach (var sample in drySamples)
            {
                dry.Points.Add(sample);
            }

            var dryData = drySamples.Select(s => new Complex(s.YValues[0], 0D)).ToArray();

            var subOctave1Data = SubOctaves(dryData);

            for (var i = 0; i < subOctave1Data.Length; i++)
            {
                oc1.Points.AddXY(drySamples[i].XValue, subOctave1Data[i].Real);
            }

            //var subOctave2Data = SubOctaves(subOctave1Data);

            //for (var i = 0; i < subOctave2Data.Length; i++)
            //{
            //    oc2.Points.AddXY(drySamples[i].XValue, subOctave2Data[i].Real);
            //}
        }

        private static Complex[] SubOctaves(Complex[] data)
        {
            var output = new Complex[data.Length];
            Array.Copy(data, output, data.Length);
            Fourier.Forward(output);
            output = PitchShift(output);
            Fourier.Inverse(output);
            return output;
        }

        private static IEnumerable<DataPoint> GenerateSamples(double[] signalFrequencies, double sampleFrequency, int numberOfSamples)
        {
            var samplingPeriod = 1D / sampleFrequency;

            for (var i = 0; i < numberOfSamples; i++)
            {
                var xValue = samplingPeriod * i;
                var yValue = signalFrequencies.Sum(x => Math.Sin(xValue * RadiansPerPhase * x) * 0.5);
                yield return new DataPoint(xValue, yValue);
            }
        }

        private static Complex[] PitchShift(Complex[] samples)
        {
            var output = new Complex[samples.Length];

            Array.Copy(samples, output, output.Length);

            //for (int i = 0, il = (samples.Count - 1) / 2; i < il; i++)
            //{
            //    var a = samples[i * 2];
            //    var b = samples[i * 2 + 1];

            //    output[i] = a.Magnitude > b.Magnitude ? a : b;
            //}

            return output;
        }
    }
}
