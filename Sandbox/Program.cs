using System;

namespace Sandbox
{
    class Program
    {
        static void Main(string[] args)
        {
            const int fftFrameSize = 32;

            for (var k = 0; k < fftFrameSize; k++)
            {
                var sin = Math.Sin(Math.PI * k / fftFrameSize);
                var cos = -0.5 * Math.Cos(2.0 * Math.PI * k / fftFrameSize) + 0.5;

                Console.WriteLine(cos);
            }

            Console.ReadKey();
        }
    }
}
