using Jacobi.Vst.Core;
using Jacobi.Vst.Framework.Plugin;

namespace Lymer.Octaver
{
    class AudioProcessor : VstPluginAudioProcessorBase
    {
        private float[] _buffer;
        private int _index;
        private int _bufferLength;

        public AudioProcessor() 
            : base(1, 1, 0)
        {

        }

        public override float SampleRate
        {
            get
            {
                return base.SampleRate;
            }
            set
            {
                _buffer = new float[(int) value / 4];
                _bufferLength = _buffer.Length;
                base.SampleRate = value;
            }
        }

        public override void Process(VstAudioBuffer[] inChannels, VstAudioBuffer[] outChannels)
        {
            var inChannel = inChannels[0];
            var outChannel = outChannels[0];
            
            for (var i = 0; i < inChannel.SampleCount; i++)
            {
                outChannel[i] = ProcessSample(inChannel[i]);
            }
        }

        private float ProcessSample(float sample)
        {
            if (_buffer == null)
            {
                return sample;
            }

            var output = sample + _buffer[_index] * 0.5f;

            _buffer[_index] = sample;

            _index++;

            if (_index >= _bufferLength)
            {
                _index = 0;
            }

            return output;
        }
    }
}
