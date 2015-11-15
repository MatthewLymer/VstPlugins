using Jacobi.Vst.Core;
using Jacobi.Vst.Framework.Plugin;

namespace Lymer.Octaver
{
    class AudioProcessor : VstPluginAudioProcessorBase
    {
        public AudioProcessor() 
            : base(1, 1, 0)
        {

        }

        public override void Process(VstAudioBuffer[] inChannels, VstAudioBuffer[] outChannels)
        {
            var inChannel = inChannels[0];
            var outChannel = outChannels[0];
            
            for (var i = 0; i < inChannel.SampleCount; i++)
            {
                outChannel[i] = inChannel[i];
            }
        }
    }
}
