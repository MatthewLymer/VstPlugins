using Jacobi.Vst.Framework;
using Jacobi.Vst.Framework.Plugin;

namespace Lymer.Octaver
{
    class FxPluginInterfaceManager : PluginInterfaceManagerBase
    {
        protected override IVstPluginAudioProcessor CreateAudioProcessor(IVstPluginAudioProcessor instance)
        {
            return instance ?? new AudioProcessor();
        }
    }
}
