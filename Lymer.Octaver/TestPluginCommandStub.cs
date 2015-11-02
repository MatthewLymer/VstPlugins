using Jacobi.Vst.Framework;
using Jacobi.Vst.Framework.Plugin;

namespace Lymer.Octaver
{
    public class TestPluginCommandStub : StdPluginDeprecatedCommandStub
    {
        protected override IVstPlugin CreatePluginInstance()
        {
            return new FxTestPlugin();
        }
    }
}
