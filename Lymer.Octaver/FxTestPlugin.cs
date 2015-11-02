using Jacobi.Vst.Core;
using Jacobi.Vst.Framework;
using Jacobi.Vst.Framework.Plugin;

namespace Lymer.Octaver
{
    internal class FxTestPlugin : VstPluginBase
    {
        private readonly FxPluginInterfaceManager _interfaceManager;

        public FxTestPlugin() : base(
            "Octaver Plugin", 
            new VstProductInfo("Octaver", "Matthew Lymer", 1000),
            VstPluginCategory.Effect,
            VstPluginCapabilities.None,
            0,
            1)
        {
            _interfaceManager = new FxPluginInterfaceManager();
        }

        public override bool Supports<T>()
        {
            return _interfaceManager.Supports<T>();
        }

        public override T GetInstance<T>()
        {
            return _interfaceManager.GetInstance<T>();
        }

        public override void Dispose()
        {
            _interfaceManager.Dispose();

            base.Dispose();
        }
    }
}