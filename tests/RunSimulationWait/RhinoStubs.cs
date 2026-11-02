// Minimal Rhino stand-ins for the headless RunSimulation wait proof (ticket 08). Only the members
// RunSimulation.cs touches are stubbed; counts let the scenarios observe UI feedback behavior.
namespace Rhino
{
    public static class RhinoApp
    {
        public static int PromptCalls;
        public static string LastPrompt = "";
        public static void SetCommandPrompt(string prompt) { PromptCalls++; LastPrompt = prompt; }
        // Source_Stage_Timing (ticket 02, compiled in from the same linked file) logs through this.
        public static void WriteLine(string message) { }
    }
}

namespace Rhino.ApplicationSettings
{
    public static class FileSettings
    {
        public static bool AutoSaveEnabled;
    }
}

namespace Pachyderm_Acoustic.UI
{
    public class PachydermAc_PlugIn
    {
        public static readonly PachydermAc_PlugIn Instance = new PachydermAc_PlugIn();
        // 2 => Normal priority, so the test process never elevates itself.
        public int TaskPriority = 2;
    }
}
