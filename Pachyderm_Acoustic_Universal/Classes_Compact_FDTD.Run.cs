using System.Threading;

namespace Pachyderm_Acoustic.Numeric.TimeDomain
{
    public partial class Acoustic_Compact_FDTD
    {
        public void RuntoCompletion(CancellationToken cancellationToken = default)
        {
            cancellationToken.ThrowIfCancellationRequested();
            for (time_ms = 0; time_ms < tmax / 1000; time_ms += dt)
            {
                cancellationToken.ThrowIfCancellationRequested();
                Increment();
            }
            cancellationToken.ThrowIfCancellationRequested();
        }
    }
}
