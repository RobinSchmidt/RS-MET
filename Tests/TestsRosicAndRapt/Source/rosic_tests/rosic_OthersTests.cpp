//#include "rosic_OthersTests.h"
using namespace rotes;

//#include "rosic/rosic.h"
//#include "../Shared/Plotting/rosic_Plotter.h"
using namespace rosic;

void rotes::testSlewRateLimiterLinear()
{
  double fs      = 10;   // samplerate in Hz
  double attack  = 1000; // attack time in ms
  double release = 2000; // release time in ms

  // Attack- and release times are 1 and 2 seconds respectively. With a sample-rate of 10 Hz, that
  // translates to 10 and 20 samples to ramp up and down between 0 and 1, respectively

  rosic::SlewRateLimiterLinear slewRateLimiter;

  slewRateLimiter.setSampleRate(fs);
  slewRateLimiter.setAttackTime(attack);
  slewRateLimiter.setReleaseTime(release);

  static const int N = 300;
  double t[N];  // time axis
  double x[N];  // input signal
  double y[N];  // output signal

  RAPT::rsArrayTools::fillWithIndex(t, N);
  RAPT::rsArrayTools::scale(t, t, N, 1.0/fs);

  RAPT::rsArrayTools::fillWithZeros(x, N);
  RAPT::rsArrayTools::fillWithZeros(y, N);

  int n;
  for(n = N/3; n < 2*N/3; n++)
    x[n] = 1.0;

  for(n = 0; n < N; n++)
    y[n] = slewRateLimiter.getSample(x[n]);


  plotData(N, t, x, y);
}
