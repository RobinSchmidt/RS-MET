#include "rosic_OthersTests.h"
using namespace rotes;




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

  rosic::fillWithIndex(t, N);
  rosic::scale(t, t, N, 1.0/fs);

  rosic::fillWithZeros(x, N);
  rosic::fillWithZeros(y, N);

  int n;
  for(n = N/3; n < 2*N/3; n++)
    x[n] = 1.0;

  for(n = 0; n < N; n++)
    y[n] = slewRateLimiter.getSample(x[n]);


  Plotter::plotData(N, t, x, y);
}
