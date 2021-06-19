using namespace rotes;
using namespace rosic;

void rotes::testSlewRateLimiterLinear()
{
  // We test the linear slew rate limiter
  // Attack- and release times are 2 and 3 seconds respectively. With a sample-rate of 10 Hz, that
  // translates to 20 and 30 samples to ramp up and down between 0 and 1, respectively.

  double fs      = 10;   // samplerate in Hz
  double attack  = 2000; // attack time in ms
  double release = 3000; // release time in ms


  rosic::SlewRateLimiterLinear slewRateLimiter;

  slewRateLimiter.setSampleRate(fs);
  slewRateLimiter.setAttackTime(attack);
  slewRateLimiter.setReleaseTime(release);

  static const int N = 300;
  double t[N];  // time axis
  double x[N];  // input signal
  double y[N];  // output signal

  using AT = RAPT::rsArrayTools;
  AT::fillWithIndex(t, N);
  AT::scale(t, t, N, 1.0/fs);
  AT::fillWithZeros(x, N);
  AT::fillWithZeros(y, N);

  int n;
  for(n = N/3; n < 2*N/3; n++)
    x[n] = 1.0;

  for(n = 0; n < N; n++)
    y[n] = slewRateLimiter.getSample(x[n]);

  plotData(N, t, x, y);

  // Observations:
  // -The input switches from 0 to 1 at sample 100 (t = 10s) and back to 0 at sample 200 (t = 20s)
  // -The output reaches 1 at sample 120 (t = 12s) and reaches 0 at sample 230 (t = 23s) as it
  //  should be with attack = 2000ms = 2s and release = 3000ms = 3s
}
