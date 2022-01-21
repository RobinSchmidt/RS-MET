
/*
std::vector<double> randomizePhases(const std::vector<double>& x, int seed, double amount)
{
  int N = (int)x.size();
  RAPT::rsAssert(RAPT::rsIsPowerOfTwo(N), "Works only for power of 2 lengths");

  // Do forward transform:
  std::vector<double> y(N), a(N), p(N);
  using FT = RAPT::rsFourierTransformerRadix2<double>;   // todo: use Bluestein later
  FT ft;
  ft.setBlockSize(N);
  ft.getRealSignalMagnitudesAndPhases(&x[0], &a[0], &p[0]);

  // Randomize phases:
  RAPT::rsNoiseGenerator<double> ng;
  ng.setSeed(seed);
  ng.setRange(0.0, amount * 2*PI);
  for(int k = 1; k < N; k++) {    // DC is unaffected
    p[k] += ng.getSample();
    p[k] =  RAPT::rsWrapToInterval(p[k], -PI, +PI); }

  // Do inverse transform and return result:
  ft.getRealSignalFromMagnitudesAndPhases(&a[0], &p[0], &y[0]);
  return y;
}
*/


/*

Ideas:
-Maybe as alternative, use spline envelopes for pitch (as midi note) and volume (in dB)

*/