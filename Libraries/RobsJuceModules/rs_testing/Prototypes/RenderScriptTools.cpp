
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

void rsNormalizeJointly(std::vector<double>& x, std::vector<double>& y)
{
  double maxX = RAPT::rsMaxAbs(x);
  double maxY = RAPT::rsMaxAbs(y);
  double scaler = 1.0 / RAPT::rsMax(maxX, maxY);
  RAPT::rsScale(x, scaler);
  RAPT::rsScale(y, scaler);
}


template<class T>
T** rsAllocateMatrix(int M, int N)
{
  T* dataFlat = new T[M*N];
  T** data    = new T*[M];
  for(int i = 0; i < M; i++)
    data[i] = &dataFlat[i*N];
  return data;
  // Compare to RAPT::rsMatrixTools::allocateMatrix...merge code. I think, this here is better. It
  // calls new only twice and the data is stored in contiguous memory. But: what if not enough
  // contiguous memory is available? Maybe rosic::readFromWaveFile should use this function
  // internally and make it clear that the function below should be used for cleanup.
}

template<class T>
void rsDeAllocateMatrix(T** A, int M)
{
  delete[] A[0];
  delete[] A;
  // Compare this to what we do in the RS-MET-Tests repo for the sinusoidal analysis
}


void rsConvolveFiles(const char* fileNameX, const char* fileNameH)
{
  // Read input signal:
  int numChannelsX, numFramesX, sampleRateX;
  double** x = rosic::readFromWaveFile(fileNameX, numChannelsX, numFramesX, sampleRateX);

  // Read impulse response:
  int numChannelsH, numFramesH, sampleRateH;
  double** h = rosic::readFromWaveFile(fileNameH, numChannelsH, numFramesH, sampleRateH);

  // Allocate memory for output signal:
  int numChannelsY = rsMax(numChannelsX, numChannelsH);
  int numFramesY   = numFramesX + numFramesH - 1;
  double** y = rsAllocateMatrix<double>(numChannelsY, numFramesY);

  // Do the convolutions of the channels:




  // Clean up allocated memory:
  rsDeAllocateMatrix(x, numChannelsX);
  rsDeAllocateMatrix(h, numChannelsH);
  rsDeAllocateMatrix(y, numChannelsY);

  // ToDo:
  // -Create also the difference signals between the convolution result y and the original input
  //  signals x and h: y-x will be the purely wet output, y-h will be...possibly also interesting.
  //  The purely wet output can be useful for later mixing it in, for example in an sfz file.
  // -The output filenames should be derived from the input filenames, maybe like:
  //  Conv_InName1_InName2.wav
  // -We should use 32 bit float format an make normalization optional...maybe even 64 bit.
}

/*

Ideas:
-Maybe as alternative, use spline envelopes for pitch (as midi note) and volume (in dB)

rsSweepDrummer:
-Try to let the user specify the envelopes in terms of 6 (x,y) datapoints. The envelope has 6 
 degrees of freedom (the 3 weights and the 3 decay times). We need to solve the nonlinear 
 system: 
   y_1 = w1*exp(-d1*x1) + w2*exp(-d2*x1) + w3*exp(-d3*x1) 
   y_2 = w1*exp(-d1*x2) + w2*exp(-d2*x2) + w3*exp(-d3*x2) 
     ...
   y_6 = w1*exp(-d1*x6) + w2*exp(-d2*x6) + w3*exp(-d3*x6) 
 for w1,w2,w3,d1,d2,d3. Maybe have some constraints like d1 < d2 < d3, x1=0. Maybe we nee to solve 
 it iterativley. Maybe start by passing 1 exponential through the average of the first 3 and lst 3
 points, then compute a 2nd envelope by corrceting the error, etc.
-Try to create more interesting attack transients by synthesizing sines with higher harmonics of 
 the present (time-variant) frequencies but (much) faster decay and use the "multiplicative 
 synthesis" ideas with them. This could also be realized as a post-processing step that plays 
 samples at higher speeds. The results will be different from just scaling the frequencies and
 both variants could be useful


*/