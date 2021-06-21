using namespace rotes;
using namespace rosic;
using namespace RAPT;


template<class T>
void rsFGHT3(T* A, int N, T a, T b, T c, T d, T e, T f, T g, T H, T I)
{
  // Experimental - we try to do something similar with a 3x3 seed matrix. We may base this on 3D
  // rotations with or without reflections. We may have 3 independent parameters to control the
  // diffusion - we could just give those 3 parameters directly to the user but maybe we can 
  // somehow give more meaningful parameters like totalDiffusion, diffusionSpread, ...
  // i suppose, the case with 27 delaylines is most useful for reverb, maybe 9 is also ok
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 3*h) {
      for(int j = i; j < i+h; j++) {
        T x = A[j+0*h];
        T y = A[j+1*h];
        T z = A[j+2*h];
        A[j+0*h] = a*x + b*y + c*z;
        A[j+1*h] = d*x + e*y + f*z;
        A[j+2*h] = g*x + H*y + I*z;  }}
    h *= 3;  }
}
// ok - works - move to rapt
// can we generalize this further for computing a matrix vector product:
// y = A°B°C°... * x, where ° denotes the kronecker product? the function could take as input the 
// in/out vector and a vector of pointers to matrices..maybe the A,B,C etc matrices do not even 
// need to be square matrices as long as the final product is a square matrix with the right 
// dimensions? maybe call it fast kronecker product transform 
// rsFKPT(T* A, const std:vector<rsMatrix*>& M

// But first generalize to an NxN seed matrix that gets kroneckered with itself
// ...here it is (needs more tests):
template<class T>
void rsFastKroneckerTrafo(std::vector<T>& x, const rsMatrix<T>& M, int L)
{
  rsAssert(M.isSquare());
  int B = M.getNumRows();    // basis, size of the seed matrix M
  int N = (int) x.size();    // size of M^L where ^ means repeated Kronecker product with itself
  rsAssert(N == pow(B, L));  // vector has wrong size
  std::vector<T> t(B);       // temp vector - todo: let user pass workspace
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += B*h) {
      for(int j = i; j < i+h; j++) {
        for(int k = 0; k < B; k++)           // establish temp vector
          t[k] = x[j+k*h];
        for(int k = 0; k < B; k++) {
          x[j+k*h] = 0;
          for(int m = 0; m < B; m++)
            x[j+k*h] += M(k, m) * t[m]; }}}  // accumulate B elements of output vector
    h *= B; }
}

/** L is the order */
bool testHadamar2x2(int L)
{
  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;

  int N = pow(2, L);                         // size of matrix is NxN
  double a = 2, b = 3, c = 5, d = 7;         // coefficients of seed matrix

  Vec x = rsRandomIntVector(N, -9, +9);
  Vec y, z;
  Mat H1(2, 2, {a,b,c,d});                   // seed matrix
  Mat HN = H1;
  for(int j = 1; j < L; j++)
    HN = Mat::getKroneckerProduct(H1, HN);
  y = HN * x;                                // reference output
  z = x;
  rsFGHT(&z[0], N, a,b,c,d);                 // output of fast transform
  Vec z2 = x;
  rsFastKroneckerTrafo(z2, H1, L);           // output of fast generalized transform
  return z == y && z2 == y;
}

bool testHadamar3x3(int L)
{
  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;

  int N = pow(3, L);                         // size of matrix is NxN
  double a =  +2, b =  -7, c =   3,          // coefficients of seed matrix
         d =  -5, e = +13, f = -11,
         g = +17, h = -23, i =  19;

  Vec x = rsRandomIntVector(N, -9, +9);
  Vec y, z;
  Mat H1(3, 3, {a,b,c,d,e,f,g,h,i});          // seed matrix
  Mat HN = H1;
  for(int j = 1; j < L; j++)
    HN = Mat::getKroneckerProduct(H1, HN);
  y = HN * x;                                 // reference output
  z = x;
  rsFGHT3(&z[0], N, a,b,c,d,e,f,g,h,i);       // output of fast transform
  Vec z2 = x;
  rsFastKroneckerTrafo(z2, H1, L);           // output of fast generalized transform
  return z == y && z2 == y;
}

bool rotes::testFastGeneralizedHadamardTransform()
{
  // ToDo: move the old implementation FDN::fastGeneralizedHadamardTransform into prototypes..maybe
  // it can eventually be deleted completely when it's clear that the new implementation does
  // the smae thing and is more efficient

  bool ok = true;

  // 4D vector:
  double x4[4] = {4, -8, 12, -4};
  double y4[4];
  double work[16];  // workspace

  typedef rosic::FeedbackDelayNetwork FDN;
  typedef RAPT::rsArrayTools AT;

  // 4-point FWHT:
  AT::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work);
  ok &= y4[0] ==   4;
  ok &= y4[1] ==  28;
  ok &= y4[2] == -12;
  ok &= y4[3] == - 4;

  // 4-point FGHT:
  AT::copy(x4, y4, 4);
  FDN::fastGeneralizedHadamardTransform(y4, 4, 2, work, 2, 3, 5, 7);
  ok &= y4[0] ==  4;
  ok &= y4[1] == 24;
  ok &= y4[2] ==  4;
  ok &= y4[3] == 44;

  // New implementation:
  AT::copy(x4, y4, 4);
  //fght(y4, 4, 2., 3., 5., 7.);
  RAPT::rsFGHT(y4, 4, 2., 3., 5., 7.); // linker error - needs instantiation
  ok &= y4[0] ==  4;
  ok &= y4[1] == 24;
  ok &= y4[2] ==  4;
  ok &= y4[3] == 44;


  // 8D vector:
  double x8[8] = {1, 4, -2, 3, 0, 1, 4, -1};
  double y8[8];

  // 8-point FWHT:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work);
  ok &= y8[0] ==  10;
  ok &= y8[1] == - 4;
  ok &= y8[2] ==   2;
  ok &= y8[3] == - 4;
  ok &= y8[4] ==   2;
  ok &= y8[5] == -12;
  ok &= y8[6] ==   6;
  ok &= y8[7] ==   8;

  // 8-point FGHT:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, 7);
  ok &= y8[0] ==   149;
  ok &= y8[1] ==   357;
  ok &= y8[2] ==   360;
  ok &= y8[3] ==   862;
  ok &= y8[4] ==   362;
  ok &= y8[5] ==   866;
  ok &= y8[6] ==   875;
  ok &= y8[7] ==  2092;

  // New implementation:
  AT::copy(x8, y8, 8);
  RAPT::rsFGHT(y8, 8, 2., 3., 5., 7.);
  ok &= y8[0] ==   149;
  ok &= y8[1] ==   357;
  ok &= y8[2] ==   360;
  ok &= y8[3] ==   862;
  ok &= y8[4] ==   362;
  ok &= y8[5] ==   866;
  ok &= y8[6] ==   875;
  ok &= y8[7] ==  2092;

  double x16[16] = { 1, 4, -2, 3, 0, 1, 4, -1, -1, -2, -1, -3, 2, 5, 1, -2, };
  double y16[16];
  double z16[16];
  AT::copy(x16, y16, 16);
  FDN::fastGeneralizedHadamardTransform(y16, 16, 4, work, 2, 3, 5, 7);
  AT::copy(x16, z16, 16);
  RAPT::rsFGHT(z16, 16, 2., 3., 5., 7.);
  ok &= AT::equal(y16, z16, 16);
  // ok: y16 and z16 are the same - why does it not work within the FDN?

  // 8D Forward/backward trafo - check if input is reconstructed:
  AT::copy(x8, y8, 8);
  FDN::fastGeneralizedHadamardTransform(       y8, 8, 3, work, 2, 3, 5, -7);
  FDN::fastInverseGeneralizedHadamardTransform(y8, 8, 3, work, 2, 3, 5, -7);
  ok &= fabs(RAPT::rsArrayTools::maxDeviation(x8, y8, 8)) < 1.e-15;

  // Tests via Kronekcer products:
                            // # delaylines
  ok &= testHadamar2x2(1);  //   2
  ok &= testHadamar2x2(2);  //   4
  ok &= testHadamar2x2(3);  //   8
  ok &= testHadamar2x2(4);  //  16
  ok &= testHadamar2x2(5);  //  32
  ok &= testHadamar2x2(6);  //  64
  ok &= testHadamar2x2(7);  // 128
  ok &= testHadamar2x2(8);  // 256
  ok &= testHadamar3x3(1);  //   3
  ok &= testHadamar3x3(2);  //   9
  ok &= testHadamar3x3(3);  //  27
  ok &= testHadamar3x3(4);  //  81
  ok &= testHadamar3x3(5);  // 243

  return ok;

  // ToDo: 
  // -compare results of bigger sizes with explicit matrix multiplication of matrices created by 
  //  the Sylvester construction (use the Kronecker-product in rsMatrix for this)
  // -make a performance comparison betweenold and new implementation
}

bool rotes::testFeedbackDelayNetwork()
{
  bool result = true;  // get rid! this is not a unit test!

  //FeedbackDelayNetwork16 *fdn16 = new FeedbackDelayNetwork16;

  double amplitude = 0.5;       // amplitude of the input impulse
  double diffusion = 120.0;     // diffusion parameter in percent


  using AT = RAPT::rsArrayTools;

  FeedbackDelayNetwork fdn;
  fdn.setDiffusion(diffusion);

  static const int N = 100000;
  //double hL[N], hR[N];
  double *hL = new double[N];
  double *hR = new double[N];


  AT::fillWithImpulse(hL, N, amplitude);
  AT::fillWithImpulse(hR, N, amplitude);
  for(int n = 0; n < N; n++)
    fdn.processFrame(&hL[n], &hR[n]);

  double t[N];
  AT::fillWithIndex(t, N);
  //Plotter::plotData(N, t, hL, hR);
  //Plotter::plotData(N, t, hL);


  //rosic::writeToStereoWaveFile("d:\\TmpData\\FDNImpulseResponse.wav", hL, hR, N, 44100, 16);
  rosic::writeToStereoWaveFile("FDNImpulseResponse.wav", hL, hR, N, 44100, 16);
  printf("%s", "Rendering FDNImpulseResponse.wav done\n");

  delete[] hL;
  delete[] hR;
  return result;

  // Observations:
  // -With diffusion set to 100, we indeed get some sort of exponentially decaying white noise, as
  //  it should be
  // -with D=200, the left channel is generally louder than the right - maybe our output vector
  //  is bad for the given distribution of delaytimes? if no soultion can be found, we can also
  //  just use the two output channels as mid and side wet signal...that might be a good idea 
  //  anyway because it makes the whole thing more robust aginst such things...and then we can also
  //  adjust the gain for mid/side separately
  // -The diffusion parameter seems to need a nonlinear mapping that gives more precision towards
  //  higher values - between D=60 and D=100, there is not so much difference

  // ToDo:
  // -Figure out the distribution of the noise (we need an infinite decay-time for this). It is 
  //  written in the literature that exponentially decaying Gaussian white noise sounds best for
  //  a reverb impulse response. Figure out, if it is indeed Gaussian. Maybe it's an Irwin-Hall 
  //  distribution of order equal to the number of the delaylines? That would seem to make some 
  //  sense.
  // -Maybe we can render Gaussian white noise impulse responses with time-variant slope filters
  //  whose slope increases over time. Maybe we can make a convolution reverb that internally 
  //  renders impulse response according to that idea.
  // -Try diffusion less than zero and greater than 100
  // -modulate the diffusion, using a filtered (and maybe levelled) version of the output signal
  // -make an APE project to experiment with the algo and its parameters
  // -figure out the amplitude distribution of the white noise when decay time is infinite
  //  is it Gaussian? if not, can we make it so, maybe by adding outputs of multiple FDNs with
  //  different settings? or maybe just some low order allpass filter could do the job? -> figure
  //  out, what an allpass does to the amplitude distribution of uniform white noise and/or
  //  look at the impulse response
}

template<class Effect>
bool testInOutEqual(Effect& eff, int numSamples, double tolerance)
{
  bool result = true;
  RAPT::rsNoiseGenerator<double> ng;
  double xL, xR, yL, yR;
  for(int n = 0; n < numSamples; n++)
  {
    yL = xL = ng.getSample();
    yR = xR = ng.getSample();
    eff.getSampleFrameStereo(&yL, &yR);
    result &= RAPT::rsIsCloseTo(yL, xL, tolerance);
    result &= RAPT::rsIsCloseTo(yR, xR, tolerance);
    rsAssert(result == true);
  }
  return result;
}

bool rotes::testMultiComp()
{  
  // We check the multiband band compressor with neutral compressor settings for each band and 
  // various splitting configurations. In any case, the output signal should equal the input
  // signal.

  bool result = true;

  int N = 200;  // number of samples

  rosic::rsMultiBandCompressor mbc;
  double tol = 1.e-14;

  result &= mbc.getNumberOfBands() == 1;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 8000.0);  // maybe pass only the frequency, not the index, 
                              // maybe rename addSplit, maybe define splitters in terms of the 
                              // highpass freq (allows to have 0 freq, when numBands == 1)

  result &= mbc.getNumberOfBands() == 2;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 4000.0);
  result &= mbc.getNumberOfBands() == 3;
  result &= testInOutEqual(mbc, N, tol);


  mbc.insertBand(0, 2000.0);
  result &= mbc.getNumberOfBands() == 4;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 1000.0);
  result &= mbc.getNumberOfBands() == 5;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 500.0);
  result &= mbc.getNumberOfBands() == 6;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 250.0);
  result &= mbc.getNumberOfBands() == 7;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 125.0);
  result &= mbc.getNumberOfBands() == 8;
  result &= testInOutEqual(mbc, N, tol);

  mbc.insertBand(0, 62.5);
  result &= mbc.getNumberOfBands() == 9;
  result &= testInOutEqual(mbc, N, tol);

  // todo: 
  // -use loop for adding bands 
  // -test removing bands (also using a loop)




  return result;
}

