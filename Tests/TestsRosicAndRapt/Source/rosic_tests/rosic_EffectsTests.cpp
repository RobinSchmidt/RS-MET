using namespace rotes;
using namespace rosic;
using namespace RAPT;

// Move to prototypes or RAPT:

template<class T>
void rsFGHT3(T* v, int N, T a, T b, T c, T d, T e, T f, T g, T H, T I)
{
  rsMatrix3x3<T> A(a, b, c, d, e, f, g, H, I);
  //rsFGHT3(v, N, A);
  rsLinearTransforms::kronecker3x3(v, N, A);
  return;


  // Experimental - we try to do something similar with a 3x3 seed matrix. We may base this on 3D
  // rotations with or without reflections. We may have 3 independent parameters to control the
  // diffusion - we could just give those 3 parameters directly to the user but maybe we can 
  // somehow give more meaningful parameters like totalDiffusion, diffusionSpread, ...
  // i suppose, the case with 27 delaylines is most useful for reverb, maybe 9 is also ok
  int h = 1;
  while(h < N) {
    for(int i = 0; i < N; i += 3*h) {
      for(int j = i; j < i+h; j++) {
        T x = v[j+0*h];
        T y = v[j+1*h];
        T z = v[j+2*h];
        v[j+0*h] = a*x + b*y + c*z;
        v[j+1*h] = d*x + e*y + f*z;
        v[j+2*h] = g*x + H*y + I*z;  }}
    h *= 3;  }
}
// -ok - works - move to rapt
// -use it in the FDN, based on a rotation matrix defined in terms of Euler angles, maybe using
//  rsRotationXYZ for the computation of the coeffs
// -can we generalize this further for computing a matrix vector product:
//  y = A°B°C°... * x, where ° denotes the kronecker product? the function could take as input the 
//  in/out vector and a vector of pointers to matrices..maybe the A,B,C etc matrices do not even 
//  need to be square matrices as long as the final product is a square matrix with the right 
//  dimensions? maybe the final matrix doesn't even need to be square? maybe call it fast kronecker
//  product transform rsFKPT(T* A, const std:vector<rsMatrix*>& M
// -maybe the seed matrix should be given as a reference to rsMatrix3x3
// -ToDo: make benchmarks measuring the version taking the 3x3 matrix versus the version 
//  taking a,b,c,d,... maybe with the latter, the compiler is more inclined to use registers for 
//  the arguments? but with so many arguments, that might not be possible

// But first generalize to an NxN seed matrix that gets kroneckered with itself L times
// ...here it is - seems to work:
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
// maybe for this, the algo that uses N memory cells for temp storage is actually better, this here
// uses log(N) extra memory, but that makes the inner loops more complicated

// computes M[0] ° M[1] ° M[2] ° ... ° M[L-1] * x 
// ..or M[L-1] ° M[L-2] ° M[L-3] ° ... ° M[0] * x 
// currently works only when for all matrices dim(input) >= dim(output)
template<class T>
void rsFastKroneckerTrafo(std::vector<T>& x, const std::vector<const rsMatrix<T>*> M)
{
  int Nx = 1; for(size_t m = 0; m < M.size(); m++) Nx *= M[m]->getNumColumns();
  int Ny = 1; for(size_t m = 0; m < M.size(); m++) Ny *= M[m]->getNumRows(); 
  int Nt = 0; for(size_t m = 0; m < M.size(); m++) Nt  = rsMax(Nt, M[m]->getNumColumns());
  std::vector<T> t(Nt);     // temp storage
  int N = rsMax(Nx, Ny);
  x.resize(N);
  int hC = 1;
  int hR = 1;
  for(size_t m = 0; m < M.size(); m++) {   // m: matrix index, loop over the matrices
    int nR = M[m]->getNumRows();
    int nC = M[m]->getNumColumns();
    //int nM = rsMin(nC, nR);
    //int hMx = rsMax(hC, hR);   // 
    //int hMn = rsMin(hC, hR); 
    for(int i = 0; i < Nx; i += nC*hC) {   // i: base read index into input vector
      for(int j = i; j < i+hR; j++) {      // j: base write index into output vector
        for(int k = 0; k < nC; k++)        // we must store 2 temps
          t[k] = x[j+k*hC];
        for(int k = 0; k < nR; k++) {      // we must write 3 outputs
          x[j+k*hR] = 0;
          for(int n = 0; n < nC; n++)      // each output must accumulate 2 values
            x[j+k*hR] += M[m]->at(k, n) * t[n]; }}}
    hC *= nC;     // old, not sure
    //hC *= nM;
    hR *= nR; }
  x.resize(Ny);
  // The comments "must store 2 temps" etc. refer to an attempt to use a product of 2 3x2 matrices,
  // for which the function does not yet work. Something is still wrong when 
  // dim(input) < dim(output) for any of the factor matrices, i guess. If the output vector is 
  // larger then the input vector, it may be plausibe, that the j-loop must iterate over more 
  // elements...maybe j = i; j < i+max(hR,hC); j++

  // hM = rsMax(hC,hR):
  // hR,hR: 3x2 fails, others ok
  // hR,hM: all fail
  // hM,hR: all fail
  // hM,hM: all fail
  //
  // hM = rsMin(hC,hR):
  // hR,hR: 3x2 fails, others ok
  // hR,hM: 3x2 fails, others ok
  // hM,hR: 3x2 fails, others ok
  // hM,hM: 3x2 fails, others ok

  // hMx,hMn: all fail
  // hMn,hMx: all fail

  // -maybe the j-loop should use max(hC,hR) but the offset should use min(hC,hR)
  // -maybe it cannot possibly work, because we store just 2 temps but overwrite 3 - we don't store 
  //  enough temp data -> implement the algo with swap/workspace-vector - this is also simpler

  // ...hM doesn't seem to be useful
  // maybe we need to use hM, hM, together with running the i-loop to N
  //
  // ...this is still under construction. other cases seem to work, though
}
// ToDo:
// -move to prototypes
// -rename to rsKroneckerVectorProduct
// -implement production version that takes a raw pointer to the x/y vector (and a workspace)
// -maybe take a larger workspace and avoid the internal copying that is needed for the "in-place"
//  operation
// -maybe the inverse trafo can be computed by using the individual inverse matrices in reverse 
//  order?
//
// See:
// http://www.mathcs.emory.edu/~nagy/courses/fall10/515/KroneckerIntro.pdf
// https://math.stackexchange.com/questions/1879933/vector-multiplication-with-multiple-kronecker-products
// https://gist.github.com/ahwillia/f65bc70cb30206d4eadec857b98c4065

template<class T>
void rsFastKroneckerTrafo2(std::vector<T>& v, const std::vector<const rsMatrix<T>*> M)
{
  int Nx = 1; for(size_t m = 0; m < M.size(); m++) Nx *= M[m]->getNumColumns();
  int Ny = 1; for(size_t m = 0; m < M.size(); m++) Ny *= M[m]->getNumRows();
  int N = rsMax(Nx, Ny);
  std::vector<T> t(N);     // temp storage
  v.resize(N);
  T* x = &v[0];
  T* y = &t[0];
  int hC = 1;
  int hR = 1;
  bool xContainsResult = true;
  for(size_t m = 0; m < M.size(); m++) 
  {
    int nR = M[m]->getNumRows();
    int nC = M[m]->getNumColumns();
    int NC = Nx / nC;  // or N / nR?  rename to NC ..maybe it should be Nx/nR
    int NR = Ny / nR;  // should be Ny/nR

    for(int i = 0; i < Nx; i += nC*hC) 
    {
      for(int j = 0; j < NC; j++) 
      {
        for(int k = 0; k < nR; k++)  
        {
          y[j+k*NC] = 0;
          for(int n = 0; n < nC; n++)
            y[j+k*NC] += M[m]->at(k, n) * x[nC*j+n];  // or nC*j+n?  
        }

        // y[j+k*NC]: wrong order, a 0 in between
        // y[j+k*NR]: wrong numbers
        // y[j+k*nR]: wrong numbers (only the 0th matches)
        // y[j+k*nC]: wrong order, a 0 in between
        // y[j+k]:    wrong numbers

        // from FDN:
        //y[j]    = a*x[2*j] + b*x[2*j+1];
        //y[j+N2] = c*x[2*j] + d*x[2*j+1];
      }
    }
    hC *= nC;
    hR *= nR; 
    RAPT::rsSwap(x, y);
    xContainsResult = !xContainsResult;
  }

  if( !xContainsResult )
    memcpy(y, x, N*sizeof(T)); 
  v.resize(Ny);
}




// Move to TransformsTests.cpp:

/** L is the order */
bool testKronecker2x2(int L)
{
  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;

  int N = (int)pow(2, L);                    // size of matrix is NxN
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
// maybe rename to testKronecker

bool testKronecker3x3(int L)
{
  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;

  int N = (int) pow(3, L);                   // size of matrix is NxN
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
  rsFastKroneckerTrafo(z2, H1, L);            // output of fast generalized transform
  return z == y && z2 == y;
}

bool testKroneckerProductTrafo()
{
  bool ok = true;

  using Vec = std::vector<double>;
  using Mat = rsMatrix<double>;

  Mat M23(2, 3, {1,2,3,4,5,6});
  Mat M24(2, 4, {1,2,3,4,5,6,7,8});
  Mat M32(3, 2, {1,7,3,4,5,6});
  Mat M33(3, 3, {1,2,3,4,5,6,7,8,9});

  std::vector<const Mat*> M(2);


  int Nx = 9;  // input dimensionality, product of all numbers of columns
  int Ny = 9;  // output dimensionality, product of all numbers of rows
  int N  = rsMax(Nx, Ny);  // space needed for in place trafo

  Vec x,y,z;

  x = rsRandomIntVector(N, -9, +9);
  Mat M_33_33 = Mat::getKroneckerProduct(M33, M33);
  y = M_33_33 * x;
  M[0] = &M33;
  M[1] = &M33;
  z = x; rsFastKroneckerTrafo( z, M); ok &= z == y;
  //z = x; rsFastKroneckerTrafo2(z, M); ok &= z == y;

  // x has dim 9, y has dim 4:
  Mat M_23_23 = Mat::getKroneckerProduct(M23, M23);
  y = M_23_23 * x; M[0] = &M23; M[1] = &M23;
  z = x; rsFastKroneckerTrafo( z, M); ok &= z == y; 
  //z = x; rsFastKroneckerTrafo2(z, M); ok &= z == y; // FAILS!!!


  // x has dim 16, y has dim 4:
  Mat M_24_24 = Mat::getKroneckerProduct(M24, M24);
  x = rsRandomIntVector(16, -9, +9);
  y = M_24_24 * x;
  M[0] = &M24;
  M[1] = &M24;
  z = x;
  rsFastKroneckerTrafo(z, M);
  ok &= z == y;

  // x has dim 2, y has dim 3:
  x = rsRandomIntVector(2, -9, +9);
  y = M32 * x; M.resize(1); M[0] = &M32;
  z = x; rsFastKroneckerTrafo(z, M);
  ok &= z == y;
  // ok, works

  // todo: try 2x2 ° 3x2 and 3x2 ° 2x2

  // x has dim 4, y has dim 9:
  Mat M_32_32 = Mat::getKroneckerProduct(M32, M32);
  x = rsRandomIntVector(4, -9, +9);
  y = M_32_32 * x; M.resize(2); M[0] = &M32; M[1] = &M32;
  //z = x; rsFastKroneckerTrafo(z, M); ok &= z == y;  // FAILS!!!
  //z = x; rsFastKroneckerTrafo2(z, M); ok &= z == y;  // FAILS!!!
  // FAILS! maybe we need to zero out the y extra values before running the algo - nope
  // they are zero - maybe the upper limit for i must be max(Nx,Ny) -> nope, access violation
  // ...maybe one of the inner loops must be run up to max(nC,nR) ..but no: the innermost loops
  // must run over exactly as many values as they are - anything else would make no sense
  // or: Nx must be replaced by N in the outer loop and nC by nM = min(nC,nR)
  // in this case, nR < nC, so nR is the minimum of both (in the previous tests it was always 
  // >= nC and these cases work). We also have Nx < Ny here for the first time in the tests
  // ..in each inner iteration, we must read 2 x-values and write 3 y-values - check that
  // maybe the loop limits are ok, but i must add offsets to the write locations, i.e. use
  // x[j+k*hM] =  instead of  x[j+k*hR] = where hM = max(hR,hC)

  // 3 3x3 matrices:
  Mat M_33_33_33 = Mat::getKroneckerProduct(M_33_33, M33);
  x = rsRandomIntVector(27, -9, +9);
  y = M_33_33_33 * x;
  M.resize(3); M[0] = &M33; M[1] = &M33; M[2] = &M33;
  z = x; rsFastKroneckerTrafo(z, M);
  ok &= z == y;

  // There are interesting identities for Kronecker products, for example:
  //   vec(B*Q*A^T) = (A°B)*vec(Q)
  // where A,B,Q are matrices, vec() turns a matrix into a vector by stacking its columns, * is
  // matrix multiplication and ° is the Kronecker product. See:
  //   https://ieeexplore.ieee.org/document/7864371
  // or:
  //  (A°B)^T  = A^T  ° B^T
  //  (A°B)^-1 = A^-1 ° B^-1
  //  (A°B)*(C°D) = A*C°B*D  ...i don't know, which product takes precedence here -> figure out!
  // ToDo: dig out more such identities and verify them numerically

  // see also: HADAMARD, KHATRI-RAO, KRONECKER AND OTHER MATRIX PRODUCTS
  // https://www.math.ualberta.ca/ijiss/SS-Volume-4-2008/No-1-08/SS-08-01-17.pdf

  // todo: Use a product of a 2x5 and 4x3 matrix

  return ok;
}

// move to rapt math tests:
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
  RAPT::rsFGHT(y4, 4, 2., 3., 5., 7.);
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
  ok &= testKronecker2x2(1);  //   2
  ok &= testKronecker2x2(2);  //   4
  ok &= testKronecker2x2(3);  //   8
  ok &= testKronecker2x2(4);  //  16
  ok &= testKronecker2x2(5);  //  32
  ok &= testKronecker2x2(6);  //  64
  ok &= testKronecker2x2(7);  // 128
  ok &= testKronecker2x2(8);  // 256
  ok &= testKronecker3x3(1);  //   3
  ok &= testKronecker3x3(2);  //   9
  ok &= testKronecker3x3(3);  //  27
  ok &= testKronecker3x3(4);  //  81
  ok &= testKronecker3x3(5);  // 243

  // Test the more general Kronecker product transform:
  ok &= testKroneckerProductTrafo();

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

