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
//  y = A�B�C�... * x, where � denotes the kronecker product? the function could take as input the 
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

// computes M[0] � M[1] � M[2] � ... � M[L-1] * x 
// ..or M[L-1] � M[L-2] � M[L-3] � ... � M[0] * x 
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
  using LT  = rsLinearTransforms;

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
  //rsFGHT(&z[0], N, a,b,c,d);                 // output of fast transform
  LT::kronecker2x2(&z[0], N, a,b,c,d);       // output of fast transform
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

  // todo: try 2x2 � 3x2 and 3x2 � 2x2

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
  //   vec(B*Q*A^T) = (A�B)*vec(Q)
  // where A,B,Q are matrices, vec() turns a matrix into a vector by stacking its columns, * is
  // matrix multiplication and � is the Kronecker product. See:
  //   https://ieeexplore.ieee.org/document/7864371
  // or:
  //  (A�B)^T  = A^T  � B^T
  //  (A�B)^-1 = A^-1 � B^-1
  //  (A�B)*(C�D) = A*C�B*D  ...i don't know, which product takes precedence here -> figure out!
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

  //typedef rosic::FeedbackDelayNetwork FDN;
  //typedef RAPT::rsArrayTools AT;
  using FDN = rosic::FeedbackDelayNetwork;
  using AT  = RAPT::rsArrayTools;
  using LT  = RAPT::rsLinearTransforms;

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
  //RAPT::rsFGHT(y4, 4, 2., 3., 5., 7.);
  LT::kronecker2x2(y4, 4, 2., 3., 5., 7.);
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
  //RAPT::rsFGHT(y8, 8, 2., 3., 5., 7.);
  LT::kronecker2x2(y8, 8, 2., 3., 5., 7.);
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
  ok &= fabs(AT::maxDeviation(x8, y8, 8)) < 1.e-15;

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
  // -Try running two networks in parallel with exact same settings except for the reverb time. 
  //  That allows us to mix two different exponential decays for early and late stage and mix them
  //  to taste, similar to the way it is done withe modal filters. If done via SIMD, it doesn't 
  //  even need to make the processing much more expensive. the ratio of both outputs can be used 
  //  to measure time which in turn can be used to further shape the envelope
  // -Try to use a different set of a,b,c,d parameters in the kronekcer trafe for each level, i.e.
  //  a different diffusion coeff for each level. Figure out, what difference it makes if high
  //  coeffs are at lower or higher levels. Maybe define 2 2x2 matrices M1 = a1,b1,c1,d1 and 
  //  M2 = a2,b2,c2,d2 and visualize their Kronecker products kron(M1,M2) and kron(M2,M1) as 
  //  heatmaps
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

void rotes::spectralFilter()
{
  // Demonstrates usage of rosic::SpectralFilter by creating a sawtooth wave as test input signal
  // and passing it through the filter. The results will be written to wavefiles.

  using Flt = rosic::SpectralFilter;
  using Vec = std::vector<double>;

  // Setup:

  // Technical parameters:
  double sampleRate    = 44100;
  int    maxBlockSize  = 4096;
  int    blockSize     = 1024;
  int    numSamples    = 2*sampleRate;

  // Input signal parameters:
  double sawFreq       = 100;

  // Filter parameters:
  double lowerCutoff   = 500;
  double upperCutoff   = 5000;
  Flt::modes mode      = Flt::modes::BANDPASS;



  // Create sawtooth test input signal:
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 1, sawFreq, sampleRate, 0.0, true);
  x = 0.5 * x; // Reduce volume to avoid clipping when the filter overshoots

  // Create and set up filter object:
  Flt filter(maxBlockSize);
  filter.setSampleRate(sampleRate);
  filter.setInputBlockSize(blockSize);
  filter.setMode(mode);
  filter.setLowerCutoff(lowerCutoff);
  filter.setUpperCutoff(upperCutoff);

  // Create output by applying the filter to the input:
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = filter.getSample(x[n]);

  // Write input and output into wave files:
  rosic::writeToMonoWaveFile("SpectralFilterInput.wav",  &x[0], N, sampleRate, 16);
  rosic::writeToMonoWaveFile("SpectralFilterOutput.wav", &y[0], N, sampleRate, 16);


  // ToDo:
  // -Use also white noise as test input
  // -Experiment with overlap and paddingFactor parameters
  // -Make a similar test for rosic::FormantShifter
}

void rotes::formantShifter()
{
  // Tests the class rosic::FormantShifter. We create as test input a sawtooth wave that has some 
  // formants imposed on it which are created via rosic::VowelFilter. Then we pass that 
  // sawtooth-with-formants into the formant shifter. Then we write input and output signals into 
  // wavefiles for inspection and audition.


  // Setup:

  // Output file parameters:
  double sampleRate    = 44100;         // Sample rate for the signals in Hz
  int    numSamples    = 2*sampleRate;  // We create a 2 seconds long signal.

  // Input signal parameters:
  double sawFreq       = 80;            // Fundamental frequency of the sawtooth
  double vowel         = 0.5;           // The vowels are on a scale form 0..1. The middle is "ah".
  double formantAmount = 1;             // 0: none, 1: normal, >1: overpronounced

  // Formant shifter parameters:
  double formantScale  = 1.5;           // Scaling factor for the formant frequencies
  int    maxBlockSize  = 4096;          // Maximum block size - must be passed to constructor
  int    blockSize     = 512;           // Block size, must be <= maxBlockSize


  // Create raw sawtooth signal:
  using Vec = std::vector<double>;
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 1, sawFreq, sampleRate, 0.0, true);
  x = 0.1 * x;  // To reduce the volume to avoid clipping when the filter overshoots.

  // Apply formants:
  rosic::VowelFilterStereo vowelFilter;
  vowelFilter.setSampleRate(sampleRate);
  vowelFilter.setVowel(vowel);
  vowelFilter.setAmount(formantAmount);
  Vec y(N);
  for(int n = 0; n < N; n++)
    vowelFilter.getSampleFrameStereo(&x[n], &x[n], &y[n], &y[n]);

  // Shift formants:
  rosic::FormantShifter formantShifter(maxBlockSize);
  formantShifter.setSampleRate(sampleRate);
  formantShifter.setFormantScale(formantScale);
  formantShifter.setInputBlockSize(blockSize);
  Vec z(N);
  for(int n = 0; n < N; n++)
    z[n] = formantShifter.getSample(y[n]);

  // Write input and output into wave files:
  rosic::writeToMonoWaveFile("FormantShifterInput.wav",  &y[0], N, sampleRate, 16);
  rosic::writeToMonoWaveFile("FormantShifterOutput.wav", &z[0], N, sampleRate, 16);

  // ToDo:
  // -formantAmount seems to have no effect. Must be a bug in VowelFilterStereo. Fix it!
  // -Add setFormantEmphasis to FormantShifter. Should make it possible to de/emphasize the 
  //  formants.
  // -Test setting up different block sizes, overlap factors, etc.
  // -Change the formants in the input dynamically over time.
  // -Maybe change the formant-shift also dynamically over time.
  // -Test, how different block sizes affect the sound when the formants are moving around.
  // -Experiment with different window functions. This stuff is implemented in the baseclass
  //  functions setInputBlockSize, setOverlapFactor, setPaddingFactor in 
  //  rosic::OverlapAddProcessor
  // -Implement a performance test and test how these parameters affect the CPU load.
}


void testSpectralShiftViaJH()
{
  // Tests our implementation of the spectral pitch shifting algorithm of  Nicolas Juillerat and 
  // Beat Hirsbrunner described in this paper:
  // https://www.researchgate.net/publication/261078164_Low_latency_audio_pitch_shifting_in_the_frequency_domain
  //
  // Notation used in the paper and in this implementation:
  //   a        : Source bin index
  //   b        : Destination bin index
  //   Om_x     : Complex STFT value "Omega_x" at bin with index x, x is placeholder for a or b
  //   O        : Overlap factor (typically 2,4,8)
  //   N        : FFT size (typically 512..2048)
  //   k        : Frequency scaling factor (typically 0.5..2.0)
  //   m        : Multiplier for synthesis FFT size (typically 2 or 4)
  //
  // Other notation from paper not used or needed here:
  //   f_s      : Sample rate (typically 44100)
  //   f_a      : Frequency of input sine
  //   B        : Bandwidth of an FFT bin in Hz (== f_s / N)
  //   phi      : Phase of input sine in first STFT frame
  //   p        : STFT frame index
  //   E:       : Error ratio between actually synthesized freq and desired freq of output sine
  //   s1,s2,s3 : 3 input sinusoids
  //   eps      : freq difference between the 3 test input sines
  //
  // My spectral processing framework actually has not (yet?) any concept of using a multiplier for 
  // the synthesis FFT size with respect to the analysis FFT size. But I think we can achieve a 
  // similar effect using the zero padding feature. This is actually supposed to give even higher 
  // quality results because it increases the FFT size at the analysis *and* synthesis side. Maybe 
  // later we can introduce this additional multiplier for an optimization of the algorithm.





  // Setup:

  // Output file parameters:
  double sampleRate  = 44100;         // Sample rate for the signals in Hz
  int    numSamples  = sampleRate/10; // We create a 1/10 seconds long signal.

  // Input signal parameters:
  double inputPeriod = 128;           // Length of one cycle in samples
  double inputPhase  = 90;            // Phase in degrees

  // Spectral shifter parameters:
  double freqScale   = 1.2;           // Scaling factor for the frequencies
  int    blockSize   = 1024;          // Block size. Must be power of 2
  int    overlap     = 2;             // Overlap factor. Must be power of 2
  int    zeroPad     = 1;             // Zero padding factor. Must be power of 2


  // Create sinusoidal test signal:
  using Vec = std::vector<double>;
  double inputFreq = sampleRate / inputPeriod;
  //int N = numSamples;
  Vec x(numSamples);
  createWaveform(&x[0], numSamples, 0, inputFreq, sampleRate, 
    RAPT::rsDegreeToRadiant(inputPhase), true);
  x = 0.5 * x;

  // Apply pitch shifting:
  using SS = rosic::SpectralShifter;
  rosic::SpectralShifter pitchShifter(blockSize, overlap, zeroPad);
  pitchShifter.setFrequencyScale(freqScale);
  pitchShifter.setInputBlockSize(blockSize);
  pitchShifter.setOverlapFactor(overlap);
  pitchShifter.setPaddingFactor(zeroPad);
  pitchShifter.setPhaseFormula(SS::PhaseFormula::useMultiplier);
  Vec y1(numSamples);
  for(int n = 0; n < numSamples; n++)
    y1[n] = pitchShifter.getSample(x[n]);

  // For reference, generate output without without the phase formula:
  pitchShifter.setPhaseFormula(SS::PhaseFormula::keepOriginal);
  pitchShifter.reset();
  Vec y2(numSamples);
  for(int n = 0; n < numSamples; n++)
    y2[n] = pitchShifter.getSample(x[n]);


  // Plot input and output signals:
  //rsPlotVectors(x, y1, y2, y2-y1);
  rsPlotVectors(x, y1, y2);

  // Observations:
  // -With freqScale = 1.2, The first few buffers look like a mess but then it settles and the only
  //  remaining artifact is a (strong) amplitude modulation. It looks less messy when we use a 
  //  start-phase of zero instead of 90, i.e. sine instead of cosine. In this case, there are only 
  //  corners in the signal whereas with cosine, there are jumps.
  // -With freqScale = 2, the discontinuities in the first frames disappear but the amp modulation 
  //  gets worse. There is also some amount of negative DC between 512 and 768 when the phase is 
  //  zero. With 90, this is not the case. The amp-modulation period is 512 samples.
  // -For a downward shift with freqScale = 0.5, there is no such amplitude modulation. Strangely, 
  //  there doesn't seem to be any difference between multiplying by w or not. Maybe that's a 
  //  special case?
  // -For input period = 128 and k = 1.25 and k = 0.8, we do not see amp-mod. For k = 1.2, there's
  //  string amp mod. For k = 1.6, there's little amp mod but the amp is too low overall. Looks 
  //  like it's half of what it should be.
  // -The difference between using and not using the phase-multiplier becomes apparent with
  //  inputPeriod = 32 - at least in the FFT spectrum. ToDo: generate output with and without the
  //  phase multiplier and plot both - and their difference.

  // Notes:
  // -We are not yet applying an output window. Try using one! But maybe this requires to use an
  //  overlap of 4 instead ot 2. I think a Hann window squared overlpas to one only with hopSize
  //  = blockSize / 4. The hann window itself overslap to one already with hopSize = blockSize / 2.
  //  Another option could be to use a sqrt(Hann) window. That should overlap perfectly when 
  //  squared at O = 2. Maybe try other windows. Maybe try also a demodulation approach.
}

void rotes::spectralShifter()
{
  // Under construction - does not yet work.
  //
  // We want to build a pitch shifter based on spectral processing. It should have a transient 
  // preservation feature.

  testSpectralShiftViaJH();


  // Setup:

  // Output file parameters:
  double sampleRate   = 44100;         // Sample rate for the signals in Hz
  int    numSamples   = sampleRate/10; // We create a 1/10 seconds long signal.

  // Input signal parameters:
  double inputPeriod  = 128;           // Length of one cycle in samples
  double inputPhase   = 90;            // Phase in degrees

  // Spectral shifter parameters:
  double freqScale    = 2.0;           // Scaling factor for the frequencies
  int    maxBlockSize = 4096;          // Maximum block size. Must be passed to constructor.
  int    maxOverlap   = 4;             // Maximum overlap factor. Can be passed to constructor.
  int    maxZeroPad   = 4;             // Max. zero padding factor. Can be passed to constructor.
  int    blockSize    = 1024;          // Block size, must be <= maxBlockSize
  int    overlap      = 2;             // Overlap factor. Must be power of 2 and <= maxOverlap.
  int    zeroPad      = 2;             // Zero padding factor. Must be power of 2 and <= maxZeroPad


  // Create sinusoidal test signal:
  using Vec = std::vector<double>;
  double inputFreq = sampleRate / inputPeriod;
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 0, inputFreq, sampleRate, RAPT::rsDegreeToRadiant(inputPhase), true);
  x = 0.5 * x;  


  rosic::SpectralShifter pitchShifter(maxBlockSize, maxOverlap, maxZeroPad);
  //pitchShifter.setSampleRate(sampleRate);
  pitchShifter.setInputBlockSize(blockSize);
  pitchShifter.setFrequencyScale(freqScale);
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = pitchShifter.getSample(x[n]);


  // Plot input and output signals:
  rsPlotVectors(x, y);

  // Write input and output into wave files:
  rosic::writeToMonoWaveFile("SpectralShifterInput.wav",  &x[0], N, sampleRate, 16);
  rosic::writeToMonoWaveFile("SpectralShifterOutput.wav", &y[0], N, sampleRate, 16);
  int dummy = 0;


  // Observations:
  // -This does not yet work. But it's not really supposed to. I'm currently interpolating the 
  //  complex spectrum. What we need to do is convert to magnitude/phase, interpolate the 
  //  magnitudes and adjust the phases according to a prediction from the previous frame.
  // -With blockSize = 1024, freqScale = 2, inputPeriod = 100, the output shows strong amplitude
  //  modulation. This remains to be the case with inputPeriod = 128 so the misalignment of the 
  //  cycles with the window is not the reason for this amp-mod. The modulation period is given by
  //  512 samples. This aligns with the hopSizes as well as with blockSize/2. ToDo: figure out,
  //  if it remains the same, if we decrease the hopSize to 256 or if it remains as is. From that
  //  we can conclude if it is indeed related to the hopSize or always blockSize/2.
  // -With these settings, the output signal starts at sample 768. How does that number come about?
  //  Try to use a different start phase for the input sine and see, if this changes this. OK - 
  //  yes - it shows that the first nonzero output block starts at 512 but there's only some small
  //  scale rippling going on for the first half of the first nonzero output buffer and the real
  //  signal starts at 768.
  //
  // ToDo:
  // -Add the phase multiplication step. 
  // -Try it on a sinusoid that aligns with an FFT bin. Then try it on a  sinusoid that is in 
  //  between two bins.
  //
  // Ideas:


  // Resources:
  // https://www.reddit.com/r/DSP/comments/k6t24c/pitch_shifting_algorithm_in_frequency_domain/
  //
  // https://www.guitarpitchshifter.com/algorithm.html
  // https://www.guitarpitchshifter.com/software.html
  // -> uses blockSize 1024, hopSize 256
  // https://www.guitarpitchshifter.com/matlab.html
  //
  // https://lcav.gitbook.io/dsp-labs/dft/implementation  has python implementation
  //
  // https://pitchtech.ch/
  // https://www.researchgate.net/publication/261078164_Low_latency_audio_pitch_shifting_in_the_frequency_domain

  // https://quod.lib.umich.edu/cgi/p/pod/dod-idx/real-time-low-latency-audio-processing-in-java.pdf
  // ...there is supposed to be a java implementation of that approach but googling Decklight 4
  // as they have called the framework in the paper reveals nothing. But the paper is from 2007, so
  // that framework may not exist anymore or maybe go by another name now
  // 
  // https://github.com/JorenSix/TarsosDSP/tree/master/core/src/main/java/be/tarsos/dsp
  // ...doesn't seem to have pitch-shift but some other interesting stuff
  //
  // https://wiki.linuxaudio.org/apps/categories/pitch_effect
  //
  // https://www.ee.columbia.edu/~dpwe/papers/LaroD99-pvoc.pdf
  //
  // https://blogs.zynaptiq.com/bernsee/pitch-shifting-using-the-ft/
  //
  // https://github.com/BelaPlatform/bela-online-course/tree/master/lectures/lecture-18
  // https://www.youtube.com/watch?v=xGmRaTaBNZA
  //
  // https://www.youtube.com/watch?v=fJUmmcGKZMI  Four Ways To Write A Pitch-Shifter - Geraint Luff - ADC22
  //
  // https://www.youtube.com/watch?v=PjKlMXhxtTM  Making a Pitch Shifter
  // -says at 1:30 that a time-stretch - then resample algo for pitch shift is easier than direct 
  //  pitch shift
  // -at 14:55 says that at transients, phases in the resynthesized signal should not be modified with 
  //  respect to the input STFT at that bin
  // -at 15:30: mentions librosa -  a python library that implements PV based pitch-shift
  //
  // Phase Vocoder (Flanagan/Golden):
  // https://archive.org/details/bstj45-9-1493/mode/2up
  // The classic paper on the topic, I think.

}