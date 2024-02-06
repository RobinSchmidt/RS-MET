using namespace rotes;
using namespace rosic;
using namespace RAPT;

// rename to allpassDiffusor
void rotes::allpassDisperser()
{
  // We plot some impulse responses of allpass filters. The goal is to build some intuition for
  // designing a disperser/diffusor stage for an algorithmic reverb.

  // User parameters:
  double sampleRate = 44100;
  int    numStages  = 4;
  bool   biquads    = false;     // Switches between 1st and 2nd order stages
  double freq       = 3900;
  double quality    = 2.0;       // Quality factor "Q" for 2nd order stages
  int    N          = 4096;      // Number of samples for the plot

  // Create and set up the allpass chain:
  rosic::AllpassChain apf1;
  apf1.setSampleRate(sampleRate);
  apf1.setFrequency(freq);
  apf1.setNumStages(numStages);
  apf1.setQ(quality);
  if(biquads)
    apf1.setMode(apf1.SECOND_ORDER_ALLPASS);
  else
    apf1.setMode(apf1.FIRST_ORDER_ALLPASS);

  // Record the impulse response of the allpass:
  using Vec = std::vector<double>;
  Vec y(N);
  y[0] = apf1.getSample(1);
  for(int n = 1; n < N; n++)
    y[n] = apf1.getSample(0);

  // Plot the impulse response:
  //rsPlotVector(y);



  // Now apply the allpass delay chain diffusor:
  std::vector<int> delays = { 13, 17, 23, 29 };
  //std::vector<int> delays = { 7, 13, 17, 23, 29 };
  //std::vector<int> delays = { 7, 13, 17, 23 };
  //std::vector<int> delays = { 5, 7, 13, 17, 23 };
  //std::vector<int> delays = { 17, 23, 29, 37 };
  //std::vector<int> delays = { 3, 5, 7, 13, 17, 23 };
  double coeff  = +0.75;
  double sclAmt =  0.24;

  // Compute the coeffs for the stages:
  numStages = (int) delays.size();   // we repurpose this variable here (that's kinda a ugly!)
  Vec coeffs(numStages);
  for(int i = 0; i < numStages; i++)
  {
    double scale = double(delays[0]) / double(delays[i]);
    scale = pow(scale, sclAmt);
    coeffs[i] = scale * coeff;
  }

  // Create and set up the allpass chain:
  rsAllpassDelayChain<double, double> apdc;
  apdc.setMaxNumStages(numStages);
  apdc.setNumStages(numStages);
  for(int i = 0; i < numStages; i++)
  {
    apdc.setMaxDelayInSamples(i, delays[i]);
    apdc.setDelayInSamples(   i, delays[i]);
    apdc.setAllpassCoeff(     i, coeffs[i]);
  }

  // Apply it to y:
  Vec z(N);
  for(int n = 0; n < N; n++)
    z[n] = apdc.getSample(y[n]);

  // Write the impulse response to a wave file for listening:
  RAPT::rsArrayTools::normalize(&z[0], N);
  rosic::writeToMonoWaveFile("Diffusor.wav", &z[0], N, 44100, 16);

  // Plot the signal:
  rsPlotVector(z);

  // Observations for y:
  // -For a single 1st order allpass stage, we get an initial negative spike followed by a 
  //  positive value which then turns into an exponential decay. Adding stages seems to directly
  //  translate to adding additional local maxima to teb impulse response.
  // -For 2nd order allpass chains with 1 stage, we observe an impuse response that looks like
  //  and initial positive/negative spike pair that then turns into a a damped sine. More stages
  //  will make the initial transient look more chaotic.
  // -With sampleRate = 44100, numStages = 1, freq = 3900, quality = 2.0, we see a nice first 
  //  sinusoidal peak at n = 5.
  //
  // Observations for z:
  // -With the allpass coeff of around 0.9, we can here ha clear toanlity in the output. With 
  //  coeff = 0.7 that tonality is less pronounced - but the impulse response is shorter. I think,
  //  the sweet spot is perhaps around 0.75.
  // -delays = 3,5,7,13,17,23, coeff = +0.75, sclAmt = 0.5 gives a nice zap. ..oH - maybe my audio 
  //  interface was distorting. However - using sclAmt = 0.2..0.25 seems a good compromise.
  // -Using these vaules with only 4 delays of 13,17,23,29 seems to work also well. Experimenting 
  //  with more stages for the initial AllpassChain also helps to shape the sound.

  // ToDo:
  // -Alternatively to the single 2nd order stage, use a 4-stage 1st order chain - yeah - that's 
  //  what we currently do.
  // -The goal is to get a dense noisy impulse response with (short) exponential envelope of a few
  //  hundreds of samples length.
  // -Maybe notch out the Nyquist freq to get rid of the oscillation there
}

void rotes::allpassDelay()
{
  // User parameters:
  double coeff = +0.9;
  int    delay = 10;
  int    N     = 50*delay;     // Number of samples to plot

  // Create and set up the allpass delay:
  //rsAllpassDelayNaive<double, double> apd;
  rsAllpassDelay<double, double> apd;
  apd.setMaxDelayInSamples(delay);
  apd.setDelayInSamples(delay);
  apd.setAllpassCoeff(coeff);

  // Record the impulse response of the allpass delay:
  using Vec = std::vector<double>;
  Vec y(N);
  y[0] = apd.getSample(1);
  for(int n = 1; n < N; n++)
    y[n] = apd.getSample(0);

  // Plot the impulse response:
  rsPlotVector(y);

  // Observations:
  // -A positive coeff leads to alternating spikes. A negative coeff leads to positive spikes only.
  //
  // Conclusions:
  // -I think, for a diffuser in a reverb, a positive coeff should be used. The alternation seems 
  //  desirable. ...but that's just a hunch and needs to be verified by listening tests.

}

void rotes::allpassDelayChain()
{
  // User parameters:
  //std::vector<int> delays = { 13, 17, 23, 29, 37 };
  std::vector<int> delays = { 13, 17, 23, 29 };
  double coeff  = +0.9;
  int    N      = 4000;
  double sclAmt = 0.2;   // Amount of scaling of the coeff as function of delay

  // Compute the coeffs for the stages:
  int numStages = (int) delays.size();
  using Vec = std::vector<double>;
  Vec coeffs(numStages);
  for(int i = 0; i < numStages; i++)
  {
    double scale = double(delays[0]) / double(delays[i]);
    scale = pow(scale, sclAmt);
    coeffs[i] = scale * coeff;
  }


  // Create and set up the allpass delays:
  std::vector<rsAllpassDelayNaive<double, double>> apds(numStages);
  for(int i = 0; i < numStages; i++)
  {
    apds[i].setMaxDelayInSamples(delays[i]);
    apds[i].setDelayInSamples(delays[i]);
    apds[i].setAllpassCoeff(coeffs[i]);
  }

  // Record the impulse response of the allpass delay chain:
  Vec y(N);
  double tmp = apds[0].getSample(1);
  for(int i = 1; i < numStages; i++)
    tmp = apds[i].getSample(tmp);
  y[0] = tmp;

  for(int n = 1; n < N; n++)
  {
    tmp = apds[0].getSample(0);
    for(int i = 1; i < numStages; i++)
      tmp = apds[i].getSample(tmp);
    y[n] = tmp;
  }


  // Now do the same thing with the class rsAllpassDelayChain
  rsAllpassDelayChain<double, double> apdc;
  apdc.setMaxNumStages(numStages);
  apdc.setNumStages(numStages);
  for(int i = 0; i < numStages; i++)
  {
    apdc.setMaxDelayInSamples(i, delays[i]);
    apdc.setDelayInSamples(   i, delays[i]);
    apdc.setAllpassCoeff(     i, coeffs[i]);
  }
  Vec z(N);
  z[0] = apdc.getSample(1);
  for(int n = 1; n < N; n++)
    z[n] = apdc.getSample(0);
  //bool ok = y == z;
  bool ok = rsIsCloseTo(y, z, 1.e-16);
  rsAssert(ok);


  // Write the impulse response to a wave file for listening:
  RAPT::rsArrayTools::normalize(&y[0], N);
  rosic::writeToMonoWaveFile("AllpassDelayChain.wav", &y[0], N, 44100, 16);

  // Plot the impulse response:
  rsPlotVector(y);

  // Observations:
  // -Towards the end, it looks like there is an oscillation at the Nyquist freq when we don't 
  //  scale the coeff at all, i.e. choose sclAmt = 0, For sclAmt = 1, we see a repetitive pattern
  //  at the end.
  // -Higher values of sclAmt tend to concentrate the energy more at the start.
  // -It sounds like a kind of breath noise (like in a flute) but with percussive envelope. With 
  //  higher values of sclAmt, it begins to sound more like a mallet-on-glass sound. Negative 
  //  values (like -0.1) elongate it - it begins to sound like bowed glass or something. I think,
  //  sclAmt in the range 0.0...0.2 should be OK for a diffusor
  //
  // ToDo:
  // -Maybe try to use a different coeff for each stage. Maybe stages with longer delay should have 
  //  smaller coeffs such that the overall decay time is the same for each stage. I think, the 
  //  coeff for stage i should be multiplied by double(delays[0]) / double(delays[i]).
  //  ...hmm...doing so leads to a clear repetitive pattern at the end
  // -Write a class rsAllpasDelayChain that encapsulates the chaining of several allpass delays 
  //  (done). Then transform this test here into a unit test for that class
  // -Combine the allpass delay chain with the disperser from allpassDisperser. That should 
  //  give a nice "random noise" kind of signal.
}

bool rotes::testAllpassDelayNested()
{
  // Turn this into a unit test!

  using Real = double;

  int N = 4096;


  std::vector<int>  delays;
  std::vector<Real> coeffs;
  using Vec = std::vector<Real>;
  bool ok = true;



  // Create and set up the 0-level nested allpass delay structure:
  delays = {   7 };
  coeffs = { 0.8 };
  rsAllpassDelay<Real, Real> apdn0;
  apdn0.setMaxDelayInSamples(delays[0]);
  apdn0.setDelayInSamples(   delays[0]);
  apdn0.setAllpassCoeff(     coeffs[0]);

  Vec x(N), y0(N);
  x[0] = 1;  
  for(int n = 0; n < N; n++)
    y0[n] = apdn0.getSample(x[n]);

  // Create and set up the multi-level nested allpass delay structure and set it up to zero levels of
  // nesting such that its output should match y0.
  rsAllpassDelayNested<Real, Real> apdn;
  apdn.setMaxNumStages(4);
  apdn.setNumStages(1);                      // 1 stage means 0 levels of nesting
  apdn.setMaxDelayInSamples(0, delays[0]);
  apdn.setDelayInSamples(   0, delays[0]);
  apdn.setAllpassCoeff(     0, coeffs[0]);

  //  Create impulse response of multi-level nested allpass:
  Vec z0(N);
  for(int n = 0; n < N; n++)
    z0[n] = apdn.getSample(x[n]);
  ok &= z0 == y0;
  RAPT::rsAssert(ok);
  //rsPlotVectors(y0, z0);
  // OK - nice - we have a match!

  // Create and set up the 1-level nested allpass delay structure:
  delays = {   7,  11 };
  coeffs = { 0.8, 0.7 };
  rsAllpassDelayNestedL1<Real, Real> apdn1;
  apdn1.setMaxDelayInSamples(0, delays[0]);
  apdn1.setMaxDelayInSamples(1, delays[1]);
  apdn1.setDelayInSamples(   0, delays[0]);
  apdn1.setDelayInSamples(   1, delays[1]);
  apdn1.setAllpassCoeff(     0, coeffs[0]);
  apdn1.setAllpassCoeff(     1, coeffs[1]);

  // Create impulse response:
  Vec y1(N);
  for(int n = 0; n < N; n++)
    y1[n] = apdn1.getSample(x[n]);

  // Set up the multi-level nested allpass delay structure and set it up to one level of
  // nesting such that its output should match y1.
  apdn.reset();
  apdn.setNumStages(2);                      // 2 stages means 1 level of nesting
  apdn.setMaxDelayInSamples(0, delays[0]);
  apdn.setMaxDelayInSamples(1, delays[1]);
  apdn.setDelayInSamples(   0, delays[0]);
  apdn.setDelayInSamples(   1, delays[1]);
  apdn.setAllpassCoeff(     0, coeffs[0]);
  apdn.setAllpassCoeff(     1, coeffs[1]);

  // Create impulse response of multi-level nested allpass with the general and the unrolled
  // special case getSample-function:
  Vec z1(N);
  for(int n = 0; n < N; n++)
    z1[n] = apdn.getSample2Stages(x[n]);   // unrolled special case for numStages = 2
  ok &= z1 == y1;

  apdn.reset();
  for(int n = 0; n < N; n++)
    z1[n] = apdn.getSample(x[n]);    // general case
  ok &= z1 == y1;

  //RAPT::rsAssert(ok);

  //rosic::writeToMonoWaveFile("AllpassDelaysNested1.wav", &y1[0], N, 44100, 16);
  //rsPlotVectors(y1, z1); 
  // When looking at a spectrum in Audacity, we need to ensure to use the right settings for the 
  // spectrum analysis: window should be rectangular and FFT size should be N

  //bool ok = rsIsWhite(y1, tol);  // Check that the filter is allpass
  // maybe make sure that it has decayed sufficiently such that truncating the impulse response 
  // does not lead to a large deviation from whiteness


  // Now test a 2-level nesting:
  delays = {   7,  11,  17 };
  coeffs = { 0.7, 0.6, 0.5 };
  rsAllpassDelayNestedL2<Real, Real> apdn2;
  apdn2.setMaxDelayInSamples(0, delays[0]);
  apdn2.setMaxDelayInSamples(1, delays[1]);
  apdn2.setMaxDelayInSamples(2, delays[2]);
  apdn2.setDelayInSamples(   0, delays[0]);
  apdn2.setDelayInSamples(   1, delays[1]);
  apdn2.setDelayInSamples(   2, delays[2]);
  apdn2.setAllpassCoeff(     0, coeffs[0]);
  apdn2.setAllpassCoeff(     1, coeffs[1]);
  apdn2.setAllpassCoeff(     2, coeffs[2]);

  Vec y2(N);
  for(int n = 0; n < N; n++)
    y2[n] = apdn2.getSample(x[n]);

  // Set up the multi-level nested allpass delay structure and set it up to two levels of
  // nesting such that its output should match y2.
  apdn.reset();
  apdn.setNumStages(3);                      // 3 stages means 2 levels of nesting
  apdn.setMaxDelayInSamples(0, delays[0]);
  apdn.setMaxDelayInSamples(1, delays[1]);
  apdn.setMaxDelayInSamples(2, delays[2]);
  apdn.setDelayInSamples(   0, delays[0]);
  apdn.setDelayInSamples(   1, delays[1]);
  apdn.setDelayInSamples(   2, delays[2]);
  apdn.setAllpassCoeff(     0, coeffs[0]);
  apdn.setAllpassCoeff(     1, coeffs[1]);
  apdn.setAllpassCoeff(     2, coeffs[2]);

  Vec z2(N);
  for(int n = 0; n < N; n++)
    z2[n] = apdn.getSample3Stages(x[n]);
  ok &= z2 == y2;
  //RAPT::rsAssert(ok);
  //rsPlotVectors(y2, z2);

  apdn.reset();
  for(int n = 0; n < N; n++)
    z2[n] = apdn.getSample(x[n]);
  ok &= z2 == y2;

  //RAPT::rsAssert(ok);
  //rsPlotVectors(y2, z2);
  //rosic::writeToMonoWaveFile("AllpassDelaysNested2.wav", &y2[0], N, 44100, 16);
  //rsPlotVector(y2);


  // Now test a 3-level nesting:
  delays = {   7,  11,  17,  23 };
  coeffs = { 0.6, 0.5, 0.4, 0.3 };
  rsAllpassDelayNestedL3<Real, Real> apdn3;
  apdn3.setMaxDelayInSamples(0, delays[0]);
  apdn3.setMaxDelayInSamples(1, delays[1]);
  apdn3.setMaxDelayInSamples(2, delays[2]);
  apdn3.setMaxDelayInSamples(3, delays[3]);
  apdn3.setDelayInSamples(   0, delays[0]);
  apdn3.setDelayInSamples(   1, delays[1]);
  apdn3.setDelayInSamples(   2, delays[2]);
  apdn3.setDelayInSamples(   3, delays[3]);
  apdn3.setAllpassCoeff(     0, coeffs[0]);
  apdn3.setAllpassCoeff(     1, coeffs[1]);
  apdn3.setAllpassCoeff(     2, coeffs[2]);
  apdn3.setAllpassCoeff(     3, coeffs[3]);
  Vec y3(N);
  for(int n = 0; n < N; n++)
    y3[n] = apdn3.getSample(x[n]);

  apdn.reset();
  apdn.setNumStages(4);                      // 4 stages means 3 levels of nesting
  apdn.setMaxDelayInSamples(0, delays[0]);
  apdn.setMaxDelayInSamples(1, delays[1]);
  apdn.setMaxDelayInSamples(2, delays[2]);
  apdn.setMaxDelayInSamples(3, delays[3]);
  apdn.setDelayInSamples(   0, delays[0]);
  apdn.setDelayInSamples(   1, delays[1]);
  apdn.setDelayInSamples(   2, delays[2]);
  apdn.setDelayInSamples(   3, delays[3]);
  apdn.setAllpassCoeff(     0, coeffs[0]);
  apdn.setAllpassCoeff(     1, coeffs[1]);
  apdn.setAllpassCoeff(     2, coeffs[2]);
  apdn.setAllpassCoeff(     3, coeffs[3]);

  Vec z3(N);
  for(int n = 0; n < N; n++)
    z3[n] = apdn.getSample(x[n]);

  ok &= z3 == y3;
  //rsPlotVectors(y3, z3);






  //rosic::writeToMonoWaveFile("AllpassDelaysNested3.wav", &y3[0], N, 44100, 16);
  //rsPlotVector(y3);


  return ok;
  //RAPT::rsAssert(ok);
  //int dummy = 0;

  // ToDo:
  // -Make a lattice-based implementation for an arbitrary amount of nesting and compare the 
  //  results to the explicit implementations of the various nesting levels up to 3 (done)
  // -Implement an rsIsWhite method and use it on the impulse responses to test, if they are white.
  //  It should take a tolerance in dB.
}


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

bool rotes::feedbackDelayNetwork()
{
  bool result = true;  // get rid! this is not a unit test!

  //FeedbackDelayNetwork16 *fdn16 = new FeedbackDelayNetwork16;

  double amplitude = 0.5;       // amplitude of the input impulse
  double diffusion = 100.0;     // diffusion parameter in percent


  using AT = RAPT::rsArrayTools;

  FeedbackDelayNetwork fdn;
  fdn.setDiffusion(diffusion);
  //fdn.setReverbTime(5.0); // Should be scaled in seconds

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

bool rotes::testAllpassDelay()
{
  bool  ok    = true;
  using Real  = double;
  using Vec   = std::vector<Real>;

  Real  coeff = 0.9;  // The allpass coefficient
  int   N     = 100;  // Number of samples to render

  // Create first order allpass filter:
  RAPT::rsOnePoleFilter<Real, Real> apf1;
  apf1.setCoefficients(coeff, 1.0, -coeff);

  // Create and set up M-th order allpass filter with M=1:
  rsAllpassDelayNaive<Real, Real> apd1;
  apd1.setMaxDelayInSamples(1);
  apd1.setDelayInSamples(1);
  apd1.setAllpassCoeff(coeff);

  // Create input and target output:
  Vec x(N), t(N);
  x[0] = 1;
  for(int n = 0; n < N; n++)
    t[n] = apf1.getSample(x[n]);

  // Create output of our M-th order allpass:
  Vec y(N);
  for(int n = 0; n < N; n++)
    y[n] = apd1.getSample(x[n]);

  // Compare both outputs - they should be equal:
  ok &= y == t;
  //rsPlotVectors(x, y, t);

  // Now use a higher delay:
  int M = 5;
  apd1.setMaxDelayInSamples(M);
  apd1.setDelayInSamples(M);
  for(int n = 0; n < N; n++)
    y[n] = apd1.getSample(x[n]);

  // The expected result is that y equals t with M-1 zeros inserted between adjacent samples of t, 
  // i.e. y[M*n] = t[n]  or  y[n] = t[n/M]  when n is divisible by M and zero otherwise.
  for(int n = 0; n < n/M; n++) {
    if( n % M == 0 )
      ok &= y[n] == t[n/M];
    else
      ok &= y[n] == 0; }
  //rsPlotVectors(x, y, t);

  // At this point, the naive implementaion of the allpass delay seems to work. Now test if the 
  // optimized variant produces the same result:
  rsAllpassDelay<Real, Real> apd2;
  apd2.setMaxDelayInSamples(M);
  apd2.setDelayInSamples(M);
  apd2.setAllpassCoeff(coeff);
  Vec y2(N);
  for(int n = 0; n < N; n++)
    y2[n] = apd2.getSample(x[n]);
  Real tol = 1.e-16;
  ok &= rsIsCloseTo(y2, y, tol);
  //rsPlotVectors(y2 - y);

  return ok;

  // ToDo:
  // -Write an optimized version of the class rsAllpassDelay (using DF2 or TDF2) and compare it to
  //  the current DF1 implementation.
  // -Maybe write a function rsIsWhite and feed it the impulse response of our supposed allpass
  //  filter. It should take an FFT and check if the spectrum is white up to some tolerance which
  //  we need because we work with a truncated impulse response. That function can then be used
  //  in statements like ok &= rsIsWhite(y, tol); etc.
}



// Move to rs_testing:
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

bool rotes::testFreqShifter()
{
  bool ok = true;


  // Try to create a rosic::FrequencyShifter object. We had an access violation in its
  // constructor, specifically in the line:
  //   halfbandFilter2.setApproximationMethod(rsPrototypeDesignerD::ELLIPTIC)
  rosic::FrequencyShifter freqShifter;
  // It happened in  rsEngineersFilter<TSig, TPar>::updateCoefficients(bool resetState) in the line
  // rsBiquadCascade<TSig, TPar>::initBiquadCoeffs(); and I think TSig=TPar=double. Apparently, 
  // when trying to write a1[i] = 0.0;  where i=0  is the violation. 
  //
  // OK - but this bug was fixed. Now we can repurpose this test here to actually produce some 
  // output of the freq-shifter and check if it looks as we expect it. Maybe produce a combination 
  // of sine-waves, shift it and programmatically check the magnitude spectrum of the output signal.


  return ok;
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





std::vector<double> getSpectralShifterInput(int numSamples, 
  int inputWaveform, double inputPeriod, double inputPhase)
{
  int sampleRate = 44100;            // Doesn't really matter
  using Vec = std::vector<double>;
  double inputFreq = sampleRate / inputPeriod;
  Vec x(numSamples);
  createWaveform(&x[0], numSamples, inputWaveform, inputFreq, (double)sampleRate, 
    RAPT::rsDegreeToRadiant(inputPhase), true);
  x = 0.5 * x;
  return x;
}
std::vector<double> getSpectralShifterOutput(const std::vector<double> x, double freqScale,
  rosic::SpectralShifter::Algorithm algo, int blockSize, int overlap, int zeroPad,
  bool useAnalysisWindow, bool useSynthesisWindow, int windowPower,
  rosic::SpectralShifter::PhaseFormula phaseFormula)
{
  int numSamples = (int) x.size();
  //int sampleRate = 44100;          // Needed for output file

  // Set up pitch shifter:
  using SS = rosic::SpectralShifter;
  rosic::SpectralShifter ps(blockSize, overlap, zeroPad);
  ps.setAlgorithm(algo);
  ps.setFrequencyScale(freqScale);
  ps.setInputBlockSize(blockSize);
  ps.setOverlapFactor(overlap);
  ps.setPaddingFactor(zeroPad);
  ps.setUseInputWindow(useAnalysisWindow);
  ps.setUseOutputWindow(useSynthesisWindow);
  ps.setWindowPower(windowPower);
  ps.setPhaseFormula(phaseFormula);

  // Set up plotting:
  bool plot = true;  // maybe make it a function parameter
  if(plot)
  {
    ps.blocksToPlot ={ 1,2,3,4,5,6,7,8,9,10 };            // The blocks for which plots are to be produced
    //ps.plotRawInputBlock       = true;
    //ps.plotWindowedInputBlock  = true;
    //ps.plotPaddedInputBlock    = true;
    //ps.plotInputSpectrum       = true;
    //ps.plotOutputSpectrum      = true;
    //ps.plotRawOutputBlock      = true;
    //ps.plotWindowedOutputBlock = true;
    // ToDo: Maybe let the plotter plot to files
  }

  // Apply the pitch shifting:
  std::vector<double> y(numSamples);
  for(int n = 0; n < numSamples; n++)
    y[n] = ps.getSample(x[n]);
  return y;
}
// The function parameters form 3 groups: (1) creative parameters (currently only 1: freqScale),
// (2) technical parameters (algo, blocksize, ...), (3) input signal parameters.
// inputWaveform: 0: sine, 1: saw, 2: square, 3: triangle
void testSpectralShifter(double freqScale, 
  rosic::SpectralShifter::Algorithm algo, int blockSize, int overlap, int zeroPad, 
  bool useAnalysisWindow, bool useSynthesisWindow, int windowPower, 
  rosic::SpectralShifter::PhaseFormula phaseFormula,
  int inputWaveform, double inputPeriod, double inputPhase)
{
  int numSamples = 6 * blockSize;  // We should produce enough blocks to pass the transient phase
  using Vec = std::vector<double>;
  Vec x = getSpectralShifterInput(numSamples, inputWaveform, inputPeriod, inputPhase);
  Vec y = getSpectralShifterOutput(x, freqScale, algo, blockSize, overlap, zeroPad,
    useAnalysisWindow, useSynthesisWindow, windowPower, phaseFormula);
  rsPlotVectors(x, y);

  // Maybe write wavefile:
  // -generate a filename from the parameters like 
  //  SpectralShift_FrqScl=0.8_Alg=LD_BlkSz=1024_OvLp=2_ZrPd=2_Phs=Mul__Wv=Sin_Cyc=128__Output.wav
  //  and write an output file with that name to disk, maybe optionally also write in corresponding
  //  inputs signal...or maybe write input into left channel and output into right. Thenwe can drop
  //  the _Output part.
}

void testSpectralShift()
{
  using SS = rosic::SpectralShifter;
  SS::Algorithm    JH   = SS::Algorithm::JuilHirs;
  SS::Algorithm    RS1  = SS::Algorithm::RobSchm1;
  SS::Algorithm    RS2  = SS::Algorithm::RobSchm2;
  SS::PhaseFormula Mul  = SS::PhaseFormula::useMultiplier;
  SS::PhaseFormula Keep = SS::PhaseFormula::keepOriginal;
  std::vector<double> x, y1, y2;

  /*
  //-----------------------------------------------------------------------------------------------
  // Experiments with the Juillerat/Hirsbrunner (JH) algorithm:

  //-------------------
  //| JH Downshifting |
  //-------------------

  // Copied from below - for debug plotting purposes:
  //testSpectralShifter(0.80, JH, 1024, 2, 1, true, true,   1, Mul,  0, 128, 90.0);

  testSpectralShifter(0.30, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Produces a good 2/5=0.4 shift from sample 1536 onwards. 
  // -Amplitude looks good.

  testSpectralShifter(0.40, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Produces a short transient and the almost silence. Its actually a very low amplitude signal 
  //  at amplitude of around 0.0027 that looks like a sine segment repeated over and over, i.e. 
  //  with phase resets.

  testSpectralShifter(0.50, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Looks good from sample 1535 onwards, i.e. after the transient/warm-up phase.


  testSpectralShifter(0.55, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -3 input peaks align with 2 output peaks at around 1536 and 1795. That's a ratio of 2/3, not
  //  the requested 3/5
  // -The amplitude looks good.

  // OK - now we try the same scale factor 0.55 with some amount of zero padding:
  testSpectralShifter(0.55, JH, 1024, 2, 2, true, false,  2, Mul,  0, 128, 90.0);
  // -Produces a freq shift of 1/2
  // -Amplitude is too low
  // -Phase is not aligned - shifted by half a cycle.

  // ...let's try even more zero padding:
  testSpectralShifter(0.55, JH, 1024, 2,  4, true, false,  2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.55, JH, 1024, 2,  8, true, false,  2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.55, JH, 1024, 2, 16, true, false,  2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.55, JH, 1024, 2, 32, true, false,  2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.55, JH, 1024, 2, 64, true, false,  2, Mul,  0, 128, 90.0);
  // -Produces a freq shift of 1/2
  // -Amplitude is too low
  // -Phase is not aligned - shifted by a quarter of a cycle.
  // -Modifying the phase multiplier formula to remove the zero-padding factor "m" from it does
  //  not make a big difference - the difference is just in some tiny details
  // -Increasing the zero padding up to 16 gives quality improvements (flatter amplitude, less 
  //  transient artifacts). Goig even higher than 16 has diminishing returns.


  testSpectralShifter(0.60, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -3 input peaks align with 2 output peaks at around 1666 and 2180. That's a ratio of 2/3, not
  //  the requested 3/5
  // -Output amplitude is half of what is should be.

  testSpectralShifter(0.65, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -After 5 peaks of input and 4 peaks of output, we get phase aslignment again. This happens
  //  for example at samples around 1800 and 2300. But that's a freq-ration of 4/5 = 0.8 not the
  //  desired 0.65
  // -Also, the output amplitude is too low. Roughly half of what it should be

  testSpectralShifter(0.75, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Produces a rate of 4/5=0.7 instead of 3/4=0.75
  // -Has weird initial transient


  testSpectralShifter(0.80, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Looks pretty good!

  // Use cosine window (not-squared) window for analysis and synthesis:
  testSpectralShifter(0.80, JH, 1024, 2, 1, true, true,   1, Mul,  0, 128, 90.0);
  // -Freq ratio looks right
  // -Has some amplitude modulation

  // OK - back to cosine-squared but with higher overlap factor of 4 to compensate for the
  // synthesis window:
  testSpectralShifter(0.80, JH, 1024, 4, 1, true, true,   2, Mul,  0, 128, 90.0);
  // -output has same freq as input and amplitude is far too low

  // ..and now also with a bit of zero-padding:
  testSpectralShifter(0.80, JH, 1024, 4, 2, true, true,   2, Mul,  0, 128, 90.0);
  // -Does no freq scale at all
  // -Amplitude is far too low and modulated


  // Now with an input freq that doesn't fit so well into the block:
  testSpectralShifter(0.80, JH, 1024, 2, 1, true, false,  2, Mul,  0, 150, 90.0);
  // -Gives really bad artifacts (discontinuities). I think these may be mitigated by a synthesis 
  //  window.

  // ...OK - so now with synthesis window:
  testSpectralShifter(0.80, JH, 1024, 2, 1, true, true,  2, Mul,  0, 150, 90.0);
  // -Produces a freq-ratio of about 8/11 and heavy amplitude and shape modulations.

  testSpectralShifter(1.25, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Looks pretty good!
  testSpectralShifter(1.5,  JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);

  // Observations
  // -For some settings, it seems to work well, for other it doesn't work at all. It may produce
  //  -wrong frequency rations and the amplitude may also be wrong. 
  // -It seems like it wants to "lock in" to certain freq ratios, for example 4/5. 
  //
  // Explanations:
  // -I think, that we sometimes get perfect alignment and sometimes total misalignment at 2048 but
  //  nothing really in between can be explained as effect of the rounding. Sometimes it seems to 
  //  be rounded to the right bin and sometimes to the wrong one. The target bin for the source bin
  //  where the sinusoid is gets always quantized to the same bin, no matter whether 
  //  k = 0.79, 0.8 or 0.81 - so all these settings lead to *exactly* the same results
  // -Whe zero padding is used, the resulting time-domain signal after IFFT may not be in the 
  //  correct part of the buffer. Maybe the fact that we get weird gains has to do with the energy
  //  of the buffer concentrated at another point that what we cut out. Maybe we should use a 
  //  shifted input and output buffer, i.e. zero-ohase windows.
  //
  // ToDo:
  // -We should really inject some polt calls to take a look at all the buffers: input, windowed
  //  input, zero-padded windowed input, spectrum, modified spectrum, IFFT result, windowed IFFT
  //  result.
  //
  // Conclusion:
  // -The JH algorithm doesn't really look too promising but maybe I'm doing something wrong. 
  //  this here is supposed to be an implementaion of the algo:
  //  https://gist.github.com/jconst/dfded80e037490d0d27fe821f18a8dee
  //  Compare the out puts of this to my implementation!


  //-----------------
  //| JH Upshifting |
  //-----------------

  testSpectralShifter(1.25, JH, 1024, 2, 1, true, false,  2, Keep, 0, 128, 90.0);
  // -After the transient has passed, i.e. from sample 1536 onwards, this result looks really good!
  //  I'd say perfect!
  // -Freq ratio is correct, amplitude is corrrect, no amp-mod and phases align periodically.
  //  ..but wait ...NOooo! the freq ration is 6/5 but it is supposed to be 5/4!

  // Now with the phase multiplier formula
  testSpectralShifter(1.25, JH, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -This looks almost perfect but the cosine peaks seem  to be shifted by one sample to the 
  //  right. This is weird!

  // Now without phase formula but with synthesis window, using cosine window for analysis and 
  // synthesis:
  testSpectralShifter(1.25, JH, 1024, 2, 1, true, true,   1, Keep, 0, 128, 90.0);
  // -Looks quite good
  // -The phases of input and output are aligned at sample indices n*512 when n >= 3. That is: the 
  //  first such alignment occurs at sample 1536, the next at 2048. Before that, we are still in 
  //  the messy transient phase. That is a strong indicator, that the formula for computing the 
  //  twiddle factor is actually correct. Without the twiddle formula, we never get a perfect phase 
  //  alignment between input and output.
  // -Has some amp-mod but not too bad.

  // Now with phase formula and synthesis window:
  testSpectralShifter(1.25, JH, 1024, 2, 1, true, true,   1, Mul,  0, 128, 90.0);
  // -some amp-mod but not too much. freq-ratio is as expected and phases align periodcially

  // To investigate the effect of the phase formula, directly compare the outputs with and without
  // using the formula:
  x  = getSpectralShifterInput(5000, 0, 128, 90);
  y1 = getSpectralShifterOutput(x, 1.25, JH, 1024, 2, 1, true, true, 1, Keep);
  y2 = getSpectralShifterOutput(x, 1.25, JH, 1024, 2, 1, true, true, 1, Mul);
  rsPlotVectors(x, y1, y2);
  // -They look almost the same - just slightly phase shifted - by just one sample.
  // -They both end abruptly at 4095 - why don't we see a smooth fade out via the synthesis window?

  testSpectralShifter(1.25, JH, 1024, 2, 1, true, true,   2, Mul,  0, 128, 90.0);
  // -Has amp-mod

  // Try to fix it with more overlap:
  testSpectralShifter(1.25, JH, 1024, 4, 1, true, true,   2, Mul,  0, 128, 90.0);
  // -This is really weird! It makes it actually much worse!

  x  = getSpectralShifterInput(5000, 0, 128, 90);
  y1 = getSpectralShifterOutput(x, 1.25, JH, 1024, 2, 1, true, false, 2, Keep);
  y2 = getSpectralShifterOutput(x, 1.25, JH, 1024, 2, 1, true, false, 2, Mul);
  rsPlotVectors(x, y1, y2);
  // -Has freq ratio of 6/5 rather than 5/4

  */

  //-----------------------------------------------------------------------------------------------
  // Experiments with my first attempt for an algorithm:

  //testSpectralShifter(0.6, RS1, 1024, 2, 8, true, false,  2, Mul,   0, 128, 90.0);


  // Now let's be brave and don't use an input window:
  //testSpectralShifter(0.80, RS1, 1024, 2, 16, false, false,  2, Mul,   0, 128, 90.0);
  // -Output has discontinuities but has the right frequency
  // -Amplitude is too low
  // -Parasitic oscillation at Nyquist freq. I guess, it comes from a sidelobe hat is not correctly
  //  reflected around zero...maybe...or maybe not. Should we actually do the reflection? Isn't 
  //  that a kind of aliasing that we should better try to avoid?




  //testSpectralShifter(0.80, RS1, 1024, 4, 8, true, true,  2, Mul,   0, 128, 90.0);
  // -Has no pitch shift at all. Output has same freq as input but smaller amplitude.
  // -Why does increasing the overlap have such an effect?
  //  -> Inspect the blocks!
  // -Maybe the pahse is wrong such that we get (soft) phase resets or some sort of osc-sync?

  // Let's ty it Without the output window:
  //testSpectralShifter(0.80, RS1, 1024, 4, 16, true, false,  2, Mul,   0, 128, 90.0);
  // -Now the output is almost silent. The blocks look good, though. The output length is 1280 
  //  samples and the input 1024. 1024 * 1.25 = 1024 / 0.8 = 1250. That checks out exactly.
  // -Must be a phase-cancellation between the blocks or something? This test clearly exposes this 
  //  behavior - so let's keep it.
  // -Multplying the phase by and ad-hoc twiddle of expC(-i * 1.0 * p); does indeed increase the
  //  output - so: yes - it has indeed to do with phase cancellation. using expC(-i * 2.0 * p);
  //  the amplitude looks (almost) good but the freq-ratio is around 6/7
  // -
  // 

  // Let's see what happens without the phase formula:
  //testSpectralShifter(0.80, RS1, 1024, 4, 16, true, false,  2, Keep,   0, 128, 90.0);
  // -Looks similar as with the formula.
  // -Using overlap of 8 makes it even more quiet
  // -With a freqScale of 0.9 we seem to get half-cancellation?




  //testSpectralShifter(0.80, RS1, 1024, 2, 8, true, false,  2, Mul,   0, 128, 90.0);

  //testSpectralShifter(0.80, RS1, 1024, 2, 16, true, false,  2, Mul,   0, 128, 90.0);

  //testSpectralShifter(0.80, RS1, 1024, 2, 1, false, false,  2, Mul,  0, 128, 90.0);



  //testSpectralShifter(0.80, RS1, 1024, 2, 1, true, false,  2, Mul,   0, 128, 90.0);
  // -The input spectrum is centerd at FFT bin 8, the output spectrum is centered at bin 6.
  //  Shouldn't it be at 0.8*8 = 6.4? But the output signal actually does have the right frequency.
  //  This is strange.
  // -Out has too low amplitude. Apparently, the problem is the linear interpolations of the 
  //  complex values. Adjacent bins have negative values. Maybe we should interpolate magnitude
  //  and phase instead. But we must think about how to unwrap the phase...or do we? Will that be a 
  //  problem? Hmm. I think, maybe taking the center bin's phase as reference and the left and 
  //  right neighbours should somehow be "unwrapped" with respect to that. Maybe zero-padding will 
  //  help to ensure some coherence of phases of neighboring bins? Yrs - that seems to be the case 
  //  indeed. The actual spectrum looks like a kind of sinusoidal blip. Interpolating magnitudes
  //  should indeed work better because the peaks are much wider than the re and im parts.
  // -Maybe let's plot the phase, too


  // Experiment to figure out a formula for the mainlobe width empirically:
  //testSpectralShifter(0.80, RS1, 1024, 2, 1, true, false,  7, Mul,   0, 128, 90.0);
  // - ZP: zero-padding, WP: window power, W: mainlobe width in bins

  // - ZP = 1, WP = 1  ->  W =  4  ...roughly (there are no sidelobes)
  // - ZP = 1, WP = 2  ->  W =  4
  // - ZP = 1, WP = 3  ->  W =  6
  // - ZP = 1, WP = 4  ->  W =  6 
  // - ZP = 1, WP = 5  ->  W =  8
  // - ZP = 1, WP = 6  ->  W =  8
  // - ZP = 1, WP = 7  ->  W =  8
  // - ZP = 1, WP = 8  ->  W = 10

  // - ZP = 2, WP = 1  ->  W =  6
  // - ZP = 2, WP = 2  ->  W =  8
  // - ZP = 2, WP = 3  ->  W = 10
  // - ZP = 2, WP = 4  ->  W = 12    D = 2

  // - ZP = 4, WP = 1  ->  W = 12
  // - ZP = 4, WP = 2  ->  W = 16
  // - ZP = 4, WP = 3  ->  W = 20
  // - ZP = 4, WP = 4  ->  W = 24    D = 4

  // - ZP = 8, WP = 1  ->  W = 24
  // - ZP = 8, WP = 2  ->  W = 32
  // - ZP = 8, WP = 3  ->  W = 40
  // - ZP = 8, WP = 4  ->  W = 48    D = 8

  // -Generally: W = 2*ZP + ZP*WP = ZP * (WP + 2)  ...formula found empirically. But it seems to 
  //  fail for ZP = 1




  //testSpectralShifter(0.80, RS1, 1024, 2, 1, true, false,  2, Keep,  0, 128, 90.0);
  // -The amplitude is too low.






  //testSpectralShifter(0.80, RS1, 1024, 2, 4, true, false,  2, Mul,  0, 128, 90.0);
  // -The phase of the output periodically aligns with the phase of the input (at peak) at samples: 
  //  1152, 1664, 2176, ... in general at: 1152 + n*512. The difference between these alignment 
  //  instants is 512 = 4*inCyc. 
  // -The output is a little bit too quiet, though  

  
  //testSpectralShifter(0.80, RS1, 1024, 2, 2, true, false,  2, Keep,  0, 128, 90.0);
  // -Without the phase formula, no phase alignment occurs but the non-aligned pitch-shifted signal 
  //  looks actually similar, too. Has also the right frequency and amplitude (??). It's just a bit 
  //  phase-shifted. ..right amplitude?...nope that's not true. that must ahve been a different
  //  setting! But I remember having seen really good outputs with correct amp and phase. Check
  //  older versions of the code for the settings that have achieved this!


  //-----------------------------------------------------------------------------------------------
  // Experiments with my second attempt for an algorithm:


  //testSpectralShifter(1.0, RS2, 1024, 2, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -All buffers are circularly shifted by half the buffer length. It doesn't seem to matter
  //  if we use the phase formula with the plus or minus in the epxonent

  //testSpectralShifter(1.0, RS2, 1024, 2, 2, true, false,  2, Mul,  0, 128, 90.0);
  // -Padded buffer is 2048 samples long
  // -Required shifts: 0: ?, 1: +256, 2: -512, 3: +512, 4: -512, 5: +512, 6: -512 (with the 
  //  positive or negative sign in the exponent)


  //testSpectralShifter(1.0, RS2, 1024, 2, 4, true, false,  2, Mul,  0, 128, 90.0);
  // -Padded buffer is 4096 samples long
  // -Required shifts with minus in the formula:
  //  0: ?, 1: ?, 2: -2560, 3: +512, 4: -512,5: -1536, 6: -2560
  // -Required shifts with plus in the formula:
  //  0: ?, 1: ?, 2: -512, 3: +512, 4: -2560, 5: -1536, 6: -512
  // -Maybe it's easier to just shift the input and output buffers instead of trying to achieve
  //  this effect via phase twiddling. Maybe to this pre- and post-shifting before calling 
  //  processSpectrum. ProcessSpectrum should receive a spectrum whose phase reference time instant
  //  is in the middle...maybe 
  // -I think, for the formula using the plus, the required phase-shift just increases by 1024 
  //  samples between the blocks and wraps around at the padded buffer length. For the formula with
  //  the minus, it decreases instead of increasing..or wait...is that true?

  //testSpectralShifter(1.0, RS2, 1024, 2, 8, true, false,  2, Mul,  0, 128, 90.0);
  // -Padded buffer is 8192 samples long
  // -Required shifts with minus in the formula:
  //  0: ?, 1: ?, 2: -2560, 3: -3584, 4: -4608, 5: -5632, 6: -6656,


  //testSpectralShifter(0.50, RS2, 1024, 8, 4, true, true,  2, Mul,  0, 128, 0.0);
  // -Very quiet output


  // Obsolete:
  // For figuring out the sampleShift formula:
  //testSpectralShifter(1.0, RS2, 1024,  4, 2, true, false,  2, Mul,  0, 128, 90.0);
  // -Let M = frameIndex
  //
  // -With O=2, P=1, every block is correct
  // -With O=2, P=2, odd blocks are correct, even block shifted by half a padded block
  // -With O=2, P=4, blocks where M % 4 = 3 are correct
  // -With O=2, P=8, blocks where M % 8 = 7 are correct
  //  -> generally, with O=2 and any P, blocks with M % P = P-1 will be correct?
  //
  // -With O=4, P=1, 

  // Different overlaps, no zero padding, cos^2 input window, no output window:
  //testSpectralShifter(1.0, RS2, 1024,  2, 1, true, false,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4, 1, true, false,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  8, 1, true, false,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024, 16, 1, true, false,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024, 32, 1, true, false,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024, 64, 1, true, false,  2, Mul,  0, 128, 90.0);
  // -Overlap 16 and higher fails. 1,2,4,8 looks good.
  // -Check the phase-shift for higher overlaps. But maybe it's a phase cancellation in the sine
  //  itself? I've adjusted the phase shift so as to put the envelope in the right position without
  //  regard to the phase of the sine within the envelope.

  // Different overlaps, no zero padding, cos^2 input window, cos^2 output window:
  //testSpectralShifter(1.0, RS2, 1024,  2, 1, true, true,  2, Mul,  0, 128, 90.0); // amp-mod
  //testSpectralShifter(1.0, RS2, 1024,  4, 1, true, true,  2, Mul,  0, 128, 90.0); // good
  //testSpectralShifter(1.0, RS2, 1024,  8, 1, true, true,  2, Mul,  0, 128, 90.0); // good
  //testSpectralShifter(1.0, RS2, 1024, 16, 1, true, true,  2, Mul,  0, 128, 90.0); // zero
  //testSpectralShifter(1.0, RS2, 1024, 32, 1, true, true,  2, Mul,  0, 128, 90.0); // zero
  //testSpectralShifter(1.0, RS2, 1024, 64, 1, true, true,  2, Mul,  0, 128, 90.0); // zero
  // -With overlap = 2, we see ampltude modulation - but that is exactly as expected
  // -Overlap = 4,8 work fine
  // -Overlap >= 16 show the phase cancellation issue again - output is close to zero.

  // Overlap = 4, cos^2 window for input and output, different zero-paddings
  //testSpectralShifter(1.0, RS2, 1024,  4,  1, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4,  2, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4,  4, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4,  8, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4, 16, true, true,  2, Mul,  0, 128, 90.0);
  // -They look all fine.

  // Test a cos^4 window. We expect to need more overlap than for cos^2. we use a fixed 
  // zero-padding factor of 4::
  //testSpectralShifter(1.0, RS2, 1024,  4,  4, true, true,  4, Mul,  0, 128, 90.0); // too quiet
  //testSpectralShifter(1.0, RS2, 1024,  8,  4, true, true,  4, Mul,  0, 128, 90.0); // too quiet
  //testSpectralShifter(1.0, RS2, 1024, 16,  4, true, true,  4, Mul,  0, 128, 90.0); // zero
  // -I expected to see amp-mod for overlap = 4 which should go away when increasing overlap to 8.
  //  But instead, both 4 and 8 give a non-modulated signal which is a bit too quiet. Going up to
  //  16, we get the zero signal again.

  // I think, using overlap = 4, cos^2 for input and output should be the best choice. Zero-padding
  // can be adjusted to taste. Maybe 1,2,4 are reasonable.

  // Try some other blockSizes
  //testSpectralShifter(1.0, RS2, 2048,  4,  8, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2, 1024,  4,  8, true, true,  2, Mul,  0, 128, 90.0);
  //testSpectralShifter(1.0, RS2,  512,  4,  8, true, true,  2, Mul,  0, 128, 90.0);  // OK
  //testSpectralShifter(1.0, RS2,  256,  4,  8, true, true,  2, Mul,  0, 128, 90.0);  // Nope!
  //testSpectralShifter(1.0, RS2,  256,  2,  8, true, true,  2, Mul,  0, 128, 90.0);  // fast amp-mod?
  // 

  // With B=512, O=4, P=2, W=cos^2 (I/O), try different input periods:
  //testSpectralShifter(1.0, RS2,  512, 4, 2, true, true,  2, Mul,  0,  64, 90.0);  // OK
  //testSpectralShifter(1.0, RS2,  512, 4, 2, true, true,  2, Mul,  0, 128, 90.0);  // OK
  //testSpectralShifter(1.0, RS2,  512, 4, 2, true, true,  2, Mul,  0, 256, 90.0);  // Nope!
  // -I think, we get problems as soon as the input period is as long as 2*H?
  //  -> check the phase


  // Let's try to confirm this by trying a cycle length of 256 and different ways of achievinh 
  // H=256 (as 512/2 and 1024/4):
  //testSpectralShifter(1.0, RS2,  512, 2, 2, true, false, 2, Mul,  0, 256, 90.0);  // H=256 - OK
  //testSpectralShifter(1.0, RS2, 1024, 4, 2, true, true,  2, Mul,  0, 256, 90.0);  // H=256 - OK
  //testSpectralShifter(1.0, RS2,  512, 2, 2, true, true,  2, Mul,  0, 256, 90.0);  // H=256 - AmpMod
  //testSpectralShifter(1.0, RS2,  512, 2, 2, true, true,  1, Mul,  0, 256, 90.0);  // H=256 - mmhh




  // ...

  // Check, if this phase-cancellation phenomenon happens later or earlier with lower or higher
  // input freqs


  // Now let's try some actual shift:
  //testSpectralShifter(1.0, RS2, 1024,  4,  8, true, true,  2, Mul,  0, 128, 90.0);
  // ...aah - it's not yet implemented in the experimental RS2 algo


  //testSpectralShifter(1.0, RS2, 1024, 2, 4, true, false,  2, Mul,  0, 128, 90.0);
  // -Input spectrum's real and imag part wiggle with an oscillation period of 2*ZP bins under 
  //  the magnitude envelope


  //testSpectralShifter(1.0, RS2, 1024, 4, 4, true, true,  2, Mul,  0, 128, 0.0);


  testSpectralShifter(0.5, RS2, 1024, 2, 4, true, false,  2, Mul,  0, 128, 0.0);
  // -When tweakking the freq-scale between 0.5 and 1.0 it seems, it want to "lock in" only to 
  //  ratios of 0.5, 0.8 and 1.0. Anything in between looks like one of these ratios but with 
  //  artifacts.





  testSpectralShifter(0.75, RS2, 1024, 4, 4, true, true,  2, Mul,  0, 128, 0.0);
  // -Looks totally wrong!


  testSpectralShifter(0.50, RS2, 1024, 4, 4, true, true,  2, Mul,  0, 128, 0.0);
  // -Is a little bit too loud but looks pretty good otherwise.
  // -Maybe try a spectral energy normalization

  //testSpectralShifter(0.50, RS2, 1024, 2, 4, true, true,  1, Mul,  0, 128, 0.0);
  // -Has amplitude modulation artifact and is a bit too loud

  //testSpectralShifter(0.80, RS2, 1024, 4, 4, true, true,  2, Mul,  0, 128, 0.0);
  // -Totally wrong: output freq equlas input freq, amplitude is too low

  //testSpectralShifter(0.80, RS2, 1024, 2, 4, true, false,  2, Mul,  0, 128, 0.0);

  //testSpectralShifter(0.80, RS2, 1024, 2, 2, true, false,  2, Mul,  0, 128, 0.0);
  // -The whole output buffer needs a circular shift. The freq looks about right.
  // -It seems liek the amount of circular shift is not the same in each buffer.
  // -...OK - with a constant phase-shift, we get closer. But the output has amp-mod and some 
  //  discontinuities


  //testSpectralShifter(0.80, RS2, 1024, 2, 4, true, false,  2, Mul,  0, 128, 0.0);
  // -Hmm - but with zero padding of 4, the phase shift is totally different. Now the center of
  //  gravity is in the center of the padded block

  //testSpectralShifter(0.80, RS2, 1024, 2, 8, true, false,  2, Mul,  0, 128, 0.0);

  //testSpectralShifter(0.80, RS2, 1024, 2, 1, true, false,  2, Mul,  0, 128, 0.0);


  //-----------------------------------------------------------------------------------------------
  // Older:
  // -We use the Juillerat/Hirsbrunner algorithm with blockSize of 1024 and a sinusoidal input with
  //  a cycle length of 128 samples ...TBC..
  testSpectralShifter(0.80, JH, 1024, 2, 4, true, true,  2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.80, JH, 1024, 2, 4, true, false, 2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.80, JH, 1024, 2, 2, true, false, 2, Mul,  0, 128, 90.0);
  testSpectralShifter(0.80, JH, 1024, 2, 2, true, false, 2, Keep, 0, 128, 90.0);
  testSpectralShifter(1.25, JH, 1024, 2, 2, true, false, 2, Mul,  0, 128, 90.0);


  int dummy = 0;

  // -Shouldn't we accumulate the phase shifts/twiddles? Maybe if we do, we can get rid of that 
  //  weird formula with the modulo operation? It seems to be unnatural to have the frameIndex
  //  explicitly entering the computations.
  // -Maybe test it with an impulse-like input signal.
  // -We need to reflect frequencies below zero. Maybe with complex conjugation or something.
  // -Try to interpolate magnitude (squared?) and phase instead. Or maybe let the phases run freely
  //  and onle reset it on transients to the input phase
  // -Try energy normalization - or RMS try it in time and freq-domain to see what works better
  // -I think, the amount of zero padding also limits by hwo much we can down-shift because it 
  //  limits by how much a block can be lengenthed. Zero-padding of 2 should allow a downshift by
  //  1 octave, zero-padding by 4 a downshift of two octaves. ...Or maybe not? 
  // -Implement unit tests and/or experiments with the following test signals: 
  //  -sines whose period fits nicely in the window
  //  -sine whose freq fits badly (half a cycle overhang or something)
  //  -mix of sines
  //  -DC
  //  -Nyquist freq
  //  -Impulses
  // -On these signals: 
  //  -try identity resynthesis (with delay)
  //  -shift octave down
  //  -shift ovtave up
  //  -shifts by 0.5, 0.75, 0.8, 1.0, 1.2, 1.25, 1.5, 2.0  ...maybe also try something very 
  //   irrational such as the golden ratio and its reciprocal

  // This function should eventually replace testSpectralShiftViaJH/RS. We need to go through the 
  // comments there and for each setting make a correspoding function call here and the copy the 
  // observations´directly under the call to the experiment. That makes it much cleaner and more
  // orderly to conduct the experiments and document the results.
}



// Soon to be obsolete:
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
  // BUT: figure out if this should have an effect on the phase formula! Is it really correct to 
  // use the padding factor in the phase-shift?



  // Setup:

  // Output file parameters:
  double sampleRate  = 44100;         // Sample rate for the signals in Hz
  int    numSamples  = sampleRate/10; // We create a 1/10 seconds long signal.

  // Input signal parameters:
  double inputPeriod = 128;           // Length of one cycle in samples
  double inputPhase  = 90;            // Phase in degrees

  // Spectral shifter parameters:
  double freqScale   = 0.80;          // Scaling factor for the frequencies
  int    blockSize   = 1024;          // Block size. Must be power of 2
  int    overlap     = 2;             // Overlap factor. Must be power of 2
  int    zeroPad     = 4;             // Zero padding factor. Must be power of 2
  bool   anaWindow   = true;          // Use analysis window or not
  bool   synWindow   = false;          // Use synthesis window or not



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
  rosic::SpectralShifter ps(blockSize, overlap, zeroPad);
  ps.setAlgorithm(SS::Algorithm::JuilHirs);
  ps.setFrequencyScale(freqScale);
  ps.setInputBlockSize(blockSize);
  ps.setOverlapFactor(overlap);
  ps.setPaddingFactor(zeroPad);
  ps.setUseInputWindow(anaWindow);
  ps.setUseOutputWindow(synWindow);
  ps.setPhaseFormula(SS::PhaseFormula::useMultiplier);
  Vec y1(numSamples);
  for(int n = 0; n < numSamples; n++)
    y1[n] = ps.getSample(x[n]);

  // For reference, generate output without without the phase formula:
  ps.setPhaseFormula(SS::PhaseFormula::keepOriginal);
  ps.reset();
  Vec y2(numSamples);
  for(int n = 0; n < numSamples; n++)
    y2[n] = ps.getSample(x[n]);


  // Plot input and output signals:
  rsPlotVectors(x, y1);
  //rsPlotVectors(x, y2);

  // Observations:
  // -With freqScale = 1.2, The first few buffers look like a mess but then it settles and the only
  //  remaining artifact is a (strong) amplitude modulation. It looks less messy when we use a 
  //  start-phase of zero instead of 90, i.e. sine instead of cosine. In this case, there are only 
  //  corners in the signal whereas with cosine, there are jumps.
  // -With freqScale = 2, the discontinuities in the first frames disappear but the amp modulation 
  //  gets worse. There is also some amount of negative DC between 512 and 768 when the phase is 
  //  zero. With 90, this is not the case. The amp-modulation period is 512 samples.


  // -freqScale = 0.9; blockSize 1024; overlap = 2; zeroPad = 1;
  //  inputFreq = 128; inputPhase = 90;
  //  -> very strong amplitude modulation
  //  -> at 2048 samples, the shifted signal has the opposite phase to the input. This was one of 
  //     the instants where where k = 0.8 and k = 1.25 gave a perfect phase match.
  //  -> it happens also for k = 0.85. 
  //  -> at k = 0.825, there's no modulation, opposite alignment at 2048, perfect alignment at 1792
  //     and the signal is too quiet by a factor of 2.
}

// Soon to be obsolete:
void testSpectralShiftViaRS()
{
  // Tests an algorithm that was my own idea. ...TBC...

  // Setup:

  // Output file parameters:
  double sampleRate   = 44100;         // Sample rate for the signals in Hz
  int    numSamples   = sampleRate/10; // We create a 1/10 seconds long signal.

  // Input signal parameters:
  double inputPeriod  = 128;           // Length of one cycle in samples
  double inputPhase   = 90;            // Phase in degrees

  // Spectral shifter parameters:
  double freqScale    = 0.8;           // Scaling factor for the frequencies
  int    blockSize    = 1024;          // Block size. Must be power of 2
  int    overlap      = 2;             // Overlap factor. Must be power of 2
  int    zeroPad      = 4;             // Zero padding factor. Must be power of 2
  bool   anaWindow    = true;          // Use analysis window or not
  bool   synWindow    = false;         // Use synthesis window or not
  int    winPower     = 2;             // power/exponent for the cos^n window


  // Create sinusoidal test signal:
  using Vec = std::vector<double>;
  double inputFreq = sampleRate / inputPeriod;
  int N = numSamples;
  Vec x(N);
  createWaveform(&x[0], N, 0, inputFreq, sampleRate, RAPT::rsDegreeToRadiant(inputPhase), true);
  x = 0.5 * x;

  using PS = rosic::SpectralShifter;
  PS ps(blockSize, overlap, zeroPad);
  ps.setAlgorithm(PS::Algorithm::RobSchm1);
  ps.setFrequencyScale(freqScale);
  ps.setInputBlockSize(blockSize);
  ps.setOverlapFactor(overlap);
  ps.setPaddingFactor(zeroPad);
  ps.setUseInputWindow(anaWindow);
  ps.setUseOutputWindow(synWindow);
  ps.setWindowPower(winPower);
  ps.setPhaseFormula(PS::PhaseFormula::useMultiplier);
  Vec y1(N);
  for(int n = 0; n < N; n++)
    y1[n] = ps.getSample(x[n]);

  // For reference, generate output with different settings: 
  //ps.setPhaseFormula(PS::PhaseFormula::keepOriginal); // without the phase formula
  ps.setAlgorithm(PS::Algorithm::JuilHirs);             // with other algo
  ps.reset();
  Vec y2(numSamples);
  for(int n = 0; n < numSamples; n++)
    y2[n] = ps.getSample(x[n]);


  // Plot input and output signals:
  rsPlotVectors(x, y1);
  //rsPlotVectors(x, y2);
  //rsPlotVectors(x, y1, y2);


  // Write input and output into wave files:
  //rosic::writeToMonoWaveFile("SpectralShifterInput.wav",  &x[0], N, sampleRate, 16);
  //rosic::writeToMonoWaveFile("SpectralShifterOutput.wav", &y[0], N, sampleRate, 16);
  int dummy = 0;

  // We use the following abbrevations:
  // inCyc: inputPeriod (cycle), inPhs: inputPhase
  // frqScl: freqScale, blkSz: blockSize, ovrLp: overlap, zrPd: zeroPad,  
  // anaWn: anaWindow, synWn: synWindow, WnPw: window power

  // Observations:
  //
  /// -inCyc=128, inPhs=90, frqScale=0.8, blkSz=1024, ovrLp=2, zrPd=4, anaWn=yes, synWn=no, WnPw=2:
  ///  -Using the phase formula causes the phase of the output periodically align with the phase
  ///   of the input (at peak) at samples: 1152, 1664, 2176, ... in general at: 1152 + n*512
  ///   The difference between these alignment instants is 512 = 4*inCyc. Without the formula, no
  ///   alignment occurs but the non-aligned pitch-shifted signal looks actually good, too. Has also
  ///   the right frequency and amplitude. It's just a bit phase-shifted.
  ///  -The output is a little bit too quiet, though

  //  -Trying to lower the frqScl to 0.75 or 2./3 does not really seem to lower the pitch of the
  //   output sine. WTF? It's always 5 input peaks and 4 output peaks between the insteants of
  //   alignments. It's as if the freqScale factor has locked in to a quantized value of 0.8=4/5
  //  -Lowering it further to 0.6=3/5 seems to produce the right frequency - but the amplitude is 
  //   too low and there are no instants of phase alignment. Increasing the zrPd does not help.
  //  -Setting frqScl=1/3 produces garbage - just one blip
  //  -Setting frqScl=1.1 produces strong amp-mod plus some sort of phase-mod and the output 
  //   doesn't seem to have the right freq. There is no real single freq anyway
  //  -setting frqScl=1.2 works better again. The output does have the right freq. But there's some
  //   amp-mod going on.
  //
  // -inCyc=128, inPhs=90, frqScale=1.0, blkSz=1024, ovrLp=2, zrPd=2, anaWn=yes, synWn=yes, WnPw=1:
  //  -perfect reconstruction as expected. when using both windows with a power of 1, the combined
  //   ana/syn windows will give cos^2 which is suitable for overlap of 2
  //  -frqScl=0.5 is too quiet and has amp-mod
  //  -frqScl=0.8 is okayish, 0.75 produces roughly the same frequency

  // -When using the synthesis window, the output is too quiet.
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
  // -With zeroPad = 1, the output is too quite. I guess we need at least padding of tow because
  //  the output my become longer than the input due to these time aliasing effects
  // -When the input period is 500 with blockSize = 1024 and freqScale = 1.1, it looks like the 
  //  freq has not really been changed at all
  // -Without analysis window, it doesn't seem to work at all
  //
  // Conclusions:
  // -zeroPad needs to be at least 2. 4 seems to work well enough for a nice signal and shift 
  //  factor
  //
  // ToDo:
  // -[Done] Add the phase multiplication step. 
  // -[Done] Use output window
  // -Try it on a sinusoid that aligns with an FFT bin. Then try it on a  sinusoid that is in 
  //  between two bins.
  //
  // Ideas:
  // -Try to preserve energy of spectrum. That could perhaps fix the problem of the wrong 
  //  amplitudes. Compute energy of input and output spectrum and scale output spectrum 
  //  accordingly. But maybe it's better to do that in the time domain because then it will also
  //  take all the windows into account. Maybe that can be done in class OverlapAddProcessor. It 
  //  could have a member normalizeEnergy. Or it could have a member normalizationMode with 
  //  settings none, peak, energy, rms (I think, energy and rms will give the same results). Maybe
  //  the setting can be set up continuously from none to full (by raising the resulting 
  //  normalization factor to a power)
  // -experiment with setWindowPower. The default is 2 which is the right choice when overlap=2
  //  and only an analysis window is used. When both analysis and synthesis windows are used,
  //  we could use 1 for the window power or 4 for the overlap.
}



void rotes::spectralShifter()
{
  // Under construction - does not yet work.
  //
  // We want to build a pitch shifter based on spectral processing. It should have a transient 
  // preservation feature.

  testSpectralShift();
  testSpectralShiftViaJH();
  testSpectralShiftViaRS();


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
  // https://github.com/Signalsmith-Audio/pitch-time-example-code
  // https://signalsmith-audio.co.uk/writing/2023/stretch-design/
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


  // Paper: Low latency audio pitch shifting in the time domain
  // https://ieeexplore.ieee.org/document/4590019

  

}