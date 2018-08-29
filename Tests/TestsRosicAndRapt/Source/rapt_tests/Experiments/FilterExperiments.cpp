#include "FilterExperiments.h"
using namespace RAPT;

void bandSplittingTwoWay()
{
  // user parameters:
  float sampleRate = 44100;
  float splitFreq  = sampleRate/4; // 11025

  // set up the splitter:
  float omega = 2*float(PI)*splitFreq/sampleRate;
  rsTwoBandSplitterFF splitter;
  splitter.setOmega(omega);

  // plot poles/zeros and magnitude responses:
  FilterPlotter<float> plt;
  // ...
}

void bandSplittingMultiWay()
{
  // user parameters:
  int numSamples = 100;
  float sampleRate = 44100;
  //vector<float> splitFreqs = { 100, 300, 1000, 3000, 10000 };
  //vector<float> splitFreqs = { 80, 150, 250, 500, 1000, 2000, 4000, 8000 };
  vector<float> splitFreqs = { 125, 250, 500, 1000, 2000, 4000, 8000 };  // 8 bands
  //vector<float> splitFreqs = { 60, 90, 135, 200, 300, 450, 675, 1000, 1500, 2200, 3000, 4500, 6000, 9000, 13000 };  // 16 bands

  // set up dsp objects and signal arrays:
  int n, k;
  rsMultiBandSplitterFF splitter;
  splitter.setSplitFrequencies(splitFreqs);
  splitter.setSplitMode(rsMultiBandSplitterFF::ACCUMULATE_INTO_LOWPASS);
  //splitter.setSplitMode(rsMultiBandSplitterFF::ACCUMULATE_INTO_HIGHPASS);
  //splitter.setSplitMode(rsMultiBandSplitterFF::BINARY_TREE);
  int numBands = splitter.getNumActiveBands();
  std::vector<float> x(numSamples);
  std::vector<std::vector<float>> y(numBands);
  for(k = 0; k < numBands; k++)
    y[k].resize(numSamples);
  RAPT::rsNoiseGenerator<float> ng;
  std::vector<float> tmp(numBands);  // holds the array of bandpass output samples

  // produce input and output signals:
  x[0] = 1;
  for(n = 0; n < numSamples; n++)
  {
    //x[n] = ng.getSample();
    splitter.processSampleFrame(x[n], &tmp[0]);
    for(k = 0; k < numBands; k++)
      y[k][n] = tmp[k];

    // test, if the band outputs recombine correctly:
    float xr = 0.f;
    for(k = 0; k < numBands; k++)
      xr += y[k][n];
    float error = x[n] - xr;
    rsAssert(fabs(error) < 1.e-7);
  }

  // write outputs to files:
  rosic::writeToMonoWaveFile("Input.wav", &x[0], numSamples, 44100, 16);
  for(k = 0; k < numBands; k++)
    writeToMonoWaveFile("Output"+to_string(k)+".wav", &y[k][0], numSamples, (int) sampleRate, 16);

  //FilterPlotter<float> plt;
  GNUPlotter plt;
  for(k = 0; k < numBands; k++)
    plt.addDataArrays(numSamples, &y[k][0]);
  plt.plot();
}

//void append(std::string& str, std::string& appendix


int updateBandSplitStrings(std::vector<std::string>& str, int start, int length)
{
  static int numCalls = 0;
  numCalls++;

  if(length == 1)
    return numCalls;

  int L = length/2;
  str[start+L] = str[start] + "H";
  str[start]   = str[start] + "L";  

  // recursion:
  updateBandSplitStrings(str, start+L, L); // upper branch (highpass)
  updateBandSplitStrings(str, start,   L); // lower branch (lowpass)

  return numCalls;
}
void bandSplittingTreeAlgo()
{
  // We test the strategy to build a binary tree of splits symbolically by using strings that
  // represent the path of partial signal along the splitter tree, for example, "LHL" means
  // the signal went through a lowpass then a highpass then another lowpass

  int numBands = 8; // algo currently assumes a power of 2
  std::vector<std::string> str1(numBands);
  int numCalls = updateBandSplitStrings(str1, 0, numBands); 
    // numCalls should be 2*numBands-1 (for numBands a power of two)

  // ok, the recursive function seems to work - now let's try to convert it into an iterative algo
  std::vector<std::string> str2(numBands);
  int inc = numBands; // nextPowerOfTwo(numBands) for non-powers-of-2?
  int numIts = 0;
  while(inc > 1)
  {
    int pos = 0;
    while(pos < numBands)
    {
      str2[pos+inc/2] = str2[pos] + "H"; // maybe call inc/2 offset
      str2[pos]       = str2[pos] + "L";
      pos += inc;
      numIts++;
    }
    inc /= 2;
  }
  // ok - this seems to work, too. In the actual splitter, we need to initialize y[0] with the 
  // input and then run exactly that type of iteration using y[pos] as input to the two-way 
  // splitters, y[pos] (also) as lowpass output and y[pos+inc/2] as highpass output

  bool works = (str2 == str1);
}

void bandSplitFreqResponses()
{
  // Computes and plots frequency response curves of the multiband bandsplitter outputs. The 
  // response is calculated in two ways: (1) FFT of the impulse-response and (2) using the built-in
  // getMagnitude... functions. The prupose is mainly to compare the two to verify that 
  // implementation of the latter is correct.

  // maybe turn into unit test

  typedef std::vector<float> Vec;

  // user parameters:
  int N = 2048;     // should be power of two for easy FFT
  float sampleRate = 44100;
  //Vec splitFreqs = { };
  //Vec splitFreqs = { 10000 };
  //Vec splitFreqs = { 7000, 14000 };
  Vec splitFreqs = { 5000, 10000, 15000 };
  //Vec splitFreqs = { 4000, 8000, 12000, 16000 };
  //Vec splitFreqs = { 3000, 6000, 9000, 12000, 15000, 18000 };

  // create and set up the splitter:
  rsMultiBandSplitterFF splitter;
  splitter.setSampleRate(sampleRate);
  splitter.setSplitFrequencies(splitFreqs);
  splitter.setSplitMode(rsMultiBandSplitterFF::ACCUMULATE_INTO_LOWPASS);
  //splitter.setSplitMode(rsMultiBandSplitterFF::ACCUMULATE_INTO_HIGHPASS);
  //splitter.setSplitMode(rsMultiBandSplitterFF::BINARY_TREE);

  // create and set up the fourier transformer:
  typedef rsFourierTransformerRadix2<float> FT;
  rsFourierTransformerRadix2<float> ft;
  ft.setBlockSize(N);
  //ft.setNormalizationMode(FT::NORMALIZE_ON_FORWARD_TRAFO);
  ft.setNormalizationMode(FT::NORMALIZE_ON_INVERSE_TRAFO);

  //
  GNUPlotter plt;
  int numBands = splitter.getNumActiveBands();
  Vec tmp(numBands);  // array of bandpass output samples
  Vec fftFreqs(N);    // fft bin frequencies
  Vec impResp(N);     // impulse response of current band
  Vec magFFT(N);      // magnitude response computed by FFT
  Vec magTF(N);       // magnitude response computed by transfer function
  int i, k, n;        // indices for band, fft-bin and sample
  for(k = 0; k < N; k++)
    fftFreqs[k] = ft.binIndexToFrequency(k, N, sampleRate);
  for(i = 0; i < numBands; i++)
  {
    splitter.reset();

    // compute impulse response of i-th band:
    splitter.processSampleFrame(1, &tmp[0]);
    impResp[0] = tmp[i];
    for(n = 1; n < N; n++) {
      splitter.processSampleFrame(0, &tmp[0]);
      impResp[n] = tmp[i];
    }

    // get FFT magnitudes from impulse response:
    ft.getRealSignalMagnitudes(&impResp[0], &magFFT[0]);

    // get band magnitude reponse from transfer function:
    for(k = 0; k < N; k++)
      magTF[k] = splitter.getBandMagnitudeAt(i, fftFreqs[k]);

    // now, the contents of magFFT and magTF should be the same up to errors due to 
    // truncation of impulse-reponse before FFT and roundoff - plot them:
    //if(i == 3) // just to pick out one for test
    {
      plt.initialize();
      plt.plotFunctionTables(N, &fftFreqs[0], &magFFT[0], &magTF[0]);
    }
    // the last 2 are correct, the others only similar
  }
}

//-------------------------------------------------------------------------------------------------
// experiments for perfect reconstruction IIR filters:

typedef std::complex<double> Complex;

// analog 1-pole/0-zero
void splitterPrototypeA_1_0(double* k, Complex* p, Complex* z)
{
  *k   =  1.0;
  p[0] = -1.0;
}

// analog 2-pole/0-zero
void splitterPrototypeA_2_0(double* k, Complex* p, Complex* z)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = Complex(-s, s);
  p[1] = conj(p[0]);
}

// analog 2-pole/1-zero
void splitterPrototypeA_2_1(double* k, Complex* p, Complex* z)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = Complex(-s, s);
  p[1] = conj(p[0]);
  z[0] = -1;
  // when we place a zero on the real axis, the ultimate slope is again only 6 dB/oct - seen from a 
  // distance (i.e. being far away on the imaginary axis) the zero cancels the effect of one of the
  // poles
}

// analog 2-pole/2-zero
void splitterPrototypeA_2_2(double* k, Complex* p, Complex* z)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = Complex(-s, s);
  p[1] = conj(p[0]);
  z[0] = Complex(0, 2);
  z[1] = conj(z[0]);
}

/*
// move to RAPT:
//Complex getDigitalTransferFunctionAt(Complex* zeros, int numZeros, Complex* poles, int numPoles,
Complex digitalTransferFunctionZPK(const Complex* zeros, int numZeros, 
  const Complex* poles, int numPoles, Complex k, Complex z)
{
  Complex zr  = 1.0/z;       // z^-1
  Complex num = 1, den = 1;  // numerator and denominator
  for(int i = 0; i < numZeros; i++)  num *= (1.0 - zeros[i] * zr);
  for(int i = 0; i < numPoles; i++)  den *= (1.0 - poles[i] * zr);
  return k * num/den;
  // see https://ccrma.stanford.edu/~jos/filters/Factored_Form.html for formula
} // moved to rsFilterSpecificationZPK<T>::transferFunctionAt

Complex transferFunctionZPK(const FilterSpecificationZPK<double>& zpk, Complex z)
{
  if(zpk.sampleRate != RS_INF(double))
    return digitalTransferFunctionZPK(&zpk.zeros[0], (int)zpk.zeros.size(), 
      &zpk.poles[0], (int)zpk.poles.size(), zpk.gain, z);
  else
  {
    //return analogTransferFunctionZPK(&zpk.zeros[0], (int)zpk.zeros.size(), 
    //  &zpk.poles[0], (int)zpk.poles.size(), zpk.gain, z); // z should be named s here ..maybe sz? s_or_z?
    rsAssertFalse; // not yet implemented
    return 0;
  }
} // moved to rsFilterSpecificationZPK<T>::transferFunctionAt

Complex digitalTransferFunctionBA(const Complex* b, int Nb, const Complex* a, int Na, Complex z)
{
  Complex num = 0, den = 0;
  for(int i = 0; i < Nb; i++)  num += b[i] * pow(z, -i); // can be optimized
  for(int i = 0; i < Na; i++)  den += a[i] * pow(z, -i);
  return num/den;
}

Complex transferFunctionBA(const FilterSpecificationBA<double>& ba, Complex z)
{
  if(ba.sampleRate != RS_INF(double))
    return digitalTransferFunctionBA(&ba.b[0], (int)ba.b.size(), &ba.a[0], (int)ba.a.size(), z);
  else
  {
    rsAssertFalse; // not yet implemented
    return 0;
  }
}
*/





double dcGainNormalizer(Complex* zeros, int numZeros, Complex* poles, int numPoles)
{
  // set gain factor k to normalize DC gain to 1:
  Complex H1 = digitalTransferFunctionZPK(
    zeros, numZeros, poles, numPoles, Complex(1.0, 0.0), 
    Complex(1.0, 0.0)); // H(z) at z=1
  return 1.0 / abs(H1);  
}

// digital 1-pole/1-zero - works
void splitterPrototypeD_1_1(double* k, Complex* p, Complex* z)
{
  //typedef rsInfiniteImpulseResponseDesigner<double> DSGNR;
  //DSGNR dsgnr;
  //dsgnr.setApproximationMethod(rsPrototypeDesigner<double>::BUTTERWORTH);
  //dsgnr.setPrototypeOrder(1);
  //dsgnr.setSampleRate(1);
  //dsgnr.setFrequency(0.25);      // halfband filter, 0.5 is Nyquist freq
  //dsgnr.getPolesAndZeros(p, z);


  p[0] =  0;
  z[0] = -1;

  *k = dcGainNormalizer(z, 1, p, 1);
}

// digital 2-pole/2-zero (i think, this is a 2nd order digital Butterworth halfband filter via 
// bilinear transform):
void splitterPrototypeD_2_2(double* k, Complex* p, Complex* z)
{
  double s = sqrt(2)-1;  // ad-hoc, i have no derivation for this
  p[0] = Complex(0, s);
  p[1] = conj(p[0]);
  z[0] = -1.0;
  z[1] = conj(z[0]);
  *k = dcGainNormalizer(z, 2, p, 2);

  //Complex H1 = digitalTransferFunctionZPK(z, 2, p, 2, 1, Complex(1.0, 0.0)); // H(z) at z=1
  //*k = 1 / abs(H1);  
}
// this doesn't work

// digital 2-pole/3-zero - old - worked when we didn't really treat the reversed order of digital
// filter coeff arrays right in conversions between zpk/ba form
void splitterPrototypeD_2_3(double* k, Complex* p, Complex* z)
{
  double  s = sqrt(2)-1;
  //s = 0.5;  // test
  Complex j = Complex(0, 1);
  p[0] =  j*s;   // p1
  p[1] = -j*s;   // p2
  z[0] = -1;     // q1
  z[1] = -1;     // q2
  z[2] = -s;     // q3
  //z[2] = -0.5;   // test
  *k = dcGainNormalizer(z, 3, p, 2);  // gain factor k to normalize DC gain to 1

  // I arrived at these poles and zeros by just starting with a (bilinear transform based) 
  // Butterworth halfband lowpass which fixed p1,p2 and q1,q2 to the values above and manually 
  // added the 3rd zero q3 by trial and error. It turned out that it had to be placed along the 
  // negative real axis the same distance as the two poles sit along the (positive and negative)
  // imaginary axes. Unfortunately, such a simple strategy doens't seem to generalize to higher 
  // order Butterworth filters. Maybe, instead, we have to consider the formula:
  // https://ccrma.stanford.edu/~jos/filters/Factored_Form.html and write it down for the 
  // particular N,M (here N=2, M=3):
  //
  //             (1-q1/z)*(1-q2/z)*(1-q3/z)                        (1-r1/z)*(1-r2/z)*(1-r3/z)
  // H(z) = k * ----------------------------, G(z) = 1-H(z) = c * ----------------------------
  //                 (1-p1/z)*(1-p2/z)                                 (1-p1/z)*(1-p2/z)
  //  
  // define some constraints and solve for the poles, zeros and k (or, alternatively, for the 
  // polynomial coefficients). We have: G(z) = 1-H(z) = 1 - (B(z)/A(z)) = (A(z)-B(z))/A(z) where 
  // G(z) is the complementary filter to H(z), i.e. the highpass that is complementary to the 
  // lowpass H(z) in the sense that is obtained by subtracting the lowpass output signal from the 
  // unfiltered input. Maybe in this case, we should have required q1=q2=-1, r1=r2=+1, r3=-q3, 
  // p1=-p2,|H(1)|=1, |H(-1)|=0, |G(1)|=0, |G(-1)|=1, maybe |H(z)| = |G(-z)| in general
  // (here: r1,r2,r3 are the zeros of the highpass, q1,q2,q3 are the zeros of the lowpass - the 
  // poles are the same in lowpass and highpass).

  // Try, if the poles/zeros for the (2,3)-case can be recovered with this method and then try to 
  // generalize the method to higher order filters (maybe (3,4), (4,5), etc?). Maybe try first the
  // (1,1)-case. If no general formula can be derived, maybe it's possible to devise a numerical 
  // algorithm (based on multi-dimensional root-finding or gradient-descent) to find suitable 
  // poles/zeros and tabulate them for various orders.

  // I think, A(z)-B(z) should have its roots at positions opposite to those of B(z), i.e. 
  // reflected about the imaginary axis. Maybe we need the condition A(z)-B(z) = B(-z)?
  // Or A(z)-B(z) = B(-conj(z))? Maybe, we should take the poles already as given (may have to be
  // determined by other conditions such as monotonicity of magnitude response?). Let's define
  // C(z) = A(z)-B(z). It seems, the condition C(z)=B(-z) leads to a1=0, a2=2*b2 - which holds for 
  // the 2,3 filter. It also seems liek C(z) has the coeffs of B(z) but the the odd-numbered coeffs 
  // sign-inverted - which makes sense because a sign-change in the input translates to a sign 
  // change in odd-numbered coefficients. Maybe when the poles and some of the zeros are fixed, 
  // equations for the remaining coeffs/zeros can be obtained? Actually, maybe we should require
  // G(z)=H(-z) instead of C(z)=B(-z), but if we place all poles on the imaginary axis, we will 
  // have A(z)=A(-z), so G(z)=H(-z) = C(z)/A(z) = B(-z)/A(-z) -> C(z)=B(-z). Let's consider the
  // transfer functions in sum form - here the 2,3 case:
  // 
  //         b0 + b1/z + b2/z^2 + b3/z^3           c0 + c1/z + c2/z^2 + c3/z^3
  // H(z) = -----------------------------, G(z) = -----------------------------
  //         a0 + a1/z + a2/z^2                    a0 + a1/z + a2/z^2
  //
  // by:  C(z)=B(-z)   C(z)=A(z)-B(z)
  // c0 =    b0      = a0-b0             // we can't assume a0=1 at this point
  // c1 =   -b1      = a1-b1
  // c2 =    b2      = a2-b2
  // c3 =   -b3      = a3-b3 = 0-b3      // ok, -b3 = 0-b3 is not very useful
  //
  // For each even i, we get two equations: ci =  bi, ci = ai-bi -> bi = ai-bi -> bi = 2*ai and
  // for each odd i, we get:                ci = -bi, ci = ai-bi -> ai = 0, so the general rules
  // seems to be:
  // -odd-numbered a-coeffs must be 0 (yes, that also holds for the 1,1 case)
  // -even numbered b-coeffs must be half of the corresponding a-coeffs
  // -poles must be on the imaginary axis (this places further constraints on the a-coeffs)
  // -N (= numPoles) zeros must be at z=-1 (to get a proper lowpass response) 
  // -maybe the other M-N zeros should also be in the left half-plane?

  //
  // hmmm...with these equations together with the poles and some of the zeros already fixed, we 
  // may be able to solve for the remaining coeffs/zeros.

  // is the condition H(z)=G(-z) actually correct? or should it be |H(z)|=|G(-z)|?
  // -verify numerically, if H(z)=G(-z) for the 2,3 case
  // -check the value of a1 in the 1,1 case - it doens't follow the pattern ai=0 for odd i
}

void splitterPrototypeD_2_3_new(double* k, Complex* p, Complex* z)
{
  double  s = sqrt(2)-1;   // like in 2nd order butterworth
  double  t = s;
  //s = 0.5; t = 0.5;
  Complex j = Complex(0, 1);
  p[0] =  j*s;   // p1
  p[1] = -j*s;   // p2
  z[0] = -1;     // q1
  z[1] = -1;     // q2
  z[2] = -t;     // q3
  *k = dcGainNormalizer(z, 3, p, 2); 
}

rsFilterSpecificationBA<double> splitterPrototype_2_3_new()
{
  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a.resize(3);
  ba.b.resize(4);

  // fixed coeffs (by constraints):
  ba.a[0] = 1;              // 0th a coeff is always 1
  ba.a[1] = 0;              // odd a-coeffs must be zero
  ba.b[0] = ba.a[0] / 2.0;  // b0 = 0.5

  // freely choosable coeffs:
  ba.a[2] =  0.5;           // must be > 0 for a complex conjugate pair
  ba.b[1] =  ba.b[0];       // b1 == b0 puts one zero at z=-1 (i think)
  ba.b[3] = 0.7;


  // even b-coeffs must be half of their corresponding a-coeffs:
  ba.b[2] = ba.a[2] / 2.0;  // depends on choice for a[2]
  //ba.b[3] = ba.b[2];        // b3 == b2 puts another zero at z=-1 (i think..or maybe not?)


  // hmm - with such a pole/zero placement, two zeros end up on the two poles


  // i think, we get two zeros at z = -1 if b1==b0 and b3==b2 ..nooo - that doesn't work
  // let H(z) = B(z)/A(z) = b0*((1-q1/z)*(1-q2/z)*(1-q3/z))/((1-p1/z)*(1-p2/z))
  // and fix q1 = q2 = -1 - that gives:
  // B(z) = b0 * ( (1+1/z)*(1+1/z)*(1-q3/z) )
  // multiply out to obtain b0,..,b3 in terms of q3 (our additional zero to be freely placed)
  // let r = 1/z = z^-1 for convenience, then:
  // B(z) = b0 * (1 + (2-q3)*r + (1-2*q3)*r^2 - q3*r^3)
  // so, with b0 = 1/2: b1 = (1-q3/2), b2 = (1/2 - q3), b3 = -q3/2
  // we also need: b2 = a2/2 giving: a2 = 1 - 2*q3



  // hmmm...we have 7 degrees of freedom all in all. constraints uniquely fix
  // a0 = 1, a1 = 0, b0 = a0/2 = 1/2. once a2 is chosen, we must have b2 = a2/2.
  // maybe with the equations above, we can reduce the number of degrees of freedom to 1 and then
  // tweak that remaining variable?





  // choose a2 such that we get a monotonic response and b1,b3 such that two zeros are at z=-1
  // ...maybe the other way around

  return ba;
}

// digital 3-pole/3-zero - doesn't work:
void splitterPrototypeD_3_3(double* k, Complex* p, Complex* z)
{
  Complex j = Complex(0, 1);
  double  a = 0.4;
  p[0] =  j*a;
  p[1] = -j*a;
  p[3] =  0;

  double  b = 0.2;
  z[0] = -1;
  z[1] = -1;
  z[2] = -b;

  *k = dcGainNormalizer(z, 3, p, 3);
}

// digital 4-pole/6-zero (does not work)
void splitterPrototypeD_4_6(double* k, Complex* p, Complex* z)
{
  typedef rsInfiniteImpulseResponseDesigner<double> DSGNR;
  DSGNR dsgnr;
  dsgnr.setApproximationMethod(rsPrototypeDesigner<double>::BUTTERWORTH);
  dsgnr.setPrototypeOrder(4);
  dsgnr.setSampleRate(1);
  dsgnr.setFrequency(0.25);      // halfband filter, 0.5 is Nyquist freq
  dsgnr.getPolesAndZeros(p, z);

  z[4] = -p[0].imag(); // test
  z[5] = -p[2].imag(); // test

  *k = dcGainNormalizer(z, 6, p, 4);
}

rsFilterSpecificationBA<double> complementaryFilter(const rsFilterSpecificationBA<double>& baSpec)
{
  rsFilterSpecificationBA<double> ba = baSpec, r;

  r.sampleRate = ba.sampleRate;
  int Na = (int)ba.a.size()-1;
  int Nb = (int)ba.b.size()-1;

  r.b.resize(std::max(Na,Nb)+1);
  r.a = ba.a;                // denominator is the same

  //if(ba.isDigital()) {
  //  rsReverse(ba.b);
  //  rsReverse(ba.a);
  //}

  rsPolynomial<complex<double>>::subtractPolynomials(&ba.a[0], Na, &ba.b[0], Nb, &r.b[0]);

  //if(r.isDigital())
  //  rsReverse(r.b);

  return r;
} // move to FilterPlotter or rapt rsFilterSpecificationBA

template<class T>
void plotMagnitudesBA(int numFreqs, T lowFreq, T highFreq,
  bool logFreqAxis, bool decibels,
  const std::vector<rsFilterSpecificationBA<T>>& filterSpecs)
{
  FilterPlotter<T> plt;
  for(size_t i = 0; i < filterSpecs.size(); i++)
    plt.addFilterSpecificationBA(filterSpecs[i]);
  plt.plotMagnitude(numFreqs, lowFreq, highFreq, logFreqAxis, decibels);
  //plt.plotPolesAndZeros(); // test
} // maybe move as static function to class FilterPlotter


bool testSplitConditions(const rsFilterSpecificationBA<double>& lpfBA)
{
  // Give the filter-prototype specifications for a lowpass filter, this function checks, if the 
  // filter satisfies the conditions for a perfect reconstruction crossover, assuming the highpass
  // signal is obtained by subtracting the lowpass signal fro the original input - that in itself 
  // ensures perfect reconstruction, but we check here for additional conditions such as symmetry 
  // of the responses.

  bool result = true;
  rsFilterSpecificationBA<double> hpfBA = complementaryFilter(lpfBA);
  rsFilterSpecificationZPK<double> lpfZPK = lpfBA.toZPK();
  rsFilterSpecificationZPK<double> hpfZPK = hpfBA.toZPK();

  // check, if the poles are equal:
  // ...


  // check if zeros are mirrored along the imaginary axis:
  // ...


  // check symmetry: H(z) = G(-z), use some random values for z for that
  int numValues = 100;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(-2.0, +2.0);
  double tol = 1.e-13;
  for(int i = 0; i < numValues; i++)
  {
    Complex z   = Complex(prng.getSample(), prng.getSample());
    Complex Hz  = lpfBA.transferFunctionAt( z);   // H(z)
    Complex Gz  = hpfBA.transferFunctionAt( z);   // G(z)
    Complex Gzm = hpfBA.transferFunctionAt(-z);   // G(-z)
    Complex sum = Hz + Gz;   // Complement: H(z) + G(z) = 1, ensured by complementaryFilter
    Complex dif = Hz - Gzm;  // Symmetry:   H(z) = G(-z) -> H(z)-G(-z) = 0
    result &= abs(1.0-sum) < tol;
    result &= abs(dif)     < tol;
    // for the 2p3z filter, the sum is totally off and Hz+Gzm is zero instead of Hz-Gzm
    // how can the sum be so totally off? is there a bug in complementaryFilter?
    // ok - the sum is correct now - we should not reverse arrays in complementaryFilter
    // but if we don't reverse, the plots are messed up again...wtf?
    // if we do reverse: sum is wrong, if we don't reverse: plots are wrong
    // try to not reverse (the sum *must* be correct) and choose a different additional zero - 
    // maybe it's wrongly placed, after all

    // is our symmetry constraint tto restrictive? maybe it should be |H(z)| = |G(-z)|:
    double absHz  = abs(Hz);
    double absGzm = abs(Gzm);
    // nope - that doesn't work either


    Complex test1 = Hz+Gzm;
    Complex test2 = Hz-Gz;

    //rsAssert(result);
    int dummy = 0;
  }

  return result;
}

void bandSplitHighOrderIIR()
{
  // Experiment to figure out pole/zero placements in the s-domain to obtain a high/low IIR 
  // splitter with perfect reconstruction...

  // ...under construction...

  int N, M;                     // number of poles and zeros
  Complex p[10], z[10];         // poles and zeros
  double k;                     // gain
  double fs;                    // sample-rate
  double inf = RS_INF(double);


  //splitterPrototypeA_1_0(&k, p, z); N = 1; M = 0; fs = inf;  // analog 1-pole
  //splitterPrototypeA_2_0(&k, p, z); N = 2; M = 0; fs = inf;  // analog 2-pole (2nd order Butterworth)
  //splitterPrototypeA_2_1(&k, p, z); N = 2; M = 1; fs = inf;  // analog 2-pole/1-zero - 
  //splitterPrototypeA_2_2(&k, p, z);   N = 2; M = 2; fs = inf;  // analog 2-pole/2-zero
  // ...

  double fsd = 0.5/PI;  // sample-rate for digital filters
  //splitterPrototypeD_1_1(&k, p, z); N = 1; M = 1; fs = fsd;  // digital 1-pole/1-zero - works
  splitterPrototypeD_2_3(&k, p, z); N = 2; M = 3; fs = fsd;  // digital 2-pole/3-zero - old
  splitterPrototypeD_2_3_new(&k, p, z); N = 2; M = 3; fs = fsd;



  //splitterPrototypeD_2_2(&k, p, z); N = 2; M = 2; fs = fsd;  // digital 2-pole/2-zero
  //splitterPrototypeD_3_3(&k, p, z); N = 3; M = 3; fs = fsd;  // test - not yet working
  //splitterPrototypeD_4_6(&k, p, z); N = 4; M = 6; fs = fsd;    // nope - that doesn't work

  // create filter specification objects for lowpass and highpass filter:
  //rsFilterSpecificationZPK<double> lowpassZPK(toVector(z, M), toVector(p, N), k, fs);
  //rsFilterSpecificationBA<double>  lowpassBA  = lowpassZPK.toBA();


  rsFilterSpecificationBA<double>  lowpassBA  = splitterPrototype_2_3_new();
  rsFilterSpecificationBA<double>  highpassBA = complementaryFilter(lowpassBA);

  bool splitConditionsMet = testSplitConditions(lowpassBA);
  // for the 2,3 filter a[2] is > 5 -> unstable filter...i think in converting between zpk/ba, i 
  // need to reverse the polynomial coefficient arrays in case of digital filters...
  // it works

  //FilterPlotter<double> testPlt;
  //testPlt.addFilterSpecificationZPK(lowpassZPK); // these two specs should lead to the same plots
  //testPlt.addFilterSpecificationBA( lowpassBA);  // this is too large!!
  //testPlt.plotMagnitude(1000, 0.0, 0.5, false, false);
  //int dummy = 0;
  // ok, this seems to work


  // plot frequency response:
  FilterPlotter<double> plt;
  //plt.addFilterSpecificationZPK(N, p, M, z, k, fs);
  //plt.addFilterSpecificationZPK(lowpassZPK);
  plt.addFilterSpecificationBA(lowpassBA);
  //plt.addFilterSpecificationBA(highpassBA);
  plt.plotPolesAndZeros();
  plotMagnitudesBA(1000, 0.0, 0.5, false, false, { lowpassBA, highpassBA });
  //plt.plotMagnitude(1000, 0.0, 0.5, false, false);

  //plt.plotMagnitude(1000, 0.01, 100, true, true);  // suitable for analog filters
  //plt.plotMagnitude(1000, 0.0, 2*PI, false, false);    // todo: rescale the freq-axis such that PI maps to 0.5 or 1.0


  // putting additional finite zeros into the s-plane is not a good idea - it makes the final slope
  // shallower - the lowpass should have all of its zeros at infinity

  // ...soooo that means we have to use an allpole lowpass filter and therefore the highpass should
  // have all its zeros at s=0. the only wiggle room is the exact placement of the poles

  // OR: we design a halfband lowpass prototype in the digital domain, leave the zeros at z = -1 
  // and add *additional* zeros. this seems to work for the 2nd order case at least

  // try "contracted Butterworth" - all pole angles are scaled by a factor < 1
  // try to place poles on a ellipse instead of a circle - let the user select a width/height 
  // ratio or eccentricity

  // maybe experiment with multiplicities

  // maybe start with z-plane prototype poles (and zeros), maybe aligned along the imaginary axis
  // and spread in various ways

  // place N poles along the imaginary axis and N zeros at z = -1. then try to place additional 
  // zeros into the z-plane such that we get a nice crossover...maybe we need a GUI for freely
  // placing poles and zeros into the z-plane

}

// end of experiments for perfect reconstruction IIR filters
//-------------------------------------------------------------------------------------------------

void ladderResonanceManipulation()
{
  // We create two ladder filters with the same cutoff frequency, one with and one without 
  // resonance and plot the difference of the states of the two ladder filters when we pass
  // a sawtooth wave into them.

  // parameters:
  int   N   =  2000;      // number of samples
  float fs  = 44100;      // sample rate
  float fc  =  500;       // cutoff frequency
  float res =    0.99f;   // resonance
  float fIn =   80;       // input frequency

  // create and set up the filters:
  typedef RAPT::rsLadderFilter<float, float> LDR;  // for convenience

  LDR withReso;
  withReso.setSampleRate(fs);
  withReso.setCutoff(fc);
  withReso.setResonance(res);
  withReso.setMode(LDR::LP_24);

  LDR noReso;
  noReso.setSampleRate(fs);
  noReso.setCutoff(fc);
  noReso.setResonance(0);
  noReso.setMode(LDR::LP_24);

  // create time axis and input signal;
  vector<float> t(N), x(N);
  createTimeAxis(N, &t[0], fs);
  createWaveform(&x[0], N, 1, fIn, fs, 0.f, false);

  // filter input signal and record the difference between the state variables of teh two filters:
  vector<float> y0(N), y1(N), y2(N), y3(N), y4(N);
  vector<float> z0(N), z1(N), z2(N), z3(N), z4(N);
  vector<float> r0(N), r1(N), r2(N), r3(N), r4(N);
  vector<float> offset(N);
  float y[5], z[5];
  float dummy;
  for(int n = 0; n < N; n++)
  {
    // let the filter compute outputs (we don't actuualy need them):
    dummy = withReso.getSample(x[n]);
    dummy = noReso.getSample(  x[n]);

    // obtain the filter states:
    withReso.getState(y);
    noReso.getState(z);

    // record the states and state-differences:
    y0[n] = y[0];
    y1[n] = y[1];
    y2[n] = y[2];
    y3[n] = y[3];
    y4[n] = y[4];

    z0[n] = z[0];
    z1[n] = z[1];
    z2[n] = z[2];
    z3[n] = z[3];
    z4[n] = z[4];

    r0[n] = y[0] - z[0];
    r1[n] = y[1] - z[1];
    r2[n] = y[2] - z[2];
    r3[n] = y[3] - z[3];
    r4[n] = y[4] - z[4];

    // scale difference-states with appropriate factor to make the sinusoids have the same 
    // amplitude:
    r1[n] *= (float)SQRT2;      // s2^1, s2: sqrt(2)
    r2[n] *= 2;                 // s2^2
    r3[n] *= 2*(float)SQRT2;    // s2^3
    r4[n] *= 4;                 // s2^4

    // According to the model, this value would be zero because r0 and r4 are 180° out of phase. 
    // Any deviation from 0 is a shortcoming of the model and considered an offset to the idealized
    // value:
    offset[n] = r0[n] + r4[n];
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &x[0]);
  //plt.addDataArrays(N, &t[0], &y0[0], &y1[0], &y2[0], &y3[0], &y4[0]);
  //plt.addDataArrays(N, &t[0], &z0[0], &z1[0], &z2[0], &z3[0], &z4[0]);
  //plt.addDataArrays(N, &t[0], &r0[0], &r1[0], &r2[0], &r3[0], &r4[0]);
  plt.addDataArrays(N, &t[0], &r0[0], &r4[0]); // should be 180° out of phase
  plt.addDataArrays(N, &t[0], &offset[0]);
  //plt.addDataArrays(N, &t[0], &y0[0], &z0[0], &r0[0]);
  plt.plot();

  // Observations:
  // The states of the imaginary "pure-resonance" difference filters (resonant - nonResonant)
  // represent (decaying) sinusoids, where each of the successive states has the same sinusoid
  // multiplied by 1/sqrt(2) and phase-shifted by -45 degrees with respect to the stage before. 
  // We assume that at any time instant the resonance waveform has instantaneous amplitude a and 
  // phase p. We must then have:
  //
  // r4 = a             * sin(p)
  // r3 = a * sqrt(2)^1 * sin(p + 1*pi/4)
  // r2 = a * sqrt(2)^2 * sin(p + 2*pi/4) = a * 2 * sin(p + pi/2)
  // r1 = a * sqrt(2)^3 * sin(p + 3*pi/4)
  // r0 = a * sqrt(2)^4 * sin(p + 4*pi/4) = a * 4 * sin(p + pi) 
  //
  // That means, given the states r4,r2 we should be able to figure out the instantaneous amplitude
  // a and phase p. We just divide the state r2 by two and together with r4 we can take it as the 
  // sine and cosine part of our sine from which amplitude and phase can be computed. Having a and
  // p, we can do whatever we want to them (for example sync the phase to an input signal) and then 
  // compute our new states r0,...,r4 using the formulas above. Having done that, we may update the 
  // states of the resonant filter according to y0 = z0 + r0, ..., y4 = z4 + r4.  This should give 
  // us a ladder filter in which we can freely mess with the instantaneous phase of the resonance 
  // without needing to post-process the resonance signal.
  // hmm...well, all of this holds only approximately. Maybe we should estimate the sine parameters
  // not from r2,r4 but from r3,r4 (so we use the 2 most filtered outputs) and maybe we should 
  // compute the difference between the actual states and the idealized states according to the 
  // model and add this difference back after modification of the states, such that when we don't
  // apply any modification, we really don't apply any modification (when recomuting the states
  // according to the formulas)

  // Note: it works only when the compensation gain is applied at the input and when the filter
  // is linear.
}


/** Moving average filter for non-uniformly spaced samples. N: num samples, x: abscissa-values
(mostly time), y: ordinate values, avg: average - the output (must be distinct from y), width: 
length/width/range of the support of the filter, weightFunc: normalized weighting function, 
should have a support in the range -1..+1, this is the filter kernel as continuous function. 
\todo: maybe let the user pass a functor for the weighting function isntead of a function 
pointer (maybe it's possible to write a functor-wrapper for ordinary functions with automatic
type casting such that ordinary functions can still be passed as usual?) */
template<class Tx, class Ty>
void movingAverage(int N, Tx* x, Ty* y, Ty* avg, Tx width, Ty (*weightFunc)(Tx))
{
  Tx w2  = 0.5f * width;                       // half width
  Tx w2r = 1 / w2;                             // reciprocal of half width
  Tx dist;                                     // distance
  int k;                                       // inner loop index
  for(int n = 0; n < N; n++){                  // outer loop over all points
    Ty wgt = weightFunc(0);                    // weight
    Ty sw  = wgt;                              // sum of weights
    Ty swv = wgt * y[n];                       // sum of weighted values
    k = n-1;                                   // immediate left neighbour
    while(k >= 0 && (dist = x[n]-x[k]) <= w2){ // left side loop
      wgt  = weightFunc(dist * w2r);           // compute weight for distance
      sw  += wgt;                              // accumulate weight sum
      swv += wgt * y[k];                       // accumulate weighted values
      k--; }                                   // jump to next neighbour
    k = n+1;
    while(k < N  && (dist = x[k]-x[n]) <= w2){ // right side loop
      wgt  = weightFunc(dist * w2r);
      sw  += wgt;
      swv += wgt * y[k];
      k++; }
    avg[n] = swv / sw; }
}
// Weighting functions for nonuniform MA. Maybe implement more, see here:
// https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
// actually, for use in movingAverage, the checks against > 1 are superfluous bcs the loops there
// ensure already that x <= 1. maybe we can wrap all this stuff into a class 
// NonUniformMovingAverage. We could have functions like setNominalWidth (the width used for 
// uniform weighting), setWeightingFunction(T (*func)(T), T area) and for functions with unknown
// area just setWeightingFunction(T (*func)(T)) which finds the area by numerical integration
// and then calls setWeightingFunction(T (*func)(T), T area) with the numerically computed area
// ...finally an opportunity to implement and use numerical integration algorithms :-)
// other functions to try: (1-x^a)/(1+x^a), a = 1, 1/2, 1/3, 2, 
float uniform(float x)       // rename to uniform, maybe scale by 0.5
{
  if(fabs(x) > 1)
    return 0;
  return 1;
}
float triangular(float x)      // half-area: 0.5 (integral from 0 to 1)
{
  x = fabs(x);
  if(x > 1)
    return 0;
  return 1 - x;
}
float rationalTent(float x)   // half-area: log(4)-1
{
  x = fabs(x);
  if(x > 1)
    return 0;
  return (1-x) / (1+x);
}
float parabolic(float x)      // half-area: 2/3  (Epanechikov)
{
  if(fabs(x) > 1)
    return 0;
  return 1 - x*x;
}
float cubicBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::cubic(fabs(x));
}
float quinticBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::quintic(fabs(x));
}
float hepticBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::heptic(fabs(x));
}
void nonUniformMovingAverage()
{
  // user parameters:
  static const int N = 200;  // number of samples
  float minDist = 1.f;       // minimum distance between successive samples
  float maxDist = 10.f;      // maximum ...
  float minY    = -10.f;     // minimum y value
  float maxY    = +10.f;     // maximum
  float width   =  30.f;     // filter support width for box, others are scaled by their area
  int   seed    = 0;

  // create input data:
  float x[N], y[N];
  float t = 0.f; // time
  float dt;      // time delta between samples
  x[0] = t;
  y[0] = (float)round(RAPT::rsRandomUniform(minY, maxY, seed)); 
  for(int n = 1; n < N; n++){
    //dt = (float)round(rsRandomUniform(minDist, maxDist));
    dt = (float)RAPT::rsRandomUniform(minDist, maxDist);
    t += dt;
    x[n] = t;
    y[n] = (float)RAPT::rsRandomUniform(minY, maxY); 
    //y[n] = (float)round(rsRandomUniform(minY, maxY)); 
  }

  // create filtered versions:
  float a = 1.f / float(log(4.f)-1); // reciprocal of area under rationalTent weighting function
                                     // ...the definite integral from 0 to 1
  float y0[N], y1[N], y2[N], y3[N], y4[N],
    yR[N],  // rational
    yP[N];  // parabolic

  // 0,1,2,3,4 is the smoothness of the weight function

  movingAverage(N, x, y, y0,  width,      uniform);
  movingAverage(N, x, y, y1,  width*2,    triangular);
  movingAverage(N, x, y, yR,  width*a,    rationalTent); 
  movingAverage(N, x, y, yP,  width*1.5f, parabolic);
  movingAverage(N, x, y, y2,  width*2,    cubicBell);
  movingAverage(N, x, y, y3,  width*2,    quinticBell);
  movingAverage(N, x, y, y4,  width*2,    hepticBell);

  // For all weighting functions other than the box, we scale the support width by factor to
  // equal to the reciprocal of the half of the area under the function to make the plots more 
  // easily comparable. 

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, x, y, y0, y1, y2, y3, y4);
  //plt.addDataArrays(N, x, y0, y1, yR);
  //plt.addDataArrays(N, x, y0, y1, y2, y3, y4, yR);
  //plt.addDataArrays(N, x, y0, y2, y4);
  plt.addDataArrays(N, x, y, y0, y1, yR, yP);
  plt.plot();

  // Observations: The rationalTent weighting seems best in terms of smoothing. It shows the 
  // smallest deviations from the average, i.e. the curve is the most "inside" of all (stays
  // closer to the mean than the others). 
  // ToDo: maybe try to find even better weighting functions. Maybe the smoothness is related to 
  // the area? smaller area -> more width widening -> better smoothing? maybe because the jitter
  // gets averaged out better?
}

void smoothingFilterOrders()
{
  // We plot the step responses of the rsSmoothingFilter for various orders.

  static const int numOrders = 4;   // number of filters with different orders
  bool expSpacing = false;          // if true, orders are 1,2,4,8,.. else 1,2,3,4,..
  int orderIncrement = 2;           // or 1,3,5,.. or 1,4,7,...
  static const int N = 300;         // number of samples
  float fs  = 1.0f;                  // sample rate
  float tau = 101.0f;                  // time constant

  // create and set up the smoother:
  rsSmoothingFilterFF smoother;
  //rsSmoothingFilterDD smoother;
  smoother.setTimeConstantAndSampleRate(tau, fs);
  //smoother.setNumSamplesToReachHalf(tau*fs);
  //smoother.setShape(rsSmoothingFilterFF::FAST_ATTACK);
  smoother.setShapeParameter(0.5f);

  // compute step responses:
  int order = 1;
  float y[numOrders][N];
  for(int i = 0; i < numOrders; i++)
  {
    smoother.setOrder(order);
    if(expSpacing) // update order for next iteration
      order *= 2;
    else
      order += orderIncrement;

    smoother.reset();
    y[i][0] = smoother.getSample(1.f);
    for(int n = 1; n < N; n++)
    {
      //y[i][n] = smoother.getSample(0.f);   // impulse response
      y[i][n] = smoother.getSample(1.f); // step response
    }
  }

  // plot:
  float t[N];
  createTimeAxis(N, t, fs);
  GNUPlotter plt;
  for(int i = 0; i < numOrders; i++)
    plt.addDataArrays(N, t, y[i]); 
  plt.plot();

  // Observations:
  // The step responses of the different orders are comparable in terms of overall transition time.
  // However, they do not meet in a common point. From a user's perspective, it would perhaps be
  // most intuitive, if he could just set up the time-instant, where the step-response goes through
  // 0.5. No matter what the order is, the time instant where it passes through 0.5 should remain 
  // fixed. To achieve that, i think, we need a closed form expression for the step-response and 
  // then set that expression equal to 0.5 and solve for the time-constant. If it's too hard to 
  // find such a formula, we could also create a table by just reading off the time-instants where
  // the various step responses go through 0.5 and then use that value as divider.

  // ToDo:
  // Measure the sample instants where the the step responses pass through 0.5 for some normalized
  // filter setting and create tables from that. These tables will eventually be used for scaling. 
  // It should be a 2D table, 1st dimension is the order (say 1..32) and 2nd dimension is the value
  // of the shape parameter (maybe from 0 to 4 in 32 steps?) the 2nd dimenstion can use linear
  // interpolation, the 1st dimension doesn't need that because inputs are integers anyway.
  // The values of the table should be overall scalers for all the time constants.


  // try higher orders - like 100 - see what kind of shaped is approached for order -> inf
}

void smoothingFilterTransitionTimes()
{
  // We plot the step response for smoothing filters of a given order (and shape/asymmetry) for
  // different transition time constants. What we should see is time-stretched/compressed versions
  // of the same response. But the nonworking tuning tables suggest that somthing might be wrong
  // with that assumption, so we check.

  // user parameters:
  static const int N = 400;           // number of samples
  int order = 5;                      // order of the filters
  float asymmetry = 1.f;              // asymmetry parameter
  static const int numFilters = 5;    // number of transition times
  float transitionTimes[numFilters] = { 51.f, 101.f, 151.f, 200.1f, 251.f }; // transition times in samples


  // create and set up the smoother:                                                             
  rsSmoothingFilterFF smoother;
  smoother.setShapeParameter(asymmetry);
  smoother.setOrder(order);

  // create the step responses:
  float y[numFilters][N];
  for(int i = 0; i < numFilters; i++)
  {
    smoother.setNumSamplesToReachHalf(transitionTimes[i]);
    smoother.reset();
    for(int n = 0; n < N; n++)
      y[i][n] = smoother.getSample(1.f);
  }

  // plot:
  GNUPlotter plt;
  for(int i = 0; i < numFilters; i++)
    plt.addDataArrays(N, y[i]); 
  plt.plot();

  // Observations:
  // The intuitive assumption that scaling all the time constants by the same factor has the effect
  // of stretching/compressing the step-response by that amount turns out to be false. It seems to
  // be true for 1st order filters but for higher orders, especially when there's asymmetry, the 
  // actually observed step repsonses deviate from that simple rule more and more. 
  // Why is that? Is it because the design formula uses impulse-invariant transform and should use
  // step-invariant? ...or maybe even bilinear?

  // see here: http://web.cecs.pdx.edu/~tymerski/ece452/Chapter4.pdf
  // http://homes.esat.kuleuven.be/~maapc/Sofia/slides_chapter8.pdf
  // http://dspcan.homestead.com/files/IIRFilt/zfiltsii.htm
}




template<class T>
void reflectRoots(complex<T>* roots, int N) // rename to mirrorFirstHalf, make generic, add to library
{
  for(int i = 0; i <= N/2; i++)
    roots[N-1-i] = conj(roots[i]);
}
template<class T>
void removeInfiniteValues(vector<complex<T>>& z)
{
  for(size_t i = 0; i < z.size(); i++) {
    if(RAPT::isInfinite(z[i])) {   // as rs prefix
      z.erase(z.begin() + i);
      i--; }}
}
template<class T>
rsFilterSpecificationZPK<T> getFilterSpecificationZPK(RAPT::rsPrototypeDesigner<T>& pd)
{
  int nz = pd.getNumFiniteZeros();
  int np = pd.getNumFinitePoles();
  vector<complex<T>> z(np);          // np is correct, we may get infinite zeros which will be..
  vector<complex<T>> p(np);          // ..removed later
  pd.getPolesAndZeros(&p[0], &z[0]); // returns only the non-redundant upper halfplane poles/zeros
  reflectRoots(&p[0], np);           // create full pole/zero set by reflection
  reflectRoots(&z[0], np);
  removeInfiniteValues(z);
  rsFilterSpecificationZPK<T> spec;
  spec.p = p;
  spec.z = z;
  spec.sampleRate = RS_INF(T);
  //spec.gain  = pd.getGain();
  return spec;
}
void prototypeDesign()
{
  //rsEngineersFilterFF ef;

  /*
  // something is wrong with odd order bessel filters (exept for order = 1), looks like it never 
  // worked in rapt - ive gone back to the commit of 21.11.2017, 00:06:28 with message:
  // implemented FilterPlotter<T>::plotMagnitude
  // code for debugging :
  // compare float/double versions:
  rsPrototypeDesignerD pdD;
  complex<double> polesD[5], zerosD[5];
  pdD.setApproximationMethod(pdD.BESSEL);
  pdD.setOrder(3);
  pdD.getPolesAndZeros(polesD, zerosD);

  rsPrototypeDesignerF pdF;
  complex<float> polesF[5], zerosF[5];
  pdF.setApproximationMethod(pdF.BESSEL);
  pdF.setOrder(3);
  pdF.getPolesAndZeros(polesF, zerosF);

  int dummy = 0;
  // ok, yes - that's the problem: the double precision poles are correct (they match with what the 
  // rosic version computes) but the single precision float poles are wrong. It's probably 
  // something about the root finder
  // in rsPrototypeDesigner<T>::makeBesselLowShelv the pickNonRedundantPolesAndZeros fills the p 
  // member array differently for double and float versions, the pTmp array is the same in both
  // cases

  // in rsZeroNegligibleImaginaryParts the "real" pole has some small negative imaginary part which 
  // doesn't get zeroed out and in the subsequent call to rsOnlyUpperHalfPlane, the real pole
  // gets thrown away - we need a different threshold for float than we use for double
  // ...ok - seems to be fixed.
  */

  typedef float Real;
  typedef RAPT::rsPrototypeDesigner<Real> PD;


  // min and max filter order to plot:
  int minOrder = 1;
  int maxOrder = 10;

  // create and set up prototype designer:
  PD pd;
  pd.setApproximationMethod(PD::GAUSSIAN);
  //pd.setApproximationMethod(PD::BESSEL);
  //pd.setApproximationMethod(PD::BUTTERWORTH);
  //pd.setApproximationMethod(PD::PAPOULIS);
  //pd.setApproximationMethod(PD::HALPERN);
  //pd.setApproximationMethod(PD::CHEBYCHEV);
  //pd.setApproximationMethod(PD::INVERSE_CHEBYCHEV);
  //pd.setApproximationMethod(PD::ELLIPTIC);

  //pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE); // comment fo lowpass
  pd.setGain(+6.02f); // needs to be nonzero for plots

  pd.setPassbandRipple(1); 
  pd.setStopbandRejection(20);

  // create plotter, add filter specs for the desired orders to it and plot:
  FilterPlotter<Real> plt;
  for(int i = minOrder; i <= maxOrder; i++)
  {
    pd.setOrder(i);
    rsFilterSpecificationZPK<Real> spec = getFilterSpecificationZPK(pd);
    plt.addFilterSpecificationZPK(spec);
  }
  //plt.plotPolesAndZeros(600);
  plt.plotMagnitude(1000, 0, 3, false, false);

  // issues:

  // it seems, even order elliptic prototypes have an overall gain equal to the reciprocal of the 
  // linear stopband rejection (passband ripple seems to have no effect), odd order elliptics have
  // a different gain and it seems to depend on the ripple (might be computed by evaluating DC 
  // gain). Papoulis design has also a wrong DC gain -> add overall gain to the prototype designer
  // in EngineersFilter, it's not a problem, because i have a gain normalization step at the very
  // end of the design pipeline - but it would be desirable to have correct gains at all stages
  // of the design process rather than renormalizing at the end (which is sort of dirty)

  // todo: test the gaussian filter design
}


// this is for testing the rewrite attempt - which i then decided to give up:

// helper function to convert form raw arrays of poles and zeros to the FilterSpecificationZPK 
// structure used by the plotter
template<class T>
rsFilterSpecificationZPK<T> analogPrototypeSpecZPK(
  int N, std::complex<T>* za, std::complex<T>* pa, T k)
{
  vector<complex<T>> z = RAPT::toVector(za, N);
  vector<complex<T>> p = RAPT::toVector(pa, N);
  reflectRoots(&p[0], N);           // create full pole/zero set by reflection
  reflectRoots(&z[0], N);
  removeInfiniteValues(z);
  rsFilterSpecificationZPK<T> spec;
  spec.p = p;
  spec.z = z;
  spec.sampleRate = RS_INF(T); // indicates analog (s-plane) poles and zeros
  spec.k = k;
  return spec;
}

void poleZeroPrototype()
{
  // tests the new implementation of analog filter prototype design

  typedef float Real;
  typedef PoleZeroPrototype<Real> PZP;
  static const int minOrder = 4;   // minimum filter order to plot
  static const int maxOrder = 4;  // maximum filter order to plot
  Real G0 = 0.f;    // use G0=0, G=1 for lowpass and G0=2, G=8 for low-shelf
  Real G  = 1.f;     

  PZP pzp;
  pzp.setApproximationMethod(PZP::BUTTERWORTH);

  // maybe the new version should use optimized versions of halpernT2, papoulisL2 from 
  // Prototypes.h/cpp (avoids dynamic memory allocation)

  Real k; 
  std::complex<Real> p[maxOrder], z[maxOrder]; // (maxOrder+1)/2 should be enough
  FilterPlotter<Real> plt;
  for(int n = minOrder; n <= maxOrder; n++)
  {
    pzp.setOrder(n);
    pzp.getPolesZerosAndGain(p, z, &k);
    rsFilterSpecificationZPK<Real> spec = analogPrototypeSpecZPK(n, z, p, k);
    plt.addFilterSpecificationZPK(spec);
  }
  plt.plotPolesAndZeros(600);
  //plt.plotMagnitude(500, 0, 3, false, false);

  // something is wrong...

  int dummy = 0;
}

