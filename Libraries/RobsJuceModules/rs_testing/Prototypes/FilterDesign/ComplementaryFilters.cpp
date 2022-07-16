//typedef std::complex<double> Complex;


rsFilterSpecificationBA<double> weightedSumWithOne(
  const rsFilterSpecificationBA<double>& baSpec, double w1, double wH)
{
  rsFilterSpecificationBA<double> ba = baSpec, r;
  r.sampleRate = ba.sampleRate;
  int Na = (int)ba.a.size()-1;
  int Nb = (int)ba.b.size()-1;
  r.b.resize(rsMax(Na,Nb)+1);
  r.a = ba.a;                // denominator is the same
  rsPolynomial<std::complex<double>>::weightedSum(&ba.a[0], Na, w1, &ba.b[0], Nb, wH, &r.b[0]);
  return r;
}
// Verify this formula!


rsFilterSpecificationBA<double> complementaryFilter(const rsFilterSpecificationBA<double>& baSpec)
{
  return weightedSumWithOne(baSpec, 1, -1);
} 
// -Move to FilterPlotter or rapt rsFilterSpecificationBA
// -Maybe implement an even more general function that adds two arbitrary transfer functions. But 
//  actually it would make more sense to just use rsRationalFunction for that

bool isComplementary(const rsFilterSpecificationBA<double>& lpfBA)
{
  // Given the filter-prototype specifications for a lowpass filter, this function checks, if the 
  // filter satisfies the conditions for a perfect reconstruction crossover, assuming the highpass
  // signal is obtained by subtracting the lowpass signal from the original input. That in itself 
  // ensures perfect reconstruction, but we check here for additional conditions such as symmetry 
  // of the responses.

  using Complex = std::complex<double>; 

  bool result = true;
  rsFilterSpecificationBA<double> hpfBA = complementaryFilter(lpfBA);
  rsFilterSpecificationZPK<double> lpfZPK = lpfBA.toZPK();
  rsFilterSpecificationZPK<double> hpfZPK = hpfBA.toZPK();

  // Check, if the poles are equal:
  // ...


  // Check if zeros are mirrored along the imaginary axis:
  // ...


  // Check symmetry H(z) = G(-z) numerically using some random values for z:
  int numValues = 100;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(-2.0, +2.0);
  double tol = 1.e-13;
  for(int i = 0; i < numValues; i++)
  {
    Complex z   = std::complex<double>(prng.getSample(), prng.getSample());
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


    //rsAssert(result);
    //int dummy = 0;
  }

  return result;
}

template<class T>
void plotMagnitudesBA(int numFreqs, T lowFreq, T highFreq, bool logFreqAxis, bool decibels,
  const std::vector<rsFilterSpecificationBA<T>>& filterSpecs)
{
  FilterPlotter<T> plt;
  for(size_t i = 0; i < filterSpecs.size(); i++)
    plt.addFilterSpecificationBA(filterSpecs[i]);
  plt.plotMagnitude(numFreqs, lowFreq, highFreq, logFreqAxis, decibels);
  //plt.plotPolesAndZeros(); // test
} // maybe move as static function to class FilterPlotter
// rename to analyzeFiltersBA, put more tests here...or: write a function that just takes the BA
// array
// move to plotting functions

template void plotMagnitudesBA(
  int numFreqs, double lowFreq, double highFreq, 
  bool logFreqAxis, bool decibels,
  const std::vector<rsFilterSpecificationBA<double>>& filterSpecs);



std::vector<double> impulseResponse(const RAPT::rsFilterSpecificationBA<double>& specBA, int N)
{
  typedef std::vector<double> Vec;

  // to make the filtering routine work, we have to convert the complex coeffs to reals:
  Vec a, b;
  a.resize(specBA.a.size());
  b.resize(specBA.b.size());
  for(size_t i = 0; i < a.size(); i++) a[i] = specBA.a[i].real();
  for(size_t i = 0; i < b.size(); i++) b[i] = specBA.b[i].real();

  // comput impulse response by filtering a unit impulse:
  Vec x = createImpulse(N);
  Vec y = createSilence(N);  // direct filter output
  RAPT::rsArrayTools::filter(&x[0], N, &y[0], N, &b[0], (int)b.size()-1, &a[0], (int)a.size()-1);
  return y;
} // move somewhere else

bool analyzeComplementaryFilter(const RAPT::rsFilterSpecificationBA<double>& specBA)
{
  bool result = isComplementary(specBA);

  // get the b,a specification of the complementary filter, i.e. the filter that is obtained by
  // subtracting the given filter's output from the original input:
  RAPT::rsFilterSpecificationBA<double> compBA = complementaryFilter(specBA);


  // obtain impulse-response and write to wave-file:
  typedef std::vector<double> Vec;
  int N = 1024;
  Vec x  = createImpulse(N);
  Vec yd = impulseResponse(specBA, N);  // direct filter output
  Vec yc = x-yd;                        // complementary filter output
  //rosic::writeToMonoWaveFile("ComplementaryFilterOut1.wav", &yd[0], N, 44100, 16);
  //rosic::writeToMonoWaveFile("ComplementaryFilterOut2.wav", &yc[0], N, 44100, 16);
  // for the files, we should really choose a lower cutoff frequency to actually see something


  // make some plots:
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(specBA);
  plt.addFilterSpecificationBA(compBA);
  //plt.usePiAxisTics();    // should make the axis tics multiples of pi


  //plt.initialize();  // needed?
  plt.addCommand("set xtics ('0' 0, 'sr/4' pi/2, 'sr/2' pi)");
  plt.plotMagnitude(1000, 0.0, PI, false, false); // todo: write pi/2, pi, 2pi etc on the w-axis

  //plt.initialize();  // needed?
  //plt.plotPolesAndZeros();

  return result;

  // Notes:
  // -Something seems to be wrong about the frequency axis. The responses look repeated several
  //  times. Actually, it's a bit more than 3 times....like...pi times maybe? Maybe I have confused
  //  normalized freq from 0..0.5 with normalized *radian* freq from 0..pi? ...seems fixed now
  // -Plotting the poles and zeros first and then the magnitude produces a messed up plot. -> Debug
  //  FilterPlotter (or maybe GNUPlotCPP). It may have something to do with a not properly 
  //  (re)initialized commandFile and/or dataFile
}

bool analyzeComplementaryAllpass(const RAPT::rsFilterSpecificationBA<double>& apfBA)
{
  // Obtain the derived lowpass and highpass responses by forming half the sum and difference of
  // allpass transfer function and unity:
  RAPT::rsFilterSpecificationBA<double> lpfBA, hpfBA;
  lpfBA = weightedSumWithOne(apfBA, 0.5,  0.5);
  hpfBA = weightedSumWithOne(apfBA, 0.5, -0.5);

  // Plot lowpass and highpass responses:
  FilterPlotter<double> plt;
  plt.addFilterSpecificationBA(lpfBA);
  plt.addFilterSpecificationBA(hpfBA);
  plt.addCommand("set xtics ('0' 0, 'sr/4' pi/2, 'sr/2' pi)");
  plt.plotMagnitude(1000, 0.0, PI, false, false); 
  // This shows a bandpass and bandreject response

  //analyzeComplementaryFilter(lpfBA);
  // This shows the responses, as it should

  return isComplementary(lpfBA) && isComplementary(hpfBA);
}

//-------------------------------------------------------------------------------------------------
// Rules for designing filters H(z) = B(z)/A(z) that satisfy the constraint that subtracting the 
// filter output from the input G(z) = 1-H(z) = (A(z)-B(z))/A(z) = C(z)/A(z) leads to a frequency 
// response G(z) that is a mirror-image of the response of the original filter H(z), i.e.
// G(z) = H(-z):
// -(1) odd a-coeffs are zero
// -(2) even b-coeffs are half of corresponding a-coeffs (even a-coeffs are twice the b-coeffs)
// -poles are symmetrical with respect to imaginary axis (check this)
//  ...(implying A(z)=A(-z) - right?) -> G(z) = H(-z) reduces to C(z) = B(-z)
// -from normalization, we have a0 = 1, so with even-b rule, b0 = 0.5 always
// -when the number of poles is even, we may have either the same number of poles and zeros or one 
//  zero more than poles - otherwise constraint (2) is violated
// -an odd number of poles is not possible (but what about 1st order?) because the odd a-coeffs are
//  constrained to be zero
// After a prototype halfband filter is designed, it can be tuned to any frequency by applying the 
// Constantinides frequency warping formulas to the poles and zeros.


// the 1-pole,1-zero case is equivalent to a first order Butterworth filter via bilinear transform
RAPT::rsFilterSpecificationBA<double> complementaryLowpass1p1z()
{
  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a.resize(2);
  ba.b.resize(2);

  ba.a[0] = 1;              // 0th a-coeff is always 1 (normalization)
  ba.a[1] = 0;              // odd a-coeffs must be zero
  ba.b[0] = ba.a[0] / 2.0;  // = 0.5, even b-coeffs are half of corresponding a-coeffs
  ba.b[1] = 0.5;            // odd b-coeffs are determined by zeros, 
                            // b1 = b0 = 0.5 places the zero at z=-1

  return ba;
}

RAPT::rsFilterSpecificationBA<double> complementaryLowpass2p2z()
{
  // We start with the numerator B(z) of the transfer function H(z) = B(z)/A(z) and write it as:
  //
  //   B(z) = b0*(1-q1/z)*(1-q2/z)      defining r = 1/z = z^-1, this becomes:
  //   B(z) = b0*(1-q1*r)*(1-q2*r)      we fix the first zero q1 at z = -1, i.e. q1 = -1, so:
  //   B(z) = b0*(1+r)*(1-q2*r)         multiplying out and using r=1/z:
  //   B(z) = b0 + b0*(1-q2)/z - b0*q2/z^2
  //
  // With our constraints, we get:
  //
  //   b0 = 0.5, b1 = 0.5*(1-q2), b2 = -0.5*q2, a0 = 1, a1 = 0, a2 = -q2
  //
  // This leaves us one tweakable zero q2 to adjust the filter response to taste.

  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a.resize(3);
  ba.b.resize(3);

  //double q2 = -0.101;
  double q2 = -0.1010454;
  // Our tweakable. q2 = -0.101 seems to be (near) the value where there's no overshoot. With
  // -0.102, we can seen some slight overshoot/resonant peak when zooming in.
  // ToDo: numerically optimize that value. So far, I have only eyeballed it. The goal is to get
  // to the point just before an overshoot peak forms. Maybe find the peak of the freq-response and
  // then adjust q2 to the maximum value possible such that the peak is still at DC.

  ba.a[0] = 1;
  ba.b[0] = 0.5;

  ba.a[1] = 0;
  ba.b[1] = 0.5*(1-q2);

  ba.a[2] = -q2;
  ba.b[2] = -0.5*q2;

  return ba;
}
// ToDo: maybe factor out a function that designs the filter directly as zpk specification

rsFilterSpecificationBA<double> complementaryLowpass2p3z()
{
  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a.resize(3);
  ba.b.resize(4);

  // fixed coeffs (by constraints):
  ba.a[0] = 1;              // 0th a-coeff is always 1
  ba.a[1] = 0;              // odd a-coeffs must be zero
  ba.b[0] = ba.a[0] / 2.0;  // b0 = 0.5

  // Let H(z) = B(z)/A(z) = b0*((1-q1/z)*(1-q2/z)*(1-q3/z))/((1-p1/z)*(1-p2/z))
  // and fix q1 = q2 = -1. That gives: B(z) = b0 * ( (1+1/z)*(1+1/z)*(1-q3/z) )
  // Multiply out to obtain b0,..,b3 in terms of q3 (our additional zero to be freely placed)
  // Let r = 1/z = z^-1 for convenience, then:
  // B(z) = b0 * (1 + (2-q3)*r + (1-2*q3)*r^2 - q3*r^3)
  // So, with b0 = 1/2: b1 = (1-q3/2), b2 = (1/2 - q3), b3 = -q3/2
  // We also need: b2 = a2/2 giving: a2 = 1 - 2*q3

  // ..ok - let's try it:
  double q3 = 0.2;  // 0.5 leads to a2=0 - does this mean, there are no poles?
  ba.a[2] = 1 - 2*q3;
  ba.b[1] = 1 - q3/2;
  ba.b[2] = ba.a[2] / 2.0; // == 1/2 - q3
  ba.b[3] = -q3/2;
  // ok... the pole-zero-plot looks good - also for the complementary filter - but the magnitude 
  // response plot is messed up...and then we need to optimize q3

  // hmmm...we have 7 degrees of freedom all in all. constraints uniquely fix
  // a0 = 1, a1 = 0, b0 = a0/2 = 1/2. once a2 is chosen, we must have b2 = a2/2.
  // maybe with the equations above, we can reduce the number of degrees of freedom to 1 and then
  // tweak that remaining variable?

  // it seems that however we choose q3, the response shape is nonmonotonic
  // ...maybe we need two additional zeros to get a good response...maybe first try to use a pair
  // of complex conjugate zeros (giving only one adjustable parameter - the imaginray part of both 
  // zeros)
  // ...but maybe we should first verify that our plots are actually correct by takeing an FFT of
  // an impulse response and compare to our plot - maybe try this also for the old filter
  // with this sqrt(2)-1 stuff...
  // at least, our symmetry constraints seem to work

  // choose a2 such that we get a monotonic response and b1,b3 such that two zeros are at z=-1
  // ...maybe the other way around

  // maybe we should generally consider the filter's zeros as out tweakable degrees of freedom and
  // let the poles fall whereever they may due to the constraints

  return ba;
}

RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p4z1t()
{
  // sage:
  // var("r b0 q3 q4")
  // B(r) = b0*(1+r)*(1+r)*(1+r)*(1-q4*r)
  // (expand(B)).collect(r)
  // gives:
  // -b0*q4*r^4 - (3*b0*q4 - b0)*r^3 - 3*(b0*q4 - b0)*r^2 - (b0*q4 - 3*b0)*r + b0

  typedef std::complex<double> Complex;
  Complex q4, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4;
  //q4 = 1.5;
  //q4 = -0.5;
  q4 = -1.0;
  //q4 = -0.3;  // unstable
  // q4 < 0: unstable
  // gives nice magnitude responses - but is unstable :'-(
  // can we reflect the poles or will this destroy meeting the constraints? probably..
  // ...maybe try 5 zeros...

  b0 =  0.5;    // always, because a0 = 1 due to normalization and ai = 2*bi for even i
  b1 = -b0*(q4 - 3.0);
  b2 = -3.0*b0*(q4 - 1.0);
  b3 = -b0*(3.0*q4 - 1.0);
  b4 = -b0*q4;

  // a-coeffs are determined by b-coeffs (maybe factor that out):
  a0 =  1;
  a1 =  0;      // as always
  a2 =  2.0*b2; // as always
  a3 =  0;
  a4 =  2.0*b4;

  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a = { a0, a1, a2, a3, a4 };
  ba.b = { b0, b1, b2, b3, b4 };
  return ba;
}

RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p4z()
{
  // again, let r = 1/z = z^-1, to get the expressions for b0,b1,b2,b3,b4, this time, we use sage:
  // http://sagecell.sagemath.org/?z=eJwrSyzSUCpSSDJQKDRWKDRR0uTlctIo0lSwBQppaRhqF2nCSd1CYy0owwTI4OXSSK0oSMxL0XDS1NRLzs_JSU0uAWoFADI3FM4=&lang=sage
  // var("r b0 q3 q4")
  // B(r) = b0*(1+r)*(1+r)*(1-q3*r)*(1-q4*r)
  // (expand(B)).collect(r)
  // gives:
  // b0*q3*q4*r^4 + (2*b0*q3*q4 - b0*q3 - b0*q4)*r^3 + (b0*q3*q4 - 2*b0*q3 - 2*b0*q4 + b0)*r^2 
  // - (b0*q3 + b0*q4 - 2*b0)*r + b0

  typedef std::complex<double> Complex;
  Complex q3, q4, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4;

  // tweakable:
  //q3 = Complex(0.3, 0.2); q4 = conj(q3);
  q3 = Complex(0.5, 0.25); q4 = conj(q3);
  //q3 =  0.8; q4 = -q3;

  // coeffs can be computed by constraint equations, once q3,q4 are chosen:

  b0 =  0.5;    // as always
  b1 = -b0*(q3 + q4 - 2.0);
  b2 =  b0*(q3*q4 - 2.0*q3 - 2.0*q4 + 1.0);
  b3 =  b0*(2.0*q3*q4 - q3 - q4);
  b4 =  b0*q3*q4;

  a0 =  1;      // as always
  a1 =  0;      // as always
  a2 =  2.0*b2; // as always
  a3 =  0;
  a4 =  2.0*b4;

  // tweak q3,q4 and maybe try to put 3 zeros at z = -1 leaving only q4 tweakable

  // b1 =
  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a = { a0, a1, a2, a3, a4 };
  ba.b = { b0, b1, b2, b3, b4 };
  return ba;
}

RAPT::rsFilterSpecificationBA<double> complementaryLowpass4p5z()
{
  // 4 zeros at z = -1, q5 tweakable
  // var("r b0 q5")
  // B(r) = b0*(1+r)*(1+r)*(1+r)*(1+r)*(1-q5*r)
  // (expand(B)).collect(r)
  // -b0*q5*r^5 - (4*b0*q5 - b0)*r^4 - 2*(3*b0*q5 - 2*b0)*r^3 - 2*(2*b0*q5 - 3*b0)*r^2 
  // - (b0*q5 - 4*b0)*r + b0


  typedef std::complex<double> Complex;
  Complex q5, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4, b5;

  q5 = -0.5;

  // 
  b0 = 0.5;
  b1 = -b0*(q5 - 4.0);
  b2 = -2.0*b0*(2.0*q5 - 3.0);
  b3 = -2.0*b0*(3.0*q5 - 2.0);
  b4 = -b0*(4.0*q5 - 1.0);
  b5 = -b0*q5;

  // same as always - maybe make a function setFeedbackCoeffsForSymmetry:
  a0 =  1;
  a1 =  0;
  a2 =  2.0*b2;
  a3 =  0;
  a4 =  2.0*b4;

  // unstable - maybe try to fix only 3 zeros at z = -1

  // todo: maybe make an explorer plugin where i can manually place all zeros - select numPoles,
  // select, if there should be an additional zero and then move them around in the z-plane

  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a = { a0, a1, a2, a3, a4 };
  ba.b = { b0, b1, b2, b3, b4, b5 };
  return ba;
}

RAPT::rsFilterSpecificationBA<double> complementaryLowpass3p3z()
{
  // preliminary:
  rsFilterSpecificationBA<double> ba = complementaryLowpass2p3z();
  ba.a.resize(4); 
  return ba;
}



RAPT::rsFilterSpecificationBA<double> complementaryAllpass1p1z()
{
  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;

  ba.a.resize(2);
  ba.a[0] = 1;
  ba.a[1] = 0;

  ba.b = ba.a;  
  rsReverse(ba.b);

  return ba;
}


RAPT::rsFilterSpecificationBA<double> complementaryAllpass2p2z()
{
  // We replicate the poles (i.e. the feedback coeffs) of the complementaryLowpass2p2z design and
  // obtain the feedforward coeffs by array reversal of the feedback coeffs.

  rsFilterSpecificationBA<double> ba;
  ba.sampleRate = 1;
  ba.a.resize(3);

  double q2 = -0.101; // Our tweakable. See comments in the 2p2z lowpass design.

  ba.a[0] =  1;
  ba.a[1] =  0;
  ba.a[2] = -q2;

  ba.b = ba.a;
  rsReverse(ba.b);

  return ba;
}
// Hmmm...the responses of the derived lowpass and highpass do not look right. Maybe the idea that
// we can just use a lowpass design, turn it into an allpass by using the same poles and invert 
// them to obtain the zeros and then getting the lowpass back by doing (1 + H(z)) / 2 is 
// fundamentally flawed? I kinda assumed that it may work without much justification other than 
// having it seen working for the 1st order case. The responses we get here look more like the
// bandpass/bandstop responses that we would get by starting from a 1st order prototype.



// Maybe requiring C(z) = B(-z) is too restrictive - we assumed that A(z) = A(-z) in the step going 
// from requiring G(z) = H(-z) - maybe by not imposing such a symmetry constraint on A(z), more 
// degrees of freedom can be unlocked?
// G(z) = H(-z) -> C(z)/A(z) =  (A(z)-B(z))/A(z) = B(-z)/A(-z) ...this is less restrictive than
// C(z) = B(-z)

// let's start with (A(z)-B(z))/A(z) = B(-z)/A(-z) and define A=A(z),A'=A(-z),B=B(z),B'=B(-z), so:
// (A-B)/A = B'/A' -> A'*(A-B) = A*B' -> A' * A - A' * B = A * B' 
// -> A'A - A'B - AB' = 0 for all z...maybe by strategically selecting some particular z, this 
// leads to useful equations? and/or maybe we can throw a multidimensinal root-finder at this 
// problem? or maybe we can cast this as a constrained optimization problem - the target function
// is some function of the magnitude response - for example, how far it is away from the 
// butterworth lowpass of the same order and the constraint equation is A'A - A'B - AB' = 0?
// maybe start the optimization with a butterworth filter as initial guess?


//-------------------------------------------------------------------------------------------------

void zMapFirstOrder(rsFilterSpecificationZPK<double>& zpk, double g, double c, 
  std::complex<double> zNorm)
{
  using Complex = std::complex<double>;  // Maybe use rsComplex
  Complex one(1);
  Complex k = zpk.k;
  zpk.k = one;
  double kp = rsAbs(zpk.transferFunctionAt(one));
  auto mapRoot = [&](Complex r) { return (g*r - c) / (one - g*r*c); };
  for(int i = 0; i < zpk.z.size(); i++) zpk.z[i] = mapRoot(zpk.z[i]);
  for(int i = 0; i < zpk.p.size(); i++) zpk.p[i] = mapRoot(zpk.p[i]);
  double kt = rsAbs(zpk.transferFunctionAt(zNorm));
  zpk.k = k * (kp/kt);

  // See:
  // http://www.rs-met.com/documents/dsp/TwoInterpretationsOfFrequencyWarpedTransferFunctions.pdf
}

void zMapSecondOrder(rsFilterSpecificationZPK<double>& zpk, double g, double c, double d,
  std::complex<double> zNorm)
{
  using Complex = std::complex<double>;
  Complex one(1);
  Complex k = zpk.k;
  zpk.k = one;
  double kp = rsAbs(zpk.transferFunctionAt(one));
  auto mapRoot = [&](Complex r, Complex* r1, Complex* r2) 
  { 
    Complex k1, k2, k3, k4;
    k1 = 1.0 / (2.*d*g*r - 2.0);
    k2 = k1 * c * (1.0 - g*r);
    k3 = c*c - 4.0*d;
    k4 = k1 * rsSqrt(k3*(1.0+g*g*r*r) + (4.*d*d-2.*c*c+4.)*g*r);
    *r1 = k2 - k4;
    *r2 = k2 + k4;
  };
  int M = (int) zpk.z.size(); zpk.z.resize(2*M);
  int N = (int) zpk.p.size(); zpk.p.resize(2*N);
  for(int i = M-1; i >= 0; --i) mapRoot(zpk.z[i], &zpk.z[2*i], &zpk.z[2*i+1]);
  for(int i = N-1; i >= 0; --i) mapRoot(zpk.p[i], &zpk.p[2*i], &zpk.p[2*i+1]);
  double kt = rsAbs(zpk.transferFunctionAt(zNorm));
  zpk.k = k * (kp/kt);
}

RAPT::rsFilterSpecificationBA<double> zLowpassToLowpass(
  const RAPT::rsFilterSpecificationBA<double>& baProto, double wp, double wt)
{
  rsFilterSpecificationZPK<double> zpk = baProto.toZPK();
  zMapFirstOrder(zpk, 1.0, -sin(0.5*(wp-wt)) / sin(0.5*(wp+wt)), 1.0);
  return zpk.toBA();
}

RAPT::rsFilterSpecificationBA<double> zLowpassToHighpass(
  const RAPT::rsFilterSpecificationBA<double>& baProto, double wp, double wt)
{
  rsFilterSpecificationZPK<double> zpk = baProto.toZPK();
  zMapFirstOrder(zpk, -1.0, -cos(0.5*(wp+wt)) / cos(0.5*(wp-wt)), -1.0);
  return zpk.toBA();
}

RAPT::rsFilterSpecificationBA<double> zLowpassToBandpass(
  const RAPT::rsFilterSpecificationBA<double>& baProto, double wp, double wl, double wu)
{
  using Complex = std::complex<double>;
  rsFilterSpecificationZPK<double> zpk = baProto.toZPK();

  // Maybe factor out (is used also for LP -> BR):
  double a  = - cos(0.5*(wu+wl)) / cos(0.5*(wu-wl));
  double t1 = tan(0.5*wp);
  double t2 = tan(0.5*(wu-wl));

  double k  = t1/t2;
  double wc = 2*atan(sqrt(tan(wl/2.)*tan(wu/2.))); // Normalize at prewarped center frequency
  Complex j(0,1);
  Complex zNorm = exp(j*wc); // Maybe use std::polar instead or maybe write an rsPolarToCartesian
  zMapSecondOrder(zpk, -1.0, 2*a*k/(k+1.0), (k-1.0)/(k+1.0), zNorm);
  return zpk.toBA();
}

RAPT::rsFilterSpecificationBA<double> zLowpassToBandreject(
  const RAPT::rsFilterSpecificationBA<double>& baProto, double wp, double wl, double wu)
{
  using Complex = std::complex<double>;
  rsFilterSpecificationZPK<double> zpk = baProto.toZPK();

  // Maybe factor out (is used also for LP -> BP):
  double a  = - cos(0.5*(wu+wl)) / cos(0.5*(wu-wl));
  double t1 = tan(0.5*wp);
  double t2 = tan(0.5*(wu-wl));

  double k = t1*t2;
  zMapSecondOrder(zpk, 1.0, 2.0*a/(k+1.0), (1.0-k)/(1.0+k), 1.0);
  return zpk.toBA();
}

// -Maybe these functions should go into class rsPoleZeroMapper. The functionality fits there well.
//  But there, they should operate on raw arrays for production use. We can the keep these 
//  functions here as convenience functions that call the other ones.
// -These formulas are quite expensive to implement. Maybe it would be cheaper to do a bilinear
//  z-to-s transform, then do an s-domain frequency transformation, then transform back to the
//  z-domain by a bilinear s-to-z transform? How should be deal with wl, wu? How did RBJ do it? He
//  says something about prewarping the bandwidth - how does that work? Maybe it would make sense
//  to match the center frequency and lower bandedge and just ignore the upper bandedge and let it 
//  fall whereever it falls?
// -What about a lowpass-to-lowshelf trafo in the z-domain?
// -Maybe implement the same functions for rsFilterSpecificationZPK. the BA versions may then
//  just call return zLowpassToLowpass(baProto.toZPK(), wp, wt).toBA(); etc.

//-------------------------------------------------------------------------------------------------
// code below doesn't give rise to working complementary filters but contains a lot of info in the
// comments:

// analog 1-pole/0-zero
void splitterPrototypeA_1_0(double* k, std::complex<double>* p, std::complex<double>* /*z*/)
{
  *k   =  1.0;
  p[0] = -1.0;
}

// analog 2-pole/0-zero
void splitterPrototypeA_2_0(double* k, std::complex<double>* p, std::complex<double>* /*z*/)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = std::complex<double>(-s, s);
  p[1] = conj(p[0]);
}

// analog 2-pole/1-zero
void splitterPrototypeA_2_1(double* k, std::complex<double>* p, std::complex<double>* z)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = std::complex<double>(-s, s);
  p[1] = conj(p[0]);
  z[0] = -1;
  // when we place a zero on the real axis, the ultimate slope is again only 6 dB/oct - seen from a 
  // distance (i.e. being far away on the imaginary axis) the zero cancels the effect of one of the
  // poles
}

// analog 2-pole/2-zero
void splitterPrototypeA_2_2(double* k, std::complex<double>* p, std::complex<double>* z)
{
  double s = SQRT2_INV;
  *k   =  1.0;
  p[0] = std::complex<double>(-s, s);
  p[1] = conj(p[0]);
  z[0] = std::complex<double>(0, 2);
  z[1] = conj(z[0]);
}

// digital 1-pole/1-zero - works
void splitterPrototypeD_1_1(double* k, std::complex<double>* p, std::complex<double>* z)
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

  *k = abs(dcGainNormalizer(z, 1, p, 1));
}
// todo: return a BA specification

// digital 2-pole/2-zero (i think, this is a 2nd order digital Butterworth halfband filter via 
// bilinear transform):
void splitterPrototypeD_2_2(double* k, std::complex<double>* p, std::complex<double>* z)
{
  double s = sqrt(2)-1;  // ad-hoc, i have no derivation for this
  p[0] = std::complex<double>(0, s);
  p[1] = conj(p[0]);
  z[0] = -1.0;
  z[1] = conj(z[0]);
  *k = abs(dcGainNormalizer(z, 2, p, 2));
}
// this doesn't work

// digital 2-pole/3-zero - old - worked when we didn't really treat the reversed order of digital
// filter coeff arrays right in conversions between zpk/ba form - check, why this seemed to work, 
// and if it actually does
void splitterPrototypeD_2_3(double* k, std::complex<double>* p, std::complex<double>* z)
{
  double  s = sqrt(2)-1;
  //s = 0.5;  // test
  std::complex<double> j = std::complex<double>(0, 1);
  p[0] =  j*s;   // p1
  p[1] = -j*s;   // p2
  z[0] = -1;     // q1
  z[1] = -1;     // q2
  z[2] = -s;     // q3
                 //z[2] = -0.5;   // test
  *k = abs(dcGainNormalizer(z, 3, p, 2));  
  // gain factor k to normalize DC gain to 1

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
  //  ...or maybe not - maybe it's sufficient if they are symmetric with respect to the imaginary
  //  axis - being *on* the axis is a special case of that
  // -N (= numPoles) zeros must be at z=-1 (to get a proper lowpass response - not needed for 
  //  complementariness itself) 
  // -maybe the other M-N zeros should also be in the left half-plane?
  //
  // hmmm...with these equations together with the poles and some of the zeros already fixed, we 
  // may be able to solve for the remaining coeffs/zeros.
  // is the condition H(z)=G(-z) actually correct? or should it be |H(z)|=|G(-z)|?
  // -verify numerically, if H(z)=G(-z) for the 2,3 case
  // -check the value of a1 in the 1,1 case - it doens't follow the pattern ai=0 for odd i
}

void splitterPrototypeD_2_3_new(double* k, std::complex<double>* p, std::complex<double>* z)
{
  double  s = sqrt(2)-1;   // like in 2nd order butterworth
  double  t = s;
  //s = 0.5; t = 0.5;
  std::complex<double> j = std::complex<double>(0, 1);
  p[0] =  j*s;   // p1
  p[1] = -j*s;   // p2
  z[0] = -1;     // q1
  z[1] = -1;     // q2
  z[2] = -t;     // q3
  *k = abs(dcGainNormalizer(z, 3, p, 2));
}

// digital 3-pole/3-zero - doesn't work:
void splitterPrototypeD_3_3(double* k, std::complex<double>* p, std::complex<double>* z)
{
  std::complex<double> j = std::complex<double>(0, 1);
  double  a = 0.4;
  p[0] =  j*a;
  p[1] = -j*a;
  p[3] =  0;

  double  b = 0.2;
  z[0] = -1;
  z[1] = -1;
  z[2] = -b;

  *k = abs(dcGainNormalizer(z, 3, p, 3));
}

// digital 4-pole/6-zero (does not work)
void splitterPrototypeD_4_6(double* k, std::complex<double>* p, std::complex<double>* z)
{
  typedef RAPT::rsInfiniteImpulseResponseDesigner<double> DSGNR;
  DSGNR dsgnr;
  dsgnr.setApproximationMethod(RAPT::rsPrototypeDesigner<double>::BUTTERWORTH);
  dsgnr.setPrototypeOrder(4);
  dsgnr.setSampleRate(1);
  dsgnr.setFrequency(0.25);      // halfband filter, 0.5 is Nyquist freq
  dsgnr.getPolesAndZeros(p, z);

  z[4] = -p[0].imag(); // test
  z[5] = -p[2].imag(); // test

  *k = abs(dcGainNormalizer(z, 6, p, 4));
}





/*

Notes:


-putting additional finite zeros into the s-plane is not a good idea - it makes the final slope
 shallower - the lowpass should have all of its zeros at infinity...soooo that means we have to use
 an allpole lowpass filter and therefore the highpass should have all its zeros at s=0. the only 
 wiggle room is the exact placement of the poles
-OR: we design a halfband lowpass prototype in the digital domain, leave the zeros at z = -1 
 and add *additional* zeros. this seems to work for the 2nd order case at least
-try "contracted Butterworth" - all pole angles are scaled by a factor < 1
 try to place poles on a ellipse instead of a circle - let the user select a width/height 
 ratio or eccentricity
-maybe experiment with multiplicities
-maybe start with z-plane prototype poles (and zeros), maybe aligned along the imaginary axis
 and spread in various ways
-place N poles along the imaginary axis and N zeros at z = -1. then try to place additional 
 zeros into the z-plane such that we get a nice crossover...maybe we need a GUI for freely
 placing poles and zeros into the z-plane


*/