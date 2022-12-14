#ifndef RAPT_FILTERDESIGNFORMULAS_H
#define RAPT_FILTERDESIGNFORMULAS_H

/** Given a biquad coefficient set, this function manipulates the coefficients, such that the
inverse filter results. Be careful to use this only with minimum-phase biquads (zeros inside the
unit circle) because otherwise the inverse filter will be unstable. */
template<class T>
void invertBiquad(T &b0, T &b1, T &b2, T &a1, T &a2);
// Move to some sort of FilterTransformations class, maybe rename to invertBiquadDigital or 
// inverBiquadZ, take pointers instead of references

/** Given two coefficient sets of 1st order filter stages, this function consolidates them into a
single biquad stage. b0[0] is the b0 coefficient of the 1st stage, b0[1] the b0 coefficient of
the 2nd stage, etc. The biquad coefficients are returned in B0, B1, etc. */
template<class T>
void twoOnePolesToBiquad(T b0[2], T b1[2], T a1[2], T &B0, T &B1, T &B2, T &A1, T &A2);
// Move to some sort of FilterStructureConversions class
// use pointers for output variables

/** Given either the numerator coefficiets b0, b1 or the denominator coefficients a0, a1 of a
1st order filter (in c0, c1), this function ensures that the root (zero or pole) will be inside
the unit circle. Poles inside the unit circle are required for stability and zeros inside the
unit circle are required for minimum-phase response.
Note: when doing this for a numerator and denominator, the overall gain factor in the magnitude
response may have changed after this reflection (i'm currently not sure, why) - the caller needs
to compensate this. */
template<class T>
void ensureOnePoleRootInsideUnitCircle(T &c0, T &c1);
// use pointers

/** Given 3 pairs of of a normalized radian frequency w[i] (== 2*PI*f/fs) and a desired
magnitude m[i] at this frequency, this function computes coefficients for a first order filter
that matches these contraints (if such a filter is possible - if it is impossible, it will return
a coefficient set that realizes an identity filter) */
template<class T>
void magnitudeMatchedOnePoleCoeffs(T &b0, T &b1, T &a1, T w[3], T m[3]);
// Move into some rsFilterDesignFormulas class.
// use pointers



//-------------------------------------------------------------------------------------------------

/** This class collects a bunch of filter design formulas and algortihms from various sources into 
a utility class of static functions. We generally assume a filter of the form:

  y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-1] - a1*y[n-1] - a2*y[n-2]

but for some filters, there may be less coefficients - for example, for a 2-pole-1-zero filter, the 
b2 coeff would be missing. 

The API of the functions is as follows:
-There are sometimes two template parameters, one for the user parameters (TPar), one for the 
 coefficients (TCof). The rationale is that the formulas may use double-precision but the end 
 result may be assigned to a single precision variable.


...tbc... */

class rsFilterDesignFormulas
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Simple designs

  /** Assigns the biquad coefficients so as to realize a neutral filter. */
  template<class T>
  static inline void bypassBiquad(T* b0, T* b1, T* b2, T* a1, T* a2);


  //-----------------------------------------------------------------------------------------------
  // \name Martin Vicanek's magnitude matched designs

  /** Lowpass design by Martin Vicanek, simplified (cheaper) version. Uses MZT for the poles and
  magnitude matching at DC and fs/2 for the zeros. Approximates the analog lowpass prototype 
  transfer function:

                   w0^2
    H(s) = ---------------------
            w0^2 + s*w0/Q + s^2     */
  template<class T>
  static inline void mvLowpassSimple(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2);
  // maybe rename to mvLowpass2p1zSimple - yes, it actually has only one zero, so b2 is always
  // assigned to zero - but for consistency of the API, we require the caller to pass a pointer to 
  // a b2 variable

  /** Highpass design by Martin Vicanek, simplified (cheaper) version. Uses MZT for the poles and
  requires a double zero at DC and magnitude match at fs/2 for the zeros. Approximates the analog 
  highpass prototype transfer function:

                   s^2
    H(s) = ---------------------
            w0^2 + s*w0/Q + s^2   */
  template<class T>
  static inline void mvHighpassSimple(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2);
  // maybe rename to mvLowpass2p2zSimple

  /** Bandpass design by Martin Vicanek, simplified (cheaper) version. Uses MZT for the poles and
  requires a single zero at DC, a slope match at DC and a magnitude match at fs/2. Approximates the
  analog bandpass prototype transfer function:

                   s*w0/Q
    H(s) = ---------------------
            w0^2 + s*w0/Q + s^2    */
  template<class T>
  static inline void mvBandpassSimple(T w0, T Q, bool constSkirt, 
    T* b0, T* b1, T* b2, T* a1, T* a2);


  //-----------------------------------------------------------------------------------------------
  // \name Robin Schmidt's designs:

  /** Designs a two-pole-one-zero filter in terms of its impulse response which is given as the 
  damped sinusoid:

    h[n] = A * exp(-n/d) * sin(w*n + p) * u[n]

  where u[n] is the unit step function. The filter can be implemented as:

    y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] - a2*y[n-2]

  The parameters are:
  w: Normalized radian frequency (= 2*pi*f/fs, f: frequency, fs: samplerate)
  A: Amplitude
  d: Normalized decay time constant (= tau*fs, tau: time (in s) to decay to A/e = A*0.3678...)
  p: Start phase in radians.

  You can use different datatypes for the parameters and coefficients (for example, double and 
  float). */
  template<class TPar, class TCof>
  static inline void dampedSine(
    TPar w, TPar A, TPar d, TPar p, TCof* b0, TCof* b1, TCof* a1, TCof* a2);
  // ToDo: implement inverse function that computes the parameters from the coeffs


protected:

  /** Computes feedback coeffs for the Vicanek designs. Implements the matched-z-trnasform (MZT) 
  which leads to an invariance of the impulse response when there are no zeros (verify!). */
  template<class T>
  static inline void mvFeedbackCoeffs(T w0, T Q, T* a1, T* a2);
  // Maybe make protected...but maybe it's useful for client code? Maybe it can be used in more
  // general contexts? Maybe for an impulse invariant 2-pole-0-zero resonator filter?

};

template<class T>
inline void rsFilterDesignFormulas::bypassBiquad(T* b0, T* b1, T* b2, T* a1, T* a2)
{
  *b0 = T(1);
  *b1 = T(0);
  *b2 = T(0);
  *a1 = T(0);
  *a2 = T(0);
}

//-------------------------------------------------------------------------------------------------
// Martin Vicanek's Designs:

template<class T>
inline void rsFilterDesignFormulas::mvFeedbackCoeffs(T w0, T Q, T* a1, T* a2)
{
  T q = T(0.5) / Q;
  if(q <= 1)
    *a1 = -2*exp(-q*w0) * cos( sqrt(1-q*q)*w0);
  else
    *a1 = -2*exp(-q*w0) * cosh(sqrt(q*q-1)*w0);
  *a2 = exp(-2*q*w0);
}
// Maybe try to optimize away one of the exp calls:
//   t   = exp(-q*w0);
//   *a1 = -2*t * cos( sqrt(1-q*q)*w0);  etc
//   *a2 = t*t

template<class T>
inline void rsFilterDesignFormulas::mvLowpassSimple(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  mvFeedbackCoeffs(w0, Q, a1, a2);

  T f0  = w0 * T(1.0/PI);
  T f02 = f0*f0;
  T k   = (1-f02);
  T r0  = 1 + *a1 + *a2;
  T r1  = (1 - *a1 + *a2)*f02 / sqrt(k*k + f02/(Q*Q));

  *b0 = T(0.5) * (r0 + r1);
  *b1 = r0 - *b0;
  *b2 = T(0);
}

template<class T>
inline void rsFilterDesignFormulas::mvHighpassSimple(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  mvFeedbackCoeffs(w0, Q, a1, a2);

  T f0  = w0 * T(1.0/PI);
  T f02 = f0*f0;
  T k   = (1-f02);
  T r1  = (1 - *a1 + *a2) / sqrt(k*k + f02/(Q*Q));

  *b0 = T(0.25) * r1;
  *b1 = T(-2)   * *b0;
  *b2 = *b0;
}

template<class T>
inline void rsFilterDesignFormulas::mvBandpassSimple(
  T w0, T Q, bool constSkirt, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  mvFeedbackCoeffs(w0, Q, a1, a2);

  T f0  = w0 * T(1.0/PI);
  T f02 = f0*f0;
  T k   = (1-f02);
  T r0  = (1 + *a1 + *a2) / (PI * f0 * Q);
  T r1  = ((1 - *a1 + *a2)*f0/Q) / sqrt(k*k + f02/(Q*Q));

  *b1 = T(-0.5) * r1;
  *b0 = T( 0.5) * (r0 - *b1);
  *b2 = -*b0 - *b1;

  if(constSkirt) {
    *b0 *= Q; *b1 *= Q; *b2 *= Q; } // Formula needs to be verified
}

//-------------------------------------------------------------------------------------------------
// Robin Schmidt's Designs:

template<class TPar, class TCof>
inline void rsFilterDesignFormulas::dampedSine(
  TPar w, TPar A, TPar d, TPar p, TCof* b0, TCof* b1, TCof* a1, TCof* a2)
{
  TPar cw, sw, cp, sp, P;
  rsSinCos(w, &sw, &cw);
  rsSinCos(p, &sp, &cp);
  P   = exp(TPar(-1)/d);          // = exp(-alpha), pole radius
  *a1 = TCof(-2*P*cw);            // = -2*P*cos(w)
  *a2 = TCof(P*P);                // = P^2
  *b0 = TCof(A*sp);               // = A*sin(p)
  *b1 = TCof(A*P*(sw*cp-cw*sp));  // = A*P*sin(w-p) via addition theorem
}



/*

ToDo:
-Drag the RBJ cookbook biquad designs into this class, maybe deprecate the old BiquadDesigner
 class, they should all get the same API as the Vicanek designs: in terms of w0, Q, G and be 
 named rbjLowpass, rbjHighpass, etc.
 -the old BiquadDesigner class may then be deprecated for good - the API is a mess anyway.
-Drag also the Orfanidis "prescribed Nyquist frequency gain" peak EQ design in, name it 
 soPeak (for Sophocles Orfanidis)
-Implement the formulas from the DAFX book - maybe prefix them by uz for Udo Zoelzer
-Move the first order filter designs that are currently implemented in rsFirstOrderFilterBase
 here
-Move the modal and time-domain-biquad filter designs here, i.e. rsDampedSineFilterCoeffs in 
 ModalFilterBank.h/cpp
-Implement the 5-point matched designs
-Implement the slop/tilt filter here
-Maybe rename to rsSmallFilterDesigner
-Try to get rid of the code duplication in the mv-designs. The computations of f0, f02, k, r1 are
 always the same so maybe factor them out into a (protected or private) function.
-Add the non-simplified Vicanek designs for lowpass, highpass, bandpass, peak. The formulas are
 already implemented and ready for copy-and-paste in biquadDesignVicanek in FilterExperiments.cpp
 -Check, if peak-design needs some fudging when G > 1 such that the bandwidth always matches 
  such that cuts and boost by inverse amount cancel each other out

*/



#endif
