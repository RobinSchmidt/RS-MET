#ifndef RAPT_FILTERDESIGNFORMULAS_H
#define RAPT_FILTERDESIGNFORMULAS_H

/** Given a biquad coefficient set, this function manipulates the coefficients, such that the
inverse filter results. Be careful to use this only with minimum-phase biquads (zeros inside the
unit circle) because otherwise the inverse filter will be unstable. */
template<class T>
void invertBiquad(T &b0, T &b1, T &b2, T &a1, T &a2);
// Move to some sort of FilterTransformations class, maybe rename to invertBiquadDigital or 
// inverBiquadZ

/** Given two coefficient sets of 1st order filter stages, this function consolidates them into a
single biquad stage. b0[0] is the b0 coefficient of the 1st stage, b0[1] the b0 coefficient of
the 2nd stage, etc. The biquad coefficients are returned in B0, B1, etc. */
template<class T>
void twoOnePolesToBiquad(T b0[2], T b1[2], T a1[2], T &B0, T &B1, T &B2, T &A1, T &A2);
// Move to some sort of FilterStructureConversions class

/** Given either the numerator coefficiets b0, b1 or the denominator coefficients a0, a1 of a
1st order filter (in c0, c1), this function ensures that the root (zero or pole) will be inside
the unit circle. Poles inside the unit circle are required for stability and zeros inside the
unit circle are required for minimum-phase response.
Note: when doing this for a numerator and denominator, the overall gain factor in the magnitude
response may have changed after this reflection (i'm currently not sure, why) - the caller needs
to compensate this. */
template<class T>
void ensureOnePoleRootInsideUnitCircle(T &c0, T &c1);

/** Given 3 pairs of of a normalized radian frequency w[i] (== 2*PI*f/fs) and a desired
magnitude m[i] at this frequency, this function computes coefficients for a first order filter
that matches these contraints (if such a filter is possible - if it is impossible, it will return
a coefficient set that realizes an identity filter) */
template<class T>
void magnitudeMatchedOnePoleCoeffs(T &b0, T &b1, T &a1, T w[3], T m[3]);
// Move into some rsFilterDesignFormulas class.



//-------------------------------------------------------------------------------------------------

/** This class collects a bunch of filter design formulas and algortihms from various sources into 
a utility class of static functions. We generally assume a filter of the form:

  y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-1] - a1*y[n-1] - a2*y[n-2]

but for some filters, there may be less coefficients - for example, for a 2-pole-1-zero filter, the 
b2 coeff would be missing. ...tbc... */

class rsFilterDesignFormulas
{

public:

  /** Assigns the coefficients so as to realize a neutral filter. */
  template<class T>
  static inline void bypassBiquad(T* b0, T* b1, T* b2, T* a1, T* a2);

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

// ToDo:
// -Drag the RBJ cookbook biquad designs into this class, maybe deprecate the old BiquadDesigner
//  class, they should all get the same API as the Vicanek designs: in terms of w0, Q, G
// -Move the first order filter designs that are currently implemented in rsFirstOrderFilterBase
//  here
// -Move the modal and time-domain-biquad filter designs here
// -Implement the 5-point matched designs
// -Implement the slop/tilt filter here
// -Maybe rename to rsSmallFilterDesigner





#endif
