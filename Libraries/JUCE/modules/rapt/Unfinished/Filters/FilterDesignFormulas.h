#ifndef RS_FILTERDESIGNFORMULAS_H
#define RS_FILTERDESIGNFORMULAS_H

namespace RSLib
{

  /** Given a biquad coefficient set, this function manipulates the coefficients, such that the 
  inverse filter results. Be careful to use this only with minimum-phase biquads (zeros inside the 
  unit circle) because otherwise the inverse filter will be unstable. */
  void RSLib_API invertBiquad(double &b0, double &b1, double &b2, double &a1, double &a2);

  /** Given two coefficient sets of 1st order filter stages, this function consolidates them into a 
  single biquad stage. b0[0] is the b0 coefficient of the 1st stage, b0[1] the b0 coefficient of 
  the 2nd stage, etc. The biquad coefficients are returned in B0, B1, etc. */
  void RSLib_API twoOnePolesToBiquad(double b0[2], double b1[2], double a1[2],                                 
    double &B0, double &B1, double &B2, double &A1, double &A2);

  /** Given either the numerator coefficiets b0, b1 or the denominator coefficients a0, a1 of a 
  1st order filter (in c0, c1), this function ensures that the root (zero or pole) will be inside
  the unit circle. Poles inside the unit circle are required for stability and zeros inside the 
  unit circle are required for minimum-phase response.
  Note: when doing this for a numerator and denominator, the overall gain factor in the magnitude
  response may have changed after this reflection (i'm currently not sure, why) - the caller needs
  to compensate this. */
  void RSLib_API ensureOnePoleRootInsideUnitCircle(double &c0, double &c1);

  /** Given 3 pairs of of a normalized radian frequency w[i] (== 2*PI*f/fs) and a desired 
  magnitude m[i] at this frequency, this function computes coefficients for a first order filter 
  that matches these contraints (if such a filter is possible - if it is impossible, it will return 
  a coefficient set that realizes an identity filter) */
  void RSLib_API magnitudeMatchedOnePoleCoeffs(double &b0, double &b1, double &a1, 
    double w[3], double m[3]);

} 

#endif
