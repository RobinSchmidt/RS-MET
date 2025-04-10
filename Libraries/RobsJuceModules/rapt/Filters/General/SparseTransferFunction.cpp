
template<class T, class TTol>
void rsSparseTransferFunction<T, TTol>::addPreDelay(int amountInSamples)
{
  if(amountInSamples < 0)
  {
    rsError("Negative predelay is not allowed");
    return;
    // We could allow it if we already have some predelay and the given amount would just 
    // reduce it. Maybe we can relax the restriction to amountInSamples >= -num.getMinPower() or
    // something. If the minimum power is 3, we could allow a predelay amount of -3. This would 
    // then just reduce the predelay to zero.
  }

  num.shiftPowers(amountInSamples);

  // Maybe do num.shiftPowers(rsMax(amountInSamples, -getPreDelay()) ); and write unit tests for
  // this
}

template<class T, class TTol>
bool rsSparseTransferFunction<T, TTol>::_isCanonical() const
{
  bool ok = true;

  // Numerator and denominator polynomials should not be empty:
  ok &= num.getNumTerms() > 0;    // Maybe empty numerator is admissible to represent H(z) = 0? 
  ok &= den.getNumTerms() > 0;

  // We assume the num, den polynomials to be in canonical representation:
  ok &= num._isCanonical();
  ok &= den._isCanonical();

  // Filter should satisfy the a0 == 1 normalization property:
  ok &= den.getPower(0) == 0;
  ok &= den.getCoeff(0) == T(1);  // Should we use a tolerance? ...but maybe not.

  // The roundoff tolerances of numerator and denominator should match:
  ok &= num.getRoundoffTolerance() == den.getRoundoffTolerance();

  return ok;


  // ToDo:
  //
  // - Maybe we are too strict here. Maybe we should allow for empty numerators. I'm also not sure
  //   if an exact comparison to 1 is appropriate the const coeff of the denominator. Maybe we 
  //   should use an inexact comparison, i.e. use some tolerance. We'll see...
  //
  // - Maybe put the tests that access getPower(0), getCoeff(0) into a conditional to avoid access
  //   violations when the den actually is empty.
}

template<class T, class TTol>
void rsSparseTransferFunction<T, TTol>::invert()
{
  rsAssert(_isCanonical());

  // A filter with predelay cannot be inverted (at least not if it operates in realtime). The 
  // best thing we can do in this case is to invert the filter up to the predelay. Removing the
  // predelay ensures that the 0-th term in the numerator has power of 0, i.e. it's a  b0 * z^-0  
  // term:
  removePreDelay();
  rsAssert(num.getPower(0) == 0); // Numerator is already asserted to be non-empty in isCanoncial
  rsAssert(num.getCoeff(0) != 0); // ..so we can access the 0-th element without risk here

  // Swap numerator and denominator while maintaining the a0 = 0 normalization condition by 
  // appropriately scaling the numerator before and after the swap:
  T s = T(1) / num.getCoeff(0);  // Desired scaler s is 1/b0
  _scale(s);                     // Scale old numerator to achieve b0 = 1 before swap
  std::swap(num, den);           // Swap numerator and denominator. b0 is now 1 because a0 was.
  _scale(s);                     // Scale new numerator to achieve desired overall gain

  // I'm pretty sure it doesnt' allocate. The swap of the underyling std::vectors should use move 
  // semantics. Verify and document this. We need to be able to call this function on a realtime 
  // thread in some damped comb allpass filters, so it is important that this function is 
  // non-allocating.
}

template<class T, class TTol>
void rsSparseTransferFunction<T, TTol>::reflectZeros()
{
  int deg = num.getDegree();
  for(int i = 0; i < num.getNumTerms(); i++)
    num._setPower(i, deg - num.getPower(i));
  num._reverse();                               // Order array by ascending powers again

  // How about a reflectPoles() function? But that would turn stable filters into unstable ones,
  // so it's usefulness is questionable. For the time being, we can do without. Maybe we should 
  // also apply complex conjugation of the coeffs in case of complex coeffs? If so, maybe
  // do it in a function num.conjugateCoeffs() which just calls rsConj() on each coeff (which is
  // an empty function for real types)
}

//=================================================================================================
/*

Notes:

- The class was initially named rsSparseDigitalTransferFunction because we may possibly also
  have situations where we need to deal with analog transfer functions. But the name was too long
  and in a DSP library, we may just treat the "digital" qualifier as some default thing that goes 
  without saying. When we need a class for analog transfer functions, then maybe only that should
  have a a qualifier like rsSparseAnalogTransferFunction. But actually, for these, we can just use
  the raw class rsSparseRationalFunction so we don't even need a dedicated class for these.

*/