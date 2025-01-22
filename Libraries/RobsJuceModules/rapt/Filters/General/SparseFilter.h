#ifndef RAPT_SPARSEFILTER_H
#define RAPT_SPARSEFILTER_H


//=================================================================================================


/** A subclass of rsSparseRationalFunction that is meant to deal specifically with transfer 
functions of digital filters. A general transfer function for a digital filter looks like:

          b0 + b1 * z^-1 + b2 * z^-2 + ... + bN * z^-N
  H(z) = ----------------------------------------------
          1  + a1 * z^-1 + a2 * z^-2 + ... + aM * z^-M

Note that this is actually a rational function not in z itself but in z^-1. We override the 
function evaluation operator () to take care of this reciprocation of z. We also implement some 
additional functionality on top of the baseclass that is specific to such transfer functions. For 
example, digital filter transfer functions are usually normalized to a0 = 1, as seen above. We 
implement a check for that condition (and a few others) isCanonical(). We also provide functions to
invert the transfer function (basically, swapping numerator and denominator but maintaining that 
the a0 = 1 still holds after the swap), reflecting the zeros about the unit circle (turning 
minimum phase filters into maximum phase ones), etc. ...TBC... */


template<class T>
class rsSparseDigitalTransferFunction : public rsSparseRationalFunction<T>
{

public:

  using Base = rsSparseRationalFunction<T>;    // For convenience
  using Base::Base;                            // Inherit constructors



  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given
  amount. */
  void addPreDelay(int amountInSamples)
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

  /** Removes the predelay from this filter, if any is present. This makes sure that the lowest
  exponent of z^-1 in the numerator is zero. see getPreDelay(), addPreDelay()  */
  void removePreDelay() { num.shiftPowers(-getPreDelay()); }

  /** Turns the filter into its inverse. This basically amounts to swapping numerator and
  denominator and possibly applying some scaling of the coefficients if b0 != 1. */
  void invert()
  {
    rsAssert(isCanonical());

    // A filter with predelay cannot be inverted in realtime. The best thing we can do in this 
    // case is to invert the filter up to the predelay. Removing the predelay ensures that the
    // 0-th term in the numerator has power of 0, i.e. it's a  b0 * z^-0  term and not some crazy
    // b7 * z^-7  term:
    removePreDelay();
    rsAssert(num.getPower(0) == 0); // Numerator is already asserted to be non-empty in isCanoncial
    rsAssert(num.getCoeff(0) != 0); // so we can access the 0-th element without risk here

    // Swap numerator and denominator while maintaining the a0 = 0 normalization condition:
    T s = T(1) / num.getCoeff(0);
    scale(s);
    std::swap(num, den);
    scale(s);

    // I'm pretty sure it doesnt' allocate. Verify and document.
  }

  /** Reflects the zeros of the filter about the unit circle. This will turn a minimum phase
  filter into a maximum phase one and vice versa. For mixed phase filters, it inverts the mix.
  It doesn't affect stability or filter order. */
  void reflectZeros()
  {
    int deg = num._getDegree();                   // ToDo: use canonical getDegree
    for(int i = 0; i < num.getNumTerms(); i++)
      num._setPower(i, deg - num.getPower(i));
    num._reverse();                               // Order array by ascending powers again

    // How about a reflectPoles() function? But that would turn stable filters into unstable ones,
    // so it's usefulness is questionable. For the time being, we can do without. And maybe we
    // should also apply complex conjugation in case of complex coeffs?
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the predelay introduced by this filter. A predelay is characterized by the fact that
  the lowest exponent of z^-1 in the numerator is not zero. That means the numerator does not look
  like b0 + b1*z-^1 + b2*z^-2 + ... but rather something like b7*z^-7 + b8*z^-8 + b9*z^-9 + ...
  for a predelay of 7 samples, for example. */
  int getPreDelay() const { return num.getPower(0); }

  /** Returns the order of the filter. This is the maximum exponent of z^-1 that occurs in the
  transfer function. */
  int getFilterOrder() const { return rsMax(num._getDegree(), den._getDegree()); }
  // ToDo: use canonical getDegree

  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isCanonical() const
  {
    bool ok = true;

    // Numerator and denominator polynomials should not be empty:
    ok &= num.getNumTerms() > 0 && den.getNumTerms() > 0;

    // We assume the filter polynomials to be in canonical shape:
    ok &= num.isCanonical();
    ok &= den.isCanonical();

    // Filter should satisfy the a0 == 1 normalization property:
    ok &= den.getPower(0) == 0 && den.getCoeff(0) == T(1);

    return ok;
  }


  /** Computes the density of the numerator defined as the number of actual nonzero coeffs divided
  by the number of potentially nonzero coeffs given the degree of the numerator. */
  double getNumeratorDensity() const
  { return double(num.getNumTerms()) / double(num._getDegree()+1); }
  // ToDo: use canonical getDegree()

  /** Computes the density of the denominator defined as the number of actual nonzero coeffs 
  divided by the number of potentially nonzero coeffs given the degree of the denominator. */
  double getDenominatorDensity() const
  { return double(den.getNumTerms()-1) / double(den._getDegree()); }
  // ToDo: use canonical getDegree()

  /** Returns the "combined density" defined as the number of actual nonzero coeffs of the filter 
  divided by the number of potential nonzero coeffs for the given filter order. The a0 coeff 
  doesn't count because it's always 1. */
  double getCombinedDensity() const
  {
    int numPossibleCoeffs = 2*getFilterOrder() + 1;
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }

  /** Returns the "separated density" defined as the as number of actual nonzero coeffs in 
  numerator and denominator divided by the number of potential nozero coeffs for the given orders
  of numerator and denominator. */
  double getSeparatedDensity() const
  {
    int numPossibleCoeffs = num._getDegree()+1 + den._getDegree();  // ToDo: use canonical getDegree
    int numActualCoeffs   = num.getNumTerms() + (den.getNumTerms()-1);
    return double(numActualCoeffs) / double(numPossibleCoeffs);
  }
  // ToDo: Explain this better! We consider numerator and denominator as separate filters that can
  // have their own orders and therefore the computation of the number of possibly nonzero coeffs
  // is different. For example, in an Nth order allpole filter, we have a 0th order numerator. In 
  // getCombinedDensity, we would assume that it could potentially have N+1 coeffs. Here, we assume
  // that it can have only one because the order of the numerator is zero.

  // Maybe this could be moved into the baseclass. But I'm not sure, if the +1 fo the num and
  // no +1 for the den also applies there...well...I think, it does when we assume a canonical
  // representation with monic denomionator. Here we normalize the denominator to a0=1 - but
  // whatever way we normalize the function, the denominator has always one degree of freedom
  // less than the numerator (for the same degree). A degree 1 polynomial has 2 coeffs and a degree
  // N polynomial has N+1 coeffs. That number applies to the numerator as is. But the denominator
  // is normalized so we lose one degree of freedom and subtract 1 again.

  // ToDo: Maybe return the densities as rsFraction<int>. They are rational numbers so maybe we 
  // should treat them as such.





  T operator()(T x) const 
  { 
    T xr = T(1) / x;
    return num(xr) / den(xr); 
  }


  template<class TArg>
  TArg operator()(TArg z) const 
  { 
    TArg zr = TArg(1) / z;
    return num(zr) / den(zr);
  }
  // Reciprocation of z needed because we store the coeffs of H(z^-1)

  // We also need to override the evaluateAt functions...or maybe get rid of them in the baseclass.
  // Oh - looks like, we don't have such functions there. OK - good.

  // This boilerplate is needed to have the desired arithmetic operators available also for the 
  // derived class. They can not be inherited from the baseclass because their parameter and
  // return types are different. It may work with pointer-types but not with value-types (I guess):


  rsSparseDigitalTransferFunction<T> operator+(const rsSparseDigitalTransferFunction<T>& q) const 
  { rsSparseDigitalTransferFunction<T> r; weightedSum(*this, T(1), q, T(1), &r, T(0)); return r; }

  rsSparseDigitalTransferFunction<T> operator-(const rsSparseDigitalTransferFunction<T>& q) const 
  { rsSparseDigitalTransferFunction<T> r; weightedSum(*this, T(1), q, T(-1), &r, T(0)); return r; }

  rsSparseDigitalTransferFunction<T> operator*(const rsSparseDigitalTransferFunction<T>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.num, den * q.den); }

  rsSparseDigitalTransferFunction<T> operator/(const rsSparseDigitalTransferFunction<T>& q) const 
  { return rsSparseDigitalTransferFunction(num * q.den, den * q.num); }


};

/** Multiplies a coefficient and a sparse digital transfer function. */
template<class T>
inline rsSparseDigitalTransferFunction<T> operator*(
  const T& s, const rsSparseDigitalTransferFunction<T>& p)
{
  rsSparseDigitalTransferFunction<T> r(p);
  r.scale(s);
  return r;
}




//=================================================================================================

/** A class for representing sparse filters in direct form. We represent them using an object of
type rsSparseRationalFunction to store the coefficients of the transfer function H(z). Well, we
actually store the coeffs of H(z-^1) there because that's what's needed for implementation of the
difference equation. The filter is implemented in direct form 2 using a single delayline. */

template<class TSig, class TPar>
class rsSparseFilter
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the filter from dense arrays of numerator and denominator coeffs. When a coefficient
  in the dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<TPar>& numCoeffs, 
    const std::vector<TPar>& denCoeffs, TPar tol)
  {
    H.num.setupFromDenseCoeffs(numCoeffs, tol);
    H.den.setupFromDenseCoeffs(denCoeffs, tol);
    updateDelayLineLength();
  }
  // This may allocate!
  // ToDo: use num.setupFromDenseCoeffs(&numCoeffs[0], (int) numCoeffs.size(), tol)
  // use  H.setupFromDenseCoeffs(numCoeffs, denCoeffs, tol)


  void setup(const rsSparseRationalFunction<TPar>& newTransferFunction)
  { H.copyDataFrom(newTransferFunction); updateDelayLineLength(); }
  // The parameter should really be of type rsSparseDigitalTransferFunction


  void setNumNumeratorTerms(int newNumTerms) { H.num._setNumTerms(newNumTerms); }
  // This may allocate!

  void setNumDenominatorTerms(int newNumTerms) { H.den._setNumTerms(newNumTerms); }
  // This may allocate!

  void setNumeratorTerm(int index, TPar coeff, int delay) { H.num._setTerm(index, coeff, delay); }

  void setDenominatorTerm(int index, TPar coeff, int delay) { H.den._setTerm(index, coeff, delay); }
  // Actually, we really should call updateDelayLineLength() after setting a term because it 
  // potentially requires a change of the length. But: updateDelayLineLength() is expensive and 
  // setting terms is an operation that might be called in a loop or sequence in which case only
  // one update after the sequence of calls should be done. Maybe when calling it in a sequence,
  // we should use other functions like setNumeratorTermSuppressDelayUpdate. Or maybe give the 
  // function a boolean parameter updateDelayLength which defaults tor true

  // Maybe call them setNumeratorTermNoUpdate(). ...or maybe get rid of these function entirely and
  // instead give the user functions like
  // rsSparsePolynomial<TPar>& getNumeratorRef() for low level access. That breaks encapsulation
  // though. Maybe move them into some extra section "Low level API"



  /** Ensures that the delayline has enough memory allocated to support the desired transfer
  function. This may re-allocate heap memory in cases where the delayline does not already have
  enough capacity, so you don't want to call it on a realtime thread. */
  void ensureEnoughDelayMemory()
  { setMaxDelayInSamples(rsMax(getMaxDelayInSamples(), getFilterOrder())); }

  void setMaxDelayInSamples(int newMaxDelay)
  { delayLine.setMaxDelayInSamples(newMaxDelay); }


  /** Updates the length of the delayline according to the maximum power of z^-1 that occurs in
  numerator and denominator polynomial. */
  void updateDelayLineLength()
  { 
    rsAssert(hasEnoughDelayMemory(), "Not enough delay memory!");
    // When this happens, it means that you are trying to request a transfer function from the
    // filter that it can't support because you didn't pre-allocate enough delay memory. At some 
    // point, you need to call setMaxDelayInSamples() where you set up the maximum possible delay
    // and therefore the maximum possible filter order. If you don't allocate enough and request a
    // too high filter order later, this assertion will trigger. We do not just re-allocate here 
    // because this function is designed to be realtime safe. In a realtime context, you really
    // need to pre-allocate enough on construction or in some prepareToPlay() function or something
    // like that.

    delayLine.setDelayInSamples(getFilterOrder()); 
  }
  // Maybe we actually should re-allocate if necessary but still leave the assertion in. Maybe 
  // that's the most benign way to recover from the error condition in a release build? But nah!
  //
  // ToDo: Figure out and document what happens, when we ignore this error. I think, the delay
  // time will be wrapped around / bitmasked by the actual maxDelay in rsDelay member.
  // I say "actual" because the vaue you set up via setMaxDelayInSamples() may be "rounded up"
  // to the next power of two minus one...or something. See implementation of 
  // rsDelay::readOutputAt(). So if the actual maxDelay is 15, all delays will be 
  // interpreted modulo 16, so readOutputAt(20) would actually amount to a delay of 
  // 20 % 16 = 4 rather than 20. ...verify this!




  /** Applies a scaling factor to the filter. This basically means to scale all numerator coeffs by
  that factor. */
  void scale(TPar scaler) { H.num.scale(scaler); }

  /** Adds an overall predelay to the whole filter by shifting all exponents of z^-1 by the given 
  amount. */
  void addPreDelay(int delayInSamples) { H.addPreDelay(delayInSamples); updateDelayLineLength(); }

  /** Removes any predelay that may be present in the filter. */
  void removePreDelay() { H.removePreDelay(); updateDelayLineLength(); }

  /** Turns the filter into its inverse. */
  void invert() { H.invert(); updateDelayLineLength(); }

  /** Reflects the zeros of the filter about the unit circle. */
  void reflectZeros() { H.reflectZeros(); }
  // Doesn't change the required delayline length so we don't need to call updateDelayLineLength


  /** Copies the settings (i.e. the coefficients) from the given other filter into this object and
  possibly adjusts the delayline length, if necessary. */
  void copySettingsFrom(const rsSparseFilter<TSig, TPar>& other)
  { H.copyDataFrom(other.H); updateDelayLineLength(); }

  // Maybe also provide a copyStateFrom() method which would also copy the content of the 
  // delayline. We could then also have a function copyDataFrom that calls both. Maybe the 
  // delayline class should also have a function copyStateFrom (or maybe copyDataFrom)



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the order of the filter. This is the maximum amount of delay needed to implement 
  the filter. */
  int getFilterOrder() const 
  { 
    return H.getFilterOrder();
    //return rsMax(H.num.getDegree(), H.den.getDegree()); 
  }

  int getMaxDelayInSamples() const { return delayLine.getMaxDelayInSamples(); }
  // Maybe rename to getMaxFilterOrder


  /** Performs some sanity checks. Is meant for debug assertions. */
  bool isFilterValid() const { return H.isCanonical() && areDelaysConsistent(); }
  // ToDo: Elaborate documentation. Give some details about what it checks.





  /** Computes the transfer function H(z) of this filter at the given complex value z. */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const { return H(z); }
  // Yes - that's right! "return H(z)" is the whole implementation. Isn't that elegant? :-D


  // Reciprocation of z needed because H actually stores the coeffs of H(z^-1)
  // ToDo: factor the reciprocation out into the () operator of
  // rsSparseDigitalTransferFunction ...done!


  /** Returns a const reference to our transfer function object H(z). */
  const rsSparseDigitalTransferFunction<TPar>& getTransferFunction() const { return H; }




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  /** Computes one output sample at a time using a direct form 2 implementation. */
  TSig getSample(TSig in)
  {
    rsAssert(isFilterValid());

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < H.den.getNumTerms(); i++)
      tmp -= H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply numerator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < H.num.getNumTerms(); i++)
      tmp += H.num.getCoeff(i) * delayLine.readOutputAt(H.num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }

  /** Computes a sample at a time of the inverse filter. We apply the desired transformation to the
  filter into its inverse on the fly. */
  TSig getSampleInverse(TSig in)
  {
    rsAssert(isFilterValid());
    rsAssert(H.num.getPower(0) == 0);
    rsAssert(H.num.getCoeff(0) != 0);
    // Maybe we can relax this? If we do not expect the power of the 0-th coeff to be 0, we will 
    // just produce an inverted filter up to delay?

    // Apply scaled numerator of H as feedback part:
    TPar s = TPar(1) / H.num.getCoeff(0);
    TSig tmp = in;
    for(int i = 1; i < H.num.getNumTerms(); i++)
      tmp -= s * H.num.getCoeff(i) * delayLine.readOutputAt(H.num.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply scaled denominator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < H.den.getNumTerms(); i++)
      tmp += s * H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }


  /** Computes a sample at a time of a filter that has the numerator transformed from min-phase to
  max-phase or vice versa. For mixed phase filters, it inverts the mix. We reflect the zeros in
  the unit circle. */
  TSig getSamplePhased(TSig in)
  {
    rsAssert(isFilterValid());

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < H.den.getNumTerms(); i++)
      tmp -= H.den.getCoeff(i) * delayLine.readOutputAt(H.den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply reversed numerator of H as feedforward path:
    int deg = H.num._getDegree();  // ToDo: use canonical getDegree()
    tmp = 0;
    for(int i = 0; i < H.num.getNumTerms(); i++)
      tmp += H.num.getCoeff(i) * delayLine.readOutputAt(deg - H.num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }
  // Needs test! This is perhaps not great for realtime use because the H.num.getDegree() call must
  // iterate through the whole numerator. Maybe that value could be cached. Not sure. Although,
  // If we assume H.num to be in canonical representation (which it is, I think), then getDegree()
  // can be replaced by  H.getCoeff(H.getNumTerms()-1)  which avoids the iteration. Maybe such a 
  // call could even be encapsulated into something like H.getLastPower(). Maybe add an
  // H.isCanonical() check to isFilterValid().
  //
  // Maybe factor out the getSample() functions into free functions like 
  //
  // TSig rsGetSample(TSig in,
  //        const rsSparseDigitalTransferFunction<TPar>& H, rsDelay<TSig>* delay);
  //
  // This facilitates memory optimizations in situations where we have multiple filters with the
  // same set of coeffs but independent states, i.e. independent delaylines. We can store the
  // coeffs of H once and use them with different delaylines. As is currently is, we would have to
  // create a sparse filter object for each of the filters and therfore store the coeffs 
  // redundantly. This might become a general pattern for implementing filters: separate coeffs and
  // state and provide free functions that take a const ref to the coeffs and (mutable) a pointer
  // to the state.







  /** Resets the filter's state. This clears the delayline. */
  void reset() { delayLine.reset(); }



protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Self test */

  /** Checks if the maximum delay required by the transfer function (i.e. the highest power of 
  z^-1 that occurrs) is consistent with the length of the delayline. This is meant for internal 
  sanity checks. */
  bool areDelaysConsistent() const { return delayLine.getDelayInSamples() == getFilterOrder(); }
  // Maybe rename to areDelayAndOrderConsistent, doesDelayMatchOrder

  /** Checks, if the delayline has enough memory allocated to support the transfer function H(z).
  The maximum possible delay must be greater or equal to the order of the filter. */
  bool hasEnoughDelayMemory() const 
  { return delayLine.getMaxDelayInSamples() >= getFilterOrder(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  rsDelay<TSig> delayLine;                  // Delayline for the direct form 2 implementation.
  rsSparseDigitalTransferFunction<TPar> H;  // Transfer function H(z). Contains filter coeffs.

};






#endif