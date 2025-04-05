#ifndef RAPT_SPARSEFILTER_H
#define RAPT_SPARSEFILTER_H






//=================================================================================================

/** A class for representing sparse filters in direct form. We represent them using an object of
type rsSparseRationalFunction to store the coefficients of the transfer function H(z). Well, we
actually store the coeffs of H(z-^1) there because that's what's needed for implementation of the
difference equation. The filter is implemented in direct form 2 using a single delayline. */

template<class TSig, class TPar, class TTol>  // ToDo: have optional TTol template parameter for passing to H
class rsSparseFilter
{

public:


  using SparsePoly = rsSparsePolynomial<TPar, TTol>;  // For convenience

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the filter from dense arrays of numerator and denominator coeffs. When a coefficient
  in the dense representation is zero, we not create a term for that. */
  void setupFromDenseCoeffs(const std::vector<TPar>& numCoeffs, 
    const std::vector<TPar>& denCoeffs, TPar tol)
  {
    H.num.setRoundoffTolerance(tol);
    H.den.setRoundoffTolerance(tol);
    H.num.setupFromDenseCoeffs(numCoeffs);
    H.den.setupFromDenseCoeffs(denCoeffs);
    updateDelayLineLength();
  }
  // This may allocate!
  // ToDo: use num.setupFromDenseCoeffs(&numCoeffs[0], (int) numCoeffs.size(), tol)
  // use  H.setupFromDenseCoeffs(numCoeffs, denCoeffs, tol)


  void setup(const rsSparseRationalFunction<TPar, TTol>& newTransferFunction)
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
  void copySettingsFrom(const rsSparseFilter<TSig, TPar, TTol>& other)
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
  const rsSparseDigitalTransferFunction<TPar, TTol>& getTransferFunction() const { return H; }




  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  /** Computes one output sample at a time using a direct form 2 implementation. */
  TSig getSample(TSig in)
  {
    rsAssert(isFilterValid());

    const SparsePoly& num = H.getNumeratorConst();
    const SparsePoly& den = H.getDenominatorConst();

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < den.getNumTerms(); i++)
      tmp -= den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply numerator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < num.getNumTerms(); i++)
      tmp += num.getCoeff(i) * delayLine.readOutputAt(num.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }

  /** Computes a sample at a time of the inverse filter. We apply the desired transformation to the
  filter into its inverse on the fly. */
  TSig getSampleInverse(TSig in)
  {
    rsAssert(isFilterValid());

    const SparsePoly& num = H.getNumeratorConst();
    const SparsePoly& den = H.getDenominatorConst();

    rsAssert(num.getPower(0) == 0);
    rsAssert(num.getCoeff(0) != 0);
    // Maybe we can relax this? If we do not expect the power of the 0-th coeff to be 0, we will 
    // just produce an inverted filter up to delay?

    // Apply scaled numerator of H as feedback part:
    TPar s = TPar(1) / num.getCoeff(0);
    TSig tmp = in;
    for(int i = 1; i < num.getNumTerms(); i++)
      tmp -= s * num.getCoeff(i) * delayLine.readOutputAt(num.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply scaled denominator of H as feedforward path:
    tmp = 0;
    for(int i = 0; i < den.getNumTerms(); i++)
      tmp += s * den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;
  }
  // Check, if this has unit tests!


  /** Computes a sample at a time of a filter that has the numerator transformed from min-phase to
  max-phase or vice versa. For mixed phase filters, it inverts the mix. We reflect the zeros in
  the unit circle. */
  TSig getSamplePhased(TSig in)
  {
    rsAssert(isFilterValid());

    const SparsePoly& num = H.getNumeratorConst();
    const SparsePoly& den = H.getDenominatorConst();

    // Apply denominator of H as feedback part:
    TSig tmp = in;
    for(int i = 1; i < den.getNumTerms(); i++)
      tmp -= den.getCoeff(i) * delayLine.readOutputAt(den.getPower(i));
    delayLine.writeInputNoUpdate(tmp);

    // Apply reversed numerator of H as feedforward path:
    int deg = num.getDegree();
    tmp = 0;
    for(int i = 0; i < num.getNumTerms(); i++)
      tmp += num.getCoeff(i) * delayLine.readOutputAt(deg - num.getPower(i));

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

  rsDelay<TSig> delayLine;                        // Delayline for direct form 2 implementation.
  rsSparseDigitalTransferFunction<TPar, TTol> H;  // Transfer function H(z). Has filter coeffs.

};






#endif