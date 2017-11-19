#ifndef RAPT_QUADRATURENETWORK_H_INCLUDED
#define RAPT_QUADRATURENETWORK_H_INCLUDED

/** This class implements a pair of filters which approximates a 90 degree phase shift between 
both output signals. It accomplishes this by starting from an elliptic halfband lowpass filter and
rotating its pole/zero pattern by 90 degrees, thus leaving only the positive frequencies in the
(now complex) output signal. The real and imaginary parts of this complex signal are now
phase-shifted by 90 degrees with respect to one another. */

template<class TSig, class TPar>
class QuadratureNetwork
{

  // preliminary:
  typedef std::complex<TPar> ComplexPar; 
  typedef std::complex<TSig> ComplexSig;

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  QuadratureNetwork();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Chooses one of the approximation methods as defined in
  rsPrototypeDesigner::approximationMethods. */
  void setApproximationMethod(int newApproximationMethod);

  /** Selects the order of the filter. */
  void setOrder(int newOrder);

  /** Sets the ripple in the passband in decibels. */
  void setRipple(TPar newPassbandRipple);

  /** Sets the rejection in the stopband in decibels. */
  void setStopbandRejection(TPar newStopbandRejection);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the approximation method to be used
  @see enum rsPrototypeDesigner::approximationMethods. */
  int getApproximationMethod() { return designer.getApproximationMethod(); }

  /** Returns true if the currently selected mode supports a ripple parameter. */
  bool hasCurrentModeRippleParameter() { return designer.hasCurrentModeRippleParameter(); }

  /** Returns true if the currently selected mode supports a rejection parameter. */
  bool hasCurrentModeRejectionParameter() { return designer.hasCurrentModeRejectionParameter(); }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Returns a pair of output samples that represent the in-phase component (real part) and
  the quadrature phase component (imaginary part). */
  inline void getOutputSamplePair(TSig in, TSig *outReal, TSig *outImag)
  {
    Complex tmp = gain*in;
    for(int i = 0; i < order; i++)
    {
      y[i] = tmp - zeros[i]*x[i] + poles[i]*y[i];
      x[i] = tmp;
      tmp  = y[i];
      // can(?) be streamlined by noting that x[i] == y[i-1] for i >= 1
    }
    tmp = y[order-1];
    *outReal = tmp.re;
    *outImag = tmp.im;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the effect. */
  void reset();

protected:

  /** Triggers a re-calculation of the filter coefficients. */
  void updateCoefficients();

  static const int maxOrder = 20;
  ComplexPar poles[maxOrder];
  ComplexPar zeros[maxOrder];
  TPar gain;
  int  order;

  ComplexSig x[maxOrder];  // past inputs of the individual stages
  ComplexSig y[maxOrder];  // past outputs of the individual stages 

  rsInfiniteImpulseResponseDesigner designer;
};

#endif
