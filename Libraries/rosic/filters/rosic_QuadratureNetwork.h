#ifndef rosic_QuadratureNetwork_h
#define rosic_QuadratureNetwork_h

// rosic-indcludes:
#include "rosic_InfiniteImpulseResponseDesigner.h"

namespace rosic
{

  /**

  This class implements a pair of filters which approximates a 90 degree phase shift between both
  output signals. It accomplishes this by starting from an elliptic halfband lowpass filter and
  rotating its pole/zero pattern by 90 degrees, thus leaving only the positive frequencies in the
  (now complex) output signal. The real and imaginary parts of this complex signal are now
  phase-shifted by 90 degrees with respect to one another.

  */

  class QuadratureNetwork 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    QuadratureNetwork();   

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Chooses one of the approximation methods as defined in 
    PrototypeDesigner::approximationMethods. */
    void setApproximationMethod(int newApproximationMethod);

    /** Selects the order of the filter. */
    void setOrder(int newOrder);

    /** Sets the ripple in the passband in decibels. */
    void setRipple(double newPassbandRipple);

    /** Sets the rejection in the stopband in decibels. */
    void setStopbandRejection(double newStopbandRejection);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the approximation method to be used 
    @see enum PrototypeDesigner::approximationMethods. */
    int getApproximationMethod() { return designer.getApproximationMethod(); }

    /** Returns true if the currently selected mode supports a ripple parameter. */
    bool hasCurrentModeRippleParameter() { return designer.hasCurrentModeRippleParameter(); }

    /** Returns true if the currently selected mode supports a rejection parameter. */
    bool hasCurrentModeRejectionParameter() { return designer.hasCurrentModeRejectionParameter(); }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Returns a pair of output samples that represent the in-phase component (real part) and
    the quadrature phase component (imaginary part). */
    INLINE void getOutputSamplePair(double in, double *outRealPart, double *outImaginaryPart);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the effect. */
    void reset();

  protected:

    /** Triggers a re-calculation of the filter coefficients. */
    void updateCoefficients();

    static const int maxOrder = 20;
    Complex poles[maxOrder];
    Complex zeros[maxOrder];
    double  gain;
    int     order;

    Complex x[maxOrder];  // past inputs of the individual stages
    Complex y[maxOrder];  // past outputs of the individual stages 

    InfiniteImpulseResponseDesigner designer;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void QuadratureNetwork::getOutputSamplePair(double in, 
    double *outRealPart, double *outImaginaryPart)
  {
    Complex tmp = gain*in;
    for(int i=0; i<order; i++)
    {
      y[i] = tmp - zeros[i]*x[i] + poles[i]*y[i];
      x[i] = tmp;
      tmp  = y[i];  

      // can(?) be streamlined by noting that x[i] == y[i-1] for i >= 1
    }

    tmp = y[order-1];
    *outRealPart      = tmp.re;
    *outImaginaryPart = tmp.im;
  }

} // end namespace rosic

#endif // rosic_QuadratureNetwork_h
