#ifndef rosic_InstantaneousEnvelopeDetector_h
#define rosic_InstantaneousEnvelopeDetector_h

// rosic-indcludes:
#include "../filters/rosic_QuadratureNetwork.h"

namespace rosic
{

  /**

  This is an instantaneous envelope detector based on a QuadratureNetwork.

  */

  class InstantaneousEnvelopeDetector
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    InstantaneousEnvelopeDetector();  

    /** Destructor. */
    ~InstantaneousEnvelopeDetector();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:   

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Estimates the instantaneous amplitude envelope of the incoming signal (as raw amplitude). */
    INLINE double getInstantaneousEnvelope(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state. */
    void reset() { quadratureNetwork.reset(); }

    //=============================================================================================

  protected:

    QuadratureNetwork quadratureNetwork;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double InstantaneousEnvelopeDetector::getInstantaneousEnvelope(double in)
  {
    double re, im;
    quadratureNetwork.getOutputSamplePair(in, &re, &im);
    return sqrt(re*re + im*im);
  }

} // end namespace rosic

#endif // rosic_InstantaneousEnvelopeDetector
