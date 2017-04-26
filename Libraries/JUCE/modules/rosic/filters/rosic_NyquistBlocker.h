#ifndef rosic_NyquistBlocker_h
#define rosic_NyquistBlocker_h

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

  /**

  This is a simple first order lowpass filter tuned to shortly below the Nyquist frequency 
  (-3dB point is at 22040 Hz when sample-rate is at 44100 Hz). It is intended to be used to 
  counteract parasitic oscillations at the Nyquist frequency which may occur in certain signal
  processing algorithms.

  */

  class NyquistBlocker
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    NyquistBlocker() { reset(); }


    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE double getSample(double in)
    {
      y1 = 0.99928812771611986 * (in+x1) - 0.99857625543223971 * y1;
      x1 = in;
      return y1;
    }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state variables. */
    void reset() { x1 = 0.0; y1 = 0.0; }

    //=============================================================================================

  protected:

    // state variables:
    double x1, y1;

  };

} // end namespace rosic

#endif // rosic_NyquistBlocker_h
