#ifndef rosic_BiquadBase_h
#define rosic_BiquadBase_h

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

  /**

  This class serves as baseclass for biquad filters which realize the difference equation:

  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] + a_1 y[n-1] + a_2 y[n-2] \f]

  it contains only the coefficients as members - past inputs/output samples are supposed to be 
  added in subclasses because different topolies require different state variables.

  */

  class BiquadBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    BiquadBase();

    /** Destructor. */
    ~BiquadBase();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the filter coefficients to new values. */
    void setCoefficients(double newB0, double newB1, double newB2, double newA1, double newA2);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Initializes the biquad coefficients to b0=1.0, b1=b2=a1=a2=0.0 which is essentially a
    bypass 'filter'. */
    void initializeCoefficients();

    //=============================================================================================

  protected:

    // direct form coefficients:
    double b0, b1, b2, a1, a2;

  };

} // end namespace rosic

#endif // rosic_BiquadBase_h
