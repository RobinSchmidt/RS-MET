#ifndef rosic_BesselFilterForGainSignal_h
#define rosic_BesselFilterForGainSignal_h

#include <string.h> // for memmove

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

  /**

  This is an 8th order Bessel filter with a cutoff frequency of 2.5 kHz at 44.1 kHz sample-rate 
  intended to be used to filter the gain signal in dynamics processors to avoid aliasing. 

  */

  class BesselFilterForGainSignal
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. Initializes coefficients for an elliptic halfband filter. */
    BesselFilterForGainSignal();   

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the filter state. */
    void reset();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    //=============================================================================================

  protected:

    // state buffer:
    double w[8];

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

// filter coefficients:
#define B0 0.00000039003194546098
#define B1 0.00000312025556368787
#define B2 0.00001092089447290755
#define B3 0.00002184178894581509
#define B4 0.00002730223618226886
#define B5 0.00002184178894581509
#define B6 0.00001092089447290755
#define B7 0.00000312025556368787
#define B8 0.00000039003194546098

#define A0   1.00000000000000000000
#define A1  -5.97776588510238400000
#define A2  15.78002108385409500000
#define A3 -24.01192954605262000000
#define A4  23.02389936487794500000
#define A5 -14.23800525393077400000
#define A6   5.54298593603332530000
#define A7  -1.24156017742317170000
#define A8   0.12245470040715851000

  INLINE double BesselFilterForGainSignal::getSample(double in)
  {
    // calculate intermediate and output sample via direct form II - the parentheses facilitate 
    // out-of-order execution of the independent additions (for performance optimization):
    double tmp =   (in + TINY)
                 - ( (A1*w[0] + A2*w[1] ) + (A3*w[2]   + A4*w[3]  ) ) 
                 - ( (A5*w[4] + A6*w[5] ) + (A7*w[6]   + A8*w[7]  ) );
   
    double y =     B0*tmp 
                 + ( (B1*w[0] + B2*w[1])  +  (B3*w[2]   + B4*w[3]  ) )  
                 + ( (B5*w[4] + B6*w[5])  +  (B7*w[6]   + B8*w[7]  ) );

    // update state variables:
    memmove(&w[1], &w[0], 7*sizeof(double));
    w[0] = tmp;

    return y;
  }

#undef B0
#undef B1
#undef B2
#undef B3
#undef B4
#undef B5
#undef B6
#undef B7
#undef B8

#undef A0
#undef A1
#undef A2
#undef A3
#undef A4
#undef A5
#undef A6
#undef A7
#undef A8

} // end namespace rosic

#endif // rosic_BesselFilterForGainSignal_h
