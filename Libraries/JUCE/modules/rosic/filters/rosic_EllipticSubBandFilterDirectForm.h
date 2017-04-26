#ifndef rosic_EllipticSubBandFilterDirectForm_h
#define rosic_EllipticSubBandFilterDirectForm_h

#include <string.h> // for memmove

// rosic-indcludes:
#include "rosic_EllipticSubBandFilter.h"

namespace rosic
{

  /**

  This is an elliptic subband filter of 12th order using a Direct Form II implementation structure.

  */

  class EllipticSubBandFilterDirectForm
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. Initializes coefficients for an elliptic halfband filter. */
    EllipticSubBandFilterDirectForm();   

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the subdivision factor, for example 2 for a halfband filter (which passes everything
    below half the Nyquist frequency and stops everything above) or 4 for a quarterband filter. */
    void setSubDivision(double newSubDivision);

    /** Resets the filter state. */
    void reset();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    //=============================================================================================

  protected:

    // state buffer:
    double w[12];

    // filter coefficients:
    double a[13];
    double b[13];

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double EllipticSubBandFilterDirectForm::getSample(double in)
  {
    /*
    // this is the straightforward, non-optimized version:
    double tmp =   in + TINY
                 - a[1]*w[0] - a[2]*w[1] - a[3]*w[2] - a[4]*w[3]  - a[5]*w[4]   - a[6]*w[5]
                 - a[7]*w[6] - a[8]*w[7] - a[9]*w[8] - a[10]*w[9] - a[11]*w[10] - a[12]*w[11];
   
    double y =     b[0]*tmp 
                 + b[1]*w[0] + b[2]*w[1] + b[3]*w[2] + b[4]*w[3]  + b[5]*w[4]   + b[6]*w[5]          
                 + b[7]*w[6] + b[8]*w[7] + b[9]*w[8] + b[10]*w[9] + b[11]*w[10] + b[12]*w[11];
    */

    // calculate intermediate and output sample via direct form II - the parentheses facilitate 
    // out-of-order execution of the independent additions (for performance optimization):
    double tmp =   (in + TINY)
                 - ( (a[1]*w[0] + a[2]*w[1] ) + (a[3]*w[2]   + a[4]*w[3]  ) ) 
                 - ( (a[5]*w[4] + a[6]*w[5] ) + (a[7]*w[6]   + a[8]*w[7]  ) )
                 - ( (a[9]*w[8] + a[10]*w[9]) + (a[11]*w[10] + a[12]*w[11]) );
   
    double y =     b[0]*tmp 
                 + ( (b[1]*w[0] + b[2]*w[1])  +  (b[3]*w[2]   + b[4]*w[3]  ) )  
                 + ( (b[5]*w[4] + b[6]*w[5])  +  (b[7]*w[6]   + b[8]*w[7]  ) )
                 + ( (b[9]*w[8] + b[10]*w[9]) +  (b[11]*w[10] + b[12]*w[11]) );

    // update state variables:
    memmove(&w[1], &w[0], 11*sizeof(double));
    w[0] = tmp;

    return y;
  }

} // end namespace rosic

#endif // rosic_EllipticSubBandFilterDirectForm_h
