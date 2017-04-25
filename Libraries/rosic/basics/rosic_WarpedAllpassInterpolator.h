#ifndef rosic_WarpedAllpassInterpolator_h
#define rosic_WarpedAllpassInterpolator_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a warped allpass interpolator. Allpass interpolators provide a flat 
  magnitude response, regardless of the value of the fractional part of the sample-index. This is 
  in contrast to - for example - a linear interpolator which has a lowpass frequency response when 
  the fractional part is 0.5 (in which case it constitutes a first order moving average filter).
  This comes at the expense of introducing transient artifacts when the fractional part is 
  modulated (due to the IIR nature of the interpolator) and an abberation of the actual 
  delay / readout position in a frequency dependent manner. The coefficient warping in the 'warped' 
  allpass provides for a correction of the delay abberation at DC.

  References:
   -Jon Datorro: Effect Design Part 2 - Delay-Line Modulation and Chorus (JAES, Vol 45, No. 10)

  */

  class WarpedAllpassInterpolator
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    WarpedAllpassInterpolator();

    /** Destructor. */
    ~WarpedAllpassInterpolator();

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the fractional part of the readout position and calculates the allpass
    coefficient (with warping). Note that the value used here as input argument corresponds to 
    1-frac in the paper because the WarpedAllpassInterpolator class sees it from the viewpoint of 
    a phase-accumulating oscillator whereas the paper deals with delaylines. */
    void setFractionalPart(double newFractionalPart)
    { coeff = (1.0-newFractionalPart) / (1.0+newFractionalPart);  }

    /** Sets up the fractional part of the readout position without warping in which case the 
    coefficient is equal to the fractional part itself. */
    void setFractionalPartNonWarped(double newFractionalPart)
    { coeff = newFractionalPart;  }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calclulates one output sample by interpolating between y[0] and y[1] at position as set by
    setFractionalPart. y should point to an array of at least two values. */
    INLINE double getSample(double *y);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state variable to zero (relevant only for allpass interpolators). */
    void reset();

    //=============================================================================================

  protected:

    double coeff;           // coefficient for the interpolator
    double previousOutput;  // previous output sample

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double WarpedAllpassInterpolator::getSample(double *y)
  {
    double result  = y[0] + coeff*y[1] - 0.999*coeff*previousOutput; 
       // multiplication of feedback with 0.999... avoids self-oscillation at fs/2

    previousOutput = result;
    return result;
  }

}  // end namespace rosic

#endif // rosic_WarpedAllpassInterpolator_h
