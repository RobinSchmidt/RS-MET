#ifndef rosic_DirectFormFilter_h
#define rosic_DirectFormFilter_h

#include <string.h> // for memmove
//#include <math.h>

//#include "../basics/GlobalDefinitions.h"
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a generic filter using a Direct Form II implementation structure.

  */

  class DirectFormFilter
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - the passed parameter will determine the maximum order which this filter will be able to realize.  */
    DirectFormFilter(int maximumOrder);   

    /** Destructor. */
    ~DirectFormFilter();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Srts up the filter coefficients to use. */
    void setCoefficients(double *newCoeffsA, double *newCoeffsB, int newOrder);

    /** Applies a global gain factor (by multiplying all b-coefficients with tha factor). */
    void setGlobalGainFactor(double newFactor);

    /** Resets the filter state. */
    void reset();

    /** Initializes the filter coefficients so as to realize a 'bypass' filter. */
    void initializeCoefficients();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inqiury:

    /** Returns the magnitude response at the given normalized radian frequency 'omega' where pi corresponds to half the sample-rate. */
    double getMagnitudeResponseAt(double omega);

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array "frequencies" (in Hz) and stores them 
    in the array "magnitudes". Both arrays are assumed to be "numBins" long. "inDecibels" indicates, if the frequency response should be 
    returned in decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be multiplied with (or added to, 
    when "inDecibels" is true) to the magnitudes which are already there in the "magnitudes"-array. This is useful for calculating the 
    magnitude response of several filters in series. */
    void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, double sampleRate, bool inDecibels = false, 
      bool accumulate = false);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample. */
    INLINE double getSample(double in);

    //=====================================================================================================================================

  protected:

    double *w;       // state buffer
    double *a, *b;   // filter coefficients
    int    order;    // currently used filter order - \todo maybe have separate orders for numerator and denominator
    int    maxOrder; // maximum possible filter order

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double DirectFormFilter::getSample(double in)
  {
    // obtain intermediate signal as output from feedback part:
    int k;
    double tmp = in + TINY;
    for(k=1; k<=order; k++)
      tmp -= a[k]*w[k];

    // obtain output signal as output from feedforward part:
    double y = tmp*b[0];
    for(k=1; k<=order; k++)
      y += b[k]*w[k];

    // update state variables:
    memmove(&w[2], &w[1], (order-1)*sizeof(double));
    w[1] = tmp;

    return y;
  }

} 

#endif 
