#ifndef RAPT_DIRECTFORMFILTER_H_INCLUDED
#define RAPT_DIRECTFORMFILTER_H_INCLUDED

/** This is a generic filter using a Direct Form II implementation structure. 

This code needs some seriuos rework!

  (1) Use a fixed maxOrder - maybe 10 - that should be more than one should ever implement in 
      direct form

*/

template<class TSig, class TCoef>
class rsDirectFormFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor - the passed parameter will determine the maximum order which this filter will be
  able to realize.  */
  rsDirectFormFilter(int maximumOrder = 8);

  /** Destructor. */
  ~rsDirectFormFilter();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Srts up the filter coefficients to use. */
  void setCoefficients(TCoef* newCoeffsA, TCoef* newCoeffsB, int newOrder);
  // This is a bad API! It's easy for the caller to assume that the newOrder is the length of the 
  // arrays when their length is in fact, supposed to be order+1. In such a case we may read 
  // garbage from a memory location that we are not supposed to!

  /** Applies a global gain factor (by multiplying all b-coefficients with tha factor). */
  void setGlobalGainFactor(TCoef newFactor);

  /** Inverts the filter, i.e. swaps poles and zeros and inverts the gain. */
  void invert();


  // Under construction - needs tests:
  void reflectZeros();
  void reflectPoles();


  // ToDo: reflectPoles, reflectZeros

  /** Initializes the filter coefficients so as to realize a 'bypass' filter. */
  void initializeCoefficients();

  //-----------------------------------------------------------------------------------------------
  /** \name Inqiury */

  /** Returns the magnitude response at the given normalized radian frequency 'omega' where pi 
  corresponds to half the sample-rate. */
  TCoef getMagnitudeResponseAt(TCoef omega);

  /** Calculates the magnitudes of the frequency-response at the frequencies given in the array 
  "frequencies" (in Hz) and stores them in the array "magnitudes". Both arrays are assumed to be 
  "numBins" long. "inDecibels" indicates, if the frequency response should be returned in decibels. 
  If "accumulate" is true, the magnitude response of this biquad-cascade will be multiplied with 
  (or added to,when "inDecibels" is true) to the magnitudes which are already there in the 
  "magnitudes"-array. This is useful for calculating the magnitude response of several filters in
  series. */
  void getMagnitudeResponse(TCoef* frequencies, TCoef* magnitudes, int numBins, 
    TCoef sampleRate, bool inDecibels = false, bool accumulate = false);
  // This is also a bad API. These boolean parameters to switch behavior are annyoing and also lead
  // to a complicated implementation

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Calculates a single filtered output-sample. */
  inline TSig getSample(TSig in)
  {
    // obtain intermediate signal as output from feedback part:
    int k;
    TSig tmp = in;
    for(k = 1; k <= order; k++)
      tmp -= a[k]*w[k];

    // obtain output signal as output from feedforward part:
    TSig y = tmp*b[0];
    for(k = 1; k <= order; k++)
      y += b[k]*w[k];

    // update state variables:
    memmove(&w[2], &w[1], (order-1)*sizeof(TSig));
    w[1] = tmp;

    return y;
  }
  // ToDo: Provide getSample versions that apply only the poles or only the zeros.

  /** Resets the filter state. */
  void reset();



protected:

  TSig   *w;       // state buffer
  TCoef  *a, *b;   // filter coefficients
  int    order;    // currently used filter order - \todo maybe have separate orders for numerator and denominator
  int    maxOrder; // maximum possible filter order

  //static const int maxOrder = 8;  // For later
};

#endif 
