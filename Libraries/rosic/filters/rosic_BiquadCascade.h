#ifndef rosic_BiquadCascade_h
#define rosic_BiquadCascade_h

// rosic-indcludes:
#include "rosic_FilterAnalyzer.h"

namespace rosic
{

  /**

  This class implements a cascade of biquad-filter stages where each of the biquads implements the difference equation:

  y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

  Each stage has its own set of coefficients which has to be set from outside this class - this class does not do the filter-design. The 
  coefficients can be calculated by one of the "Designer" classes such as for example the BiquadDesigner class.

  */

  class BiquadCascade
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You may pass the maximum number of stages that this cascade will be able to realize here (default is 12). */
    BiquadCascade(int newMaxNumStages = 12);  
 
    /** Destructor. */
    ~BiquadCascade();  

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets a new sample-rate. */
    //void setSampleRate(double newSampleRate);

    /** Sets the number of biquad stages. */
    void setNumStages(int newNumStages);    

    /** Sets the order of the filter. The number of stages will we either half this value (if the order is even) or (order+1)/2 if the 
    order is odd. */
    void setOrder(int newOrder);    

    /** Sets up the global gain factor (by multiplying the feedforward coefficients of the first stage by the factor). */
    void setGlobalGainFactor(double newGainFactor);

    /** Copies the settings (numStages and the coefficients) from another instance of this class. */
    void copySettingsFrom(BiquadCascade *other);

    /** Turns this biquad-cascade into an allpass filter that has the same poles as the original filter. The zeros are moevd to positions
    that are reflections of the poles in the unit circle. */
    void turnIntoAllpass();

    /** Multiplies all feedforward ('b') coefficients by some factor. */
    //void multiplyFeedforwardCoeffsBy(double factor);

    /** Allows the user to set the filter coefficients for the individual biquad-stages. The difference-equation of each of the biquad 
    stages is: \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \f] */
    INLINE void setCoeffs(double *newB0, double *newB1, double *newB2, double *newA1, double *newA2); //, double newGain = 1.0);
                                       
    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the current number of stages. */
    int getNumStages() const { return numStages; }

    /** Returns the memory-address of the b0 array. */
    double* getAddressB0() const { return b0; }

    /** Returns the memory-address of the b1 array. */
    double* getAddressB1() const { return b1; }

    /** Returns the memory-address of the b2 array. */
    double* getAddressB2() const { return b2; }

    /** Returns the memory-address of the a1 array. */
    double* getAddressA1() const { return a1; }

    /** Returns the memory-address of the a2 array. */
    double* getAddressA2() const { return a2; }

    /** Returns the global gain factor. */
    //double getGlobalGainFactor() const { return gain; }

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array "frequencies" (in Hz) and stores them 
    in the array "magnitudes". Both arrays are assumed to be "numBins" long. "inDecibels" indicates, if the frequency response should be 
    returned in decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be multiplied with (or added to, 
    when "inDecibels" is true) to the magnitudes which are already there in the "magnitudes"-array. This is useful for calculating the 
    magnitude response of several biquad-cascades in series. */
    //void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, bool inDecibels = false, bool accumulate = false);

    /** Writes the complex frequency-response of a biquad-cascade at the normalized radian frequencies given in 'w' into the array 'H'. */
    void getFrequencyResponse(double *w, Complex *H, int numBins, int accumulationMode = FilterAnalyzer::NO_ACCUMULATION); 

    /** Writes the magnitdue-response of a biquad-cascade at the normalized radian frequencies given in 'w' into the array 
    'magnitudes'. */
    void getMagnitudeResponse(double *w, double *magnitudes, int numBins, bool inDecibels = false, bool accumulate = false);

    /** Writes the magnitdue-response of a biquad-cascade at the physical frequencies given in 'frequencies' into the array 
    'magnitudes'. */
    void getMagnitudeResponse(double *frequencies, double sampleRate, double *magnitudes, int numBins, bool inDecibels = false, 
      bool accumulate = false);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample via a cascade of biquads in Direct-Form 1 */
    INLINE double getSampleDirect1(double in); 

    /** Calculates a single filtered output-sample via a cascade of biquads in Direct-Form 2 */
    INLINE double getSampleDirect2(double in); 

    /** Calculates a single filtered output-sample via a cascade of biquads. */
    INLINE double getSample(double in) { return getSampleDirect1(in); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Initializes the biquad coefficients to neutral values. */
    void initBiquadCoeffs(); 

    /** Sets the buffers for the previous input and output doubles of all biquad stages to zero. */
    void reset(); 

    //=====================================================================================================================================

  protected:

    // filter coefficients:
    doubleA *a1;   // a1-coefficients for the individual stages
    doubleA *a2;   // a2-coefficients for the individual stages
    doubleA *b0;   // b0-coefficients for the individual stages
    doubleA *b1;   // b1-coefficients for the individual stages
    doubleA *b2;   // b2-coefficients for the individual stages

    // buffering:
    doubleA *x1;   // arrays of previous biquad input samples, array index indicates the biquad-stage
    doubleA *x2;   
    doubleA *y1;   // arrays of previous output samples
    doubleA *y2;

    // current number of biquad-stages:
    intA numStages;

    // maximum number of biquad-stages:
    intA maxNumStages;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void BiquadCascade::setCoeffs(double *newB0, double *newB1, double *newB2, double *newA1, double *newA2) //, double  newGain)
  {
    for(int i=0; i<numStages; i++)
    {
      b0[i] = newB0[i];
      b1[i] = newB1[i];
      b2[i] = newB2[i];
      a1[i] = newA1[i];
      a2[i] = newA2[i];
    }
  }

  INLINE double BiquadCascade::getSampleDirect1(double in)
  {
    doubleA tmp, tmp2;
    intA i;  // for the loop through the stages 

    tmp = in + TINY; // + TINY to avoid denorm problems

    // calculate current output-sample (y[n]) of all the BiQuad-stages (the output of one stage is the input for the next stage):
    for (i=0; i<numStages; i++)  
    {
      tmp2 = tmp; // for x1[i]

      // calculate current output-sample (y[n]) of BiQuad-stage i:
      tmp = b0[i]*tmp + (b1[i]*x1[i] + b2[i]*x2[i]) - (a1[i]*y1[i] + a2[i]*y2[i]);                      

      // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next iteration:
      x2[i] = x1[i];
      x1[i] = tmp2;
      y2[i] = y1[i];
      y1[i] = tmp;  
    }

    return tmp;
  }

  INLINE double BiquadCascade::getSampleDirect2(double in)
  {
    static doubleA x, y, g;
    static intA    i;  // for the loop through the stages

    x = in + TINY;

    // calculate current output-sample (y[n]) of all the BiQuad-stages (the output of one stage is the input for the next stage):
    for (i=0; i<numStages; i++)  
    {
      // calculate current output-sample (y[n]) of BiQuad-stage i:
      g = x - a1[i]*y1[i] - a2[i]*y2[i];
      y = b0[i]*g + b1[i]*y1[i] + b2[i]*y2[i];

      // set g[n-1], g[n-2] for the next iteration:
      y2[i] = y1[i];
      y1[i] = g;

      x = y; // output of one stage is input to the next
    }

    return y;
  }


  //=======================================================================================================================================

  /**

  This is the stereo version of the BiquadCascade class.

  */

  class BiquadCascadeStereo : public BiquadCascade
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You may pass the maximum number of stages that this cascade will be able to realize here (default is 12). */
    BiquadCascadeStereo(int newMaxNumStages = 12);  
 
    /** Destructor. */
    ~BiquadCascadeStereo();  

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a single filtered output-sample frame via a cascade of biquads in Direct-Form 1 */
    INLINE void getSampleFrameDirect1(double *inOutL, double *inOutR); 

    /** Calculates a single filtered output-sample frame via a cascade of biquads in Direct-Form 2 */
    INLINE void getSampleFrameDirect2(double *inOutL, double *inOutR); 

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Sets the buffers for the previous input and output doubles of all biquad stages to zero. */
    void reset(); 

  protected:

    // additional buffering variables for the right channel:
    doubleA *x1R;  
    doubleA *x2R;   
    doubleA *y1R;   
    doubleA *y2R;

  };

  INLINE void BiquadCascadeStereo::getSampleFrameDirect1(double *inOutL, double *inOutR)
  {
    doubleA tmpL = *inOutL + TINY;
    doubleA tmpR = *inOutR + TINY;

    doubleA tmpL2, tmpR2;

    for(int i=0; i<numStages; i++)  
    {
      tmpL2  = tmpL;
      tmpR2  = tmpR;

      tmpL   = b0[i]*tmpL + (b1[i]*x1[i]  + b2[i]*x2[i])  - (a1[i]*y1[i]  + a2[i]*y2[i]);                      
      tmpR   = b0[i]*tmpR + (b1[i]*x1R[i] + b2[i]*x2R[i]) - (a1[i]*y1R[i] + a2[i]*y2R[i]);  

      x2[i]  = x1[i];
      x1[i]  = tmpL2;
      y2[i]  = y1[i];
      y1[i]  = tmpL;  

      x2R[i] = x1R[i];
      x1R[i] = tmpR2;
      y2R[i] = y1R[i];
      y1R[i] = tmpR;  
    }

    *inOutL = tmpL;
    *inOutR = tmpR;
  }

  INLINE void BiquadCascadeStereo::getSampleFrameDirect2(double *inOutL, double *inOutR)
  {
    doubleA xL, yL, gL, xR, yR, gR;

    xL = *inOutL + TINY;
    xR = *inOutR + TINY;

    for(int i=0; i<numStages; i++)  
    {
      gL     = xL - a1[i]*y1[i]  - a2[i]*y2[i];
      gR     = xR - a1[i]*y1R[i] - a2[i]*y2R[i];

      yL     = b0[i]*gL + b1[i]*y1[i]  + b2[i]*y2[i];
      yR     = b0[i]*gR + b1[i]*y1R[i] + b2[i]*y2R[i];

      y2[i]  = y1[i];
      y1[i]  = gL;

      y2R[i] = y1R[i];
      y1R[i] = gR;

      xL     = yL; 
      xR     = yR; 
    }

    *inOutL = yL;
    *inOutR = yR;
  }

}  

#endif
