#ifndef rosic_AllpassChain_h
#define rosic_AllpassChain_h

//// rosic-indcludes:
//#include "rosic_BiquadDesigner.h"

namespace rosic
{

  /**

  This class implements a cascade connection of a number of up to 24 first or second order allpass 
  filters to be used for a phaser effect (for example). 

  */

  class AllpassChain 
  {

  public:

    /** Enumeration of the available filter modes. */
    enum modes
    {
      //BYPASS = 0,
      FIRST_ORDER_ALLPASS = 0,
      SECOND_ORDER_ALLPASS,

      NUM_FILTER_MODES
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AllpassChain();

    /** Destructor. */
    ~AllpassChain();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate (in Hz) at which the filter runs. */
    void setSampleRate(double newSampleRate, bool updateCoefficients = true);

    /** Sets the filter mode as one of the values in enum modes. */
    void setMode(int newMode, bool updateCoefficients = true);

    /** Sets the center frequency in Hz. */
    void setFrequency(double newFrequency, bool updateCoefficients = true);

    /** Sets the Q-value for second order allpasses. */
    void setQ(double newBandwidth, bool updateCoefficients = true);

    /** Sets the number of allpass stages to use. */
    void setNumStages(int newNumStages);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns true if the current mode supports a Q parameter, false otherwise. */
    bool doesModeSupportQ() const { return mode == SECOND_ORDER_ALLPASS; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal state buffers. */
    void reset();

    /** Calculates the biquad coefficients from the specification. */
    void updateCoeffs();

    //=============================================================================================

  protected:

    static const int maxNumStages = 24;

    // filter coefficients:
    double b0, b1, b2, a1, a2;

    // buffering:
    double x1[maxNumStages];
    double x2[maxNumStages];
    double y1[maxNumStages];
    double y2[maxNumStages];

    // user parameters:
    double frequency, q, oneOverSampleRate;
    double sampleRate;
    int    mode, numStages;

  };

  INLINE double AllpassChain::getSample(double in)
  {
    double tmpY = in;
    double tmpX = in;
    for(int i=0; i<numStages; i++)
    {
      //tmpY   = (b0*tmpX+TINY) + (b1*x1[i] + b2*x2[i]) + (a1*y1[i] + a2*y2[i]); // old
      tmpY   = (b0*tmpX) + (b1*x1[i] + b2*x2[i]) + (a1*y1[i] + a2*y2[i]);  // new
      x2[i]  = x1[i];
      x1[i]  = tmpX;
      y2[i]  = y1[i];
      y1[i]  = tmpY;
      tmpX   = tmpY;
    }
    return tmpY;
  }

} // end namespace rosic

#endif // rosic_AllpassChain_h
