#ifndef rosic_ScopeScreenScanner_h
#define rosic_ScopeScreenScanner_h

//// rosic-includes:
//#include "rosic_PitchDetector.h"

namespace rosic
{

  /** This is a class for generating the sawtooth-shaped waveform used for scanning over the screen
  horizontally in an oscilloscope in 1D mode. It provides synchronization with the incoming 
  waveform by using a pitch-detector.  */

  class ScopeScreenScanner 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ScopeScreenScanner();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the frequency that shoould be used when we are not in sync mode. */
    void setScanFreqNoSync(double newFrequency);

    /** Sets the minimum expected frequency for the pitch detector. */
    void setMinFrequency(double newMinFreq);

    /** Sets the maximum expected frequency for the pitch detector. */
    void setMaxFrequency(double newMaxFreq);

    /** Switches synchronization of the sawtooth to the input signal on/off */
    void setSync(bool shouldSync);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Generates one sawtooth output sample at the time. You must pass the input signal value that
    is used for the pitch analysis. The value is between 0 and 1. */
    INLINE double getSample(double in);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Updates our phase increment variable for the sawtooth generator according to the 
    parameters and (possibly) the imput signal. */
    void updateSawIncrement(double in);

    /** Resets the internal state. */
    void reset();

    //=============================================================================================

  protected:

    doubleA sampleRate, scanFreq, numCyclesShown;
    doubleA sawPhase, sawInc;  // sawtooth phase and increment

    bool sync;

    PitchDetector pitchDetector;

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double ScopeScreenScanner::getSample(double in)
  {
    updateSawIncrement(in);
    double result = sawPhase;
    sawPhase += sawInc;
    if(sawPhase > 1)
      sawPhase -= 1;
    return result;
    //return 0.0; // preliminary
  }

} // end namespace rosic

#endif // ScopeScreenScanner
