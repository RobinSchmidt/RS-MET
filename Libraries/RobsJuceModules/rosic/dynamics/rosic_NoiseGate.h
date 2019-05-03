#ifndef rosic_NoiseGate_h
#define rosic_NoiseGate_h

//// rosic-indcludes:
//#include "rosic_DynamicsProcessorBase.h"

namespace rosic
{

  /**

  This class implements a noise-gate.

  */

  class NoiseGate : public DynamicsProcessorBase
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - constructs a dynamics processor with a given maximum number of samples 
    lookahead. */
    NoiseGate(int newLookAheadBufferSize = DEFAULT_LOOKAHEAD_SIZE);

    /** Destructor */
    ~NoiseGate();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets up the samplerate. */
    void setSampleRate(double newSampleRate);

    /** Sets the attack time of the envelope-follower in milliseconds. */
    void setAttackTime(double newAttackTime)
    { attackReleaseEnveloper.setAttackTime(newAttackTime); }

    /** Sets the threshold above which the gate opens (in decibels). */
    void setThreshold(double newThreshold);

    /** Sets the level difference between the opening- and closing-threshold in decibels. */
    void setHysteresis(double newHysteresis);

    /** Sets the minimum time (in ms) for the gate to stay open after the opening threshold was
    exceeded. */
    void setHoldTime(double newHoldTime);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

  protected:

    /** Updates the opening- and closing-threshold (which are internal parameters) according to the 
    user-parameters 'threshold' and 'hysteresis'. */
    void updateThresholds();

    double openingThreshold;  // threshold for opening the gate in linear amplitude units
    double closingThreshold;  // threshold for closing the gate in linear amplitude units
    double threshold;         // threshold for opening the gate in dB
    double hysteresis;        // difference between opening- and closing-threshold in dB
    double holdTime;          // minimum time for keeping the gate open in ms
    int    numHoldSamples;    // number of samples to hold the gate open
    int    holdCounter;       // sample counter for measuring the time since the last 
                              // opening-threshold exceeding
    bool   gateClosed;        // flag to indicate that the gate is closed (or closing)

  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void NoiseGate::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // establish the mono sum of the (non-delayed) input signal and obtain its envelope:
    double tmpL = inputGainFactor * (*inOutL);
    double tmpR = inputGainFactor * (*inOutR);  
    double inM  = RAPT::rsMax(fabs(tmpL), fabs(tmpR));
    double env  = levelDetector.attackReleaseFollower.getSample(inM);

    // evaluate the opening-/closing-logic:
    double gateGain;
    if( env >= openingThreshold )
    {
      gateGain    = 1.0;
      gateClosed  = false;
      holdCounter = 0;
    }
    else if( env < closingThreshold )
    {
      if( holdCounter >= numHoldSamples || gateClosed )
      {
        gateGain    = 0.0;
        gateClosed  = true;
        holdCounter = 0;
      }
      else
      {
        gateGain    = 1.0;
        gateClosed  = false;
        holdCounter = RAPT::rsMin(++holdCounter, numHoldSamples); // don't let it grow indefinitely
      }
    }
    else  // env is in between opening- and closing threshold
    {
      gateGain = (double) (!gateClosed);  // 0.0 if closed, 1.0 if open
      holdCounter = 0;
    }

    // apply the attack/release smoothing:
    gateGain  = attackReleaseEnveloper.getSample(gateGain);
    gateGain *= outputGainFactor;

    // apply delay, gain an mix dry/wet (this can be almalgameted into a single operation):
    applyLookAheadDelay(&tmpL, &tmpR, &tmpL, &tmpR);
    tmpL    *= (dry + wet*gateGain);
    tmpR    *= (dry + wet*gateGain);
    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_NoiseGate_h
