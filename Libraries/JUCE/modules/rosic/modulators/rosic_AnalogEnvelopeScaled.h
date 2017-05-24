#ifndef rosic_AnalogEnvelopeScaled_h
#define rosic_AnalogEnvelopeScaled_h

//// rosic-indcludes:
//#include "rosic_AnalogEnvelope.h"

namespace rosic
{

  /**

  This class is a modified AnalogEnvelope that reaches its target levels in the prescribed
  times exactly by appropriately scaling and offsetting the exponential decay functions that make 
  up the segments of the envelope. The implementation is not based on the digital RC equivalent 
  anymore. Instead, it generates normalized exponential decay functions by means of a 
  multiplicative accumulator and then scales and offsets the result.

  */

  class AnalogEnvelopeScaled : public AnalogEnvelope
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    AnalogEnvelopeScaled();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the length of the attack phase (in milliseconds). */
    void setAttack(double newAttackTime);    

    /** Sets the hold time (in milliseconds). */
    void setHold(double newHoldTime);      

    /** Sets the length of the decay phase (in milliseconds). */
    void setDecay(double newDecayTime);     
 
    /** Sets the length of the release phase (in milliseconds). */
    void setRelease(double newReleaseTime);  

    /** Scales the A,D,H and R times by adjusting the increment. It is 1 if not used - a timescale 
    of 2 means the envelope is twice as fast, 0.5 means half as fast -> useful for implementing a 
    key/velocity-tracking feature for the overall length for the envelope. */
    void setTimeScale(double newTimeScale); 


    void setShape(double newShape);

    //---------------------------------------------------------------------------------------------
    // inquiry:
    /** True, if the end point is reached. */
    bool endIsReached();  

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output sample at a time. */
    INLINE double getSample();    

    //---------------------------------------------------------------------------------------------
    // others:

    /** Causes the envelope to start with its attack-phase. When the parameter 
    'startFromCurrentValue' is true, the internal state will not be reset to startLevel, such that 
    the curve begins at the level, where the envelope currently is. */  
    void noteOn(bool startFromCurrentLevel = false, int newKey = 64, int newVel = 64);

    /** Causes the envelope to start with its release-phase. */
    void noteOff();  

    /** Resets the time variable. */
    void reset();   

  protected:

    double attackAccu,    decayAccu,    releaseAccu;
    double attackAccuMin, decayAccuMin, releaseAccuMin;  // redundant? they always are equal to mu?
    double attackAccuMax, decayAccuMax, releaseAccuMax;
    double mu, shape;
    int    attackSamples, holdSamples, decaySamples, releaseSamples, sampleCounter;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double AnalogEnvelopeScaled::getSample()
  {
    double out;

    // attack phase:
    if(sampleCounter < attackSamples)   // noteIsOn has not to be checked, because, time is advanced to the 
                                         // beginning of the release phase in noteOff()
    {
      out         =  peakScale*peakLevel - attackCoeff*(attackAccu-attackAccuMin);
      attackAccu *= attackAccuMax;
      sampleCounter++;
    }
    // hold phase:
    else if(sampleCounter < attackSamples+holdSamples)
    {
      out = peakScale*peakLevel;
      sampleCounter++;
    }
    // decay phase:
    else if(sampleCounter < attackSamples+holdSamples+decaySamples)
    {
      out        =  sustainLevel - decayCoeff*(decayAccu-decayAccuMin);
      decayAccu *= decayAccuMax;
      sampleCounter++;
    }
    // sustain phase:
    else if(noteIsOn)
    {
      out = sustainLevel;
      // time is not incremented in sustain
    }
    // release phase:
    else if(sampleCounter < attackSamples+holdSamples+decaySamples+releaseSamples)
    {
      out          =  endLevel - releaseCoeff*(releaseAccu-releaseAccuMin);
      releaseAccu *= releaseAccuMax;
      sampleCounter++;
    }
    else
      out = endLevel;

    previousOutput = out; // needed to start attack or release from the current level
    return out;
  }

} // end namespace rosic

#endif // rosic_AnalogEnvelopeScaled_h
