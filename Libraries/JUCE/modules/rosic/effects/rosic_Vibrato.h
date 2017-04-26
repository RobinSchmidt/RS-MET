#ifndef rosic_Vibrato_h
#define rosic_Vibrato_h

// rosic-indcludes:
#include "rosic_ModulationEffect.h"
#include "../delaylines/rosic_DelayLineStereo.h"

namespace rosic
{

  /**

  This class implements a vibrato effect with sinusoidal modulators for left and right channel
  separately.

  */

  class Vibrato : public ModulationEffect, public DelayLineStereo
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor - will allocate a delay-buffer with a given maximum number of samples delay. */
    Vibrato(int bufferLengthToAllocate = 65536);

    /** Destructor */
    ~Vibrato();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the cycle-length in seconds or beats (depending on whether sync is active). */
    void setCycleLength(double newCycleLength);

    /** Sets the depth of the effect in semitones. */
    void setDepth(double newDepth);

    /** Switches the tempo-sync on or off. */
    void setTempoSync(bool shouldTempoSync);

    /** Sets up the tempo in  beats per minute. */
    void setTempoInBPM(double newTempoInBPM);

    /** Sets the average delay-time in milliseconds. */
    void setAverageDelayTime(double newAverageDelayTime);

    /** Sets the ratio between dry and wet between 0...1. */
    void setDryWetRatio(double newDryWet) { dryWetRatio = newDryWet; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the ratio between dry and wet between 0...1. */
    double getDryWetRatio() const { return dryWetRatio; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal buffers to zero. */
    void reset();

    //=============================================================================================

  protected:

    /** Calculates the depth of the modulation of the delaytime in samples (from the desired 
    vibrato-depth in semitones and the frequency of modulating sinusoid. */
    void calculateDepthInSamples();

    double  d;                   // modulation depth in samples
    int     dA;                  // average delay in samples
    double  dryWetRatio;         // dry/wet ratio 0....1
    double  averageDelay;        // average delay in ms      

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Vibrato::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    feedIn(*inOutL, *inOutR);

    double left, right;
    lfo.getSampleFrameStereo(&left, &right);
    double dL = dA + d*left;    // instantaneous delay for the left tap
    double dR = dA + d*right;   // instantaneous delay for the right tap

    double yL = getLeftOutputHermiteAt(dL);
    double yR = getRightOutputHermiteAt(dR);

    // read out the dry signal from the delaylines for mixing dry/wet without additional 
    // comb-filter artifacts:
    int readPos = wrapAround(tapIn-dA+1, length); 
    double dryL = bufferL[readPos];
    double dryR = bufferR[readPos];

    *inOutL = (1.0-dryWetRatio)*dryL + dryWetRatio*yL;
    *inOutR = (1.0-dryWetRatio)*dryR + dryWetRatio*yR;

    incrementWritePointer();
  }

} // end namespace rosic

#endif // #ifndef rosic_Vibrato_h
