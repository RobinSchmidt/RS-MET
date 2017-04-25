#ifndef rosic_KeyShotVoice_h
#define rosic_KeyShotVoice_h

// rosic-indcludes:
#include "../infrastructure/rosic_PolyphonicInstrumentVoice.h"
//#include "rosic_KeyShotOscSection.h"
#include "../generators/rosic_SamplePlayer.h"
#include "../filters/rosic_MultiModeFilter.h"
#include "../modulators/rosic_BreakpointModulator.h"

namespace rosic
{
  /**

  This class is the base class a single voice for a simple sampler.

  */

  class KeyShotVoice : public PolyphonicInstrumentVoice
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    KeyShotVoice();

    /** Destructor. */
    virtual ~KeyShotVoice();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this voice. */
    virtual void setSampleRate(double newSampleRate); 

    /** Sets up the waveform for one of the oscillator (1-3). This is expected to be a stereo 
    sample with the first index indicating the channel and the second the sample number. */
    //virtual void setWaveform(double **newWaveform, int newLength, int whichOscillator);

    /** Overrides PolyphonicInstrumentVoice::setMasterTuneA4. */
    virtual void setMasterTuneA4(double newTuneA4);

    /** Overrides PolyphonicInstrumentVoice::setBeatsPerMinute() in order to set up the modulation 
    generators. */
    virtual void setBeatsPerMinute(double newBpm);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one otput sample frame at a time */
    virtual void getSampleFrameStereo(double *outL, double *outR, double *voiceAmplitude);  

    //---------------------------------------------------------------------------------------------
    // event processing:

    //virtual void noteOn(int newNoteNumber, int newVelocity, int newDetune = 0);

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    SamplePlayer        samplePlayer;
    BreakpointModulator ampEnv;

  protected:

    /** Overrides triggerNote() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerNote(int newKey, int newVelocity, int newDetune = 0);

    /** Overrides triggerRelease() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerRelease(int noteToBeReleased, int noteToBeReleasedVel);

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void KeyShotVoice::getSampleFrameStereo(double *outL, double *outR, 
    double *voiceAmplitude)
  {
    // this call just sets up our currentFrequencyWithPitchBend and currentNoteAge members:
    PolyphonicInstrumentVoice::getSampleFrameStereo(outL, outR, voiceAmplitude);

    // set up the playback frequency of the sample player (nominal means without detuning here):
    samplePlayer.setPlaybackFrequencyNominal(currentFrequencyWithPitchBend);
      // may be scrapped?

    double tmpL, tmpR;
    double ampEnvOut = ampEnv.getSample();
    samplePlayer.getSampleFrameStereo(&tmpL, &tmpR);
    tmpL *= ampEnvOut;
    tmpR *= ampEnvOut;

    // increment the note age:
    currentNoteAge++;

    // set the inherited isSilent flag to true, when the amp-env has reached its end:
    if( ampEnv.endIsReached )
    {
      isSilent          = true;
      isReleasing       = false;
      noteBeingReleased = -1;
    }
  }

} // end namespace rosic

#endif //  rosic_KeyShotVoice_h