#ifndef rosic_QuadrigaVoice_h
#define rosic_QuadrigaVoice_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_PolyphonicInstrumentVoice.h"
//#include "../modulators/rosic_EnvelopeGenerator.h"
//#include "../generators/rosic_Quadrigen.h"

namespace rosic
{
  /**

  This class is the base class a single voice for a simple sampler.

  */

  class QuadrigaVoice : public PolyphonicInstrumentVoice
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    QuadrigaVoice();

    /** Destructor. */
    virtual ~QuadrigaVoice();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this voice. */
    virtual void setSampleRate(double newSampleRate);

    /** Sets up the waveform for one of the oscillator (1-3). This is expected to be a stereo
    sample with the first index indicating the channel and the second the sample number. */
    //virtual void setWaveform(double **newWaveform, int newLength, int whichOscillator);

    /** Overrides PolyphonicInstrumentVoice::setMasterTuneA4. */
    virtual void setMasterTuneA4(double newTuneA4);

    /** Overrides PolyphonicInstrumentVoice::setBeatsPerMinute() in order to set up the
    modulation generators. */
    virtual void setBeatsPerMinute(double newBpm);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** This is not to be used - we need to formally implement it in order to make this class
    non-abtract. */
    virtual void getSampleFrameStereo(double *outL, double *outR, double *voiceAmplitude);

    /** Calculates the output-samples for both channels at an oversampled rate by factor 2. The
    *outL1 and outL2 slots represent two subsequent (oversampled) samples for the left channel,
    likewise the *outR1, *outR2 slots for the right channel. The decimation to the target
    sample-rate is to be done by the outlying class in order to do that only once after all voices
    have been added (instead of doing it per voice). Note that the outputs do not overwrite the
    slots, but add the voice's output to what is already there. The fifth output tells the current
    amplitude of the voice (i.e. the output of the amp-envelope.) to enable the
    PolyphonicInstrument class to apply automatic volume scaling acording to the number of  playing
    voices and their loudnesses. */
    virtual void getSampleFrameStereo(double *outL1, double *outR1, double *outL2, double *outR2,
      double *voiceAmplitude);

    //---------------------------------------------------------------------------------------------
    // event processing:

    //virtual void noteOn(int newNoteNumber, int newVelocity, int newDetune = 0);

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    Quadrigen          quadrigen;
    EnvelopeGenerator  ampEnv;

  protected:

    /** Overrides triggerNote() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerNote(int newKey, int newVelocity, int newDetune = 0);

    /** Overrides glideToNote() inherited from class PolyphonicInstrumentVoice. */
    virtual void glideToNote(int newKey, int newVelocity, int newDetune = 0);

    /** Overrides triggerRelease() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerRelease(int noteToBeReleased, int noteToBeReleasedVel);

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void QuadrigaVoice::getSampleFrameStereo(double* outL, double* outR,
    double* /*voiceAmplitude*/)
  {
    *outL = 0.0;
    *outR = 0.0;
  }

  INLINE void QuadrigaVoice::getSampleFrameStereo(double *outL1, double *outR1,
    double *outL2, double *outR2, double *voiceAmplitude)
  {
    // this call just sets up our currentFrequencyWithPitchBend and currentNoteAge members:
    PolyphonicInstrumentVoice::getSampleFrameStereo(outL1, outR1, voiceAmplitude);

    double ampEnvOut = ampEnv.getSample();

    double tmpL, tmpR;
    quadrigen.getSampleFrameStereo(&tmpL, &tmpR);

    // accumulate the ouput of this voice to what's already there (from other voices):
    *outL1 += currentAmplitude * ampEnvOut * tmpL;
    *outR1 += currentAmplitude * ampEnvOut * tmpR;

    // again due to 2x oversampling:
    quadrigen.getSampleFrameStereo(&tmpL, &tmpR);
    *outL2 += currentAmplitude * ampEnvOut * tmpL;
    *outR2 += currentAmplitude * ampEnvOut * tmpR;

    // add the amplitude of this to the voiceAmplitude slot:
    *voiceAmplitude += ampEnvOut*ampEnvOut * currentAmplitude*currentAmplitude;

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

#endif //  rosic_QuadrigaVoice_h
