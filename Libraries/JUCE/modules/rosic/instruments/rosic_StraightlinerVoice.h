#ifndef rosic_StraightlinerVoice_h
#define rosic_StraightlinerVoice_h

//// rosic-indcludes:
//#include "../infrastructure/rosic_PolyphonicInstrumentVoice.h"
//#include "../generators/rosic_FourOscSection.h"
//#include "../filters/rosic_MultiModeFilter.h"
//#include "../filters/rosic_LeakyIntegrator.h"
//#include "../modulators/rosic_BreakpointModulator.h"

namespace rosic
{
  /**

  This class is a single voice for The Straightliner synthesizer.

  */

  class StraightlinerVoice : public PolyphonicInstrumentVoice
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    StraightlinerVoice();

    /** Destructor. */
    ~StraightlinerVoice();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this voice. */
    virtual void setSampleRate(double newSampleRate); 

    /** Sets up the waveform for one of the oscillator (1-3). This is expected to be a stereo sample with the first index indicating the 
    channel and the second the sample number. */
    //virtual void setWaveform(double **newWaveform, int newLength, int whichOscillator);

    /** Overrides PolyphonicInstrumentVoice::setMasterTuneA4. */
    virtual void setMasterTuneA4(double newTuneA4);

    /** Overrides PolyphonicInstrumentVoice::setBeatsPerMinute() in order to set up the modulation generators. */
    virtual void setBeatsPerMinute(double newBpm);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // synchronization of different voices:

    /** Implements the purely virtual method inherited from PolyphonicInstrumentVoice. */
    //virtual void copyCommonVoiceData(PolyphonicInstrumentVoice* voiceToCopyFrom);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Sets up the frequency for one of the oscillators, including the pitch envelope. */
    INLINE void setOscillatorFrequencyWithEnv(OscillatorStereo& osc, const double& envOut);

    /** This is not to be used - we need to formally implement it in order to make this class non-abtract. */
    virtual void getSampleFrameStereo(double *outL, double *outR, double *voiceAmplitude);  

    /** Calculates the output-samples for both channels at an oversampled rate by factor 2. The *outL1 and outL2 slots represent two 
    subsequent (oversampled) samples for the left channel, likewise the *outR1, *outR2 slots for the right channel. The decimation to the 
    target sample-rate is to be done by the outlying class in order to do that only once after all voices have been added (instead of doing 
    it per voice). Note that the outputs do not overwrite the slots, but add the voice's output to what is already there. The fifth output 
    tells the current amplitude of the voice (i.e. the output of the amp-envelope.) to enable the PolyphonicInstrument class to apply 
    automatic volume scaling acording to the number of  playing voices and their loudnesses. */
    virtual void getSampleFrameStereo(double *outL1, double *outR1, double *outL2, double *outR2, double *voiceAmplitude);   

    //-------------------------------------------------------------------------------------------------------------------------------------
    // event processing:

    //virtual void noteOn(int newNoteNumber, int newVelocity, int newDetune = 0);



    //-------------------------------------------------------------------------------------------------------------------------------------
    // embedded objects:

    //OscillatorStereo    osc1, osc2, osc3, osc4;
    FourOscSection      oscSection;
    MultiModeFilter     filter;
    BreakpointModulator pitchEnv, filterEnv, ampEnv;
    LeakyIntegrator     cutoffSmoother;

  protected:

    /** Overrides triggerNote() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerNote(int newKey, int newVelocity, int newDetune = 0);

    /** Overrides glideToNote() inherited from class PolyphonicInstrumentVoice. */
    virtual void glideToNote(int newKey, int newVelocity, int newDetune = 0);

    /** Overrides triggerRelease() inherited from class PolyphonicInstrumentVoice. */
    virtual void triggerRelease(int noteToBeReleased, int noteToBeReleasedVel);

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void StraightlinerVoice::setOscillatorFrequencyWithEnv(OscillatorStereo& osc, const double& envOut)
  {
    double depth = osc.getPitchEnvelopeDepth();
    if( envOut > 0.0 )
    {
      if( depth == 0.0 )
        osc.setFrequencyNominal(currentFrequencyWithPitchBend);
      else if( depth == 1.0 )
        osc.setFrequencyNominal(currentFrequencyWithPitchBend*envOut);
      else if( depth == -1.0 )
        osc.setFrequencyNominal(currentFrequencyWithPitchBend/envOut);
      else 
        osc.setFrequencyNominal(currentFrequencyWithPitchBend*pow(envOut, depth));
    }
    else
    {
      if( depth != 0.0 )
        osc.setFrequencyNominal(0.0);
      else
        osc.setFrequencyNominal(currentFrequencyWithPitchBend);
    }
  }

  INLINE void StraightlinerVoice::getSampleFrameStereo(double *outL, double *outR, double *voiceAmplitude)
  {
    *outL = 0.0;
    *outR = 0.0;
  }

  INLINE void StraightlinerVoice::getSampleFrameStereo(double *outL1, double *outR1, double *outL2, double *outR2, double *voiceAmplitude)
  {
    PolyphonicInstrumentVoice::getSampleFrameStereo(outL1, outR1, voiceAmplitude);

    // calculate modulators:
    double pitchEnvOut = pitchEnv.getSample();
    double ampEnvOut   = ampEnv.getSample();
    double fltEnvOut   = filterEnv.getSample();

    // set up the oscillator frequencies:
    setOscillatorFrequencyWithEnv(oscSection.osc1, pitchEnvOut);
    setOscillatorFrequencyWithEnv(oscSection.osc2, pitchEnvOut);
    setOscillatorFrequencyWithEnv(oscSection.osc3, pitchEnvOut);
    setOscillatorFrequencyWithEnv(oscSection.osc4, pitchEnvOut);

    // set up the filter:
    if( glideIsActive && remainingGlideSamples >= 0)
      filter.setKey(currentPitchWithPitchBend);
    double freq = cutoffSmoother.getSample(filter.getFrequencyWithKeyAndVel());  
    freq        = freq * fltEnvOut;
    filter.setFrequencyInstantaneous(freq);
    
    double osc1L, osc1R, osc2L, osc2R, osc3L, osc3R, osc4L, osc4R, tmpL, tmpR;

    //-----------------------------------------------------
    // oversampled calculations:

    // calculate the outputs signals of the 3 oscillators:
    oscSection.osc1.getSampleFrameStereo(&osc1L, &osc1R);
    oscSection.osc2.getSampleFrameStereo(&osc2L, &osc2R);
    oscSection.osc3.getSampleFrameStereo(&osc3L, &osc3R);
    oscSection.osc4.getSampleFrameStereo(&osc4L, &osc4R);

    // add the oscillator signals:
    tmpL = osc1L + osc2L + osc3L + osc4L;
    tmpR = osc1R + osc2R + osc3R + osc4R;

    // apply the filter:
    filter.getSampleFrameStereo(&tmpL, &tmpR, &tmpL, &tmpR);

    // accumulate the ouput of this voice to what's already there (from other voices):
    *outL1 += currentAmplitude * ampEnvOut * tmpL;
    *outR1 += currentAmplitude * ampEnvOut * tmpR;

    // (almost) the same again because of 2x oversampling:
    oscSection.osc1.getSampleFrameStereo(&osc1L, &osc1R);
    oscSection.osc2.getSampleFrameStereo(&osc2L, &osc2R);
    oscSection.osc3.getSampleFrameStereo(&osc3L, &osc3R);
    oscSection.osc4.getSampleFrameStereo(&osc4L, &osc4R);
    tmpL = osc1L + osc2L + osc3L + osc4L;
    tmpR = osc1R + osc2R + osc3R + osc4R;
    filter.getSampleFrameStereo(&tmpL, &tmpR, &tmpL, &tmpR);
    *outL2 += currentAmplitude * ampEnvOut * tmpL;
    *outR2 += currentAmplitude * ampEnvOut * tmpR;

    // end of oversampled calculations
    //-----------------------------------------------------


    // add the amplitude of this to the voiceAmplitude slot:
    *voiceAmplitude += ampEnvOut*ampEnvOut * currentAmplitude*currentAmplitude;

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

#endif //  rosic_StraightlinerVoice_h
