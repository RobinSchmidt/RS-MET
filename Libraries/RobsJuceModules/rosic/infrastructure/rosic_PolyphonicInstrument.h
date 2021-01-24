#ifndef rosic_PolyphonicInstrument_h
#define rosic_PolyphonicInstrument_h


namespace rosic
{

/** PolyphonicInstrument is a base class from which polyphonic instruments can be derived. One 
should derive a single voice for an instrument from the PolyphonicIntrumentVoice class and the 
whole polyphonic instrument then from this class here. These two classes will then handle the 
voice management.

This is used in Straightliner, but i actually think, it should be used anymore for new instruments.
There are better ways of doing this...

todo: introduce public flag allVoicesAreSilent
-> poll it in PolyphonicInstrumentsAudioModule process
-> return directly if the flag is true and ther are no further midi-events for the current buffer
-> removes overhead for instances that just sit there */

class PolyphonicInstrument //: public PresetRememberer
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. Sets the variable numAllocatedVoices - however: the actual allocation of the 
  voices MUST be done in your subclass! */
  PolyphonicInstrument(int numVoicesToAllocate = 17);

  /** Destructor. */
  virtual ~PolyphonicInstrument();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  /** Sets the sample-rate. */
  virtual void setSampleRate(double newSampleRate);

  /** Sets the maximum number of playable voices (polyphony) - the zeroth voice is supposed to 
  serve as master/template voice and is not playable. */
  virtual void setNumPlayableVoices(int newNumVoices);

   /** Sets the master level decibels. */
  virtual void setMasterLevel(double newMasterLevel);

  /** Sets the key dependence of the voice level in dB at 127. */
  virtual void setVoiceLevelByKey(double newVoiceLevelByKey);

  /** Sets the velocity dependence of the voice level in dB at 127. */
  virtual void setVoiceLevelByVel(double newVoiceLevelByVel);

  /** Sets up an automatic scaling of the master level by the number of playing voices in order to
  decouple the overall loudness from that. */
  virtual void setMasterLevelByVoices(double newMasterLevelByVoices);

  /** Sets the ratio between mid- and side signal (range: 0...1). */
  virtual void setMidSideRatio(double newMidSideRatio);

  /** Sets the tuning frequency for the A4 reference tone (usually 440 Hz). */
  virtual void setMasterTuneA4(double newTuneA4);

  /** Sets the range for the pitch wheel in semitones. */
  virtual void setPitchWheelRange(double newRange);

  /** Sets the time for the glide in milliseconds. */
  virtual void setGlideTime(double newGlideTime);

  /** Switches the note-glide function on or off - when on, this object will call glideTo() 
  instead of noteOn() for the the stolen voice. The glideTo()-function (virtual memeber of the 
  PolyphonicInstrumentVoice class) is supposed to cause a glide to a new pitch instead of
  re-triggering the voice. When number of voices is 1, this works as the usual glide in 
  monophonic synths. */
  virtual void setGlideMode(bool shouldUseGlide);

  /** Calls setBeatsPerMinute for all the voices - in PolyphonicInstrumentVoice, this function is 
  empty, you will need to override it in your subclass if you need sync functionality. */
  virtual void setBeatsPerMinute(double newBpm);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Sets the maximum number of playable voices (polyphony). */
  virtual int getNumPlayableVoices();

  /** Returns the level of the voice in decibels. */
  virtual double getMasterLevel();

  /** Returns the key dependence of the voice level in dB at 127. */
  virtual double getVoiceLevelByKey();

  /** Returns the velocity dependence of the voice level in dB at 127. */
  virtual double getVoiceLevelByVel();

  /** Returns the automatic scaling of the master level by the number of playing voices. */
  virtual double getMasterLevelByVoices();

  /** Returns the ratio between mid- and side signal (range: 0...1). */
  virtual double getMidSideRatio();

  /** Returns the tuning frequency for the A4 reference tone (usually 440 Hz). */
  virtual double getMasterTuneA4();

  /** Returns the range for the pitch wheel in semitones. */
  virtual double getPitchWheelRange();

  /** Returns the glide time in milliseconds. */
  virtual double getGlideTime();

  /** Informs, whether glide-mode is on or off. */
  virtual bool isInGlideMode();

  /** Returns the number of currently playing voices (polyphony). */
  virtual int getNumActiveVoices();

  /** Returns true if the whole instrument is silent (i.e. all voices are silent). */
  virtual bool isSilent();

  //-----------------------------------------------------------------------------------------------
  // event processing:

  /** Triggers a note. Depending on the current setting of polyphony, voice assignment, etc. it 
  will choose a voice and call the noteOn()-function of the chosen voice. */
  virtual void noteOn(int newNoteNumber, int newVelocity, int newDetune = 0);

  /** Triggers a note-off. This will also be called from noteOn in case the velocity was zero. */
  virtual void noteOff(int newNoteNumber);

  /** Calls allNotesOff() for all voices. */
  virtual void allNotesOff();

  /** Overrides the setMidiController function inherited from AutomatableModule in order to treat 
  controller 64 (sustain) differently. */
  virtual void setMidiController(int controllerNumber, int controllerValue);

  /** Sets the value of the pitch-wheel in normalized units between -1...+1. */
  virtual void setPitchBend(double newPitchBend);

  /** Switches sustain on or off. */
  virtual void setSustain(bool shouldBeSustained);

  /** Calls reset() for all voices. */
  virtual void resetAllVoices();

  //-----------------------------------------------------------------------------------------------
  // automation:

  /** Overrides the purely virtual parameterChanged() method of the AutomatableModule base
  class with an empty method - you may want to override this in your subclasses. */
  //virtual void parameterChanged(AutomatableParameter* parameterThatHasChanged);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates the output-samples for both channels and stores them at the adresses of *outL and 
  *outR. */
  virtual INLINE void getSampleFrameStereo(double* outL, double* outR);

  /** Calculates the output-samples for both channels and stores them at the adresses of *outL and
  *outR. */
  virtual INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);

  //===============================================================================================

  TuningTable tuningTable; // a table mapping MIDI notes to frequencies for microtuning support

protected:

  /** Applies all the volume scalings such as masterVolume, the cumulative power compensation and 
  mid/side adjustment. */
  INLINE void applyVolumeScalers(double* left, double* right, double cumulativePower);

  /** Calculates and returns a scaling factor for the output signal which compensates (to an user 
  adjustable extent) the gain in ouput power when several voices are playing together - the input
  argument should be the cumulative power of all playing voices. A reaonable way to obtain this 
  value could be to sum the squares of the amplitude envelopes of the individual voices. */
  INLINE double getCumulativePowerCompensator(double cumulativePower);

  /** This function can be used to fill the pendingNoteOff array either with trues or falses. */
  void fillPendingNoteOffArray(bool valueToFillWith);

  /** This function is called when the sustain switch is turned off - it sends out all the note
  offs which have been collected during sustain was on and resets the array of pending note-offs 
  to all false. */
  void dischargePendingNoteOffs();

  /** This is an array of pointers to the actual voices which will be in general of a subclass of 
  PolyphonicInstrumentVoice - but we use only the base class pointer here. The respective subclass
  of PolyphonicInstrument must take care to properly assign these pointers in its constructor such
  that they point to the actual objects which realize those voices. The voice with index 0 plays a 
  special role - it will serve as master or template for all the others and should be accessed by 
  the GUI. It will not be used for actual audio processing. */
  PolyphonicInstrumentVoice** voicePointers;

  /** Number of allocated voices - this must be set in the constructor. The number of playable 
  voices will be one less - the zeroth voice serves only as template for the others. */
  int numAllocatedVoices;

  /** The current number of playable voices - this can be set by setNumPlayableVoices to at most 
  numAllocatedVoices-1 . */
  int numPlayableVoices;

  /** A master level and the corresponding amplitude factor. */
  double masterLevel, masterAmplitude;

  /** An automatic scaling of the master amplitude by the number of playing voices. */
  double masterLevelByVoices;

  /** Mid/side ratio and scaling factors for mid and side signal. You may or may not use them in 
  getSampleFrameStereo in your subclass. */
  double midSideRatio, midScale, sideScale;

  /** The most recent note parameters. */
  int mostRecentNote, mostRecentNoteVel, mostRecentNoteDetune;

  /** Number of currently playing voices. */
  int numActiveVoices;

  /** An array to indicate if there are note off events pending due to sustain. */
  bool pendingNoteOffs[128];

  /** Flag to indicate whether sustain is on. */
  bool sustainIsActive;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void PolyphonicInstrument::applyVolumeScalers(double* left, double* right, 
  double cumulativePower)
{
  // apply power compensation and master amplitude:
  double compensator = getCumulativePowerCompensator(cumulativePower);
  *left  *= masterAmplitude*compensator;
  *right *= masterAmplitude*compensator;

  // apply mid/side scaling:
  double mid  = midScale  * SQRT2_INV * (*left + *right);
  double side = sideScale * SQRT2_INV * (*left - *right);
  *left       =             SQRT2_INV * (mid  +  side);
  *right      =             SQRT2_INV * (mid  -  side);
}

INLINE double PolyphonicInstrument::getCumulativePowerCompensator(double cumulativePower)
{
  double psq = RAPT::rsMax(1.0, cumulativePower);
  double a   = 0.005 * masterLevelByVoices;
  return 1.0 / pow(psq, a); // \todo maybe we can optimize the pow to an exp?
}

INLINE void PolyphonicInstrument::getSampleFrameStereo(double* outL, double* outR)
{
  if(voicePointers == NULL)
    return;

  double accuL, accuR;
  double cumulativePower;
  int    i;

  // initialize the output slots with zero in order to accumulate over the voices and then loop 
  // through the (active) voices:
  accuL             = 0.0;
  accuR             = 0.0;
  numActiveVoices   = 0;
  cumulativePower   = 0.0;

  // loop through the voices, voicePointers[0] servers only as template for the others and is not
  // invoked to process audio:
  for(i=1; i<=numPlayableVoices; i++)
  {
    if(voicePointers[i]->isSilent != true)
    {
      // let voice i generate its output-sample and accumulate it:
      voicePointers[i]->getSampleFrameStereo(&accuL, &accuR, &cumulativePower);

      // increment the variable which keeps track of the active voices:
      numActiveVoices++;
    }
  }

  // apply the amplitude-scaling by the number of voices (weigthed by their respective amplitudes):
  //double compensator = (1.0 + 0.01*masterLevelByVoices) / (1.0 + 0.01*masterLevelByVoices * sqrt(cumulativePower) );
  double compensator = getCumulativePowerCompensator(cumulativePower);
  accuL *= compensator;
  accuR *= compensator;

  // write the final output samples into the memory-slots:
  *outL = accuL;
  *outR = accuR;
}

INLINE void PolyphonicInstrument::getSampleFrameStereo(double* /*inL*/, double* /*inR*/,
  double* outL, double* outR)
{
  getSampleFrameStereo(outL, outR);
}

}

#endif
