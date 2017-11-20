#ifndef VST_CHAOSGENERATOR_H
#define VST_CHAOSGENERATOR_H

#include "../Common/Utilities/vstPlugIn.h"

/*
ToDo list:

-check the clipping function - it seems it goes above the identity function
-add smoothing parameter for phase jumping
-envelope shapes
-hard sync whenever midi is received (on/off)

// new sync behavior:
-generate an input saw
-forget about edge detection
1. at saw edge, if above threshold, oscillate until 180 degrees
when it enters threshold again, continue oscillating starting at 180 degrees
when it leaves threshold again, oscillate from 180 to 360 (back to 0) degrees
i mean you just need a few if statements checking phase at any given moment
if (above_threshold && phase != 180) continue_oscillating();

*/


/** Baseclass for monophonic synthesizers. It handles stuff like keeping track of the incoming MIDI
events and how they affect the frequency that should be played. This includes the handling of pitch
modifiers such as a global detune factor and portamento/glide, etc. 

The subclass is supposed to call getCurrentFrequency(); on a per sample basis to inquire the actual
frequency to be played

todo: 
 -add pitchwheel support
 -add microtuning
 -implement a "slide-back" behavior

*/

class MonoSynth
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  MonoSynth();


  /** \name Setup */

  /** Sets the sample rate in Hz. */
  void setSampleRate(double newSampleRate);

  /** Sets the time to glide from one note to another in seconds */
  void setGlideTime(double newTime);

  /** Sets a global detuning in semitones. */
  void setDetune(double newDetune);


  /** \name Event Handling */

  /** Triggers a note with given key and velocity. */
  void noteOn(int key, int velocity);


protected:

  /** Supposed to be overriden by subclasses to inform this baseclass object if the synth is 
  currently in a silent state (typically, when the amplitude envelope has reached its endpoint).
  This information will be used in noteOn in order to decide whether or not triggerAttack will
  be called. */
  virtual bool isSilent() = 0;

  /** Override this in your subclass to retrigger your envelopes. */
  virtual void triggerAttack(int key, int velocity) = 0;

  /** Override this in your subclass trigger the release phase of your envelopes. */
  virtual void triggerRelease() = 0;


  // data:
  double sampleRate; 
  double glideTime;             // glide time in seconds
  double targetFrequency;
  double currentFrequency;
  double freqFactorPerSample;
  double detuneFactor;
  int currentNote;
  int glideSamples;             // glide time in samples
  int remainingGlideSamples; 


  /** Starts to glide to the desired pitch */
  void glideTo(double pitch);

  /** Immediately jumps to the new pitch */
  void jumpTo(double pitch);

  /** Returns the current frequency according to glide. This is supposed to be called for each
  sample. It keeps track of which frequency should be played according to the glide setting. */
  double getCurrentFrequency();
    // supposed to be inquired by subclass on a per-sample basis (maybe inline)
    // maybe rename to getAndUpdateCurrentFrequency

};

//=================================================================================================

/** This is a sound generator based on a feedback frequency modulation setup using two oscillators
and a highpass and lowpass filter in a structure like:

-----> Osc1 ---> HPF ---> LPF ---> Osc2 ------->
  |                                      |
  \--------------------------------------/

\todo:
-separate this class into a separate pair of .h/.cpp files at some point
-have an adjustable frequency offset (or scaler) for the modulator's frequency
-maybe have individually adjustable modulation depths for osc1-->osc2 and osc2-->osc1
 (or one master-depth and a scaler for the feedback modulation)
-maybe set up the modulation depth in terms of a modulation index
-maybe have also self-modulation, i.e. not just 1->2 and 2->1 but also 1->1 and 2->2

*/

class ChaosGenerator : public MonoSynth
{

public:

  /** This is an enumeration of the available waveforms for the oscillators. */
  enum waveforms
  {
    SINE = 0,
    TRIANGLE,
    SQUARE,
    SAW_UP,
    SAW_DOWN,

    NUM_WAVEFORMS
  };

  /** Enumeration of different synchronization modes. */
  enum syncModes
  {
    SYNC_RESET,          // reset to start phase
    SYNC_JUMP,           // jump phase 
    SYNC_REVERSE,        // reverse direction

    NUM_SYNCMODES
  };


  /** \name Construction/Destruction */

  /** Constructor. */
  ChaosGenerator();


  /** \name Setup */

  /** Sets the sample rate in Hz. */
  void setSampleRate(double newSampleRate);

  /** Sets the basic frequency of the carrier oscillator (oscillator 2). */
  void setFrequency(double newFrequency);
    // make this function protected - it's a low-level function that's not supposed to be called by 
    // client code anymore (instead, use the inherited setDetune method and pass a noteOn event
    // to determine frequency)
   
  /** Sets the factor by which the filtered input signal is mixed into the output. */
  void setInputMix(double newMixFactor);

  /** Sets the scale factor for the cutoff frequency of the input filter with respect to the master 
  frequency. */
  void setInputFreqScaler(double newScaler);

  /** Sets the amount of detuning for the input filter with respect to the master frequency (in
  semitones). */
  void setInputFilterDetune(double newDetune);

  /** Sets the ratio of the frequencies of the two oscillators. It acts as a multiplier for the
  frequency of oscillator 1 (seen as modulator in the most simple FM/PM patch). Oscillator 2 (the
  carrier) will use the master frequency as is. I think, using it this way means the value is 
  actually the reciprocal of what is called C/M ratio in FM synthesis. (check this) */
  void setFrequencyRatio(double newRatio);

  /** Sets a scale factor for the depth of both modulations 1->2 and 2->1. */
  void setModulationDepthScaler(double newScaler);

  /** Sets the depth by which oscillator 1 modulates oscillator 2. */
  void setModDepth12(double newDepth);

  /** Sets the depth by which oscillator 2 modulates oscillator 1. */
  void setModDepth21(double newDepth);

  /** Adjusts the split between frequency-modulation and phase-modulation for 1->2. It's a value
  between 0 and 1 where 0 means pure frequency-modulation and 1 mean pure phase modulation. */
  void setFreqVsPhaseMod12(double newValue);
  void setFreqVsPhaseMod21(double newValue); // same for 2->1

  /** Depth of the modulation envelope as normalized value between 0..1. */
  void setModEnvDepth(double newDepth);

  /** Attack time for modulation envelope in milliseconds. */
  void setModEnvAttack(double newAttack);

  /** Decay time for modulation envelope in milliseconds. */
  void setModEnvDecay(double newDecay);

  /** Sustain level for modulation envelope as raw multiplier. */
  void setModEnvSustain(double newSustain);

  /** Release time for modulation envelope in milliseconds. */
  void setModEnvRelease(double newRelease);

  /** Sets the 1st filter's cutoff frequency as scale factor with respect to the  master 
  frequency. */
  void setFilter1FreqScaler(double newScaler);
  void setFilter2FreqScaler(double newScaler);  // same for 2nd filter

  /** Sets up the 1st filter's frequency scale factor in terms of a detune value in
  semitones.   */
  void setFilter1Detune(double newDetune);
  void setFilter2Detune(double newDetune);  // same for 2nd filter

  /** Sets the resonance for the 1st filter (normalized between 0 and 1) */
  void setFilter1Resonance(double newReso);
  void setFilter2Resonance(double newReso);  // same for 2nd filter

  /** Selects the mode for the 1st filter. this should be one of rsLadderFilter::modes */
  void setFilter1Mode(int newMode);
  void setFilter2Mode(int newMode);  // same for 2nd filter

  // \todo maybe make both filters multimode ladders, call them Filter1, Filter2 instead of 
  // Lowpass/Highpass and make their position in the signal flow adjustable

  /** Sets the order of the lowpass filter */
  //void setLowpassOrder(int newOrder);

  /** Sets the center frequency (as multiplier for the master frequency) for the feedback EQ. */
  void setFeedbackFreqScaler(double newScaler);

  /** Sets up the feedback filter's frequency scale factor in terms of a detune value in
  semitones.   */
  void setFeedbackDetune(double newDetune);

  /** Sets the amount of boost (or cut) in decibels for the feedback EQ. */
  void setFeedbackBoostAmount(double newAmount);

  /** Sets the bandwidth in octaves for the feedback EQ. */
  void setFeedbackBoostWidth(double newWidth);

  /** Sets a linear gain multiplier for final output. */
  void setGain(double newGain);

  /** Sets the detuning with respect to the master frequency for oscillator 1 (in semitones) */
  //void setDetune1(double newDetune);

  /** Sets the detuning with respect to the master frequency for oscillator 2 (in semitones) */
  //void setDetune2(double newDetune);
    // maybe it's redundant to set 2 detunes, maybe we should just set a detune for the
    // modulator

  /** Sets the absolute frequency offset for oscillator 1 with respect to the master frequency in 
  Hz */
  void setFreqOffset1(double newOffset);
  void setFreqOffset2(double newOffset); // same for osc1
    // experimental


  /** Sets the waveform for oscillator 1 (see enum waveforms). */
  void setWaveform1(int newWaveform);

  /** Sets the waveform for oscillator 2 (see enum waveforms). */
  void setWaveform2(int newWaveform);

  /** Sets the clipping level for oscillator 1. */
  void setClipLevel1(double newLevel);

  /** Sets the clipping level for oscillator 2. */
  void setClipLevel2(double newLevel);

  /** Sets the phase offset for oscillator 1 in radians. */
  void setPhase1(double newPhase);

  /** Sets the phase offset for oscillator 2 in radians. */
  void setPhase2(double newPhase);

  /** Attack time for amplitude envelope in milliseconds. */
  void setAmpEnvAttack(double newAttack);

  /** Decay time for amplitude envelope in milliseconds. */
  void setAmpEnvDecay(double newDecay);

  /** Sustain level for amplitude envelope as raw multiplier. */
  void setAmpEnvSustain(double newSustain);

  /** Release time for amplitude envelope in milliseconds. */
  void setAmpEnvRelease(double newRelease);


  /** Selects the synchronization mode - this determines, what happens when a sync-trigger is
  encountered. Possible values are: 0: phase reset, 1: phase jump, 2: direction reversal. */
  void setSyncMode(int newMode);

  /** Sets the threshold for the input edge detector. When an edge in the input is detected that is
  above this threshold, the sound generator will be reset. */
  void setSyncThreshold(double newThreshold);

  /** Sets the amount of the synchronization effect as value between 0 and 1. 0: no sync, 1: full 
  sync, between: partial sync (reset phase in between current value and start value) */
  void setSyncAmount(double newAmount);

  /** Decide, where the filters should be reset on sync events (as opposed to resetting only the
  oscillators). */
  void setFilterResetOnSync(bool shouldReset);

  // "Ducking" parameters: the instantaneous amplitude of our sound generator output will be 
  // attenuated according to the instaneous amplitude of the input signal. The louder the input,
  // the more attenuation is applied to the output. This can be used to emulate the voltage 
  // limiting effects in analog filters. It applies an input-dependent microenvelope to our output
  // signal.

  /** Sets the range of amplitudes for the input signal for which the output signal of the
  sound generator will be not ducked away completetly. This will affect the total width of the
  microenvelope. */
  void setDuckingRange(double newRange);

  /** Sets the relative (with respect to the total range) width of the flat top zone between 0 
  and 1. At 0, we see a nice bell shape for the microenvelope and towards 1 it becomes more and
  more squarish. */
  void setDuckingFlatness(double newFlatness);

  /** You may shift the center value for the ducking from 0 to somewhere else. For a sawtooth 
  input, this shifts the microenvelope back and forth within the cycle. */
  void setDuckingCenter(double newCenter);

  /** Sets the transition shape for the ducking: 0: linear, 1: cubic, 2: quintic, 3: heptic. */
  void setDuckingShape(int newShape);

  /** Sets the level at which the chaos generator's output signal is clipped. */
  void setOutputClipLevel(double newLevel);


  /** \name Audio Processing */

  /** Computes a stereo output sample frame at a time */
  void processSampleFrame(double *inL, double *inR, double *outL, double *outR);
    // todo: use the input signal for something


  /** \name Misc */

  /** Resets the internal state-variables. */
  void reset();

  /** Resets the oscillators only. */
  void resetOscillators();

  /** Resets the filters only. */
  void resetFilters();


protected:

  /** \name Internal Functions */

  /** Given a current input sample "in" and a previous input sample "inOld", this function figures 
  out,  */
  //bool didMicroEnvelopeStart(double in, double inOld);

  /** Given a current input sample "in" and a previous input sample "inOld", this function figures 
  out, if we need to trigger a sync reset event. */
  bool needsSyncTrigger(double in, double inOld);
    // todo: include a return parameter that is used for subsample precision sync (figures out,
    // when exactly between the current and previous sample, the sync event occurred

  /** Triggers oscillator reset (and possibly filter reset, too). Called from syncIfNeccessary. */
  void triggerSync();

  /** Resets the oscillators. To be called in case of occurrence of a sync event. According to our
  syncAmount parameter, it may do only a partial reset. */
  void syncOscillators();
    // maybe include an "advance" parameter later that advances the state for a subsample time
    // difference (which should be obtained by determining the subsample time instant opf the sync
    // event)

  /** Resets the filter in case of a sync event. @see syncOscillators(). */
  void syncFilters();

  /** Applies soft clipping to the input value x.  The level parameter determines the ceiling for
  the absoulte value of the signal. Distortion of the waveshape begins already at half of this 
  value, i.e. x-values for which |x| <= 0.5*level, will be passed unchanged, 
  0.5*level < |x| < level will be attenuated and |x| > level will be mapped to the saturation 
  level. */
  double clip(double x, double level);

  /** Returns the value of the waveform indexed by "shape" (see enum waveforms) given a value
  of the instantaneous phase. */
  double getWaveform(double phase, int shape);


  // mandatory overrides for the MonoSynth baseclass:
  virtual bool isSilent();
  virtual void triggerAttack(int key, int velocity);
  virtual void triggerRelease();

  // user parameters:
  double inMix;          // factor for mixing filtered input to output
  double frequency;      // the master frequency (applies to both oscillators)
  double freqRatio;
  double modScale;       // master modulation depth scaler (scales mod12, mod21)
  double mod12, mod21;   // modulation depths  (maybe add mod11, mod22 later)
  double fmpm12, fmpm21; // frequency modulation (FM) vs phase modulation (PM)
  double clip1, clip2;   // clipping levels for osc 1 and 2
  double inScl;          // cutoff scaler for input filter
  double hpScl, lpScl;   // cutoff scalers for filters
  double fbScl;          // frequency scaler for feedback boost
  double phase1, phase2; // phase offsets
  double clipOut;        // clipping level for output
  double gain;           // linear gain for output
  int    shape1, shape2; // shape of the waveforms
  double modEnvDepth;    // depth of the modulation envelope
  double freqOffset1;
  double freqOffset2;

  int    syncMode;       // 0: reset, 1: jump, 2: reverse
  double syncThresh;     // synchronization threshold for input edge
  double syncAmount;     // amount of sync, 0: hardsync, 1: no sync, between: partial sync
  bool   syncFlt;        // switch to turn on/off filter reset on sync

  double duckCenter;     // center value for ducker
  double duckWidth;      // width for ducker
  double duckLo;         // lower threshold for ducker/micro-envelope sync
  double duckHi;         // upper threshold for ducker/micro-envelope sync


  // state variables:
  double inOld;          // remembered input signal
  double pos1, pos2;     // instantaneous phases of both oscillators
  double out1, out2;     // output signal of both oscillators
  double direction;      // +1 or -1, multiplies phase increment

  // embedded objects:
  //rsOnePoleFilter           hpf;            // highpass filter -> maybe replace with ladder, too
  //rsEngineersFilter         lpf;            // lowpass filter -> replace with ladder

  rsLadderFilter            hpf;            // rename to flt1
  rsLadderFilter            lpf;            // rename to flt2
  rsStateVariableFilter     fbFilter;       // filter in feedback path (replace by biquad)
  rsBreakpointModulator     ampEnv, modEnv; // ADSR envelopes for amplitude and modulation
  rsLadderFilter            inFilter;       // filter for input signal
  rsParametricBellFunction  ducker;         // voltage limit simulation

  // todo:
  // -maybe refactor the envelope setup functions in a way, such that the user calls functions 
  // like setAttack(int whichEnvelope, double newAttack). This avoids some code duplication
  // especially when even more enevlopes are introduced)
  // -use a square rule for the FM/PM parameters - it seems the output bandwidth is proportional
  // to the modulation index squared
};


//=================================================================================================

//class vstChaosGenerator : public vstPlugIn
class vstChaosGenerator : public vstInstrument
{

public:	
  vstChaosGenerator(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void  updateCoreParameter(VstInt32 index, float value);
  virtual void  getParameterLabel  (VstInt32 index, char* label);
  virtual void  getParameterDisplay(VstInt32 index, char* text);
  virtual void  getParameterName   (VstInt32 index, char* text);

  virtual void filterModeName(int index, char* text);

  virtual void onNoteOn(int note, int velocity, int detune);
  virtual void onControlChange(int index, int value);
  virtual void onPitchWheel(int value);

protected:
  // parameter indices:
  enum parameters
  {
    TUNE,              // master tuning
    GLIDE_TIME,        // time for portamento
    IN_MIX,            // mix factor filtered input
    IN_DETUNE,         // input filter detuning
    FREQ_RATIO,        // ratio of carrier (osc2) and modulator (osc1) frequency

    // modulation:
    MOD_SCALE,         // modulation amount/scaler
    MOD_12,            // modulation depth 1->2
    FMPM_12,           // adjustment between frequency- and phase-modulation
    MOD_21,            // ...same for 2->1
    FMPM_21,           // ...
    MOD_ENV,           // modulation envelope depth
    MOD_ATT,           // attack for modulation envelope
    MOD_DEC,           // decay
    MOD_SUS,           // sustain
    MOD_REL,           // release

    // filters:
    //LOWPASS_DETUNE,    // lowpass filter detune
    //LP_ORDER,
    //HIGHPASS_SCALE,    // highpass cutoff scale factor
    FLT1_SCALE,
    FLT1_RESO,
    FLT1_MODE,
    FLT2_SCALE,
    FLT2_RESO,
    FLT2_MODE,

    FB_DETUNE,         // feedback filter detune
    FB_BOOST,
    FB_WIDTH,

    // oscillators:
    FREQ_OFFSET1,     // maybe replace by detune1/2 or get rid of that
    FREQ_OFFSET2,
    WAVE1,
    WAVE2,
    CLIP1,
    CLIP2,
    PHASE1,           // the phase parameters don't seem to be useful
    PHASE2, 

    // amplitude:
    VOLUME,            // rename to AMPLITUDE
    AMP_ATT,           // attack for amp envelope
    AMP_DEC,           // decay
    AMP_SUS,           // sustain
    AMP_REL,           // release
    CLIP_OUT,

    // sync to input:
    SYNC_MODE,
    SYNC_THRESH,       // threshold for edge detector that triggers sync
    SYNC_AMOUNT,       // amount of sync: 1: full, 0: none
    SYNC_FILTERS,      // choose to reset filters or not on sync events

    // ducking:
    DUCK_RANGE,        // input range, for which there is nonzero output gain
    DUCK_FLATNESS,     // 0: round, 1: flat/square
    DUCK_CENTER,       // shifts the center
    DUCK_SHAPE,        // 0: linear, 1: cubic, 2: quintic, 3: heptic - smoothness order

    // clipping:
    //CLIP_LEVEL,
    //CLIP_POSITION, 

    NUM_PARAMETERS
  };

  // mapped parameters:
  double tune, glide_time, in_mix, in_det, freq_ratio,
    mod_scale, mod_12, fmpm12, mod_21, fmpm21, mod_env, mod_att, mod_dec, mod_sus, mod_rel,
    fltScl1, fltRes1, fltScl2, fltRes2, fb_det, fb_boost, fb_width,
    freq_offset1, freq_offset2, phase1, phase2, clip1, clip2,
    volume, amp_att, amp_dec, amp_sus, amp_rel, clipOut,
    sync_thresh, syncAmount, duckRange, duckFlat, duckCenter;
  int fltMode1, fltMode2, wave1, wave2, syncMode, duckShape;
  bool sync_filters;

  // the core DSP object:
  ChaosGenerator chaosGen;


  // Ideas:
  // -let the user choose the position of both filters in the signal flow, make both multimode
  //  ladders
  // -have a macro-parameter for brightness that controls several low-level parameters at once
  // -


};

#endif
