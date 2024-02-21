#ifndef rosic_OscillatorStereo_h
#define rosic_OscillatorStereo_h

namespace rosic
{

/** This class defines the user parameters for the OscillatorStereo class.  */

class OscillatorStereoParameters
{

public:

  // user parameters:
  doubleA level;                  // level in decibels
  doubleA levelByKey;             // key dependence of the level in dB@127
  doubleA levelByVel;             // velocity dependence of the level in dB@127
  doubleA pan;                    // pan position between -1...+1
  doubleA panByKey;               // key dependence of the panorama
  doubleA panByVel;               // velocity dependence of the panorama
  doubleA midSide;                // mid/side adjustment bewteen 0...1
  doubleA startPhase;             // start-phase between 0...360
  doubleA startPhaseByKey;        // key dependence of the start-phase
  doubleA startPhaseByVel;        // velocity dependence of the start-phase
  doubleA detuneSemitones;        // detuning in semitones
  doubleA detuneHz;               // additional detuning in Hz
  doubleA stereoDetuneSemitones;  // detuning between left and right in semitones
  doubleA stereoDetuneHz;         // additional detuning between left and right in Hz
  doubleA pitchEnvDepth;          // depth of the pitch-enevlope

  // derived internal parameters:
  doubleA amplitude;              // linear output amplitude
  doubleA panFactorL, panFactorR; // gain factors for left and right
  doubleA midScale, sideScale;    // scale factors for mid and side signal
  doubleA startPosition;          // start position in the table in samples
  doubleA detuneFactor;           // a detune-factor derived for the semitone-detuning
  doubleA stereoDetuneFactorL;    // additional detuning factor for left channel
  doubleA stereoDetuneFactorR;    // additional detuning factor for right channel
  doubleA detuneFactorL;          // resulting overall detuning factor for left channel
  doubleA detuneFactorR;          // resulting overall detuning factor for right channel
  doubleA freqOffsetL;            // resulting overall frequency offset for left channel
  doubleA freqOffsetR;            // resulting overall frequency offset for right channel

  //double transpositionFactor;    // transpose from pitch-wheel as factor
  double sampleRate;
  double sampleRateRec;
  int    tableLength;
  bool   mute;

  OscillatorStereoParameters()
  {
    level                  = 0.0;
    levelByKey             = 0.0;
    levelByVel             = 0.0;
    pan                    = 0.0;
    panByKey               = 0.0;
    panByVel               = 0.0;
    midSide                = 0.0;
    startPhase             = 0.0;
    startPhaseByKey        = 0.0;
    startPhaseByVel        = 0.0;
    detuneSemitones        = 0.0;
    detuneHz               = 0.0;
    stereoDetuneSemitones  = 0.0;
    stereoDetuneHz         = 0.0;
    pitchEnvDepth          = 0.0;

    amplitude              = 1.0;
    panFactorL             = SQRT2_INV;
    panFactorR             = SQRT2_INV;
    midScale               = SQRT2_INV;
    sideScale              = SQRT2_INV;
    startPosition          = 0.0;
    detuneFactor           = 1.0;
    stereoDetuneFactorL    = 1.0;
    stereoDetuneFactorR    = 1.0;
    detuneFactorL          = 1.0;
    detuneFactorR          = 1.0;
    freqOffsetL            = 0.0;
    freqOffsetR            = 0.0;

    //transpositionFactor = 1.0;
    sampleRate          = 44100.0;
    sampleRateRec       = 1.0 / sampleRate;
    tableLength         = 0;
    mute                = false;
  }

};

//=================================================================================================

/** This class implements an oscillator which can produce stereo waveforms. */

class OscillatorStereo
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  OscillatorStereo();

  /** Destructor. */
  ~OscillatorStereo();

  //-----------------------------------------------------------------------------------------------
  // parameter settings:
  // ToDo: sort the setters by whether they apply to the per-sample synthesis or modify the 
  // wavetable

  /** Sets the sample-rate. */
  void setSampleRate(double newSampleRate);

  /** Sets the wavetable object which should be used by this oscillator. */
  void setWaveTableToUse(MipMappedWaveTableStereo *newTableToUse);

  /** Time-reverses the waveform. */
  void setTimeReverse(bool shouldBeReversed);

  /** Inverts the polarity of the waveform. */
  void setPolarityInversion(bool shouldBeInversed);

  /** Sets a time warping coeffcient. @see: MipMappedWaveTableStereo::setTimeWarp() */
  void setFullWavePhaseWarp(double newWarpCoefficient);

  void setHalfWavePhaseWarp(double newWarpCoefficient);

  /** Sets the current key and velocity, so the oscillator can set up it's key and velocity 
  dependent parameters. */
  void setKeyAndVel(int newKey, int newVel);

  /** Sets the nominal frequency to be played. The actual playback-frequency will additionally take 
  into account the detune-factor as set by setDetune(). */
  INLINE void setFrequencyNominal(double newFrequency);

  /** Set this to "true" when this instance should be muted. */
  void setMute(bool shouldBeMuted);

  /** Sets the output level in dB. */
  void setLevel(double newLevel);

  /** Sets the key dependence of the oscillator level in dB/octave. */
  void setLevelByKey(double newLevelByKey);

  /** Sets the velocity dependence of the oscillator level in dB@127. */
  void setLevelByVel(double newLevelByVel);

  /** Sets the panorama (range: -1...+1). */
  void setPan(double newPan);

  /** Sets the key dependence of the panorama postion. */
  void setPanByKey(double newPanByKey);

  /** Sets the velocity dependence of the panorama postion. */
  void setPanByVel(double newPanByVel);

  /** Sets the ratio between mid- and side signal (range: 0...1). */
  void setMidSide(double newMidSide);

  /** Sets the start-phase in degrees (range: 0...360). */
  void setStartPhase(double newStartPhase);

  /** Sets the key dependence of the start-phase. */
  void setStartPhaseByKey(double newStartPhaseByKey);

  /** Sets the velocity dependence of the start-phase. */
  void setStartPhaseByVel(double newStartPhaseByVel);

  /** Sets the multiplicative detuning (in semitones). */
  void setDetuneSemitones(double newDetuneSemitones);

  /** Sets the additive detuning (in Hz). */
  void setDetuneHz(double newDetuneHz);

  /** Sets the multiplicative detuning between left and right channel (in semitones). */
  void setStereoDetuneSemitones(double newStereoDetuneSemitones);

  /** Sets the additive detuning between left and right channel (in Hz). */
  void setStereoDetuneHz(double newStereoDetuneHz);

  /** Sets the modulation depth of the pitch envelope. */
  void setPitchEnvelopeDepth(double newDepth);

  /** Sets the spectral contrast as an exponent acting on the magnitude of each harmonic. */
  void setSpectralContrast(double newContrast);

  /** Sets the spectral  slope in dB/oct for this oscillator. */
  void setSpectralSlope(double newSlope);

  /** Applies a sharp (brickwall) lowpass-filter to the harmonics of the waveform. */
  void setHighestHarmonicToKeep(int newHighestHarmonicToKeep);

  /** Applies a sharp (brickwall) highpass-filter to the harmonics of the waveform. */
  void setLowestHarmonicToKeep(int newLowestHarmonicToKeep);

  /** Adjusts the ratio between even and odd harmonics. */
  void setEvenOddRatio(double newRatio);

  /*
  void setSpectralBumpFrequency(double newFrequency);

  void setSpectralBumpBandwidth(double newBandwidth);

  void setSpectralBumpAmount(double newAmount);

  void setSpectralBumpShape(double newShape);
  */

  /** Sets a phase offset between even and odd harmonics */
  void setEvenOddPhaseShift(double newPhaseShift);

  /** Sets the phase-shift of the the harmonics between left and right channel 8in degrees. The 
  left channel's harmoincs will be shifted half this value into the negative direction and the 
  right channel's harmonics by the same amount into the positive direction. */
  void setStereoPhaseShift(double newPhaseShift);

  /** Sets an additional phase-shift between even and odd harmonics, half of which is applied with
  negative sign to the left channel and with positive sign to the right channel. */
  void setEvenOddStereoPhaseShift(double newPhaseShift);

  /** Set's the transposition-factor for this oscillator which comes from the
  pitch-wheel. */
  //void setTranspositionFactor(double newTranspositionFactor);

  /** Sets a factor for the final amplitude that is supposed to come from some modulator. */
  void setModulatedAmplitude(double newModulatedAmplitude)
  {
    modulatedAmplitude = newModulatedAmplitude;
  }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the mute-status. */
  bool isMuted();

  /** Returns the levele in dB. */
  double getLevel();

  /** Returns the key dependence of the oscillator level in dB/octave. */
  double getLevelByKey();

  /** Returns the velocity dependence of the oscillator level in dB@127. */
  double getLevelByVel();

  /** Returns the panorama position. */
  double getPan();

  /** Returns the key dependence of the panorama position. */
  double getPanByKey();

  /** Returns the velocity dependence of the panorama position. */
  double getPanByVel();

  /** Returns the mid/side ratio. */
  double getMidSide();

  /** Returns the start-phase in degrees. */
  double getStartPhase();

  /** Returns the key dependence of the start-phase. */
  double getStartPhaseByKey();

  /** Returns the velocity dependence of the start-phase. */
  double getStartPhaseByVel();

  /** Returns the multiplicative detuning in semitones. */
  double getDetuneSemitones();

  /** Returns the additive detuning in Hz. */
  double getDetuneHz();

  /** Returns the multiplicative detuning in semitones. */
  double getStereoDetuneSemitones();

  /** Returns the additive detuning in Hz. */
  double getStereoDetuneHz();

  /** Returns the modulation depth of the pitch envelope. */
  double getPitchEnvelopeDepth() { return parameters->pitchEnvDepth; }

  /** Fills the targetBuffer with values suitable for displaying the current waveform in a display. 
  The targetBuffer is assumed to be of size [2][numSamplesToShow] where the first index represents 
  the channel and the second represents the sample-number on the display (which is the x-coordinate 
  in pixels). */
  void getWaveformForDisplay(double** targetBuffer, int numSamplesToShow);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates one stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double* outL, double* outR);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Calculates the phase increment from the nominal frequency and the detune-factor. */
  INLINE void calculateIncrement();

  /** Triggers a calculation of the key and/or velocity dependent parameters. */
  void calculateKeyAndVelocityDependentParameters();

  /** Calculates the phase increment from the nominal frequency and the detune-factor for this 
  instance and also triggers the increment calculation for all slaves. */
  void calculateIncrementForAllSlaves();

  /** Resets the table-pointer to its start-position. */
  void reset();

  /** Initializes all members. */
  //void initialize();

  //-----------------------------------------------------------------------------------------------
  // matser/slave configuration:

  /** Adds a  slave instance to the vector of slaves - this will also cause the 'isMaster' flag of 
  the slave to be set to false, redirect the slaves parameters-pointer to the one of this instance 
  and delete the old (now unused) parameters pointer of the slave. */
  void addSlave(OscillatorStereo* newSlave);

  //-----------------------------------------------------------------------------------------------
  // embedded publically accessible objects:

  /** This is a pointer to the wavetable object to be used. We use a pointer instead of a direct 
  member object in order to share the same wavetable among a bunch of oscillators - the idea is, 
  that an oscillator which is part of a synthesizer-voice can share the wavetable with its siblings
  of the other voices. */
  MipMappedWaveTableStereo *waveTable; // may this also be a subclass object like MultiWaveTable?

  //===============================================================================================

protected:

  doubleA positionFracL, positionFracR;
  doubleA incrementFracL, incrementFracR;
  int     positionIntL, positionIntR;
  int     incrementIntL, incrementIntR;

  /** The version of the table which is to be used (will depend on the increment) .*/
  int     mipMapIndex;

  /** A pointer to the parameters which are potentially shared by among instances. */
  OscillatorStereoParameters* parameters;

  double frequencyNominal;                      // nominal frequency (without detune factor)
  double frequencyDetunedL, frequencyDetunedR;  // final frequencies including detuning

  // the key and/or velocity dependent parameters scaled by the current key and velocity:
  double scaledAmplitude, panFactorL, panFactorR, scaledStartPhase;
  double modulatedAmplitude; // takes also modulation into account

  // the current key and velocity:
  int key, vel;

  /** A vector of pointers to other instances of this class which shall be kept in sync to this
  instance with regard to their parameters. */
  std::vector<OscillatorStereo*> slaves;

  /** A flag which indicates whether or not this instance is a master which controls other
  instances of this class - this will also determine whether or not this objects will delete the
  pointer to the parameter set on destruction. By default, instances will be constructed as
  master, later they can be re-configured as slaves by adding them as slaves via addSlave to
  another instance. */
  bool isMaster;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void OscillatorStereo::setFrequencyNominal(double newFrequency)
{
  if(newFrequency != frequencyNominal && newFrequency >= 0.0 && newFrequency <= 22000.0)
  {
    frequencyNominal = newFrequency;
    calculateIncrement();
  }
}

INLINE void OscillatorStereo::calculateIncrement()
{
  // take transposition and detuning into account for the calcluation of the
  // final frequency:
  frequencyDetunedL = parameters->detuneFactorL * frequencyNominal + parameters->freqOffsetL;
  frequencyDetunedR = parameters->detuneFactorR * frequencyNominal + parameters->freqOffsetR;

  // calculate the new phase increment:
  double phaseIncrementL =
    frequencyDetunedL * (double)parameters->tableLength * parameters->sampleRateRec;
  double phaseIncrementR =
    frequencyDetunedR * (double)parameters->tableLength * parameters->sampleRateRec;

  // setup the member 'mipMapIndex' acccording to the exponent of the increment:
  mipMapIndex  = ((int)EXPOFDBL(phaseIncrementL));
  mipMapIndex += 1; // generates frequencies up to nyquist/2 on the highest note
  if(mipMapIndex < 0)
    mipMapIndex = 0;
  if(mipMapIndex > 9)
    mipMapIndex = 9;

  // split it into integer and fractional part:
  incrementIntL  = (int)phaseIncrementL;
  incrementIntR  = (int)phaseIncrementR;
  incrementFracL = phaseIncrementL - (double)incrementIntL;
  incrementFracR = phaseIncrementR - (double)incrementIntR;
}

INLINE void OscillatorStereo::getSampleFrameStereo(double* outL, double* outR)
{

  // catch some special conditions:
  if((parameters->mute == true) || (waveTable == NULL))
  {
    *outL = 0.0;
    *outR = 0.0;
    return;
  }

  // wraparound the integer part of the position-pointer if necesarry:
  while(positionIntL >= parameters->tableLength)
    positionIntL -= parameters->tableLength;
  while(positionIntR >= parameters->tableLength)
    positionIntR -= parameters->tableLength;

  double tmpL, tmpR, tmpM, tmpS;
  //tmpM = parameters->midScale *waveTable->getValueFast(0, mipMapIndex, positionInt,positionFrac);
  //tmpS = parameters->sideScale*waveTable->getValueFast(1, mipMapIndex, positionInt,positionFrac);
  tmpL  = waveTable->getValueFast(0, mipMapIndex, positionIntL, positionFracL);
  tmpR  = waveTable->getValueFast(1, mipMapIndex, positionIntR, positionFracR);
  tmpM  = parameters->midScale  * (tmpL+tmpR);
  tmpS  = parameters->sideScale * (tmpL-tmpR);
  //*outL = parameters->amplitude * parameters->panFactorL * (tmpM+tmpS);
  //*outR = parameters->amplitude * parameters->panFactorR * (tmpM-tmpS);
  *outL = modulatedAmplitude * scaledAmplitude * parameters->panFactorL * (tmpM+tmpS);
  *outR = modulatedAmplitude * scaledAmplitude * parameters->panFactorR * (tmpM-tmpS);

  // increment position-pointer:
  positionIntL  += incrementIntL;
  positionFracL += incrementFracL;
  if(positionFracL >= 1.0)
  {
    positionFracL -= 1.0;
    positionIntL  += 1;
  }
  positionIntR  += incrementIntR;
  positionFracR += incrementFracR;
  if(positionFracR >= 1.0)
  {
    positionFracR -= 1.0;
    positionIntR  += 1;
  }

}

}

#endif
