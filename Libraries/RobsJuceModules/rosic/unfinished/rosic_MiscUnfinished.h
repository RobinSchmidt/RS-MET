#pragma once

namespace rosic
{


//=================================================================================================

template<class TSig, class TPar> // move to RAPT
class rsVectorMixer
{

public:

  enum class Mode
  {
    linear,
    sinCosApprox,
    sinCos
  };

  void setX(TPar newX) { x = newX; }

  void setY(TPar newY) { y = newY; }

  void getGains(TSig* topLeft, TSig* topRight, TSig* bottomLeft, TSig* bottomRight)
  {
    TPar xx = TPar(0.5) * (x + TPar(1));  // map -1..+1 to 0..1
    TPar yy = TPar(0.5) * (y + TPar(1));

    // todo: switch between modes - this is for linear:
    TPar right  = xx;
    TPar left   = TPar(1) - xx;
    TPar top    = yy;
    TPar bottom = TPar(1) - yy;

    *topLeft     = TSig(top    * left);
    *topRight    = TSig(top    * right);
    *bottomLeft  = TSig(bottom * left);
    *bottomRight = TSig(bottom * right);
  }

protected:

  Mode mode = Mode::linear;

  TPar x = TPar(0), y = TPar(0);

};
// todo: maybe have separate coordinates for left and right channel xL, yL, xR, yR
// may de-templatize it -> have xL, yL, xR, yR as double and the 4 gains are computed as 
// rsFloat64x2


class rsVectorMixerPoly : public rsPolyModule, public rsVectorMixer<rsFloat64x2, double>
{

public:




  /*
  void processFrame(const double* in, int numIns, double* out, int numOuts, int voice) override
  {
    RAPT::rsAssert(numOuts == 4);
    out[0] = out[1] = out[2] = out[3] = 1.0; // preliminary
  }
  */

protected:

};

//=================================================================================================

class rsModulatorArrayPoly : public rsPolyModule
{

public:

protected:

};

//=================================================================================================

class rsTriSawOscPoly : public rsPolyModule, public RAPT::rsTriSawOscillator<rsFloat64x2>
{

public:


  rsFloat64x2 getSample(const rsFloat64x2& in, int voice) override
  {
    return rsFloat64x2(0, 0);
  }


protected:

  std::vector<RAPT::rsTriSawVoice<rsFloat64x2>> voices;
  
  //RAPT::rsTriSawVoice<rsFloat64x2> voices; // the code for this is in the wrong file - move over

};

//=================================================================================================

class rsLadderFilterPoly : public rsPolyModule
{

public:

protected:

};

//=================================================================================================

class rsAttackDecayEnvelopePoly : public rsPolyModule
{

public:

protected:

};

//=================================================================================================

/** A spectral processing based pitch- and frequency shifter. Uses the so called phase-vocoder 
approach to achieve the pitch-shifting effect. This is a (somewhat ambiguous) umbrella term for 
certain manipulations of the short-time Fourier transform (STFT) of the signal. These techniques
approach the problem by shifting the magnitudes of the spectral peaks that correspond to sinusoidal 
parts of the signal and typically involve some phase prediction step for the detected sinusoids 
from the previous FFT frame to obtain a refined frequency estimate for the sinusoid.

...TBC...


References:

(LD) NEW PHASE-VOCODER TECHNIQUES FOR PITCH-SHIFTING, HARMONIZING AND OTHER EXOTIC EFFECTS
     Jean Laroche, Mark Dolson (DAFX 1999)
     https://www.ee.columbia.edu/~dpwe/papers/LaroD99-pvoc.pdf

(JH) LOW LATENCY AUDIO PITCH SHIFTING IN THE FREQUENCY DOMAIN
     Nicolas Juillerat, Beat Hirsbrunner
     https://www.researchgate.net/publication/261078164_Low_latency_audio_pitch_shifting_in_the_frequency_domain
     https://pitchtech.ch/


*/

class SpectralShifter : public SpectralProcessor
{

public:


  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  SpectralShifter(int maxBlockSize, int maxOverlapFactor = 4, int maxPaddingFactor = 4);

  virtual ~SpectralShifter();


  //-----------------------------------------------------------------------------------------------
  // \Setup

  void setFrequencyScale(double scaleFactor) { freqScale = scaleFactor; }

  enum class Algorithm
  {
    RobSchm1,   // Robin Schmidt's 1st algo
    RobSchm2,   // ...2nd
    LaroDols,   // Laroche, Dolson
    JuilHirs    // Juillerat, Hirsbrunner

  };

  void setAlgorithm(Algorithm newAlgorithm) { algo = newAlgorithm; }
  // should switch between Laroche/Dolson, Juillerat/Hirsbrunner, ...etc. algorithms. Maybe
  // call them LD, JH, etc.

  enum class PhaseFormula
  {
    keepOriginal,
    useMultiplier  // rename to useTwiddleFactor
  };

  void setPhaseFormula(PhaseFormula newFormula) { phaseFormula = newFormula; }


  //-----------------------------------------------------------------------------------------------
  // \Processing

  void reset()
  {
    SpectralProcessor::reset();
    frameIndex = 0;
  }


protected:



  
  void processSpectrum(Complex* spectrum, int spectrumSize) override;


  // rename these: maybe doAlgoLaroDols or algoLaroDols

  void shiftViaLD(Complex* spectrum, int spectrumSize);
  // stub - Implementation of the algorithm of Laroche/Dolson

  void shiftViaJH(Complex* spectrum, int spectrumSize);
  // stub - Implementation of the algorithm of Juillerat/Hirsbrunner

  void shiftViaRS1(Complex* spectrum, int spectrumSize);
  // stub - Custom algorithm by Robin Schmidt


  void shiftViaRS2(Complex* spectrum, int spectrumSize);


  Complex *tmpSpectrum = nullptr;

  double freqScale = 1.0;  // rename to scale or freqScale

  Algorithm algo = Algorithm::JuilHirs;

  int frameIndex = 0;
  // Try to get rid of this. The JuilHirs algo needs this for its formula but maybe the formula can
  // be expressed without it. Try to figure this out! If it turns out to be so important that it
  // can't be removed, we may need to take care of letting it wrap around back to zero at 
  // appropriate instants. The paper says in 3.3 that vertical phase coherence is achieved every
  // O frames where O is the overlap factor. So maybe that means we can wrap around frameIndex at 
  // every multiple of O?

  PhaseFormula phaseFormula = PhaseFormula::useMultiplier;


  std::vector<double> mag, phs, phsOld; // Buffers for magnitude and phase
  // ToDo: 
  // -Maybe use raw array like with tmpSpectrum. Actually, We could repurpose the tmpSpectrum
  //  array to save memory. Maybe we should have a member Void* workspace; that gets casted into
  //  the appropriate buffer pointer type in shiftVia... Different algorithms need different types
  //  of buffers - but they could use the same memory space nonetheless. But that's an 
  //  optimization thing for later. But these 3 buffers need space for 3*maxSpectrumSize whereas 
  //  the single complex buffer for tmpSpectrum needs only 2*maxSpectrumSize

};

// ToDo:
// -Maybe the keyword PhaseVocoder should appear somewhere in this file to make this class easier 
//  to find. Hey - now it does! :-)

//=================================================================================================

/** A chain of allpass filters that has "zap" like sound as its impulse response, i.e. a fast 
sinusoidal downward sweep. Due to its allpass nature, the overall output has a white (i.e. flat) 
spectrum. It's meant to be used to create (raw material for) synthesized drum and percussion 
sounds. It turned out to be useful to apply a first order lowpass afterwards to convert the white 
spectrum into a brown one, i.e. one with -6 dB/oct falloff. The browning filter should probably be 
tuned somewhere below the lowest allpass tuning frequency as set by setLowFreq().

Eventually, it can be driven by sources other than an impulse generator. Noise bursts are
interesting as excitation signal, too. */

class rsFlatZapper
{

public:


  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  rsFlatZapper();


  //-----------------------------------------------------------------------------------------------
  // \Setup


  /** Sets the sample rate at which this object should operate. */
  void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; setDirty(); }

  /** Enumeration of the available filter modes. */
  enum class Mode
  {
    bypass = 0,
    onePole,       // maybe rename to allOnePole
    //twoOnePoles,   // 2 one poles lumped into a biquad per stage -> optimization
    biquad,        // maybe rename to allBiquad

    // ToDo:
    // lowFirst,   // use 1st order filters for lower k filters (and biquads for the rest)
    // lowBiquad,  // use 2nd order filters for lower k filters (and one-poles for the rest)

    numModes
  };
  // Maybe have modes that alternate between first order and second order filters - but no, a 
  // similar effect can be achieved by using two rsWhiteZappers in series - one with biquads and
  // one with first order filters, so it's not worthwhile to do internally.
  // But maybe it coul be worthwhile to have the first k filters firstorder and the remaining N-k
  // filters second order (or the other way around) for some user given k. That could perhaps 
  // introduce a knee in the phase response in the log freq domain
  // 

  /** Sets the mode of the filter chain, i.e. the type of allpass filter that will be used per 
  stage. @see Mode. */
  void setMode(Mode newMode) { mode = newMode; setDirty(); }

  /** Sets the number of allpass filter stages to be used. Using more stages generally leads to a 
  more elongated sweep/zap. The sweet spot seem to be around 50. */
  void setNumStages(int newNumStages) { allpassChain.setNumStages(newNumStages); setDirty(); }

  /** Sets the frequency to which the lowest allpass filter should be tuned. */
  void setLowFreq(double newFreq) { freqLo = newFreq; setDirty(); }

  /** Sets the frequency to which the highest allpass filter should be tuned. */
  void setHighFreq(double newFreq) { freqHi = newFreq; setDirty(); }

  /** Sets the shape parameter in the frequency domain that the tuning frequencies of the allpasses
  will follow. A value of 0 will space out the tuning pitches out linearly, i.e. it will give 
  linear spacing on a log-frequency axis. Values above 0 will bend the curve upwards ...TBC... */
  void setFreqShape(double newShape) { freqShape = newShape; setDirty(); }

  /** Sets the Q (quality factor) to which the lowest allpass will be set. The sweet spot for Q 
  seems to be around 1.0. Higher values will ...TBC... */
  void setLowQ(double newQ) { qLo = newQ; setDirty(); }

  /** Sets the Q (quality factor) to which the lowest allpass will be set. The sweet spot for Q 
  seems to be around 1.0. ...I think, the highQ setting admits higher values than lowQ without 
  sounding like crap (verify!) */
  void setHighQ(double newQ) { qHi = newQ; setDirty(); }

  /** Sets the shape that the Q-values for the filters will follow. It's similar to 
  setFreqShape. */
  void setQShape(double newShape) { qShape = newShape; setDirty(); }

  /** Initializes all parameter values to their initial/default values. The sample rate may or may
  not be re-initialized also. */
  void initSettings(bool initAlsoSampleRate = false);

  //-----------------------------------------------------------------------------------------------
  // \Processing

  /** Produces one output sample from a given input sample at a time. */
  inline double getSample(double in)
  {
    if(dirty)
      updateCoeffs();
    return allpassChain.getSample(in);
  }

  /** Resets the processing state (i.e. filter buffers). Does not affect parameter setup. */
  void reset()
  {
    allpassChain.reset();
  }



protected:

  /** Sets our dirty flag to indicate that a call to updateCoeffs is necessary before producing the
  next sample. */
  void setDirty() { dirty = true; }

  /** Updates the coefficients of all the allpass filters according to the user parameters. */
  void updateCoeffs();

  // Embedded objects:
  static const int maxNumStages = 256;
  RAPT::rsBiquadCascade<double, double> allpassChain;

  // User parameters:
  double freqLo;
  double freqHi;
  double freqShape;
  double qLo;
  double qHi;
  double qShape;
  double sampleRate;
  Mode   mode;

  // Internals:
  bool dirty;

};

// ToDo:
// -Maybe introduce global feedback - but then the output will not be white anymore. But maybe we 
//  can use the allpass nesting technique to achieve a white output with feedback? Is it possible
//  (in a practical way) to implement zero-delay feedback for the whole class rsBiquadCascade? If 
//  so, maybe do it.
// -Let the user set up the frequency curve more flexibly by having an arbitrary number of nodes 
//  distributed in 0..1 and map that curve to lowPitch...highPitch. that may be a use case for
//  rsLinearFractionalInterpolator. But maybe do that in a subclass.
// -Maybe have a function rsBiquadCascade::getSampleDF1_1p1z etc. that can be called alternatively
//  if a2 = b2 = 0. Maybe such an optimizing dispatch can be done in rsBiqadCascade itself in a 
//  processBlock function (the dispatcher overhead may make it not worthwhile in getSample but for
//  a whole block, that may be a good optimization)
// -Maybe we can build a phaser from it
// -Maybe create a class rsBrownZapper that also includes the browning filter and the cleanup 
//  highpass and perhaps more post-processing. For example, a slope/tilt filter could be useful for
//  further shaping. Or maybe an 8 band EQ.
// -Make a subclass that also has some post-processing




//=================================================================================================

/** A very simple envelope whose sole purpose is to create a fade-out. This can be used to avoid 
note-off clicks in simple midi-controlled oscillator modules that should not have a full blown
envelope generator built in but still need some way to avoid those pesky clicks. It can be used 
when the osc should neither continue running after receiving a note off nor abruptly turn off (i.e.
switch output to zero) on receiving a note-off.  */

class rsFadeOutEnvelope
{

public:

  rsFadeOutEnvelope()
  {
    reset();
  }

  void setNumFadeSamples(int newNumber)
  {
    numFadeSamples = newNumber;
  }

  bool endIsReached() const { return state == State::silent; }

  void noteOn()
  {
    state = State::open;
  }

  void noteOff()
  {
    state = State::fading;
    remainingFadeSamples = numFadeSamples;
  }

  void goToEnd()
  {
    remainingFadeSamples = 0;
    state = State::silent;
  }

  double getSample()
  {
    switch(state)
    {
    case State::open:    return 1.0;
    case State::silent:  return 0.0;
    case State::fading:
    {
      if(remainingFadeSamples <= 0)
      {
        state = State::silent;
        return 0.0;
      }
      else
      {
        double out = double(remainingFadeSamples) / double(numFadeSamples+1);
        remainingFadeSamples--;
        return out;

        // ToDo: 
        // -Optimize: 
        //  -avoid division by keeping 1.0 / numfadeSamples as member
        //  -avoid int -> double conversion by using double for all members
        //  -after implementing that, run the unit test again
        // -Maybe allow for different shapes. But maybe that should be done in a subclass

      }
    }

    }
  }


  void reset()
  {
    goToEnd();
    //remainingFadeSamples = numFadeSamples;
    //state = State::open;
  }


protected:

  int numFadeSamples = 0;
  int remainingFadeSamples = 0;

  enum class State
  {
    open,      // Volume is fully open
    fading,    // Volume is fading out
    silent     // Volume has reached zero
  };

  State state;

  // For the purposes that I have in mind, float should be good enough for the numFadeSamples and
  // remainingFadeSamples - see:
  // https://stackoverflow.com/questions/3793838/which-is-the-first-integer-that-an-ieee-754-float-is-incapable-of-representing-e
  // The largest exacly representable integer is 2^24 ...but does that mean that all the smaller 
  // ones will be produced when decrementing?
};






//=================================================================================================

/** A class that generates sounds based on frequency sweeps. The original idea was to create sounds
similar to the impulse responses of rsFlatZapper by using a sinusoidal oscillator with a specially
shaped envelope for the instantaneous frequency. The advantages of using an oscillator over an 
allpass are: more direct control, different waveforms possible, stereo phase offsets possible. We 
still want to synthesize bassdrums ...TBC...  */

class rsFreqSweeper
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  rsFreqSweeper();


  //-----------------------------------------------------------------------------------------------
  // \Setup

  /** Sets the sample rate at which this object should operate. */
  void setSampleRate(double newSampleRate) { sampleRate = newSampleRate; /*setDirty();*/ }
  // Maybe we should store the reciprocal because that's what we actually need for the DSP


  /** Sets the frequency to which the lowest allpass filter should be tuned. */
  void setLowFreq(double newFreq) { freqLo = newFreq; setDirty(); }

  /** Sets the frequency to which the highest allpass filter should be tuned. */
  void setHighFreq(double newFreq) { freqHi = newFreq; setDirty(); }

  /** Sets up the time (in seconds) that it will take to sweep down from the highest frequency at
  the start (set by setHighFreq) to some fixed, hardcoded reference frequency (which happens to be
  50 Hz above the low, asymptotic frequency (set by setLowFreq)). The rationale behind choosing 50 
  Hz as reference frequency is that this makes for a nice low-freq target in bassdrums when 
  low-freq is set to zero. */
  void setSweepTime(double newTime) { sweepTime = newTime; setDirty(); }
  // Maybe rename to setSweepTimeInSecs

  void setChirpAmount(double newAmount) { chirpAmount = newAmount; setDirty(); }

  void setChirpShape( double newShape)  { chirpShape  = newShape;  setDirty(); }

  /** Sets the start phase in degrees. */
  void setStartPhase(double newPhase) { phase = newPhase * (1.0/360.0); }

  /** Sets the phase shift between left and right channel in degrees. The left channel will half of
  the shift applied negatively and the right channel half of it positively. */
  void setStereoPhaseShift(double newShift) { phaseStereo = newShift * (1.0/360.0); }

  /** Sets the parameter for the waveshape. The range is -1..+1 where -1: fat sawDown, 0: sine,
  +1: fat saw up. The "fat saw" is like a saw that bulges out. The downward fat saw is the first
  half cycle of the cosine (i.e. from 0 to pi). The upward fat saw is the second half.  */
  //void setWaveShapeParam(double newParam) { waveShape = newParam; }

  /** Sets the function that is used to generate the waveform. If you don't set this up, it will
  use a sine by default. */
  void setWaveForm(const std::function<double(double phase01)> newWaveFunc) { wave = newWaveFunc; }

  /** Initializes all parameter values to their initial/default values. The sample rate may or may
  not be re-initialized also. */
  void initSettings(bool initAlsoSampleRate = false);


  //-----------------------------------------------------------------------------------------------
  // \Inquiry

  double getSampleRate() const { return sampleRate; }


  //-----------------------------------------------------------------------------------------------
  // \Processing

  /** Implements the formula to compute the instantaneous frequency at time t (in seconds). */
  INLINE double getInstFreq(double t)
  {
    // Coefficient update if necessary:
    if(dirty)
      updateCoeffs();

    // Compute instantaneous frequency:
    return (a + b*pow(t, p*q)) / pow(1 + c*pow(t,p), q);
  }
  // ToDo: 
  // -Document how I came up with this formula.
  // -Try to parallelize it. Maybe use float instead of double to get even more bang for the buck.


  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    // Output computation:
    double tmpPhase = instPhase + phase;
    tmpPhase += fbPhsMod * oldOutL; // Experimental
    *inOutL = wave(tmpPhase - 0.5*phaseStereo);
    *inOutR = wave(tmpPhase + 0.5*phaseStereo);
    oldOutL = *inOutL;
    // ToDo: 
    // -parallelize

    // State update:
    sampleCount++;
    double newFreq = getInstFreq(sampleCount / sampleRate);
    RAPT::rsAssert(newFreq <= 0.5 * sampleRate);               // Sanity check for debug
    instPhase += (0.5/sampleRate) * (instFreq + newFreq);      // Trapezoidal integration

    //instPhase += fbPhsMod * *inOutL; // Experimental

    if(instPhase >= 1.0)    // A while-loop would be safer but when we assume sane values for
      instPhase -= 1.0;     // newFreq and instFreq, the "if" should be good enough
    instFreq  = newFreq;

    // ToDo:
    // -Optimize: keep the reciprocal of the sample rate as member to avoid division
  }

  void reset()
  {
    sampleCount = 0;
    instPhase   = 0.0;
    instFreq    = freqHi;

    oldOutL     = 0;
  }


  // Experimental:
  double fbPhsMod = 0;
  double oldOutL  = 0;


protected:

  void setDirty() { dirty = true; }

  void updateCoeffs();

  static const double refFreq; // = 50. Hard-coded reference frequency.

  // User parameters:
  double sampleRate;

  double freqLo;         // Maybe rename to endFreq
  double freqHi;         // Maybe rename to startFreq
  // Rationale for renaming: It's really about the time. One may perhaps even start at a low freq 
  // (not sure, if that would work, though).

  // Rename into chirpAmount, chirpShape
  double chirpAmount;       // Determines (mostly) attack shape
  double chirpShape;        // Determines (mostly) decay shape


  //double refFreq = 50;   // reference frequency (is hardcoded and fixed - todo: use static const)
  double sweepTime;      // Time to sweep down to refFreq in seconds

  double phase;
  double phaseStereo;

  // Wave shape parameters:
  //double waveShape;  // -1: fat saw down, 0: sine, +1 fat saw up


  // Internals:
  bool dirty;
  double a, b, c, p, q;  // Formula parameters

  double instPhase, instFreq;
  int    sampleCount;


  std::function<double(double phase)> wave = [](double p)
  { 
    // Old:
    return RAPT::rsSin<double>(2*PI*p);

    // New:
    //double y = RAPT::rsTriSaw(p, waveShape); // TriSaw
    //return RAPT::rsSin<double>(0.5*PI * y);  // SinSaw

    // ToDo: Verify if "&" is the correct/best capture mode for this purpose. Maybe "=" or "this" 
    // is better? Or maybe we shouldn't capture the waveShape parameter but instead pass it in as
    // an additional parameter
  };
};

// ToDo:
//
// - Maybe the wavesshape could be determined by a user-defined function using std::function. Maybe
//   for flexibility, the function should receive the instantaneous phase and additionaly the 
//   absolute time to (potentially) implement time-varying waveshapes.
//
// - Maybe use a sinusoidally shaped TriSaw waveform with parameter -1..+1. See triSaw() 
//   experiment. Or maybe use a wavetable (in the Waldorf sense). 


//=================================================================================================

/** UNDER CONSTRUCTION. Just a stub at the moment.

A class for phaseshaping. Phaseshaping is like waveshaping but it is applied to a normalized
phase value, i.e. a value between 0 and 1 that is used as an osillator phase. */

class rsPhaseShaper
{

public:


  static double powerLaw(double phase, double shapeParam);
  // let s = shapeParam, 
  // s < 0: sine turns to upward saw
  // s = 0: sine is unchanged
  // s > 0: sine is squeezed
  //
  // Phase should be normalized to the range 0..1, shapeParam in -1..+1. All the different 
  // phase-shaping laws should use that API.

  //static double symmetricLinFracLaw(double phase, double shapeParam);
  // See phaseShapingLinFrac() in OscillatorExperiments.cpp


};



//=================================================================================================

/** A monophonic drum synthesizer based on rsFreqSweeper. The implementation consists mostly of
keeping track of the midi state and delegating setter calls to the embedded rsFreqSweeper object 
at the right moments (e.g. calls to frequency setters deferred to noteOn events and may take into
accoutn key and velocity as modifiers for the actual frequencies etc.).  */

class rsSweepKicker
{

public:

  //-----------------------------------------------------------------------------------------------
  // \Lifetime

  rsSweepKicker();


  //-----------------------------------------------------------------------------------------------
  // \Setup

  // This code is mostly boilerplate:
  void setSampleRate(   double newRate)  
  { 
    freqSweeper.setSampleRate(newRate);
    fadeOutEnv.setNumFadeSamples(RAPT::rsRoundToInt(fadeOutTime * newRate));
  }

  void setLowFreq(      double newFreq)  { frqLo      = newFreq;  }
  void setLowFreqByKey( double newByKey) { frqLoByKey = newByKey; }
  void setLowFreqByVel( double newByVel) { frqLoByVel = newByVel; }

  void setHighFreq(     double newFreq)  { frqHi      = newFreq;  }
  void setHighFreqByKey(double newByKey) { frqHiByKey = newByKey; }
  void setHighFreqByVel(double newByVel) { frqHiByVel = newByVel; }

  //void setSweepTime(    double newTime)  { swpTm = newTime; }  // not needed

  void setSweepTimeInMs(double newTime) { swpTm = 0.001*newTime; }
  // ...ByKy, ByVel

  //void setTransientShape
  //void setFreqDecayShape(double newShape)

  void setChirpAmount(double newAmount) {  freqSweeper.setChirpAmount(newAmount); }

  void setChirpShape( double newShape)  {  freqSweeper.setChirpShape( newShape);  }


  void setStartPhase(      double newPhase) { freqSweeper.setStartPhase(newPhase);       }
  void setStereoPhaseShift(double newShift) { freqSweeper.setStereoPhaseShift(newShift); }

  /** Enumeration of the available wave forms. */
  enum WaveShape // rename to WaveForm
  {
    Sine = 0,             // Just a sine. Ignores waveParam.
    SinFatSaw,            // TriSaw (see below) with sinusoidal waveshaping
    TriSaw,               // SawDown / Triangle / SawUp

    //PowerLaw,             // Phase-shaping with power law - rename to PhaseShapePow
    //PhaseShapeLinFrac,

    NumWaveShapes
  };
  // ToDo: provide more waveshapes - for example PhaseShapeLinFrac

  /** Sets the waveshape to be used. */
  void setWaveForm(int newShape) { waveShape = (WaveShape) newShape; }


  /** Sets the parameter that controls the shape of the waveform. 0 means sine or triangle 
  (depending on waveShape), -1 is something similar to a downward saw and +1 is similar to an 
  upward saw. But the exact shapes will depend on waveShape.  */
  void setWaveFormParameter(double newParam) { waveParam = newParam; }


  // Experimental:
  void setFeedbackPhaseMod(double amount) { freqSweeper.fbPhsMod = amount; }
  // Maybe get rid of this. It seems to make more trouble than its worth. We'll see....


  void initSettings(bool initAlsoSampleRate = false);

  void setFadeOutTimeMs(double newTime)
  {
    fadeOutTime = 0.001 * newTime;
    fadeOutEnv.setNumFadeSamples(RAPT::rsRoundToInt(fadeOutTime * getSampleRate()));
  }

  //-----------------------------------------------------------------------------------------------
  // \Inquiry


  double getSampleRate() const { return freqSweeper.getSampleRate(); }
  // It's stored only in the embedded object to avoid redundancy

  //-----------------------------------------------------------------------------------------------
  // \Processing

  INLINE void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    if(fadeOutEnv.endIsReached())
    {
      *inOutL =  *inOutR = 0.0;
      return;
    }

    freqSweeper.getSampleFrameStereo(inOutL, inOutR);
    double env = fadeOutEnv.getSample();
    *inOutL *= env;
    *inOutR *= env;
  }

  void reset()
  {
    freqSweeper.reset();
    fadeOutEnv.reset();    // Will go into silent state
    currentNote = -1;
  }

  void noteOn( int key, int vel);

  void noteOff(int key);


protected:

  // The embedded objects:
  rsFreqSweeper     freqSweeper;
  rsFadeOutEnvelope fadeOutEnv;

  // We mostly have the same parameters as the embedded rsFreqSweeper object. However - here, the
  // parameters also repsond to key and velocity.

  double frqLo, frqLoByKey, frqLoByVel;
  double frqHi, frqHiByKey, frqHiByVel;
  double swpTm, swpTmByKey, swpTmByVel;

  double fadeOutTime;

  double    waveParam;    // Parameter to control/morph the waveshape in -1..+1
  WaveShape waveShape;    // Select the type of morphable waveshape



  int currentNote = -1;

  // ...tbc... have two waveShape parameters - one that morphs between sawUp/sin/sawDown, one that
  // drives the result into saturation (to squarify it) ...maybe a third that addds an offset after
  // the drive before the saturating function. Maybe the API should allow rsFreqSweeper to take a
  // std::function for producing the waveform - and we provide one here.

  // Maybe integrate an AHDSR envelope for the amplitude

};




}