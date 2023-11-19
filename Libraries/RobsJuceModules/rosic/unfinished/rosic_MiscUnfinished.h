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

  void setFrequencyScale(double scaleFactor) { shift = scaleFactor; }

  enum class Algorithm
  {
    RobSchmt,   // Robin Schmidt
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

  void shiftViaRS(Complex* spectrum, int spectrumSize);
  // stub - Custom algorithm by Robin Schmidt



  Complex *tmpSpectrum;

  double shift = 1.0;  // rename to scale or freqScale

  Algorithm algo = Algorithm::JuilHirs;

  int frameIndex = 0;
  // Try to get rid of this. The JuilHirs algo needs this for its formula but maybe the formula can
  // be expressed without it. Try to figure this out! If it turns out to be so important that it
  // can't be removed, we may need to take care of letting it wrap around back to zero at 
  // appropriate instants. The paper says in 3.3 that vertical phase coherence is achieved every
  // O frames where O is the overlap factor. So maybe that means we can wrap around frameIndex at 
  // every multiple of O?

  PhaseFormula phaseFormula = PhaseFormula::useMultiplier;

};





}