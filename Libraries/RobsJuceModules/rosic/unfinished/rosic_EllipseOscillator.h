#ifndef rosic_EllipseOscillator_h
#define rosic_EllipseOscillator_h

namespace rosic
{


class rsEllipseOscillator : public RAPT::rsEllipseOscillator<double>
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setFrequency(double newFrequency);

  void setSampleRate(double newSampleRate);

  /** Sets a detuning in semitones. */
  inline void setDetune(double newDetune) 
  { 
    tuneFactor = RAPT::rsPitchOffsetToFreqFactor(newDetune); 
    updateOmega(); 
  }

  /** Sets an additive offset/shift for the frequency. */
  inline void setFrequencyShift(double newShift) 
  { 
    freqShift = newShift; 
    updateOmega(); 
  }

  /** Sets amplitude of this oscillator (as raw multiplier). */
  inline void setAmplitude(double newAmplitude) 
  { 
    amplitude = newAmplitude; 
  }




  //setFrequencyScaleY, setFrequencyShiftY


  //-----------------------------------------------------------------------------------------------
  // \name Audio Processing

  inline void getSampleFrameStereo(double* left, double* right)
  {
    //*left = *right = amplitude * getSample(); return; // test

    getSamplePair(left, right);
    *left  *= amplitude;
    *right *= amplitude;
  }

protected:

  INLINE void updateOmega() { setOmega(omegaFactor*(tuneFactor*frequency+freqShift)); }

  double frequency   = 1000;
  double omegaFactor = 2*PI/44100;  // to convert from frequency to normalized radian frequency
  double tuneFactor  = 1;
  double freqShift   = 0;
  double amplitude   = 1;
  double midSide     = 0;  // not yet used

};


//=================================================================================================

class rsEllipseOscillatorVoice // : public rsPolyVoice  (maybe)
{

protected:

  rsEllipseOscillator* master;

};

class rsEllipseOscillatorPoly : public rsEllipseOscillator // public rsPolyModule
{


protected:

  std::vector<rsEllipseOscillatorVoice> voices;

};

}

#endif