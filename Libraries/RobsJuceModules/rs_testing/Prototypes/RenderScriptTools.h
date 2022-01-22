#pragma once

//=================================================================================================

/** A drum synthesizer based on frequency sweeps. It provides envelopes for the frequency, 
waveshape and amplitude of an oscillator. Each of these envelopes is formed as a weighted sum of
exponential decays ...tbc... */

template<class T>
class rsSweepDrummer
{


public:

  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; dirty = true;  }

  void setFreqDecay( int index, T newDecay)  { freqDecays[index]  = newDecay;  dirty = true; }
  void setFreqWeight(int index, T newWeight) { freqWeights[index] = newWeight;  }

  void setAmpDecay(  int index, T newDecay)  { ampDecays[index]   = newDecay;  dirty = true; }
  void setAmpWeight( int index, T newAmp)    { ampWeights[index]  = newAmp;  }

  void setFreqDecays(T decay1, T decay2, T decay3)
  {
    freqDecays[0] = decay1;
    freqDecays[1] = decay2;
    freqDecays[2] = decay3;
    dirty = true; 
  }

  void setAmpDecays(T decay1, T decay2, T decay3)
  {
    ampDecays[0] = decay1;
    ampDecays[1] = decay2;
    ampDecays[2] = decay3;
    dirty = true; 
  }

  void setFreqWeights(T weight1, T weight2, T weight3)
  {
    freqWeights[0] = weight1;
    freqWeights[1] = weight2;
    freqWeights[2] = weight3;
  }

  void setAmpWeights(T weight1, T weight2, T weight3)
  {
    ampWeights[0] = weight1;
    ampWeights[1] = weight2;
    ampWeights[2] = weight3;
  }

  void setFreqDepth(T newDepth) { freqDepth = newDepth; }


  // trigger...
  void noteOn(int key, int vel)
  {
    sampleCount = 0;
    phase = 0;
    for(int i = 0; i < numEnvs; i++) {
      freqEnvs[i].noteOn(key, vel);
      ampEnvs[i].noteOn( key, vel);  }
  }


  inline void processFrame(T* L, T* R);


protected:

  void updateCoeffs();

  static const int numEnvs = 3;

  // Frequency envelope parameters:
  T freqDecays[numEnvs]  = { 10,  80, 240 };   // in ms
  T freqWeights[numEnvs] = {   6, -1,   1 };   // raw factor
  T freqFloor =    0;                          // in Hz
  T freqDepth =  150;                          // in Hz

  // Amplitude envelope parameters:
  T ampDecays[numEnvs]  = { 20,   100, 400 };  // in ms
  T ampWeights[numEnvs] = {  0.75, -1,   1 };  // raw factor

  // Waveshape envelope parameters:
  //...have two parameters: one that skews the sine into a saw via phase-shaping and one that 
  // drives the sine into a square via waveshaping/saturation. Both should be envelope controlled.
  // maybe use an envelope that is derived from freq and/or amp env?
  //

  T sampleRate = 44100;
  bool dirty = true;

  int sampleCount = 0;   // sample counter
  T   phase       = 0;   // phase in waveform in nominal range 0...2*pi but can also go beyond

  // ToDo:
  //T timeScale = 1;       // scales the envelope times
  //T freqScale = 1;       // scales the oscillator frequencies

  RAPT::rsAttackDecayEnvelope<T> freqEnvs[3];
  RAPT::rsAttackDecayEnvelope<T> ampEnvs[3];
};

template<class T>
inline void rsSweepDrummer<T>::processFrame(T* L, T* R)
{
  if(dirty)
    updateCoeffs();

  // compute outputs of envelopes:
  T freq = freqFloor;
  T amp  = 0;
  for(int i = 0; i < numEnvs; i++)
  {
    freq += freqWeights[i] * freqEnvs[i].getSample() * freqDepth;
    amp  += ampWeights[i]  * ampEnvs[i].getSample();
  }

  *L = amp * sin(phase - PI/4);
  *R = amp * sin(phase + PI/4);
  // ToDo: Use wave(phi +- PI/4, waveParam1, waveParam2) where param2 controls phase-shaping 
  // (sin-to-saw) and param2 waveshaping (sin/saw-to-square). Both parameters should also have an
  // envelope.


  // state update:
  T w = 2*PI*freq/sampleRate;
  phase += w;
  sampleCount++;
}

template<class T>
inline void rsSweepDrummer<T>::updateCoeffs()
{
  for(int i = 0; i < numEnvs; i++)
  {
    freqEnvs[i].setAttackSamples(0.0);
    ampEnvs[i].setAttackSamples(0.0);
    freqEnvs[i].setDecaySamples(freqDecays[i] * 0.001 * sampleRate);
    ampEnvs[i].setDecaySamples( ampDecays[i]  * 0.001 * sampleRate);
    // todo: update envelopes for waveform parameters
  }
  dirty = false;
}

//=================================================================================================

std::vector<double> randomizePhases(const std::vector<double>& x, int seed, double amount);

void rsNormalizeJointly(std::vector<double>& x, std::vector<double>& y);
/*
{
  double maxX = RAPT::rsMaxAbs(x);
  double maxY = RAPT::rsMaxAbs(y);
  double scaler = 1.0 / RAPT::rsMax(maxX, maxY);
  RAPT::rsScale(x, scaler);
  RAPT::rsScale(y, scaler);
}
*/


//=================================================================================================