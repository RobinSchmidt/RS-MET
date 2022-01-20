#pragma once

/** A drum synthesizer based on frequency sweeps. It provides envelopes for the frequency, 
waveshape and amplitude of an oscillator. Each of these envelopes is formed as a weighted sum of
exponential decays ...tbc... */

template<class T>
class rsSweepDrummer
{


public:

  void setSampleRate(T newSampleRate) { sampleRate = newSampleRate; dirty = true;  }

  void setFreqDecay( int index, T newDecay)  { freqDecays[index]  = newDecay;  dirty = true; }
  void setAmpDecay(  int index, T newDecay)  { ampDecays[index]   = newDecay;  dirty = true; }
  void setFreqWeight(int index, T newWeight) { freqWeights[index] = newWeight; dirty = true; }
  void setAmpWeight( int index, T newAmp)    { freqAmp[index]     = newAmp;    dirty = true; }

  inline void processFrame(T* L, T* R);


protected:

  void updateCoeffs();

  static const int numEnvs = 3;

  // Frequency envelope parameters:
  T freqDecays[numEnvs]  = { 10,  80, 240 };   // in ms
  T freqWeights[numEnvs] = {   6, -1,   1 };   // raw factor
  T freqFloor =    0;                          // in Hz
  T freqDepth =  800;                          // in Hz

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

  RAPT::rsAttackDecayEnvelope<T> freqEnvs[3];
};

template<class T>
inline void rsSweepDrummer<T>::processFrame(T* L, T* R)
{
  if(dirty)
    updateCoeffs();

}

template<class T>
inline void rsSweepDrummer<T>::updateCoeffs()
{
  for(int i = 0; i < numEnvs; i++)
  {
    freqEnvs[i].setDecaySamples(freqDecays[i] * 0.001 * sampleRate);
    ampEnvs[i].setDecaySamples( ampDecays[i]  * 0.001 * sampleRate);
  }
  dirty = false;
}