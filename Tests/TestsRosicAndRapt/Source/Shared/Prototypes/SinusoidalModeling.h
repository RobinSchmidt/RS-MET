#pragma once


template<class T>
std::vector<T> synthesizeSinusoidal(const RAPT::rsSinusoidalModel<T>& model, T sampleRate);


template<class T>
RAPT::rsSinusoidalModel<T> analyzeSinusoidal(T* sampleData, int numSamples, T sampleRate);




/** Prototype of STFT based sinusoidal analyzer */

template<class T>
class SinusoidalAnalyzer
{

public:

  inline void setBlockSize(int newBlockSize)      { sp.setBlockSize(newBlockSize); }
  inline void setHopSize(int newHopSize)          { sp.setHopSize(newHopSize); }
  inline void setZeroPaddingFactor(int newFactor) { sp.setZeroPaddingFactor(newFactor); }
  inline void setAnalysisWindowType(int newType)  { sp.setAnalysisWindowType(newType); }
  //setRootKey/setFundamentalFrequency, 

  RAPT::rsSinusoidalModel<T> analyze(T* sampleData, int numSamples, T sampleRate) const;

protected:

  RAPT::rsSpectrogram<T> sp;   // spectrogram processor

};