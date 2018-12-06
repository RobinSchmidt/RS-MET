#ifndef RAPT_SINUSOIDALMODEL_H
#define RAPT_SINUSOIDALMODEL_H

//=================================================================================================

template<class T>
class rsInstantaneousSineParams
{

public:


  /** Constructor. */
  //rsInstantaneousSineParams();


  /** \name Inquiry */

  inline T getTime() const { return time; }

  inline T getFrequency() const { return freq; }

  inline T getAmplitude() const { return gain; }

  inline T getWrappedPhase() const { return phase; }

  inline T getUnwrappedPhase() const { return phase + 2*PI*cycles; }


protected:

  T time = 0;      // in seconds
  T freq = 0;      // in Hz
  T gain = 0;      // as raw factor
  T phase = 0;     // radians in [0,2*pi) ...really - or is it [-pi, pi)? ...check atan2
  int cycles = 0;  // number of cycles passed: unwrapped phase = phase + 2*pi*cycles

};


//=================================================================================================

template<class T>
class rsSinusoidalPartial
{

public:

  inline T getStartTime() const { return instParams[0].getTime(); }

  inline T getEndTime() const { return instParams[instParams.size()-1].getTime(); }


  //rsInstantaneousSineParams<T> getInstantaneousParameters() const;

  //std::vector<T> renderPartial(); // maybe it should take a (sub-sample) time-offset?


protected:

  std::vector<rsInstantaneousSineParams<T>> instParams;

};


//=================================================================================================

template<class T>
class rsSinusoidalModel
{

public:




  /** \name Setup */

  //inline void setBlockSize(int newBlockSize) const;

  //inline void setHopSize(int newHopSize) const;

  //setRootKey/setFundamentalFrequency, 






  /** \name Processing */

  void analyze(T* sampleData, int numSamples);

  std::vector<T> synthesize() const;

protected:

  std::vector<rsSinusoidalPartial<T>> partials;


  rsPhaseVocoder<T> phsVoc;

};



#endif