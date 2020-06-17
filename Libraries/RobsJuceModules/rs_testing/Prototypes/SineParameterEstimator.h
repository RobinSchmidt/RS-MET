#ifndef RAPT_SINEPARAMETERESTIMATOR_H_INCLUDED
#define RAPT_SINEPARAMETERESTIMATOR_H_INCLUDED

/** A class for estimating the instantaneous parameters (frequency and/or phase, amplitude) of a 
sinewave. */

template<class T>
class rsSineParameterEstimator
{

public:



  //void getAmpAndPhase(const T* x, int N, T* a, T* p);

  //void getAmpFreqAndPhase(const T* x, int N, T* a, T* w, T* p);

  //-----------------------------------------------------------------------------------------------
  /** \name Static functions */

  static void sigToFreqsViaRecursion(const T* x, int N, T* w);

  static void sigToAmpsViaPeaks(const T* x, int N, T* a);

  // ToDo: sigToFreqsViaZeros

protected:

};

#endif