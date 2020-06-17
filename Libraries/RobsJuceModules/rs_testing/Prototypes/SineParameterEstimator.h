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

  /** Estimates the instantaneous normalized radian frequencies ("omega") of the signal x via the
  recursion formula for 3 successive samples of a sinewave. To estimate the omega at sample n, it 
  looks at x[n-1], x[n], x[n+1] and applies the formula. But because this formula is unreliable 
  near zero-crossings, it will also clean up the result by using a weighted average of the so 
  found omegas over 3 samples, using weights determined by a reliability measure based on how close
  a sample is to zero. */
  static void sigToOmegasViaFormula(const T* x, int N, T* w);

  static void sigToAmpsViaPeaks(const T* x, int N, T* a);

  // ToDo: sigToOmegasViaZeros

protected:

};

#endif