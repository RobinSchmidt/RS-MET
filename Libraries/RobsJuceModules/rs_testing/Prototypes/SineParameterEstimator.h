#ifndef RAPT_SINEPARAMETERESTIMATOR_H_INCLUDED
#define RAPT_SINEPARAMETERESTIMATOR_H_INCLUDED
// todo: rename files and change the #define to reflect the new class name

/** A class for modeling any signal x[n] as a single sinusoid with time-varying amplitude and 
frequency and/or phase, like:

  x[n] = a[n] * sin( p[n] )

with the instantaneous amplitude a[n] and instantaneous phase p[n]. The latter can be given in 3 
ways: (1) directly, (2) as integral (represented as cumulative sum) over an instantaneous radian
frequency w[k]: p[n] = sum_{k=0}^n w[k] or (3) a combination of (2) and an instantaneous 
phase-modulation term: p[n] = pm[n] + sum_{k=0}^n w[k]. It contains functions to synthesize the 
signal from the instantaneous amp/freq/phase arrays and - more importantly - to analyze a given 
signal to produce these arrays. To this end, various algorithms are available...tbc...  */

template<class T>
class rsSingleSineModeler
{

public:

  /** The various options for the analysis algorithm. This determines the order in which the 
  2 or 3 instantaneous parameters are estimated. ...ToDo: explain this more - what effects do 
  the different choices have etc....estimating amp first gives nictes/smoothest amp-env estimate, 
  etc... */
  enum class Algorithm
  {
    ampViaPeaks,       /**< First estimates amp-env from peaks, then phase, then freq. */
    freqViaZeros,      /**< First estimates freq via zero crossings, then phase, then amp. */
    freqViaFormula     /**< First estimates freq via formula, then phase, then amp. */
  };

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setFreqSmoothing(int medianOrder, int averageOrder)
  {
    freqMedianOrder  = medianOrder;
    freqAverageOrder = averageOrder;
  }
  // when both are zero, the phase-modulation signal will come out as zero - this sort of 
  // determines the split between what is modeled as freq-modulation and what is modeled as 
  // phase-modulation, if both are used
  // maybe allow the user to specify a custom function, what to do with the w-array before 
  // computing the phase-mod array

  /** Sets the algorithm that is used to estimate the instantaneous parameters of a given 
  signal..... */
  void setAnalysisAlgorithm(Algorithm newAlgo) { algo = newAlgo; }


  //-----------------------------------------------------------------------------------------------
  /** \name Analysis */


  void analyzeAmpAndPhase(const T* x, int N, T* a, T* p) const;

  void analyzeAmpAndFreq(const T* x, int N, T* a, T* w) const;

  void analyzeAmpFreqAndPhaseMod(const T* x, int N, T* a, T* w, T* pm) const;


  //-----------------------------------------------------------------------------------------------
  /** \name Synthesis */

  static void synthesizeFromAmpAndPhase(const T* a, const T* p, int N, T* y);

  static void synthesizeFromAmpAndFreq(const T* a, const T* w, int N, T* y);

  static void synthesizeFromAmpFreqAndPhaseMod(const T* a, const T* w, const T* pm, int N, T* y);


  //-----------------------------------------------------------------------------------------------
  /** \name Static functions */

  /** There's a recursion formula for the sine with normalized radian frequeny w: 
    y[n] = a1*y[n-1] - y[n-2] 
  where 
    a1 = 2*cos(w) 
  and the states y[n-1], y[n-2] are initialized as: 
    y[n-1] = A * sin(p - w), y[n-2] = A * sin(p - 2*w) 
  which in our notation here translates to yR = a1*yC - yL. This leads to 
    a1 = (yL+yR)/yC and
    w  = acos(a1/2). 
  This formula for w is implemented here. Note that we don't check against division by zero, so yC
  should be large enough. However, we do check, if the input to acos is in -1..+1, so the formula 
  is "half-safe". */
  static T omegaFormula(T yL, T yC, T yR) 
  { return acos(rsClip(T(0.5)*(yL+yR)/yC, T(-1), T(+1))); }
  // maybe rename to freqFormula

  /** Estimates the instantaneous normalized radian frequencies ("omega") of the signal x via the
  recursion formula for 3 successive samples of a sinewave. To estimate the omega at sample n, it 
  looks at x[n-1], x[n], x[n+1] and applies the formula. But because this formula is unreliable 
  near zero-crossings, it will also clean up the result by using a weighted average of the so 
  found omegas over 3 samples, using weights determined by a reliability measure based on how close
  a sample is to zero. */
  static void sigToOmegasViaFormula(const T* x, int N, T* w);
  // rename to sigToFreq...

  /** Estimates the amplitude envelope of the signal x via coennecting peaks with linear 
  interpolants and writes the result to a. */
  static void sigToAmpsViaPeaks(const T* x, int N, T* a, int precision = 1);
  // todo: document, if x == a is allowed (i think so - but only, if precision <= 1)

  /** Given a signal x and an array of instantaneous amplitudes a, this function computes the 
  corresponding instantaneous pahses, such that x[n] = a[n] * sin(p[n]) for each n. */
  static void sigAndAmpToPhase(const T* x, const T* a, int N, T* p);


  static void sigAndFreqToPhaseAndAmp(const T* x, const T* w, int N, T* p, T* a);

  static void sigAndFreqToAmp(const T* x, const T* w, int N, T* a);


  static void phaseToFreq(const T* p, int N, T* w);



  //static void freqToPhase(const T* w, int N, T* p, bool wrap);



  static void phaseAndFreqToPhaseMod(const T* p, const T* w, int N, T* pm);





  // ToDo: sigToOmegasViaZeros

  static void exactPeakPositionAndHeight(
    const T* x, int N, int n0, int precision, T* pos, T* height);
  // move to somewhere else



  static void smoothFreqs(T* w, int N, int medianOrder, int averageOrder);



protected:

  static void connectPeaks(const T* x, int N, T* env);
  // y == a is allowed - it can overwrite the content of a given array
  // maybe move this function to somewhere else - this could be useful in various other scenarios

  static void connectPeaks(const T* xt, const T* xi, int N, T* env, int precision);
  // under construction - uses a polynomial of order 2*precision to estimate the actual locations 
  // and heights of the peaks - using a parabolo already improves results, but there are still
  // frequency jaggies, so we may need higher accuracy for the amp-env
  // we need two inputs - one to determine the peak locations and one for using in the 
  // interpolator

 
  /** When we compute the instantaneous phase from a known signal value x[n] and its instantaneous
  amplitude a[n] via p[n] = asin(x[n] / a[n]), the returned result from asin is always in the range
  -pi/2...+pi/2. When x is a sinewave, instead of sweeping from -pi to +pi and then wrapping around 
  in one cycle, it oscillates back and forth between -pi/2...+pi/2. This function takes a raw array
  of such phase values and heuristically reflects the phases around pi/2 or -pi/2 to get rid of 
  that effect. It's sort of similar to phase-unwrapping. Used in sigAndAmpToPhase.  */
  static void unreflectPhase(const T* x, T* p, int N);
  // more research necessarry to figure out what is the best algorithm for this - this here was the 
  // first one that sort of worked for the bandpass-noise


  // todo: have enum-class members for:
  // algo: freqAmpPhase, ampPhaseFreq, ...decides what is estimated first, second, third
  // freqAlgo: zeros, formula, ...
  // ampAlgo: peaks, hilbert, ...

  int freqMedianOrder  = 1;
  int freqAverageOrder = 3;

  int ampEnvPrecision  = 1;

  Algorithm algo = Algorithm::ampViaPeaks;

};

#endif