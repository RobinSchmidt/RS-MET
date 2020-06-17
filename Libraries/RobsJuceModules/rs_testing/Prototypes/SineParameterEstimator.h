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

  /** Estimates the instantaneous normalized radian frequencies ("omega") of the signal x via the
  recursion formula for 3 successive samples of a sinewave. To estimate the omega at sample n, it 
  looks at x[n-1], x[n], x[n+1] and applies the formula. But because this formula is unreliable 
  near zero-crossings, it will also clean up the result by using a weighted average of the so 
  found omegas over 3 samples, using weights determined by a reliability measure based on how close
  a sample is to zero. */
  static void sigToOmegasViaFormula(const T* x, int N, T* w);

  static void sigToAmpsViaPeaks(const T* x, int N, T* a);
  // todo: document, if x == a is allowed (i think so)

  // ToDo: sigToOmegasViaZeros

protected:

  static void connectPeaks(const T* y, int N, T* a);
  // rename: y to x, a to env
  // y == a is allowed - it can overwrite the content of a given array
  // maybe move this function to somewhere else - this could be useful in various other scenarios

};

#endif