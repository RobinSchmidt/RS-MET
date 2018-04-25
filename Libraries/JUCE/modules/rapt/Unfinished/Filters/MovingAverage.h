#ifndef RAPT_MOVINGAVERAGE_H
#define RAPT_MOVINGAVERAGE_H

/** This is a recursive implementation of a moving average filter as decribed in "The Scientist and
Engineers Guide to DSP" page 283. To counteract accumulation of roundoff error in a floating point 
arithmetics implementation, a kind of leakage for the error was implemented by using a feedback 
coefficient slightly less than unity (and compensating the side-effects with the 
feedforward-coefficient and output gain). */

template<class TSig, class TPar>
class rsMovingAverage
{

public:

  /** Constructor. */
  rsMovingAverage();


  /** Sets the sample rate in Hz. */
  void setSampleRate(TPar newSampleRate);

  /** Sets the length of the impulse-response in seconds. */
  void setLengthInSeconds(TPar newLength);

  /** Sets the length of the impulse-response as a number of samples. */
  void setLengthInSamples(int newLength);

  /** Computes the leakage coefficient a1 from a specified maximum deviation of the
  impulse-response from the ideal flat line. The non-ideal impulse response is almost flat but
  shows a decay. The deviation is defined as the maximum relative overshoot over the flat line.
  For example, with a length of 10 samples, the ideal response would be 10 samples of 0.1. A
  non-ideal response that overshoots this 0.1 by 0.01, such that the response starts with 0.11
  is defined to have a deviation of 0.1 - or 10 percent of the ideal value. The higher the value,
  the faster the accumulated error will die away (good) but it will also drive the impulse
  response farther away from an ideal moving averager (bad). */
  void setDeviation(TPar newDeviation);

  /** Calculates a single filtered output-sample. */
  RS_INLINE TSig getSample(TSig in);

  /** Resets the internal buffers. */
  void reset();

protected:

  /** \name Internal Functions */

  /** Updates the gain g and coefficient for delayed input bN from a1 and the length of the
  impulse response. */
  void updateCoefficients();


  /** \name Data */

  rsBasicDelayLine delayLine;

  TPar a1, bN, g;  // coefficients and output gain
  TSig y1;         // y[n-1]
  TPar d;          // maximum relative deviation from the ideal constant response
  TPar length;     // length of impulse response in seconds
  TPar sampleRate; // sample rate in Hz

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsMovingAverage<TSig, TPar>::getSample(TSig in)
{
  y1 = in + a1*y1 - bN*delayLine.getSample(in);
  return g * y1;
}

#endif
