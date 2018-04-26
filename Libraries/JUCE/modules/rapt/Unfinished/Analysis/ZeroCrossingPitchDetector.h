#ifndef RAPT_PITCHDETECTOR_H
#define RAPT_PITCHDETECTOR_H

/** This class estimates the fundamental of an incoming signal (which is assumed to be pitched and
monophonic) by first suppressing the formants via an adpative filter, then lowpass-filtering it and 
measuring the distance between two successive upward zero crossings in the filtered signal.  */

template<class T>
class rsZeroCrossingPitchDetector
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsZeroCrossingPitchDetector();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(T newSampleRate);

  /** Sets the minimum expected fundamental frequency in Hz. */
  void setMinFundamental(T newMinFundamental);

  /** Sets the minimum expected fundamental frequency in Hz. */
  void setMaxFundamental(T newMaxFundamental);


  /** \name Inquiry */

  /** Returns the minimum fundamental frequnecy which the algorithm expects. */
  T getMinFundamental() const { return minFundamental; }

  /** Returns the maximum fundamental frequnecy which the algorithm expects. */
  T getMaxFundamental() const { return maxFundamental; }

  /** Returns the current reliability of the period estimate (and all related estimates) as a
  value in percent - 100% means fully reliable, 0% means fully unreliable. The reliability is
  assessed via a cross-correlation measurement of two successive cycles (or non-cycles). */
  T getReliability() const { return 100.0 * cycleCorrelation; }


  /** \name Audio Processing */

  /** Estimates and returns the fundamental frequency of the incoming signal (in seconds). */
  RS_INLINE T estimatePeriod(T inputSignal);

  /** Overloaded method for cases in which a whitened signal can be provided by the caller -
  in this case, the whitening here will be bypassed. */
  RS_INLINE T estimatePeriod(T inputSignal, T whiteSignal);

  /** Estimates and returns the fundamental frequency of the incoming signal (in Hz). */
  RS_INLINE T estimateFundamentalFrequency(T inputSignal);

  /** Estimates and returns the pitch of the incoming signal (in terms of a midi note number). */
  RS_INLINE T estimatePitch(T inputSignal);


  /** \name Misc */

  /** Resets the object back into its initial state. You can pass a value for the initial
  estimate of the fundamental frequency. This initial estimate will be returned until it
  eventually settles to an actually valid measurement (which needs "warm-up" time). */
  void reset(T initialEstimate = 1000.0);


  rsFormantRemover formantRemover; // move to the protceted area?

protected:

  /** Estimates the fractional part of the zero crossing by fitting a cubic polynomial through
  the most recent four samples (including the current one) and solving for the zero-crossing of
  this cubic polynomial. */
  T getFractionalPartOfZeroCrossing();


  T sampleRate, sampleRateRec; // sample rate and its reciprocal. 

  T minFundamental, maxFundamental, minPeriod, maxPeriod;
    // Minimum and maximum values for the fundamental frequency and period.

  T y0, y1, y2, y3;
    // 3 previous and the current value of the lowpass-filtered input signal - need to be 
    // remembered between calls to estimateFundamental in oder to detect zero-crossings. */

  T fracOld;
    // The old (from the previous iteration) value of the fractional part of the zero 
    // crossing - together with the current value of this, a non-integer period estimate can be 
    // obtained.

  T periodEstimate, frequencyEstimate; // periodEstimateSmoothed;
    // This member always contains the current estimate of the pitch period the corresponding
    // fundamental frequency - it will be updated when upward zero-crossings are detected.

  T cycleCorrelation;
    // This keeps track of the cross-correlation between two successive cycles.

  rsOnePoleFilter<T, T>  dcBlocker;
    // The dc-blocker filter for the signal in which we look for zero-crossings.

  rsFourPoleFilter<T, T> lowpass;
    // The lowpass filter for the signal in which we look for zero-crossings.

  rsEnvelopeFollower<T, T> envFollower;
    // detects the loudness of the signal

  rsCyclicAutoCorrelator<T> correlator;
    // measures the correlation between successive pitch-cycles

  int sampleCounter;
    // counts the number of samples between two succesive upward zero crossings

};


//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
RS_INLINE T rsZeroCrossingPitchDetector<T>::estimatePeriod(T inputSignal)
{
  // remove (or at least, weaken) the formants from the input signal:
  T whiteSignal = formantRemover.getSample(inputSignal);

  return estimatePeriod(inputSignal, whiteSignal);
}

template<class T>
RS_INLINE T rsZeroCrossingPitchDetector<T>::estimatePeriod(T inputSignal, T whiteSignal)
{
  T tmpCandidate, filteredSignal, env, frac;
  bool   isAboveThreshold;

  sampleCounter++;

  // feed the input signal (as is, unwhitened) into the correlation measurement device:
  correlator.acceptSample(inputSignal);

  // get the current filtered sample:
  filteredSignal = dcBlocker.getSample(lowpass.getSample(whiteSignal));

  // update the remembered previous values fo the filtered signal:
  y3 = y2;
  y2 = y1;
  y1 = y0;
  y0 = filteredSignal;

  env = envFollower.getSample(fabs(inputSignal));
  isAboveThreshold = env > rsDB2amp(-60.0);

  if(y2 < 0.0 && y1 >= 0.0)
  {
    // we have detected an upward zero crossing in the lowpass filtered signal - that means the 
    // signal has just completed a period. therefore we need to update our current estimate of 
    // it....

    // get the fractional part of the zero-crossing (where time is interpreted as having its 
    // origin two samples ago):
    frac = getFractionalPartOfZeroCrossing();

    // estimate the period (taking into account the fractional part from this iteration and the 
    // previous one):
    tmpCandidate = (sampleCounter+frac-fracOld) * sampleRateRec;
    tmpCandidate = rsClipToRange(tmpCandidate, minPeriod, maxPeriod);

    // remember the current fractional part for the next time:
    fracOld = frac;

    // get the cycle correlation as a measure of reliability of the period estimate (between 
    // 0...1) - if the estimate is reliable enough, accept it:
    cycleCorrelation = correlator.getCorrelationAndReset();
    if(isAboveThreshold && cycleCorrelation > 0.9)
      periodEstimate = tmpCandidate;

    // reset the sample-counter:
    sampleCounter = 0;

    // adapt the lowpass and dc-blocking highpass:
    frequencyEstimate = 1.0 / periodEstimate;
    //dcBlocker.setCutoff(0.25*freq);
    lowpass.setFrequency(0.5*frequencyEstimate);
    lowpass.updateFilterCoefficients();
  }

  return periodEstimate;
}

template<class T>
RS_INLINE T rsZeroCrossingPitchDetector<T>::estimateFundamentalFrequency(T inputSignal)
{
  estimatePeriod(inputSignal);
  return frequencyEstimate;
}

template<class T>
RS_INLINE T rsZeroCrossingPitchDetector<T>::estimatePitch(T inputSignal)
{
  return rsFreqToPitch(estimateFundamentalFrequency(inputSignal));
    // this is inefficient (calls freqToPitch per sample, where it could have been called
    // only when the pitch has changed) 
}

template<class T>
RS_INLINE T rsZeroCrossingPitchDetector<T>::getFractionalPartOfZeroCrossing()
{
  // compute the zero crossing of a straight line y=a*x+b: x0 = -b/a, y2 is two samples ago, y1
  // is one sample ago:
  T fracLinear = y2 / (y2-y1);
    // is this really correct? verify!

  // for a faster but less precise version, we could directly return fracLinear now
  // ->wrap this into a function getApproximateFractionalPartOfZeroCrossing

  // obtain coefficients for a cubic polynomial through the most recent 4 samples (where time is 
  // interpreted as having its origin two samples ago):
  T d = y2;
  T c = y1 - (1.0/3.0)*y3 - 0.5*y2 - (1.0/6.0)*y0;
  T b = 0.5*(y3+y1) - y2;
  T a = 0.5*(y2-y1) + (1.0/6.0)*(y0-y3);

  // with these coefficients for a cubic polynomial, we now compute the zero-crossing (root) of 
  // that polynomial via Newton-Raphson iteration, using the zero-crossing from the fitted line
  // as initial estimate ...and return that value:
  return rsPolynomial::getCubicRootNear(fracLinear, a, b, c, d, 0.0, 1.0);
}

#endif
