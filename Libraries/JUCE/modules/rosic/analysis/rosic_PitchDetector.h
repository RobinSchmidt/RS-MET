#ifndef rosic_PitchDetector_h
#define rosic_PitchDetector_h

//// rosic-indcludes:
//#include "../filters/rosic_FourPoleFilter.h"
//#include "../filters/rosic_OnePoleFilter.h"
//#include "../math/rosic_PolynomialAlgorithms.h"
//#include "rosic_CyclicAutoCorrelator.h"
//#include "rosic_EnvelopeFollower.h"
//#include "rosic_FormantRemover.h"

namespace rosic
{

  /**

  This class estimates the fundamental of an incoming signal (which is assumed to be pitched and 
  monophonic) by first suppressing the formants via an adpative filter, then lowpass-filtering it 
  and measuring the distance between two successive upward zero crossings in the filtered signal.

  */

  class PitchDetector
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    PitchDetector();   


    /** \name Setup */

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the minimum expected fundamental frequency in Hz. */
    void setMinFundamental(double newMinFundamental);

    /** Sets the minimum expected fundamental frequency in Hz. */
    void setMaxFundamental(double newMaxFundamental);

        
    /** \name Inquiry */

    /** Returns the minimum fundamental frequnecy which the algorithm expects. */
    double getMinFundamental() const { return minFundamental; }

    /** Returns the maximum fundamental frequnecy which the algorithm expects. */
    double getMaxFundamental() const { return maxFundamental; }

    /** Returns the current reliability of the period estimate (and all related estimates) as a 
    value in percent - 100% means fully reliable, 0% means fully unreliable. The reliability is 
    assessed via a cross-correlation measurement of two successive cycles (or non-cycles). */
    double getReliability() const { return 100.0 * cycleCorrelation; }


    /** \name Audio Processing */

    /** Estimates and returns the fundamental frequency of the incoming signal (in seconds). */
    INLINE double estimatePeriod(double inputSignal);

    /** Overloaded method for cases in which a whitened signal can be provided by the caller - 
    in this case, the whitening here will be bypassed. */
    INLINE double estimatePeriod(double inputSignal, double whiteSignal);

    /** Estimates and returns the fundamental frequency of the incoming signal (in Hz). */
    INLINE double estimateFundamentalFrequency(double inputSignal);

    /** Estimates and returns the pitch of the incoming signal (in terms of a midi note number). */
    INLINE double estimatePitch(double inputSignal);


    FormantRemover formantRemover; // move to the protceted area?

  protected:

    /** Estimates the fractional part of the zero crossing by fitting a cubic polynomial through
    the most recent four samples (including the current one) and solving for the zero-crossing of
    this cubic polynomial. */
    double getFractionalPartOfZeroCrossing();


    double sampleRate, sampleRateRec; // sample rate and its reciprocal. 

    double minFundamental, maxFundamental, minPeriod, maxPeriod;
      // Minimum and maximum values for the fundamental frequency and period.

    double y0, y1, y2, y3;
      // 3 previous and the current value of the lowpass-filtered input signal - need to be 
      // remembered between calls to estimateFundamental in oder to detect zero-crossings. */

    double fracOld;
      // The old (from the previous iteration) value of the fractional part of the zero 
      // crossing - together with the current value of this, a non-integer period estimate can be 
      // obtained.

    double periodEstimate, frequencyEstimate; // periodEstimateSmoothed;
      // This member always contains the current estimate of the pitch period the corresponding
      // fundamental frequency - it will be updated when upward zero-crossings are detected.

    double cycleCorrelation;
      // This keeps track of the cross-correlation between two successive cycles.

    rsOnePoleFilterDD  dcBlocker;
      // The dc-blocker filter for the signal in which we look for zero-crossings.

    FourPoleFilter lowpass;
      // The lowpass filter for the signal in which we look for zero-crossings.

    EnvelopeFollower envFollower;
      // detects the loudness of the signal

    rsCyclicAutoCorrelator correlator;
      // measures the correlation between successive pitch-cycles

    int sampleCounter;
      // counts the number of samples between two succesive upward zero crossings

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double PitchDetector::estimatePeriod(double inputSignal)
  {
    // remove (or at least, weaken) the formants from the input signal:
    double whiteSignal = formantRemover.getSample(inputSignal);

    return estimatePeriod(inputSignal, whiteSignal);
  }

  INLINE double PitchDetector::estimatePeriod(double inputSignal, double whiteSignal)
  { 
    double tmpCandidate, filteredSignal, env, frac;
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
    isAboveThreshold = env > RAPT::rsDbToAmp(-60.0);

    if( y2 < 0.0 && y1 >= 0.0 )
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
      tmpCandidate = RAPT::rsClip(tmpCandidate, minPeriod, maxPeriod);

      // remember the current fractional part for the next time:
      fracOld = frac;

      // get the cycle correlation as a measure of reliability of the period estimate (between 
      // 0...1) - if the estimate is reliable enough, accept it:
      cycleCorrelation = correlator.getCorrelationAndReset();
      if( isAboveThreshold && cycleCorrelation > 0.9 )
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

  INLINE double PitchDetector::estimateFundamentalFrequency(double inputSignal)
  {
    estimatePeriod(inputSignal);
    return frequencyEstimate;
  }

  INLINE double PitchDetector::estimatePitch(double inputSignal)
  { 
    return RAPT::rsFreqToPitch(estimateFundamentalFrequency(inputSignal));
      // this is inefficient (calls freqToPitch per sample, where it could have been called
      // only when the pitch has changed) 
  }

  INLINE double PitchDetector::getFractionalPartOfZeroCrossing()
  {
    // compute the zero crossing of a straight line y=a*x+b: x0 = -b/a, y2 is two samples ago, y1
    // is one sample ago:
    double fracLinear = y2 / (y2-y1);  

    // for a faster but less precise version, we could directly return fracLinear now
    // ->wrap this into a function getApproximateFractionalPartOfZeroCrossing

    // obtain coefficients for a cubic polynomial through the most recent 4 samples (where time is 
    // interpreted as having its origin two samples ago):
    double d = y2;  
    double c = y1 - (1.0/3.0)*y3 - 0.5*y2 - (1.0/6.0)*y0;
    double b = 0.5*(y3+y1) - y2;
    double a = 0.5*(y2-y1) + (1.0/6.0)*(y0-y3);
    
    // with these coefficients for a cubic polynomial, we now compute the zero-crossing (root) of 
    // that polynomial via Newton-Raphson iteration, using the zero-crossing from the fitted line
    // as initial estimate ...and return that value:
    //return getCubicRootNear(fracLinear, a, b, c, d, 0.0, 1.0);
    return RAPT::rsPolynomial<double>::getCubicRootNear(fracLinear, a, b, c, d, 0.0, 1.0);
  }

}

#endif
