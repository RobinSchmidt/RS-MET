#ifndef RAPT_AUTOCORRELATIONPITCHDETECTOR_H
#define RAPT_AUTOCORRELATIONPITCHDETECTOR_H

/** This class estimates the fundamental of an incoming signal (which is assumed to be pitched and
monophonic) by finding the peak-value of the autocorrelation function. 

...this class supports block-processing only...hmm...maybe extend it to sample-based processing (it
would just be a matter of buffering). ...i think, it's currently used only in some offline 
processing algorithms */

template<class T>
class rsAutoCorrelationPitchDetector
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsAutoCorrelationPitchDetector();

  /** Destructor. */
  ~rsAutoCorrelationPitchDetector();


  /** \name Setup */

  /** Sets the sample-rate. */
  void setSampleRate(T newSampleRate);

  /** Sets the maximum expected fundamental frequency in Hz. */
  void setMaxFundamental(T newMaxFundamental);

  /** Sets the size of the internal data buffer (which equals the FFT size and determines the
  lowest fundamental that can be detected). */
  void setBufferSize(int newSize);

  /** Sets the interval at which the frequency estimate is updated. */
  void setUpdateInterval(int newInterval);


  /** \name Inquiry */

  /** Returns the current reliability of the period estimate as a value normlized to 0...1. The
  reliability is defined as the autocorrelation-value at the detected peak-lag. */
  T getReliability() const { return reliability; }


  /** \name Audio Processing */

  /** Estimates and returns the fundamental frequency of the incoming signal. */
  T processBlock(T *block, int blockSize);

  /** Used internally to update the current estimate for the fundamental frequency (and its
  reliability). */
  void updateFrequencyEstimate();

  //static double estimateFundamental(double *x, int N, double fs, double fMin, double fMax,
  //  double *reliability);

  static T estimateFundamental(T *x, int L, T fs, T fMin, T fMax, T *reliability = nullptr);


  /** \name Misc */

  /** Resets the internal state. */
  void reset();

protected:

  /** \name Internal Functions */

  /** Memory re-allocation (invoked when buffer-size changes). */
  void reAllocateBuffers();

  /** Applies a nonlinear function to the autocorrleation values before the parabola is fitted.
  The function was found by trial-and-error to improve the fundamental frequency estimate.
  \todo do more research on the "optimal" function */
  static T applyNonlinearity(T x);


  /** \name Data */

  T *circularBuffer;
  T *linearBuffer;

  T sampleRate;
  T maxFundamental;
  T frequencyEstimate;
  T reliability;

  int writeIndex;
  int sampleCounter;
  int updateInterval;
  int bufferSize;

};

#endif
