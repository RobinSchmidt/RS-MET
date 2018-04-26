#ifndef RAPT_CYCLICAUTOCORRELATOR_H
#define RAPT_CYCLICAUTOCORRELATOR_H

/** This class can be used to measure the autocorrelation of the incoming signal at some specified
(but possibly time-varying) time-lag. It is intended to be used to measure the correlation
between successive cycles of a (supposedly periodic) signal. It should be used by continuously
(i.e. per sample) calling acceptSample and calling getCorrelationAndReset() after each cycle.
When the successive cycles are very similar, this correlation value returned by
getCorrelationAndReset() will be close to unity. */

template<class T>
class rsCyclicAutoCorrelator
{

public:

  /** \name Construction/Destruction */

  /** Constructor - constructs a delay-line with a given maximum number of samples delay. */
  rsCyclicAutoCorrelator(int newMaxBufferSize = 16384);

  /** Destructor - frees memory for delayline. */
  ~rsCyclicAutoCorrelator();


  /** \name Audio Processing */

  /** Acepts one input at a time. */
  RS_INLINE void acceptSample(T in);

  /** Calculates the autocorrelation between the past two pitch-cycles. */
  RS_INLINE T getCorrelationAndReset();


  /** \name Misc */

  /** Resets the accumulators and the sample-counter to zero. */
  void resetAccumulators();

  /** Fills the buffer with all zeros. */
  void resetBuffer();

protected:

  /** \name Data */

  T sumOfProducts;
    // Accumulates the product between current input samples and input samples of past 
    // pitch-period which are stored in the buffer.

  T sumOfSquares;       // Accumulates the squares of the input samples.
  T oldMeanSquare;      // mean-square value of the old pitch-period
  T *buffer;            // stores the previous pitch-period

  int oldCycleLength;   // length of the old cycle (in order to not accumulate too many samples)
  int maxBufferSize;
  int sampleCounter;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class T>
RS_INLINE void rsCyclicAutoCorrelator<T>::acceptSample(T in)
{
  if(sampleCounter < maxBufferSize)
  {
    // accumulate the product of the incoming sample with the corresponding sample from the 
    // previous pitch-cycle (but only up to the length of the old cycle):
    if(sampleCounter <= oldCycleLength)
      sumOfProducts += in * buffer[sampleCounter];

    // overwrite the sample in the buffer with the new one (for the measurement in the next 
    // cycle - within this cycle, we are done with this particular buffer-cell):
    buffer[sampleCounter] = in;

    // we also keep track of the sum of the squares of the cycle for normalization:
    sumOfSquares += in*in;

    // increment the sample-counter:
    sampleCounter++;
  }
}

template<class T>
RS_INLINE T rsCyclicAutoCorrelator<T>::getCorrelationAndReset()
{
  // determine the number of samples over which we have to average:
  T numSamplesToAverage = (T)rsMin(sampleCounter, oldCycleLength);
  if(numSamplesToAverage <= 0.0)
  {
    numSamplesToAverage = 1.0;
    RS_DEBUG_BREAK;
  }

  // get the averaged product and the mean of the squares
  T meanProduct = sumOfProducts / numSamplesToAverage;
  T meanSquare  = sumOfSquares  / numSamplesToAverage;
  T normalizer  = rsSqrt(oldMeanSquare * meanSquare);

  // get the cross-correlation value with safeguard from divide by zero:
  T result;
  if(normalizer == 0.0)
    result = 0.0;
  else
    result = meanProduct / normalizer;

  // apply an additional weighting to take into account possible different lengths of cycles - 
  // the more different the lengths, the less the weighting:
  /*
  double minLength = (double) min(sampleCounter, oldCycleLength);
  double maxLength = (double) max(sampleCounter, oldCycleLength);
  if( maxLength <= 0.0 )
  {
    result = 0.0;
    DEBUG_BREAK;
  }
  else
  {
    double weight = (minLength / maxLength);
    result       *= weight;
    //double corrector = (maxLength-minLength)/maxLength;
    //result -= corrector;
  }
  */

  // update the old mean-square value and oldCycleLength-counter:
  oldMeanSquare  = meanSquare;
  oldCycleLength = sampleCounter;

  // reset the accumulators:
  sumOfProducts = 0.0;
  sumOfSquares  = 0.0;
  sampleCounter = 0;

  return result;
}

#endif
