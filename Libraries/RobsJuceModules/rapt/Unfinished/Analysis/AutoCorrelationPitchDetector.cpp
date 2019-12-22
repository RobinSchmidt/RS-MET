// Construction/Destruction:

template<class T>
rsAutoCorrelationPitchDetector<T>::rsAutoCorrelationPitchDetector()
{
  circularBuffer    = NULL;
  linearBuffer      = NULL;
  sampleRate        = 44100.0;
  maxFundamental    = 8000.0;
  updateInterval    = 256;
  bufferSize        = 2048;
  reAllocateBuffers();
}

template<class T>
rsAutoCorrelationPitchDetector<T>::~rsAutoCorrelationPitchDetector()
{
  delete[] circularBuffer;
  delete[] linearBuffer;
}

// Setup:

template<class T>
void rsAutoCorrelationPitchDetector<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;
}

template<class T>
void rsAutoCorrelationPitchDetector<T>::setMaxFundamental(T newMaxFundamental)
{
  maxFundamental = newMaxFundamental;
}

template<class T>
void rsAutoCorrelationPitchDetector<T>::setBufferSize(int newSize)
{
  if( newSize != bufferSize  && rsIsPowerOfTwo(newSize) )
  {
    bufferSize = newSize;
    reAllocateBuffers();
  }
}

template<class T>
void rsAutoCorrelationPitchDetector<T>::setUpdateInterval(int newInterval)
{
  updateInterval = newInterval;
}

// Audio Processing:

template<class T>
T rsAutoCorrelationPitchDetector<T>::processBlock(T *block, int blockSize)
{
  int readIndex, copyLength;
  if( blockSize > bufferSize )
  {
    // fill up internal buffer with last portion of incoming block:
    readIndex = blockSize-bufferSize;
    memcpy(circularBuffer, &block[readIndex], bufferSize*sizeof(T));
      // maybe use rsCopyBuffer (rewrite such that it uses memcopy)
    writeIndex = 0;
  }
  else
  {
    if( blockSize < bufferSize-writeIndex )
    {
      // copy incoming block to current position circular buffer:
      memcpy(&circularBuffer[writeIndex], block, blockSize*sizeof(T));
      writeIndex += blockSize;
    }
    else
    {
      // copy first portion of incoming block to current position circular buffer:
      copyLength = bufferSize-writeIndex;
      memcpy(&circularBuffer[writeIndex], block, copyLength*sizeof(T));

      // copy last portion of incoming block to begin of circular buffer (wrap-around):
      readIndex  = copyLength;
      copyLength = blockSize-readIndex;
      memcpy(circularBuffer, &block[readIndex], copyLength*sizeof(T));
      writeIndex = copyLength;

      // write function copyLinearToCircular
    }
  }

  sampleCounter += blockSize;
  if( sampleCounter >= updateInterval )
  {
    // copy circular into linear buffer:
    //T linearBuffer[16];
    memcpy(linearBuffer, &circularBuffer[writeIndex], (bufferSize-writeIndex)*sizeof(T));
    memcpy(&linearBuffer[bufferSize-writeIndex], circularBuffer, writeIndex*sizeof(T));
      // write a function copyCircularToLinear

    // do processing and reset counter:
    updateFrequencyEstimate();
    sampleCounter = 0;  // shouldn't it be sampleCounter %= updateInterval
  }

  return frequencyEstimate;
}

template<class T>
void rsAutoCorrelationPitchDetector<T>::updateFrequencyEstimate()
{
  T minFundamental = 0.0; // preliminary - maybe, we don't need it at all
  frequencyEstimate = estimateFundamental(linearBuffer, bufferSize, sampleRate, minFundamental, 
    maxFundamental, &reliability);
}

template<class T>
T rsAutoCorrelationPitchDetector<T>::estimateFundamental(T *x, int L, T fs, T fMin, T fMax, 
  T *reliability)
{
  // Get autocorrelation sequence r:
  T *r = new T[L];
  rsAutoCorrelationFFT(x, L, r);
  //rsPlotArray(r, L);

  // find index in autocorrelation function where a peak with maximum value occurs:
  int startIndex  = rsMax(1, (int) (fs/fMax));
  int m           = startIndex;    // index of maximum
  //T maxValue = -rsInfDouble;
  T maxValue = -RS_INF(T);
  T value;
  for(int k = startIndex; k < L-1; k++)
  {
    value = r[k];
    if( value > r[k-1] && value >= r[k+1] && value >= maxValue )
    {
      m = k;
      maxValue  = value;
      maxValue *= 1.0+4.5/(3.0+k*k);
        // The function f(k) = 1.0+4.5/(3.0+k*k) models the ratio between the value at lag k when
        // f = fs/k (in which case the peak occurs exactly at k) and when f = fs/(k+0.5) (in which
        // case there are two equal "peak" values at lag k and k+1 which are both lower and the
        // next peak occurs exactly at 2*k+1 and might be higher than the 2 values at k and k+1).
        // If the peak at 2*k+1 is higher than the "distributed peak" at k and k+1, we would have
        // an octave error - the estimated frequency would be 1 octave too low. We assume this
        // worst case and scale the maxValue accordingly in order to avoid such errors.
    }
  }

  // after the loop, we don't need the scaling by 1.0+4.5/(3.0+k*k) anymore:
  maxValue = r[m];

  // obtain unbiased estimates of the autocorrelation function for the 3 values to be considered for
  // the parabolic fit:
  T y[3];
  y[0] = r[m-1] / (L-(m-1));
  y[1] = r[m]   / (L- m   );
  y[2] = r[m+1] / (L-(m+1));

  // it has been found empirically that a transformation of the values y <- f(y) for some nonlinear
  // monotonic function f leads to less bias towards frequencies which are dividers of the
  // sample-rate (i think, the estimated exact location of the maximum of the parabola leans less
  // to the inside):
  y[0] = applyNonlinearity(y[0]);
  y[1] = applyNonlinearity(y[1]);
  y[2] = applyNonlinearity(y[2]);

  // find exact location of maximum by fitting a parabola through 3 successive autocorrelation
  // values and using the maximum of the parabola:
  T a[3];
  rsPolynomial<T>::fitQuadratic_0_1_2(a, y);   // linker error on linux - fixed
  T offset = 0.0;
  if( a[2] != 0.0 )
    offset = -0.5*a[1]/a[2];

  // Actually, the results from this parabola-fitting are not very accurate. I tried fitting a 4th
  // order polynomial to r[m-2], r[m-1], r[m], r[m+1], r[m+2], obtain the derivative and find the 
  // derivative's root using the maximum of the parabola as initial guess, but the improvement was 
  // only miniscule. Perahps a totally different approach is necessary - try to find spectral peaks
  // of a windowed magnitude spectrum near frequencies at integer multiples of the frequency that 
  // we found here
 


  // todo: compute a correction to eliminate the (frequency-dependent) bias - this will probably
  // have to be a function of the offset:
  T correction = 0.0; // preliminary

  // compute the frequency estimate:
  T exactIndex = m - 1 + offset + correction;

  T frequencyEstimate = fs / exactIndex;

  // estimate reliability of the estimate by using the autocorrelation value at the peak,
  // normalized as if, for zero-lag, the autocorrelation-function would be unity:
  if( reliability != nullptr )
    *reliability = maxValue / r[0];

  delete[] r;

  return frequencyEstimate;
}

template<class T>
T rsAutoCorrelationPitchDetector<T>::applyNonlinearity(T x)
{
  //return x;                   // linear
  //return log(x);              // worse than linear
  //return 1/x;                 // worse than log
  return rsPowBipolar(x, 0.75); // better than linear

  // under construction - maybe try to find a better function - perhaps some rational function or 
  // something which would also be more efficient to evaluate
}

/*
template<class T>
void rsAutoCorrelationPitchDetector<T>::updateFrequencyEstimate()
{
  // obtain autocorrelation function (ACF) of the linearBuffer:
  T *tmpBuffer = new T[bufferSize];  // later: have a member or use linearBuffer itself
  autoCorrelationFFT(linearBuffer, bufferSize, tmpBuffer);

  // Apply a 2-point moving average (MA) filter to the ACF (without the scaling by 1/2). This
  // avoids octave errors when a sinusoid of f = fs/(k+1/2) (k integer) is applied. Later we will
  // need to take accout of this by moving the exact location of the maximum by -0.5:
  T x;
  T xOld = tmpBuffer[0]; // initial value for x[n-1]
  for(int n = 0; n < bufferSize; n++)
  {
    x = tmpBuffer[n];
    tmpBuffer[n] += xOld;
    xOld = x;
  }

  // find index in (MA-filtered) autocorrelation function where a peak with maximum value occurs:
  int    maxIndex = findHighestPeakIndex(tmpBuffer, bufferSize);
  T maxValue = tmpBuffer[maxIndex];

  // obtain unbiased estimates of the autocorrelation function for the 3 values to be considered for
  // the parabolic fit:
  T y[3];
  y[0] = tmpBuffer[maxIndex-1] / (bufferSize - (maxIndex-1));
  y[1] = tmpBuffer[maxIndex]   / (bufferSize -  maxIndex   );
  y[2] = tmpBuffer[maxIndex+1] / (bufferSize - (maxIndex+1));

  // it has been found empirically that a transformation of the values y <- y^(2/3) leads to less
  // bias towards frequencies which are dividers of the sample-rate (i think, the estimated exact
  // location of the maximum of the parabola leans less to the inside):
  T ex = 1.0;  // bypass nonlinearity for test purposes
  //ex = 2.0/3.0;
  ex = 0.691; // seems to be even better than 2/3 ...maybe optimize further
  //ex = 0.816; // good for 3520 Hz @ 44100
  y[0] = powBipolar(y[0], ex);
  y[1] = powBipolar(y[1], ex);
  y[2] = powBipolar(y[2], ex);
    // maybe the optimal exponent could be a function of y[0], y[1], y[2]?

  // find exact location of maximum by fitting a parabola through 3 successive autocorrelation
  // values and using the maximum of the parabola:
  T a[3];
  fitQuadratic(a, y);
  T offset = 0.0;
  if( a[2] != 0.0 )
    offset = -0.5*a[1]/a[2];

  T exactIndex = maxIndex - 1 + offset;
  exactIndex -= 0.5; // due to the 2-point MA filter on the ACF
  frequencyEstimate = sampleRate / exactIndex;

  // estimate reliability of the estimate by using the autocorrelation value at the peak,
  // normalized as if, for zero-lag, the autocorrelation-function would be unity:
  reliability = maxValue / tmpBuffer[0];

  delete[] tmpBuffer;
}
*/

// Misc:

template<class T>
void rsAutoCorrelationPitchDetector<T>::reset()
{
  writeIndex        = 0;
  sampleCounter     = 0;
  frequencyEstimate = 1000.0;
  reliability       = 0.0;
  rsArrayTools::fillWithZeros(circularBuffer, bufferSize);
  rsArrayTools::fillWithZeros(linearBuffer, bufferSize);
}

// Internal Functions;

template<class T>
void rsAutoCorrelationPitchDetector<T>::reAllocateBuffers()
{
  delete[] circularBuffer;
  delete[] linearBuffer;
  circularBuffer = new T[bufferSize];
  linearBuffer   = new T[bufferSize];
  reset();
}
