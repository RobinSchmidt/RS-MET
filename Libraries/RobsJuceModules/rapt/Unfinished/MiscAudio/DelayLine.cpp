// Construction/Destruction:

template<class T>
rsBasicDelayLine<T>::rsBasicDelayLine()
{
  maxDelay  = 3;
  delayLine = new T[maxDelay+1];
  tapIn     = 0;
  tapOut    = 0;
  reset();
}

template<class T>
rsBasicDelayLine<T>::~rsBasicDelayLine()
{
  if( delayLine != NULL )
    delete[] delayLine;
}

// Setup:

template<class T>
void rsBasicDelayLine<T>::setMaximumDelayInSamples(int newMaxDelay)
{
  if( newMaxDelay > maxDelay )
  {
    delete[] delayLine;
    maxDelay  = rsNextPowerOfTwo(newMaxDelay + 1) - 1;
    delayLine = new T[maxDelay+1];
    rsArrayTools::fillWithZeros(delayLine, maxDelay+1);
  }
}

template<class T>
void rsBasicDelayLine<T>::setDelayInSamples(int newDelay)
{
  int delay = rsMax(0, newDelay);
  if( delay > maxDelay )
    setMaximumDelayInSamples(delay);

  // adjust tapOut-pointer:
  tapOut = tapIn - delay;
  if( tapOut < 0 )
    tapOut += maxDelay+1;
}

// Misc:

template<class T>
void rsBasicDelayLine<T>::reset()
{
  for(int i = 0; i < maxDelay+1; i++)
    delayLine[i] = 0.0;
}

//-------------------------------------------------------------------------------------------------

// Construction/Destruction:

template<class TSig, class TPar>
rsDelayLine<TSig, TPar>::rsDelayLine()
{
  delayInSeconds = TPar(0.001);
  sampleRate     = TPar(44100.0);
  setDelayInSeconds(delayInSeconds);
}

template<class TSig, class TPar>
rsDelayLine<TSig, TPar>::~rsDelayLine()
{

}

// Setup:

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  setDelayInSeconds(delayInSeconds);
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInSamples(int newDelayInSamples)
{
  rsBasicDelayLine<TSig>::setDelayInSamples(newDelayInSamples);
  delayInSeconds = (TPar) newDelayInSamples / sampleRate;
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInSeconds(TPar newDelayInSeconds)
{
  delayInSeconds = rsMax(0.0, newDelayInSeconds);
  setDelayInSamples( rsRoundToInt(sampleRate*delayInSeconds) );
}

template<class TSig, class TPar>
void rsDelayLine<TSig, TPar>::setDelayInMilliseconds(TPar newDelayInMilliseconds)
{
  setDelayInSeconds(0.001*newDelayInMilliseconds);
}

//-------------------------------------------------------------------------------------------------

// construction/destruction:

template<class TSig, class TPar>
rsFractionalDelayLine<TSig, TPar>::rsFractionalDelayLine(int maximumDelayInSamples)
{
  length      = maximumDelayInSamples + 1;
  delayBuffer = new TSig[length+interpolatorMargin];

  tapIn      = 0;
  tapOut     = 0;
  delayTime  = TPar(0.25);
  sampleRate = TPar(44100.0);
  bpm        = TPar(120.0);
  tempoSync  = false;

  interpolator.setInterpolationMethod(rsInterpolator<TSig>::WARPED_ALLPASS);

  interpolator.setInterpolationMethod(rsInterpolator<TSig>::LINEAR);  // for debug

  setDelayTime(delayTime);
  clearDelayBuffer();
}

template<class TSig, class TPar>
rsFractionalDelayLine<TSig, TPar>::~rsFractionalDelayLine()
{
  if( delayBuffer != NULL )
  {
    delete[] delayBuffer;
    delayBuffer = NULL;
  }
}

// parameter settings (set-functions):

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
  }
}

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setDelayTime(TPar newDelayTime)
{
  //delayTime = rsClip(newDelayTime, 0.0, 4.25);
  delayTime = newDelayTime;
  setupDelayInSamples();
}

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setSyncMode(bool shouldTempoSync)
{
  tempoSync = shouldTempoSync;
  setupDelayInSamples();
}

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setTempoInBPM(TPar newTempoInBPM)
{
  if(newTempoInBPM >= 0.0)
  {
    bpm = newTempoInBPM;
    setupDelayInSamples();
  }
  else
    rsError("Tempo < 0");
}

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setInterpolationMethod(int newMethod)
{
  if(newMethod <= rsInterpolator<TSig>::WARPED_ALLPASS)
    interpolator.setInterpolationMethod(newMethod);
  else
    rsError("Unknown interpolation method");
}

// others:

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::clearDelayBuffer()
{
  for(int i=0; i<length+interpolatorMargin; i++)
    delayBuffer[i] = 0.0;
  interpolator.reset();
}

template<class TSig, class TPar>
void rsFractionalDelayLine<TSig, TPar>::setupDelayInSamples()
{
  double delayInSeconds;
  if( tempoSync )
    delayInSeconds = rsBeatsToSeconds(delayTime, bpm);
  else
    delayInSeconds = delayTime;
    // get rid of the temposync stuff

  delayInSamples = sampleRate*delayInSeconds;
  delayInSamples = rsClip(delayInSamples, (TPar)(interpolatorMargin-1),
    (TPar) (length-1-interpolatorMargin));
  //delayInSamples = clip(delayInSamples, (double)(interpolatorMargin-1),
  //  (double) (length-1-interpolatorMargin)); // old


  // update member delayTime to the clipped value:
  delayInSeconds = delayInSamples / sampleRate;
  if( tempoSync )
    delayTime = rsSecondsToBeats(delayInSeconds, bpm);
  else
    delayTime = delayInSeconds;

  // calculate the integer and fractional parts of the delay:
  TPar tmp   = floor(delayInSamples);
  int  dInt  = (int) tmp;
  TPar dFrac = delayInSamples - tmp;
  frac         = 1.0 - dFrac; // because we look backwards

                              // adjust tapOut-pointer:
  tapOut = tapIn - dInt - 1;
  if( frac >= 1.0 )
  {
    frac    = 0.0;
    tapOut += 1;
  }
  tapOut = wrapAround(tapOut);
}
