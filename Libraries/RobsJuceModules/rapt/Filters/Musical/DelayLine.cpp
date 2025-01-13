// Construction/Destruction:

template<class T>
rsDelayLineBasic<T>::rsDelayLineBasic()
{
  maxDelay  = 3;
  delayLine = new T[maxDelay+1];
  tapIn     = 0;
  tapOut    = 0;
  reset();

  // Maybe avoid this allocation of a buffer with 3 samples. It will most certainly be wasted 
  // because soon the client will request an allocation of a more sensible size. But then we must
  // deal with the possibility of having a nullptr. That's the tradeoff here: either always do 
  // wasted allocations or accept that pointer could initially be null. Or maybe a trick could be 
  // used: just let pointer initially point to some static member variable. We might apply a sort
  // of "null object" pattern to the pointer variable. But I'm not sure, if that's workable. And 
  // why 3. If anything, we should use 1. Or maybe even 0. But verify if the bit-maksing will work
  // or if this is a weird edge case.
}

template<class T>
rsDelayLineBasic<T>::~rsDelayLineBasic()
{
  if( delayLine != nullptr )
    delete[] delayLine;
}

// Setup:

template<class T>
void rsDelayLineBasic<T>::setMaxDelayInSamples(int newMaxDelay)
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
void rsDelayLineBasic<T>::setDelayInSamples(int delay)
{
  rsAssert(delay >= 0 && delay <= maxDelay, 
           "Delay out of range in rsBasicDelayLine::setDelayInSamples");

  // Sanitize input argument and adjust tapOut pointer:
  delay = rsClip(delay, 0, maxDelay);
  tapOut = tapIn - delay;
  if( tapOut < 0 )
    tapOut += maxDelay+1;
}

// Misc:

template<class T>
void rsDelayLineBasic<T>::reset()
{
  for(int i = 0; i < maxDelay+1; i++)
    delayLine[i] = 0.0;
}


//-------------------------------------------------------------------------------------------------

// construction/destruction:

template<class TSig, class TPar>
rsDelayLineTempoSynced<TSig, TPar>::rsDelayLineTempoSynced(int maximumDelayInSamples)
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
rsDelayLineTempoSynced<TSig, TPar>::~rsDelayLineTempoSynced()
{
  if( delayBuffer != nullptr )
  {
    delete[] delayBuffer;
    delayBuffer = nullptr;
  }
}

// parameter settings (set-functions):

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if(newSampleRate > 0.01)
  {
    sampleRate = newSampleRate;
    setupDelayInSamples();
  }
}

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::setDelayTime(TPar newDelayTime)
{
  //delayTime = rsClip(newDelayTime, 0.0, 4.25);
  delayTime = newDelayTime;
  setupDelayInSamples();
}

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::setSyncMode(bool shouldTempoSync)
{
  tempoSync = shouldTempoSync;
  setupDelayInSamples();
}

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::setTempoInBPM(TPar newTempoInBPM)
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
void rsDelayLineTempoSynced<TSig, TPar>::setInterpolationMethod(int newMethod)
{
  if(newMethod <= rsInterpolator<TSig>::WARPED_ALLPASS)
    interpolator.setInterpolationMethod(newMethod);
  else
    rsError("Unknown interpolation method");
}

// others:

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::clearDelayBuffer()
{
  for(int i=0; i<length+interpolatorMargin; i++)
    delayBuffer[i] = 0.0;
  interpolator.reset();
}

template<class TSig, class TPar>
void rsDelayLineTempoSynced<TSig, TPar>::setupDelayInSamples()
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


//=================================================================================================
/*


ToDo:

- Refactor the code in such a way that delayline classes never know about sample-rates and delay
  times in seconds. They should just know about the delay in samples. This kind of information 
  should really be held at a higher level.

- Maybe bring back the delayline implementation that is set up in terms of a time (in seconds) 
  and a sample rate and call it rsDelayLineTimeBased. It's in Misc\UnusedCode\Misc.h

- Use only one pointer for tapIn and tapOut (see Julius Smith's pasp-book)

*/