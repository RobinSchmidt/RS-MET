// construction/destruction:

template<class T>
rsScopeScreenScanner<T>::rsScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; 
  //numCyclesShown = 2;
  sync = false;
  minZeroDistance = 10;
  numZerosToReset = 3;

  //sync = true;
  reset();
}

// setup:

template<class T>
void rsScopeScreenScanner<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;
}

template<class T>
void rsScopeScreenScanner<T>::setScanFreqNoSync(T newFrequency)
{
  scanFreq = newFrequency;
}

template<class T>
void rsScopeScreenScanner<T>::setSync(bool shouldSync)
{
  sync = shouldSync;
}

// misc:

template<class T>
void rsScopeScreenScanner<T>::reset()
{
  if(sync)
  {
    sawInc =  (T)samplesSinceReset / sampleRate; // correct?

    //sawInc = scanFreq / sampleRate;  
    // preliminary - todo: use the counted samples and zero-crossings to compute a value that
    // syncs
  }
  else
  {
    sawInc = scanFreq / sampleRate;
  }

  sawPhase = 0.0;
  samplesSinceReset = 0;
  samplesSinceLastZero = 0;
  zeroCrossingCount = 0;
}