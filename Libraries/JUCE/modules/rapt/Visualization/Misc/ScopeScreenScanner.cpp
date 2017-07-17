// construction/destruction:

template<class T>
rsScopeScreenScanner<T>::rsScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; 
  //numCyclesShown = 2;
  sync = false;
  minZeroDistance = 20; // maybe have a maxFrequency parameter instead
  numZerosToReset = 5;
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
    sawInc = 1 / (T)samplesSinceReset;
  else
    sawInc = scanFreq / sampleRate;

  sawPhase = 0.0;
  samplesSinceReset = 0;
  samplesSinceLastZero = 0;
  zeroCrossingCount = 0;
}