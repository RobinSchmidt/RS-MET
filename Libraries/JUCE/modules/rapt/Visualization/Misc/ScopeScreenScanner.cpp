// construction/destruction:

template<class T>
rsScopeScreenScanner<T>::rsScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; 
  //numCyclesShown = 2;
  sync = false;
  minZeroDistance = 20; // maybe have a maxFrequency parameter instead
  numZerosToReset = 2;
  //sync = true;

  lowpass.setMode(LadderFilter<T,T>::LP_6); // maybe we can use a 1-pole filter
  lowpass.setCutoff(20.0);

  reset();
}

// setup:

template<class T>
void rsScopeScreenScanner<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;
  lowpass.setSampleRate(sampleRate);
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
  sawInc = scanFreq / sampleRate;
  if(sync)
    sawInc = rsMin(T(0.5), rsMax(sawInc, 1 / (T)samplesSinceReset));

  //if(sync)
  //  sawInc = 1 / (T)samplesSinceReset;
  //else
  //  sawInc = scanFreq / sampleRate;

  xOld = 0.0;
  sawPhase = 0.0;
  samplesSinceReset = 0;
  samplesSinceLastZero = 0;
  zeroCrossingCount = 0;
}