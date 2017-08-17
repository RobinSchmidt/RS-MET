// construction/destruction:

template<class T>
rsScopeScreenScanner<T>::rsScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; ;
  sync = false;
  minZeroDistance = 20; // maybe have a maxFrequency parameter instead
  numZerosToReset = 2;  // # zeros to be seen before reset occurs (typically number of cycles seen)
  //lowpass.setMode(LadderFilter<T,T>::LP_6); // maybe we can use a 1-pole filter - optimize
  lowpass.setMode(LadderFilter<T,T>::BP_6_6);
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
  sawInc = scanFreq / sampleRate;
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
    sawInc = rsMin(T(0.5), 1 / (T)samplesSinceReset);
  //if(sync)
  //  sawInc = rsMin(T(0.5), rsMax(sawInc, 1 / (T)samplesSinceReset));
  xOld = 0.0;
  sawPhase = 0.0;
  samplesSinceReset = 0;
  samplesSinceLastZero = 0;
  zeroCrossingCount = 0;
}