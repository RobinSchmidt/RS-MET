// construction/destruction:

template<class T>
rsScopeScreenScanner<T>::rsScopeScreenScanner()
{
  sampleRate = 44100;
  scanFreq = 5; 
  zoom = 1;
  sync = false;
  minZeroDistance = 20; // maybe have a maxFrequency parameter instead
  numZerosToReset = 2;  // # zeros to be seen before reset occurs (typically number of cycles seen)
  //lowpass.setMode(LadderFilter<T,T>::Mode::LP_6); // maybe we can use a 1-pole filter - optimize
  lowpass.setMode(rsLadderFilter<T,T>::Mode::BP_6_6);
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
  sawInc = zoom * scanFreq / sampleRate;
}

template<class T>
void rsScopeScreenScanner<T>::setSync(bool shouldSync)
{
  sync = shouldSync;
}

template<class T>
void rsScopeScreenScanner<T>::setNumCyclesShown(int numCycles)
{
  numZerosToReset = numCycles;
}

//template<class T>
//void rsScopeScreenScanner<T>::setZoom(double newZoom)
//{
//  zoom = (T)newZoom;
//  sawInc = zoom * scanFreq / sampleRate;
//}

// misc:

template<class T>
void rsScopeScreenScanner<T>::reset()
{
  if(sync)
    sawInc = rsMin(T(0.5), zoom / (T)samplesSinceReset);
  else
    sawInc = zoom * scanFreq / sampleRate;
  xOld = 0.0;
  sawPhase = 0.0;
  samplesSinceReset = 0;
  samplesSinceLastZero = 0;
  zeroCrossingCount = 0;
}