template<class TSig, class TPar>
PhasorFilter<TSig, TPar>::PhasorFilter()
{
  sampleRate = 44100.0;
  frequency  = 1000.0;
  decay      = 0.01;
  updateCoefficients();
  reset();
}

// setup:

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  updateCoefficients();
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setFrequency(TPar newFrequency)
{
  frequency = newFrequency;
  updateCoefficients();
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::setDecayTime(TPar newDecay)
{
  decay = newDecay;
  updateCoefficients();
}

// audio processing:

template<class TSig, class TPar>
inline void PhasorFilter<TSig, TPar>::processFrame(TSig *x, TSig *y)
{

}

// misc:

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::reset()
{
  xOld = yOld = 0;
}

template<class TSig, class TPar>
void PhasorFilter<TSig, TPar>::updateCoefficients()
{

}
