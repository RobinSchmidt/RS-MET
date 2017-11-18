// construction/destruction:

template<class TSig, class TPar>
rsEngineersFilter<TSig, TPar>::rsEngineersFilter() : rsBiquadCascadeStereo(25)
{
  numStages  = 1;
  sampleRate = 44100.0;
}

// parameter settings:

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  //rsBiquadCascadeStereo::setSampleRate(newSampleRate);
  sampleRate = newSampleRate;
  designer.setSampleRate(newSampleRate);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setMode(int newMode)
{
  designer.setMode(newMode);
  rsBiquadCascade<TPar>::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients(true);
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setFrequency(TPar newFrequency)
{
  designer.setFrequency(newFrequency);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setPrototypeOrder(int newOrder)
{
  designer.setPrototypeOrder(newOrder);
  rsBiquadCascade<TPar>::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setBandwidth(TPar newBandwidth)
{
  designer.setBandwidth(newBandwidth);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setGain(TPar newGain)
{
  designer.setGain(newGain);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setRipple(TPar newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::setStopbandRejection(TPar newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

// inquiry:

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::getMagnitudeResponse(TPar* frequencies, TPar* magnitudes, 
  int numBins, bool inDecibels, bool accumulate)
{
  TPar* w = new TPar[numBins];
  for(int k = 0; k < numBins; k++)
    w[k] = 2 * PI * frequencies[k] / sampleRate;
  rsBiquadCascade<TPar>::getMagnitudeResponse(w, magnitudes, numBins, inDecibels, accumulate);
  for(int k = 0; k < numBins; k++)
  {
    if( w[k] > PI )
      magnitudes[k] = -200.0;
  }
  delete[] w;
}

// internal functions:

template<class TSig, class TPar>
void rsEngineersFilter<TSig, TPar>::updateCoefficients(bool resetState)
{
  rsBiquadCascade<TPar>::initBiquadCoeffs();
  designer.getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);
  if(resetState)
    rsBiquadCascade<TPar>::reset();
}