// construction/destruction:

template<class TSig, class TPar>
rsEngineersFilter<TSig, TPar>::rsEngineersFilter() : rsBiquadCascade<TSig, TPar>(25)
{
  this->numStages  = 1;
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
  rsBiquadCascade<TSig, TPar>::setNumStages(designer.getNumBiquadStages());
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
  rsBiquadCascade<TSig, TPar>::setNumStages(designer.getNumBiquadStages());
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
    w[k] = TPar(2*PI) * frequencies[k] / sampleRate;
  rsBiquadCascade<TSig, TPar>::getMagnitudeResponse(w, magnitudes, numBins, inDecibels,
    accumulate);
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
  rsBiquadCascade<TSig, TPar>::initBiquadCoeffs();
  designer.getBiquadCascadeCoefficients(this->b0, this->b1, this->b2, this->a1, this->a2);
  if(resetState)
    rsBiquadCascade<TSig, TPar>::reset();
}


/*=================================================================================================

Projects that do similar filter designs:

https://github.com/vinniefalco/DSPFilters
https://github.com/MikeCurrington/mkfilter

*/