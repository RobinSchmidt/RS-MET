// construction/destruction:

rsEngineersFilterOld::rsEngineersFilterOld() : rsBiquadCascadeStereo(25)
{
  numStages  = 1;
  sampleRate = 44100.0;
}

// parameter settings:

void rsEngineersFilterOld::setSampleRate(double newSampleRate)
{
  //rsBiquadCascadeStereo::setSampleRate(newSampleRate);
  sampleRate = newSampleRate;
  designer.setSampleRate(newSampleRate);
  updateCoefficients();
}

void rsEngineersFilterOld::setMode(int newMode)
{
  designer.setMode(newMode);
  rsBiquadCascadeStereo::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

void rsEngineersFilterOld::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients(true);
}

void rsEngineersFilterOld::setFrequency(double newFrequency)
{
  designer.setFrequency(newFrequency);
  updateCoefficients();
}

void rsEngineersFilterOld::setPrototypeOrder(int newOrder)
{
  designer.setPrototypeOrder(newOrder);
  rsBiquadCascadeStereo::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

void rsEngineersFilterOld::setBandwidth(double newBandwidth)
{
  designer.setBandwidth(newBandwidth);
  updateCoefficients();
}

void rsEngineersFilterOld::setGain(double newGain)
{
  designer.setGain(newGain);
  updateCoefficients();
}

void rsEngineersFilterOld::setRipple(double newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

void rsEngineersFilterOld::setStopbandRejection(double newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

// inquiry:

void rsEngineersFilterOld::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
  bool inDecibels, bool accumulate)
{
  double *w = new double[numBins];
  for(int k = 0; k < numBins; k++)
    w[k] = 2 * PI * frequencies[k] / sampleRate;
  rsBiquadCascadeStereo::getMagnitudeResponse(w, magnitudes, numBins, inDecibels, accumulate);
  for(int k = 0; k < numBins; k++)
  {
    if( w[k] > PI )
      magnitudes[k] = -200.0;
  }
  delete[] w;
}

// internal functions:

void rsEngineersFilterOld::updateCoefficients(bool resetState)
{
  rsBiquadCascade::initBiquadCoeffs();
  designer.getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);
  if(resetState)
    rsBiquadCascadeStereo::reset();
}