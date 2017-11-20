using namespace RSLib;

// construction/destruction:

rsEngineersFilter::rsEngineersFilter() : rsBiquadCascadeStereo(25)
{
  numStages  = 1;
  sampleRate = 44100.0;
}

// setup:

void rsEngineersFilter::setSampleRate(double newSampleRate)
{
  //BiquadCascade::setSampleRate(newSampleRate);
  sampleRate = newSampleRate;
  designer.setSampleRate(newSampleRate);
  updateCoefficients();
}

void rsEngineersFilter::setMode(int newMode)
{
  designer.setMode(newMode);
  rsBiquadCascadeStereo::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

void rsEngineersFilter::setApproximationMethod(int newApproximationMethod)
{
  designer.setApproximationMethod(newApproximationMethod);
  updateCoefficients();
}

void rsEngineersFilter::setFrequency(double newFrequency)
{
  designer.setFrequency(newFrequency);
  updateCoefficients();
}

void rsEngineersFilter::setPrototypeOrder(int newOrder)
{
  designer.setPrototypeOrder(newOrder);
  rsBiquadCascadeStereo::setNumStages(designer.getNumBiquadStages());
  updateCoefficients();
}

void rsEngineersFilter::setBandwidth(double newBandwidth)
{
  designer.setBandwidth(newBandwidth);
  updateCoefficients();
}

void rsEngineersFilter::setGain(double newGain)
{
  designer.setGain(newGain);
  updateCoefficients();
}

void rsEngineersFilter::setRipple(double newRipple)
{
  designer.setRipple(newRipple);
  updateCoefficients();
}

void rsEngineersFilter::setStopbandRejection(double newStopbandRejection)
{
  designer.setStopbandRejection(newStopbandRejection);
  updateCoefficients();
}

// inquiry:

void rsEngineersFilter::getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
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

void rsEngineersFilter::updateCoefficients()
{
  rsBiquadCascadeStereo::initBiquadCoeffs();
  designer.getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);
  rsBiquadCascadeStereo::reset();
}
