// construction/destruction:

rsLinkwitzRileyCrossOver::rsLinkwitzRileyCrossOver(int newMaxButterworthOrder) 
: lowpass1(newMaxButterworthOrder/2)
, lowpass2(newMaxButterworthOrder/2)
, sumAllpass(newMaxButterworthOrder/2)
{
  rassert( newMaxButterworthOrder >= 1 ); // filter of zero or negative order? no such thing!

  maxButterworthOrder = newMaxButterworthOrder;
  sampleRate          = 44100.0;
  crossoverFrequency  = 1000.0;
  butterworthOrder    = RAPT::rsMin(2, maxButterworthOrder);
  updateFilterCoefficients();
}

rsLinkwitzRileyCrossOver::~rsLinkwitzRileyCrossOver()
{

}

// setup:

void rsLinkwitzRileyCrossOver::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    updateFilterCoefficients();
  }
}

void rsLinkwitzRileyCrossOver::setCrossoverFrequency(double newCrossoverFrequency)
{
  if( newCrossoverFrequency <= 20000.0 )
    crossoverFrequency = newCrossoverFrequency;
  updateFilterCoefficients();
}

void rsLinkwitzRileyCrossOver::setSlope(int newSlope)
{
  rassert( newSlope%12 == 0 && newSlope >= 12 ); // slope must be a multiple of 12 dB/oct
  setButterworthOrder(newSlope/12);
}

void rsLinkwitzRileyCrossOver::setButterworthOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder <= maxButterworthOrder )
    butterworthOrder = newOrder;
  updateFilterCoefficients();
}

// inquiry:

void rsLinkwitzRileyCrossOver::getLowpassMagnitudeResponse(double* frequencies, double* magnitudes, 
  int numBins, bool inDecibels, bool accumulate)
{
  if( accumulate == false )
  {
    if( inDecibels == true )
      RAPT::rsArrayTools::fillWithValue(magnitudes, numBins, 0.0);
    else
      RAPT::rsArrayTools::fillWithValue(magnitudes, numBins, 1.0);
  }
  lowpass1.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
  lowpass2.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
}

void rsLinkwitzRileyCrossOver::getLowpassFrequencyResponse(double* frequencies, Complex* H, 
  int numBins, bool accumulate)
{
  if( accumulate == false )  
    RAPT::rsArrayTools::fillWithValue(H, numBins, Complex(1.0));

  double *w = new double[numBins];
  RAPT::rsArrayTools::copy(frequencies, w, numBins);
  RAPT::rsArrayTools::scale(w, w, numBins, 2*PI/sampleRate);

  lowpass1.getFrequencyResponse(w, H, numBins, rsFilterAnalyzerD::MULTIPLICATIVE_ACCUMULATION);
  lowpass2.getFrequencyResponse(w, H, numBins, rsFilterAnalyzerD::MULTIPLICATIVE_ACCUMULATION);

  delete[] w;
}

void rsLinkwitzRileyCrossOver::getHighpassMagnitudeResponse(double* frequencies, 
  double* magnitudes, int numBins, bool inDecibels, bool accumulate)
{
  Complex* H = new Complex[numBins];
  getHighpassFrequencyResponse(frequencies, H, numBins, false);

  if( accumulate == true )
  {
    if( inDecibels == true )
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] += RAPT::rsAmpToDb(H[k].getRadius());
    }
    else 
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] *= H[k].getRadius();
    }
  }
  else
  {
    if( inDecibels == true )
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] = RAPT::rsAmpToDb(H[k].getRadius());
    }
    else 
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] = H[k].getRadius();
    }
  }

  delete[] H;
}

void rsLinkwitzRileyCrossOver::getHighpassFrequencyResponse(double* frequencies, Complex* H, 
  int numBins, bool accumulate)
{
  double *w = new double[numBins];
  RAPT::rsArrayTools::copy(frequencies, w, numBins);
  RAPT::rsArrayTools::scale(w, w, numBins, 2*PI/sampleRate);

  Complex *tmpLowpass = new Complex[numBins];
  getLowpassFrequencyResponse(frequencies, tmpLowpass, numBins, false);

  Complex *tmpAllpass = new Complex[numBins];
  sumAllpass.getFrequencyResponse(w, tmpAllpass, numBins);

  if( accumulate == false ) 
    RAPT::rsArrayTools::subtract(tmpAllpass, tmpLowpass, H, numBins);
  else
  {
    RAPT::rsArrayTools::subtract(tmpAllpass, tmpLowpass, tmpAllpass, numBins); // tmpAllpass is now the highpass-response
    RAPT::rsArrayTools::multiply(H, tmpAllpass, H, numBins);
  }

  delete[] tmpLowpass;
  delete[] tmpAllpass;
  delete[] w;
}

// others:

void rsLinkwitzRileyCrossOver::resetBuffers()
{
  lowpass1.reset();
  lowpass2.reset();
  sumAllpass.reset();
}

void rsLinkwitzRileyCrossOver::updateFilterCoefficients()
{
  // create and set up a filter-designer object:
  rsInfiniteImpulseResponseDesignerD designer;
  designer.setSampleRate(sampleRate);
  designer.setApproximationMethod(rsPrototypeDesignerD::BUTTERWORTH);
  designer.setPrototypeOrder(butterworthOrder);
  designer.setFrequency(crossoverFrequency);
  // \todo keep this object around as a member to avoid unnecessary re-calculations of the 
  // prototype poles

  // design the lowpasses:
  designer.setMode(rsInfiniteImpulseResponseDesignerD::LOWPASS);
  lowpass1.setOrder(butterworthOrder);
  designer.getBiquadCascadeCoefficients(lowpass1.getAddressB0(), lowpass1.getAddressB1(), 
    lowpass1.getAddressB2(), lowpass1.getAddressA1(), lowpass1.getAddressA2() );
  lowpass2.copySettingsFrom(&lowpass1);

  // obtain the allpass:
  sumAllpass.copySettingsFrom(&lowpass1);
  sumAllpass.turnIntoAllpass();
}
