// construction/destruction:

template<class TSig, class TPar>
rsLinkwitzRileyCrossOver<TSig, TPar>::rsLinkwitzRileyCrossOver(int newMaxButterworthOrder) 
: lowpass1(newMaxButterworthOrder/2)
, lowpass2(newMaxButterworthOrder/2)
, sumAllpass(newMaxButterworthOrder/2)
{
  rsAssert( newMaxButterworthOrder >= 1 ); // filter of zero or negative order? no such thing!

  maxButterworthOrder = newMaxButterworthOrder;
  sampleRate          = 44100.0;
  crossoverFrequency  = 1000.0;
  butterworthOrder    = rsMin(2, maxButterworthOrder);
  updateFilterCoefficients();
}

// setup:

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    updateFilterCoefficients();
  }
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::setCrossoverFrequency(TPar newCrossoverFrequency)
{
  if( newCrossoverFrequency <= 20000.0 )
    crossoverFrequency = newCrossoverFrequency;
  updateFilterCoefficients();
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::setSlope(int newSlope)
{
  rsAssert( newSlope%12 == 0 && newSlope >= 12 ); // slope must be a multiple of 12 dB/oct
  setButterworthOrder(newSlope/12);
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::setButterworthOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder <= maxButterworthOrder )
    butterworthOrder = newOrder;
  updateFilterCoefficients();
}

// inquiry:

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::getLowpassMagnitudeResponse(TPar* frequencies, 
  TPar* magnitudes, int numBins, bool inDecibels, bool accumulate) const
{
  if( accumulate == false )
  {
    if( inDecibels == true )
      rsArrayTools::fillWithValue(magnitudes, numBins, TPar(0));
    else
      rsArrayTools::fillWithValue(magnitudes, numBins, TPar(1));
  }
  lowpass1.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
  lowpass2.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::getLowpassFrequencyResponse(TPar* frequencies, 
  Complex* H, int numBins, bool accumulate) const
{
  if( accumulate == false )  
    rsArrayTools::fillWithValue(H, numBins, Complex(1.0));

  TPar* w = new TPar[numBins];
  rsArrayTools::copy(frequencies, w, numBins);
  rsArrayTools::scale(w, w, numBins, TPar(2*PI)/sampleRate);

  lowpass1.getFrequencyResponse(w, H, numBins, rsFilterAnalyzer<TPar>::MULTIPLICATIVE_ACCUMULATION);
  lowpass2.getFrequencyResponse(w, H, numBins, rsFilterAnalyzer<TPar>::MULTIPLICATIVE_ACCUMULATION);

  delete[] w;
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::getHighpassMagnitudeResponse(TPar* frequencies, 
  TPar* magnitudes, int numBins, bool inDecibels, bool accumulate) const
{
  Complex* H = new Complex[numBins];
  getHighpassFrequencyResponse(frequencies, H, numBins, false);
  if( accumulate == true ) {
    if( inDecibels == true ) {
      for(int k=0; k<numBins; k++)
        magnitudes[k] += rsAmpToDb(abs(H[k])); }
    else {
      for(int k=0; k<numBins; k++)
        magnitudes[k] *= abs(H[k]); }}
  else {
    if( inDecibels == true ) {
      for(int k=0; k<numBins; k++)
        magnitudes[k] = rsAmpToDb(abs(H[k])); }
    else {
      for(int k=0; k<numBins; k++)
        magnitudes[k] = abs(H[k]); }}
  delete[] H;
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::getHighpassFrequencyResponse(TPar* frequencies, 
  Complex* H, int numBins, bool accumulate) const
{
  TPar* w = new TPar[numBins];
  rsArrayTools::copy(frequencies, w, numBins);
  rsArrayTools::scale(w, w, numBins, TPar(2*PI)/sampleRate);

  Complex *tmpLowpass = new Complex[numBins];
  getLowpassFrequencyResponse(frequencies, tmpLowpass, numBins, false);

  Complex *tmpAllpass = new Complex[numBins];
  sumAllpass.getFrequencyResponse(w, tmpAllpass, numBins);

  if( accumulate == false ) 
    rsArrayTools::subtract(tmpAllpass, tmpLowpass, H, numBins);
  else
  {
    rsArrayTools::subtract(tmpAllpass, tmpLowpass, tmpAllpass, numBins); // tmpAllpass is now the highpass-response
    rsArrayTools::multiply(H, tmpAllpass, H, numBins);
  }

  delete[] tmpLowpass;
  delete[] tmpAllpass;
  delete[] w;
}

// others:

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::resetBuffers()
{
  lowpass1.reset();
  lowpass2.reset();
  sumAllpass.reset();
}

template<class TSig, class TPar>
void rsLinkwitzRileyCrossOver<TSig, TPar>::updateFilterCoefficients()
{
  // create and set up a filter-designer object:
  rsInfiniteImpulseResponseDesigner<TPar> designer;
  designer.setSampleRate(sampleRate);
  designer.setApproximationMethod(rsPrototypeDesigner<TPar>::BUTTERWORTH);
  designer.setPrototypeOrder(butterworthOrder);
  designer.setFrequency(crossoverFrequency);
  // \todo keep this object around as a member to avoid unnecessary re-calculations of the 
  // prototype poles

  // design the lowpasses:
  designer.setMode(rsInfiniteImpulseResponseDesigner<TPar>::LOWPASS);
  lowpass1.setOrder(butterworthOrder);
  designer.getBiquadCascadeCoefficients(lowpass1.getAddressB0(), lowpass1.getAddressB1(), 
    lowpass1.getAddressB2(), lowpass1.getAddressA1(), lowpass1.getAddressA2() );
  lowpass2.copySettingsFrom(&lowpass1);

  // obtain the allpass:
  sumAllpass.copySettingsFrom(&lowpass1);
  sumAllpass.turnIntoAllpass();
}
