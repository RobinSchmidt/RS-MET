//#include "rosic_LinkwitzRileyCrossOver.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

LinkwitzRileyCrossOver::LinkwitzRileyCrossOver(int newMaxButterworthOrder) 
: lowpass1(newMaxButterworthOrder/2)
, lowpass2(newMaxButterworthOrder/2)
, sumAllpass(newMaxButterworthOrder/2)
{
  rassert( newMaxButterworthOrder >= 1 ); // filter of zero or negative order? no such thing!

  maxButterworthOrder = newMaxButterworthOrder;
  sampleRate          = 44100.0;
  crossoverFrequency  = 1000.0;
  butterworthOrder    = rmin(2, maxButterworthOrder);
  updateFilterCoefficients();
}

LinkwitzRileyCrossOver::~LinkwitzRileyCrossOver()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// setup:

void LinkwitzRileyCrossOver::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 && newSampleRate != sampleRate )
  {
    sampleRate = newSampleRate;
    updateFilterCoefficients();
  }
}

void LinkwitzRileyCrossOver::setCrossoverFrequency(double newCrossoverFrequency)
{
  if( newCrossoverFrequency <= 20000.0 )
    crossoverFrequency = newCrossoverFrequency;
  updateFilterCoefficients();
}

void LinkwitzRileyCrossOver::setSlope(int newSlope)
{
  rassert( newSlope%12 == 0 && newSlope >= 12 ); // slope must be a multiple of 12 dB/oct
  setButterworthOrder(newSlope/12);
}

void LinkwitzRileyCrossOver::setButterworthOrder(int newOrder)
{
  if( newOrder >= 1 && newOrder <= maxButterworthOrder )
    butterworthOrder = newOrder;
  updateFilterCoefficients();
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

void LinkwitzRileyCrossOver::getLowpassMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
                                                         bool inDecibels, bool accumulate)
{
  if( accumulate == false )
  {
    if( inDecibels == true )
      fillWithValue(magnitudes, numBins, 0.0);
    else
      fillWithValue(magnitudes, numBins, 1.0);
  }
  lowpass1.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
  lowpass2.getMagnitudeResponse(frequencies, sampleRate, magnitudes, numBins, true, true);
}

void LinkwitzRileyCrossOver::getLowpassFrequencyResponse(double *frequencies, Complex *H, int numBins, bool accumulate)
{
  if( accumulate == false )  
    fillWithValue(H, numBins, Complex(1.0));

  double *w = new double[numBins];
  copyBuffer(frequencies, w, numBins);
  scale(w, w, numBins, 2*PI/sampleRate);

  lowpass1.getFrequencyResponse(w, H, numBins, FilterAnalyzer::MULTIPLICATIVE_ACCUMULATION);
  lowpass2.getFrequencyResponse(w, H, numBins, FilterAnalyzer::MULTIPLICATIVE_ACCUMULATION);

  delete[] w;
}

void LinkwitzRileyCrossOver::getHighpassMagnitudeResponse(double *frequencies, double *magnitudes, int numBins, 
                                                          bool inDecibels, bool accumulate)
{
  Complex *H = new Complex[numBins];
  getHighpassFrequencyResponse(frequencies, H, numBins, false);

  if( accumulate == true )
  {
    if( inDecibels == true )
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] += amp2dB(H[k].getRadius());
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
        magnitudes[k] = amp2dB(H[k].getRadius());
    }
    else 
    {
      for(int k=0; k<numBins; k++)
        magnitudes[k] = H[k].getRadius();
    }
  }

  delete[] H;
}

void LinkwitzRileyCrossOver::getHighpassFrequencyResponse(double *frequencies, Complex *H, int numBins, bool accumulate)
{
  double *w = new double[numBins];
  copyBuffer(frequencies, w, numBins);
  scale(w, w, numBins, 2*PI/sampleRate);

  Complex *tmpLowpass = new Complex[numBins];
  getLowpassFrequencyResponse(frequencies, tmpLowpass, numBins, false);

  Complex *tmpAllpass = new Complex[numBins];
  sumAllpass.getFrequencyResponse(w, tmpAllpass, numBins);

  if( accumulate == false ) 
    subtract(tmpAllpass, tmpLowpass, H, numBins);
  else
  {
    subtract(tmpAllpass, tmpLowpass, tmpAllpass, numBins); // tmpAllpass is now the highpass-response
    multiply(H, tmpAllpass, H, numBins);
  }

  delete[] tmpLowpass;
  delete[] tmpAllpass;
  delete[] w;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// others:

void LinkwitzRileyCrossOver::resetBuffers()
{
  lowpass1.reset();
  lowpass2.reset();
  sumAllpass.reset();
}

void LinkwitzRileyCrossOver::updateFilterCoefficients()
{
  // create and set up a filter-designer object:
  InfiniteImpulseResponseDesigner designer;
  designer.setSampleRate(sampleRate);
  designer.setApproximationMethod(PrototypeDesigner::BUTTERWORTH);
  designer.setPrototypeOrder(butterworthOrder);
  designer.setFrequency(crossoverFrequency);
  // \todo keep this object around as a member to avoid unnecessary re-calculations of the prototype poles

  // design the lowpasses:
  designer.setMode(InfiniteImpulseResponseDesigner::LOWPASS);
  lowpass1.setOrder(butterworthOrder);
  designer.getBiquadCascadeCoefficients(lowpass1.getAddressB0(), lowpass1.getAddressB1(), lowpass1.getAddressB2(), 
                                                                 lowpass1.getAddressA1(), lowpass1.getAddressA2() );
  lowpass2.copySettingsFrom(&lowpass1);

  // obtain the allpass:
  sumAllpass.copySettingsFrom(&lowpass1);
  sumAllpass.turnIntoAllpass();
}
