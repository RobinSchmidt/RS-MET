#include "rosic_OnePoleFilterStereo.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

OnePoleFilterStereo::OnePoleFilterStereo()
{
  shelvingGain = 1.0;

  setSampleRate(44100.0);  // sampleRate = 44100 Hz by default
  setMode      (0);        // bypass by default
  setCutoff    (20000.0);  // cutoff = 20000 Hz by default

  //calcCoeffs() is called by setMode and setCutoff

  resetBuffers();          //reset memorized samples to zero
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void OnePoleFilterStereo::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  sampleRateRec = 1.0 / sampleRate;

  calcCoeffs();
  return;
}

void OnePoleFilterStereo::setMode(int newMode)
{
  mode = newMode; // 0:bypass, 1:Low Pass, 2:High Pass
  calcCoeffs();
}

void OnePoleFilterStereo::setCutoff(double newCutoff)
{
  if( (newCutoff>0.0) && (newCutoff<=20000.0) )
    cutoff = newCutoff;
  else
    cutoff = 20000.0;

  calcCoeffs();
  return;
}

void OnePoleFilterStereo::setShelvingGain(double newGain)
{
  if( newGain > 0.0 )
  {
    shelvingGain = newGain;
    calcCoeffs();
  }
  else
    DEBUG_BREAK; // this is a linear gain factor and must be >= 0.0
}

void OnePoleFilterStereo::setShelvingGainInDb(double newGainDb)
{
  setShelvingGain(dB2amp(newGainDb));
}

void OnePoleFilterStereo::setCoeffs(double newB0, double newB1, double newA1)
{
  b0 = newB0;
  b1 = newB1;
  a1 = newA1;
}

void OnePoleFilterStereo::setState(double newX1L, double newX1R, double newY1L, double newY1R)
{
  x1L = newX1L;
  x1R = newX1R;
  y1L = newY1L;
  y1R = newY1R;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

Complex OnePoleFilterStereo::getTransferFunctionAt(rosic::Complex z)
{
  return (b0 + b1*(1/z)) / (1 - a1*(1/z));
}

double OnePoleFilterStereo::getMagnitudeAt(double frequency)
{
  double omega = 2*PI*frequency*sampleRateRec;

  Complex z = expC(Complex(0.0, omega));
  Complex H = getTransferFunctionAt(z);

  return H.getRadius();
}

//-------------------------------------------------------------------------------------------------
//others:

void OnePoleFilterStereo::calcCoeffs()
{
  switch(mode)
  {
  case LOWPASS:
    {
      // intermediate variable for calculation (x is the amount of decay
      // between adjacent samples):
      double x = exp( -2.0 * PI * cutoff * sampleRateRec);
      b0 = 1-x;
      b1 = 0.0;
      a1 = x;
    }
    break;
  case HIGHPASS:  
    {
      double x = exp( -2.0 * PI * cutoff * sampleRateRec);
      b0 =  0.5*(1+x);
      b1 = -0.5*(1+x);
      a1 = x;
    }
    break;
  case LOWSHELV:
    {
      double c = 0.5*(shelvingGain-1.0);
      double t = tan(PI*cutoff*sampleRateRec);
      double a;
      if( shelvingGain >= 1.0 )
        a = (t-1.0)/(t+1.0);
      else
        a = (t-shelvingGain)/(t+shelvingGain);

      b0 = 1.0 + c + c*a;
      b1 = c + c*a + a;
      a1 = -a;
    }
    break;
  case HIGHSHELV:
    {
      double c = 0.5*(shelvingGain-1.0);
      double t = tan(PI*cutoff*sampleRateRec);
      double a;
      if( shelvingGain >= 1.0 )
        a = (t-1.0)/(t+1.0);
      else
        a = (shelvingGain*t-1.0)/(shelvingGain*t+1.0);

      b0 = 1.0 + c - c*a;
      b1 = a + c*a - c;
      a1 = -a;
    }
    break;

  case ALLPASS:  
    {
      // calculate allpass coefficients (see DAFX):
      double t = tan(PI*cutoff*sampleRateRec);
      double x = (t-1.0) / (t+1.0);

      b0 = x;
      b1 = 1.0;
      a1 = -x;
    }
    break;


  default://bypass:
    {
      b0 = 1.0;
      b1 = 0.0;
      a1 = 0.0;
    }break;
  }//end of switch

}

void OnePoleFilterStereo::resetBuffers()
{
  x1L = 0.0;
  y1L = 0.0;
  x1R = 0.0;
  y1R = 0.0;
}
