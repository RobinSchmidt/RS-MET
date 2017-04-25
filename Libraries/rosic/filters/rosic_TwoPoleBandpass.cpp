#include "rosic_TwoPoleBandpass.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

TwoPoleBandpass::TwoPoleBandpass()
{
  setSampleRate   (44100.0);  // sampleRate = 44100 Hz by default

  setLpfCutoff    (20000.0);  // lpfCutoff = 20000 Hz by default
  setHpfCutoff    (0.0);      // lpfCutoff = 0 Hz by default

  resetBuffers();             // reset memorized samples to zero
}

TwoPoleBandpass::~TwoPoleBandpass()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void TwoPoleBandpass::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
    sampleRate = newSampleRate;

  sampleRateRec = 1.0 / sampleRate;
  calcCoeffs();
  return;
}

void TwoPoleBandpass::setLpfCutoff(double newLpfCutoff)
{
  if( (newLpfCutoff>=0.0) && (newLpfCutoff<=20000.0) )
    lpfCutoff = newLpfCutoff;
  else
    lpfCutoff = 20000.0;
  calcCoeffs();
  return;
}

void TwoPoleBandpass::setHpfCutoff(double newHpfCutoff)
{
  if( (newHpfCutoff>=0.0) && (newHpfCutoff<=20000.0) )
    hpfCutoff = newHpfCutoff;
  else
    hpfCutoff = 20000.0;
  calcCoeffs();
  return;
}


//-------------------------------------------------------------------------------------------------
// inquiry:

double TwoPoleBandpass::getLpfCutoff()
{
  return lpfCutoff;
}

double TwoPoleBandpass::getHpfCutoff()
{
  return hpfCutoff;
}

//-------------------------------------------------------------------------------------------------
// others:

void TwoPoleBandpass::calcCoeffs()
{
  //calculate lowpass coefficients:
  //intermediate variable for calculation (x is the amount of decay
  //between adjacent samples):
  double x = exp( -2.0 * PI * lpfCutoff * sampleRateRec);

  b0Lpf = 1-x;
  b1Lpf = 0.0;
  a1Lpf = x;

  //calculate highpass coefficients:
  x = exp( -2.0 * PI * hpfCutoff * sampleRateRec);

  b0Hpf =  0.5*(1+x);
  b1Hpf = -0.5*(1+x);
  a1Hpf = x;
}

void TwoPoleBandpass::resetBuffers()
{
  x_1Lpf = 0.0;
  y_1Lpf = 0.0;
  x_1Hpf = 0.0;
  y_1Hpf = 0.0;
}
