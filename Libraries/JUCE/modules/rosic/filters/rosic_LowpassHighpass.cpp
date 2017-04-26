#include "rosic_LowpassHighpass.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

LowpassHighpass::LowpassHighpass()
{
  sampleRateRec = 1.0 / 44100.0;
  lpfCutoff     = 20000.0;
  hpfCutoff     = 20.0;
  calcCoeffs();                  
}

LowpassHighpass::~LowpassHighpass()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void LowpassHighpass::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRateRec = 1.0 / newSampleRate;
  calcCoeffs();
}

void LowpassHighpass::setLowpassCutoff(double newCutoff)
{
  if( (newCutoff>=0.0) && (newCutoff<=20000.0) )
    lpfCutoff = newCutoff;
  else
    lpfCutoff = 20000.0;

  calcCoeffs();
  return;
}

void LowpassHighpass::setHighpassCutoff(double newCutoff)
{
  if( (newCutoff>=0.0) && (newCutoff<=20000.0) )
    hpfCutoff = newCutoff;
  else
    hpfCutoff = 20000.0;

  calcCoeffs();
  return;
}

//-------------------------------------------------------------------------------------------------
// others:

void LowpassHighpass::calcCoeffs()
{
  // variables for the coefficients of the 3 one-pole filters:
  double b0Lpf, b1Lpf, a1Lpf; // lowpass-coeffs
  double b0Hpf, b1Hpf, a1Hpf; // highpass-coeffs
  double x;                   // amount of decay between adjacent samples for lowpass filters

  // calculate lowpass coefficients (see dspguide for details):
  if( lpfCutoff == 0.0 )
  {
    b0Lpf = 1.0;
    b1Lpf = 0.0;
    a1Lpf = 0.0;
  }
  else
  {
    x = exp( -2.0 * PI * lpfCutoff * sampleRateRec); 

    b0Lpf = 1-x;
    b1Lpf = 0.0;
    a1Lpf = -x;
  }

  // calculate highpass coefficients (see dspguide for details):
  if( hpfCutoff == 0.0 )
  {
    b0Hpf = 1.0;
    b1Hpf = 0.0;
    a1Hpf = 0.0;
  }
  else
  {
    x = exp( -2.0 * PI * hpfCutoff * sampleRateRec);

    b0Hpf =  0.5*(1+x);
    b1Hpf = -0.5*(1+x);
    a1Hpf = -x;
  }

  // combine the 2 one-pole filters into one biquad filter by multiplying out
  // the product of the 2 transfer-functions, the resulting 2-pole coefficients
  // come out as follows:
  b0 = b0Lpf*b0Hpf;
  b1 = b0Lpf*b1Hpf + b1Lpf*b0Hpf;
  b2 = b1Lpf*b1Hpf;

  a1 = -( a1Lpf+a1Hpf );
  a2 = -( a1Lpf*a1Hpf );

  // Remark:
  // the a-coefficients have been given a minus-sign in order to implement the 
  // filters difference-equation as
  // out = b0*in + b1*x1 + b2*x2 + b3*x3
  //             + a1*y1 + a2*y2 + a3*y3;
  // which turned ot to be more efficient than subtracting weighted past 
  // output samples
}
