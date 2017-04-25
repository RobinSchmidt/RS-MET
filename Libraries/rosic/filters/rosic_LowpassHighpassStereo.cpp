#include "rosic_LowpassHighpassStereo.h"
using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

LowpassHighpassStereo::LowpassHighpassStereo()
{
	setSampleRate(44100.0);   // sampleRate = 44100 Hz by default

	setLpfCutoff (20000.0);   // lpfCutoff = 20000 Hz by default
	setHpfCutoff (0.0);       // lpfCutoff = 0 Hz by default

	resetBuffers();           // reset memorized samples to zero
}

LowpassHighpassStereo::~LowpassHighpassStereo()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void LowpassHighpassStereo::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;
 sampleRateRec = 1.0 / sampleRate;

 calcCoeffs();
 return;
}

void LowpassHighpassStereo::setLpfCutoff(double newLpfCutoff)
{
 if( (newLpfCutoff>=0.0) && (newLpfCutoff<=20000.0) )
  lpfCutoff = newLpfCutoff;
	else
		lpfCutoff = 20000.0;

 calcCoeffs();
 return;
}

void LowpassHighpassStereo::setHpfCutoff(double newHpfCutoff)
{
 if( (newHpfCutoff>=0.0) && (newHpfCutoff<=20000.0) )
  hpfCutoff = newHpfCutoff;
	else
		hpfCutoff = 20000.0;

 calcCoeffs();
 return;
}

//----------------------------------------------------------------------------
// others:

void LowpassHighpassStereo::calcCoeffs()
{
 // variables for the coefficients of the 3 one-pole filters:
 double b0Lpf, b1Lpf, a1Lpf; // lowpass-coeffs
 double b0Hpf, b1Hpf, a1Hpf; // highpass-coeffs

 double x; // intermediate variable for calculation (x is the amount of
           // decay between adjacent samples for LPFs):

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

void LowpassHighpassStereo::resetBuffers()
{
 x1L  = 0.0;
 x2L  = 0.0;
 y1L  = 0.0;
 y2L  = 0.0;
 x1R  = 0.0;
 x2R  = 0.0;
 y1R  = 0.0;
 y2R  = 0.0;
}
