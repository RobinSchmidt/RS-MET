//#include "rosic_LpfHpfApf.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

LpfHpfApf::LpfHpfApf()
{
  setSampleRate(44100.0); 
  setLowpassCutoff(20000.0);  
  setHighpassCutoff(0.0);      
  setAllpassFrequency(20000.0);  
  resetBuffers();            
}

LpfHpfApf::~LpfHpfApf()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void LpfHpfApf::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  sampleRateRec = 1.0 / sampleRate;
  calcCoeffs();
}

void LpfHpfApf::setLowpassCutoff(double newCutoff)
{
  if( (newCutoff>=0.0) && (newCutoff<=20000.0) )
    lpfCutoff = newCutoff;
  else
    lpfCutoff = 20000.0;
  calcCoeffs();
}

void LpfHpfApf::setHighpassCutoff(double newCutoff)
{
  if( (newCutoff>=0.0) && (newCutoff<=20000.0) )
    hpfCutoff = newCutoff;	
  else
    hpfCutoff = 20000.0;
  calcCoeffs();
}

void LpfHpfApf::setAllpassFrequency(double newFrequency)
{
  if( (newFrequency>=0.0) && (newFrequency<=20000.0) )
    apfFreq = newFrequency;
  else
    apfFreq = 20000.0;
  calcCoeffs();
}

//-------------------------------------------------------------------------------------------------
// others:

void LpfHpfApf::calcCoeffs()
{
  // variables for the coefficients of the 3 one-pole filters:
  double b0Lpf, b1Lpf, a1Lpf; // lowpass-coeffs
  double b0Hpf, b1Hpf, a1Hpf; // highpass-coeffs
  double b0Apf, b1Apf, a1Apf; // allpass-coefss

  double x; // intermediate variable for calculation (x is the amount of
  // decay between adjacent samples for LPFs):

  // calculate lowpass coefficients (see dspguide):
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

  // calculate highpass coefficients (see dspguide):
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

  // calculate allpass coefficients (see DAFX):
  x = (tan(PI*apfFreq*sampleRateRec)-1) / (tan(PI*apfFreq*sampleRateRec)+1);

  b0Apf = x;
  b1Apf = 1.0;
  a1Apf = x;

  // combine the 3 one-pole filters into one 3-pole filter by multiplying out
  // the product of the 3 transfer-functions, the resulting 3-pole coefficients
  // come out as follows:
  b0 = b0Lpf*b0Hpf*b0Apf;
  b1 = b0Lpf*b0Hpf*b1Apf + (b0Lpf*b1Hpf+b1Lpf*b0Hpf)*b0Apf;
  b2 = b1Lpf*b1Hpf*b0Apf + (b0Lpf*b1Hpf+b1Lpf*b0Hpf)*b1Apf;
  b3 = b1Lpf*b1Hpf*b1Apf;
  a1 = -( a1Lpf+a1Hpf+a1Apf );
  a2 = -( a1Lpf*a1Hpf + (a1Lpf+a1Hpf)*a1Apf );
  a3 = -( a1Lpf*a1Hpf*a1Apf );

  // Remark:
  // the a-coefficients have been given a minus-sign in order to implement the 
  // filters difference-equation as
  // out = b0*in + b1*x1 + b2*x2 + b3*x3
  //             + a1*y1 + a2*y2 + a3*y3;
  // which turned out to be more efficient than subtracting weighted past 
  // output samples
}

void LpfHpfApf::resetBuffers()
{
  x1 = 0.0;
  x2 = 0.0;
  x3 = 0.0;
  y1 = 0.0;
  y2 = 0.0;
  y3 = 0.0;
}
