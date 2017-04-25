#include "rosic_ToneControl.h"
using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

ToneControl::ToneControl()
{
  setSampleRate(44100.0);  //sampleRate = 44100 Hz by default

  setLowCornerFreq(100.0);
  setLowGain(0.0);
  setHighCornerFreq(8000.0);
  setHighGain(0.0);

  resetBuffers();             //reset memorized samples to zero
}

ToneControl::~ToneControl()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void ToneControl::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  sampleRateRec = 1.0 / sampleRate;

  calcLowShelfCoeffs();
  calcHighShelfCoeffs();
}

void ToneControl::setLowCornerFreq(double newLowCornerFreq)
{
  if( (newLowCornerFreq>=20.0) && (newLowCornerFreq<=20000.0) )
    lowCornerFreq = newLowCornerFreq;
  else
    lowCornerFreq = 20000.0;

  calcLowShelfCoeffs();
}

void ToneControl::setLowGain(double newLowGain)
{
  lowGain = newLowGain;

  calcLowShelfCoeffs();
}

void ToneControl::setHighCornerFreq(double newHighCornerFreq)
{
  if( (newHighCornerFreq>=20.0) && (newHighCornerFreq<=20000.0) )
    highCornerFreq = newHighCornerFreq;
  else
    highCornerFreq = 20000.0;

  calcHighShelfCoeffs();
}

void ToneControl::setHighGain(double newHighGain)
{
  highGain = newHighGain;

  calcHighShelfCoeffs();
}


//----------------------------------------------------------------------------
//others:

void ToneControl::calcBiquadCoeffs()
{
  // variables for the coefficients of the 3 one-pole filters:
  //double b0Ls, b1Ls, a1Ls; // low-shelf coeffs
  //double b0Hs, b1Hs, a1Hs; // high-shelf coeffs

  b0 = b0Ls*b0Hs;
  b1 = b0Ls*b1Hs + b1Ls*b0Hs;
  b2 = b1Ls*b1Hs;
  a1 = -(a1Ls+a1Hs);
  a2 = -(a1Ls*a1Hs);

  // Remark:
  // the a-coefficients have been given a minus-sign in order to implement the 
  // filters difference-equation as
  // out = b0*in + b1*x1 + b2*x2 + b3*x3
  //             + a1*y1 + a2*y2 + a3*y3;
  // which turned ot to be more efficient than subtracting weighted past 
  // output samples
}

void ToneControl::calcLowShelfCoeffs()
{
  double G,                 // gain as factor
    G_B,               // gain at corner-frequency
    epsilon, beta,     // intermediate variables
    x, gainFactor;     // more intermediate variables

  double omegaAna,          // pre-warped analog corner freq.
    omegaDig;          // corner freq of the digital filter

  double pProto, zProto,    // prototype pole and zero
    pAna, zAna,        // denormalized analog pole and zero
    pDig, zDig;        // digital pole and zero

  // catch the special condition of zero dB-gain:
  if( abs(lowGain) < 0.00001 ) // equality ckeck with margin
  {
    b0Ls = 1.0;
    b1Ls = 0.0;
    a1Ls = 0.0;

    // trigger re-calculation of the biquad-coefficients:
    calcBiquadCoeffs();

    return; // nothing more to do here, in this special case
  }

  // assign/calculate some intermediate variables:
  G       = dB2amp(lowGain);
  G_B     = sqrt(G);
  epsilon = sqrt( (G*G-G_B*G_B)/(G_B*G_B-1.0) );
  beta    = 1/epsilon;

  // make low-shelving prototype:
  zProto = -(G*beta);
  pProto = -beta;

  // convert the cutoff-frequency into a normalized cutoff-frequency in
  // radians/sec (a frequency between 0 and pi):
  omegaDig = (2.0*PI*lowCornerFreq)/sampleRate;

  // prewarp the desired digital cutoff radian frequency to the desired analog
  // cutoff radian frequency:
  omegaAna = 2.0*sampleRate * tan(0.5*omegaDig);

  // transform the filter to the desired selectivity and cutoff-frequency via
  // the frequency transformations:
  pAna = pProto*omegaAna;
  zAna = zProto*omegaAna;

  // transform the analog filter to a digital filter via the bilinear
  // transform:
  x    = pAna / (2.0*sampleRate);
  pDig = (1.0+x)/(1.0-x);
  x    = zAna / (2.0*sampleRate);
  zDig = (1.0+x)/(1.0-x);

  // calculate the overall gain factor for the filter:
  gainFactor = (-1.0 - pDig) / (-1.0 - zDig);

  // finally, calculate the coefficients:
  b0Ls = gainFactor;
  b1Ls = -gainFactor*zDig;
  a1Ls = -pDig;

  // trigger re-calculation of the biquad-coefficients:
  calcBiquadCoeffs();
}

void ToneControl::calcHighShelfCoeffs()
{
  double G,                 // gain as factor
    G_B,               // gain at corner-frequency
    epsilon, beta,     // intermediate variables
    x, gainFactor;     // more intermediate variables

  double omegaAna,          // pre-warped analog corner freq.
    omegaDig;          // corner freq of the digital filter

  double pProto, zProto,    // prototype pole and zero
    pAna, zAna,        // denormalized analog pole and zero
    pDig, zDig;        // digital pole and zero

  // catch the special condition of zero dB-gain:
  if( abs(highGain) < 0.00001 ) // equality ckeck with margin
  {
    b0Hs = 1.0;
    b1Hs = 0.0;
    a1Hs = 0.0;

    // trigger re-calculation of the biquad-coefficients:
    calcBiquadCoeffs();

    return; // nothing more to do here, in this special case
  }

  // assign/calculate some intermediate variables:
  G       = dB2amp(highGain);
  G_B     = sqrt(G);
  epsilon = sqrt( (G*G-G_B*G_B)/(G_B*G_B-1.0) );
  beta    = 1/epsilon;

  // make low-shelving prototype:
  zProto = -(G*beta);
  pProto = -beta;

  // convert the cutoff-frequency into a normalized cutoff-frequency in
  // radians/sec (a frequency between 0 and pi):
  omegaDig = (2.0*PI*highCornerFreq)/sampleRate;

  // prewarp the desired digital cutoff radian frequency to the desired analog
  // cutoff radian frequency:
  omegaAna = 2.0*sampleRate * tan(0.5*omegaDig);

  // transform the filter to the desired selectivity and cutoff-frequency via
  // the frequency transformations:
  pAna = zProto*omegaAna;
  zAna = pProto*omegaAna;

  // transform the analog filter to a digital filter via the bilinear
  // transform:
  x    = pAna / (2.0*sampleRate);
  pDig = (1.0+x)/(1.0-x);
  x    = zAna / (2.0*sampleRate);
  zDig = (1.0+x)/(1.0-x);

  // calculate the overall gain factor for the filter:
  gainFactor = (1.0 - pDig) / (1.0 - zDig);

  // finally, calculate the coefficients:
  b0Hs = gainFactor;
  b1Hs = -gainFactor*zDig;
  a1Hs = -pDig;

  // trigger re-calculation of the biquad-coefficients:
  calcBiquadCoeffs();
}

void ToneControl::resetBuffers()
{
  x1  = 0.0;
  x2  = 0.0;
  y1  = 0.0;
  y2  = 0.0;
  out = 0.0;
}
