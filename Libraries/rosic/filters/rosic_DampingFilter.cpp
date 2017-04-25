#include "rosic_DampingFilter.h"
using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

DampingFilter::DampingFilter()
{
  sampleRate              = 44100.0;
  sampleRateRec           = 1.0 / sampleRate;
  globalGainFactor        = 1.0;
  lowCrossoverFreq        = 250.0;
  lowCrossoverGainFactor  = 0.5;
  lowGainFactor           = 0.25;
  highCrossoverFreq       = 4000.0;
  highCrossoverGainFactor = 0.5;
  highGainFactor          = 0.25;

  /*
  setSampleRate(44100.0);       // sampleRate = 44100 Hz by default
  setGlobalGainFactor(1.0);
  setLowCrossoverFreq(100.0);
  setLowCrossoverGainFactor(0.5);
  setLowGainFactor(0.25);
  setHighCrossoverFreq(8000.0);
  setHighCrossoverGainFactor(0.75);
  setHighGainFactor(0.5);
  */

  calculateCoefficients();
  reset();               // reset memorized samples to zero
}

DampingFilter::~DampingFilter()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void DampingFilter::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
  sampleRateRec = 1.0 / sampleRate;
  calculateCoefficients();
}

void DampingFilter::setGlobalGainFactor(double newGlobalGainFactor)
{
  if( _isnan(newGlobalGainFactor) )
    DEBUG_BREAK;

  globalGainFactor = newGlobalGainFactor;
  calculateCoefficients();
}

void DampingFilter::setLowCrossoverFreq(double newLowCrossoverFreq)
{
  if( (newLowCrossoverFreq>=20.0) && (newLowCrossoverFreq<=20000.0) )
    lowCrossoverFreq = newLowCrossoverFreq;
  else
    lowCrossoverFreq = 20000.0;
  calculateCoefficients();
}

void DampingFilter::setLowCrossoverGainFactor(double newLowCrossoverGainFactor)
{
  if( newLowCrossoverGainFactor >= 0.0 )
    lowCrossoverGainFactor = newLowCrossoverGainFactor;
  calculateCoefficients();
}

void DampingFilter::setLowGainFactor(double newLowGainFactor)
{
  if( newLowGainFactor >= 0.0 )
    lowGainFactor = newLowGainFactor;
  calculateCoefficients();
}

void DampingFilter::setHighCrossoverFreq(double newHighCrossoverFreq)
{
  if( (newHighCrossoverFreq>=20.0) && (newHighCrossoverFreq<=20000.0) )
    highCrossoverFreq = newHighCrossoverFreq;
  else
    highCrossoverFreq = 20000.0;
  calculateCoefficients();
}

void DampingFilter::setHighGainFactor(double newHighGainFactor)
{
  if( newHighGainFactor >= 0.0 )
    highGainFactor = newHighGainFactor;
  calculateCoefficients();
}

void DampingFilter::setHighCrossoverGainFactor(double newHighCrossoverGainFactor)
{
  if( newHighCrossoverGainFactor >= 0.0 )
    highCrossoverGainFactor = newHighCrossoverGainFactor;
  calculateCoefficients();
}

//-------------------------------------------------------------------------------------------------
// others:

bool DampingFilter::areGainsAllowed(double G, double G_B)
{
  double margin = 0.000001; // margin for equality check

  if( fabs(G-1.0) < margin || fabs(G_B-1.0) < margin || fabs(G-G_B) < margin )
    return false;

  if( (G_B > 1.0) && (G > G_B) )
    return true;
  if( (G_B < 1.0) && (G < G_B) )
    return true;

  return false;
}

void DampingFilter::calculateCoefficients()
{
  double G,                 // gain as factor
    G_B,                    // gain at corner-frequency
    epsilon, beta,          // intermediate variables
    x, gainFactor;          // more intermediate variables

  double omegaAna,          // pre-warped analog radian corner freq.
    omegaDig;               // corner freq of the digital filter

  double pProto, zProto,    // prototype pole and zero
    pAna, zAna,             // denormalized analog pole and zero
    pDig, zDig;             // digital pole and zero

  double b0Ls, b1Ls, a1Ls;  // 1st order low-shelf coefficients
  double b0Hs, b1Hs, a1Hs;  // 1st order high-shelf coefficients

  //---------------------------------------------------------------------------
  // low-shelf coefficient calculation:

  G   = lowGainFactor;  
  G_B = lowCrossoverGainFactor;

  // catch some special (not allowed) conditions:
  if( !areGainsAllowed(G, G_B) )  
  {
    b0Ls = 1.0;
    b1Ls = 0.0;
    a1Ls = 0.0;
  }
  else
  {
    // assign/calculate some intermediate variables:
    epsilon = sqrt( (G*G-G_B*G_B)/(G_B*G_B-1.0) );
    beta    = 1/epsilon;

    // make low-shelving prototype:
    zProto = -(G*beta);
    pProto = -beta;

    // convert the cutoff-frequency into a normalized cutoff-frequency in
    // radians/sec (a frequency between 0 and pi):
    omegaDig = (2.0*PI*lowCrossoverFreq)/sampleRate;

    // prewarp the desired digital cutoff radian frequency to the desired analog
    // cutoff radian frequency:
    omegaAna = 2.0*sampleRate * tan(0.5*omegaDig);

    // transform the filter to the desired selectivity and cutoff-frequency via
    // an analog frequency transformation (LowShelf->LowShelf):
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
  }

  //---------------------------------------------------------------------------
  // high-shelf coefficient calculation:

  // we must take the reciprocal and scale by G, because the prototpye's 
  // frequency response will be inverted when transforming from low-shelf to 
  // high-shelf    
  G   = highGainFactor;
  G_B = G/highCrossoverGainFactor; 

  // catch some special (not allowed) conditions:
  if( !areGainsAllowed(G, G_B) )  
  {
    b0Hs = 1.0;
    b1Hs = 0.0;
    a1Hs = 0.0;
  }
  else
  {
    // assign/calculate some intermediate variables - 
    epsilon = sqrt( (G*G-G_B*G_B)/(G_B*G_B-1.0) );
    beta    = 1/epsilon;

    // make low-shelving prototype:
    zProto = -(G*beta);
    pProto = -beta;

    // convert the cutoff-frequency into a normalized cutoff-frequency in radians/sec (a frequency 
    // between 0 and pi):
    omegaDig = (2.0*PI*highCrossoverFreq)/sampleRate;

    // prewarp the desired digital cutoff radian frequency to the desired analog cutoff radian 
    // frequency:
    omegaAna = 2.0*sampleRate * tan(0.5*omegaDig);

    // transform the filter to the desired selectivity and cutoff-frequency via an analog frequency 
    // transformation (LowShelf->High-Shelf):
    pAna = zProto*omegaAna;
    zAna = pProto*omegaAna;

    // transform the analog filter to a digital filter via the bilinear transform:
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
  }

  //---------------------------------------------------------------------------
  // combination of the two shelving filters into a single biquad stage:

  b0 = globalGainFactor * b0Ls*b0Hs;
  b1 = globalGainFactor * (b0Ls*b1Hs + b1Ls*b0Hs);
  b2 = globalGainFactor * b1Ls*b1Hs;
  a1 = -(a1Ls+a1Hs);
  a2 = -(a1Ls*a1Hs);

  // Remark:
  // the a-coefficients have been given a minus-sign in order to implement the 
  // filters difference-equation as
  // out = b0*in + b1*x1 + b2*x2 + a1*y1 + a2*y2;
  // which turned ot to be more efficient than subtracting weighted past 
  // output samples

    
  if( _isnan(b0) || _isnan(b1) || _isnan(b2) || _isnan(a1) || _isnan(a2) )  
    DEBUG_BREAK;

  /*
  // coefficients of -0.0 seem to trigger a calculation of NAN - avoid this by 
  // converting to +0.0 - nah that seems not to solve it...:
  if( b0 == -0.0 ) b0 = 0.0;
  if( b1 == -0.0 ) b1 = 0.0;
  if( b2 == -0.0 ) b2 = 0.0;
  if( a1 == -0.0 ) a1 = 0.0;
  if( a2 == -0.0 ) a2 = 0.0;
  */

  //reset();
}

void DampingFilter::reset()
{
  x1  = 0.0;
  x2  = 0.0;
  y1  = 0.0;
  y2  = 0.0;
}
