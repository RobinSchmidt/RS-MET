#include "DampingFilter.h"

//----------------------------------------------------------------------------
//construction/destruction:

DampingFilter::DampingFilter()
{
	setSampleRate(44100.0);       // sampleRate = 44100 Hz by default

 setGlobalGainFactor(1.0);

	setLowCrossoverFreq(100.0);
	setLowCrossoverGainFactor(0.5);
 setLowGainFactor(0.25);

	setHighCrossoverFreq(8000.0);
	setHighCrossoverGainFactor(0.75);
	setHighGainFactor(0.5);

	resetBuffers();               // reset memorized samples to zero
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

 calcLowShelfCoeffs();
 calcHighShelfCoeffs();
 calcBiquadCoeffs();  
}

void DampingFilter::setGlobalGainFactor(double newGlobalGainFactor)
{
 if( newGlobalGainFactor >= 0.0 )
  globalGainFactor = newGlobalGainFactor;

 // trigger re-calculation of the biquad-coefficients:
 calcBiquadCoeffs();
}

void DampingFilter::setLowCrossoverFreq(double newLowCrossoverFreq)
{
 if( (newLowCrossoverFreq>=20.0) && (newLowCrossoverFreq<=20000.0) )
  lowCrossoverFreq = newLowCrossoverFreq;
	else
		lowCrossoverFreq = 20000.0;

 calcLowShelfCoeffs();
}

void DampingFilter::setLowCrossoverGainFactor(double newLowCrossoverGainFactor)
{
 if( newLowCrossoverGainFactor >= 0.0 )
  lowCrossoverGainFactor = newLowCrossoverGainFactor;
 calcLowShelfCoeffs();
}

void DampingFilter::setLowGainFactor(double newLowGainFactor)
{
 if( newLowGainFactor >= 0.0 )
  lowGainFactor = newLowGainFactor;
 calcLowShelfCoeffs();
}

void DampingFilter::setHighCrossoverFreq(double newHighCrossoverFreq)
{
 if( (newHighCrossoverFreq>=20.0) && (newHighCrossoverFreq<=20000.0) )
  highCrossoverFreq = newHighCrossoverFreq;
	else
		highCrossoverFreq = 20000.0;

 calcHighShelfCoeffs();
}

void DampingFilter::setHighGainFactor(double newHighGainFactor)
{
 if( newHighGainFactor >= 0.0 )
  highGainFactor = newHighGainFactor;
 calcHighShelfCoeffs();
}

void DampingFilter::setHighCrossoverGainFactor(double newHighCrossoverGainFactor)
{
 if( newHighCrossoverGainFactor >= 0.0 )
  highCrossoverGainFactor = newHighCrossoverGainFactor;
 calcHighShelfCoeffs();
}


//----------------------------------------------------------------------------
// others:

void DampingFilter::calcBiquadCoeffs()
{
 // combine the two first order filter-coefficient sets and the global gain 
 // factor into the biquad-coefficient set:
 b0 = globalGainFactor * b0Ls*b0Hs;
 b1 = globalGainFactor * (b0Ls*b1Hs + b1Ls*b0Hs);
 b2 = globalGainFactor * b1Ls*b1Hs;
 a1 = -(a1Ls+a1Hs);
 a2 = -(a1Ls*a1Hs);

 // Remark:
 // the a-coefficients have been given a minus-sign in order to implement the 
 // filters difference-equation as
 // out = b0*in + b1*x1 + b2*x2
 //             + a1*y1 + a2*y2;
 // which turned ot to be more efficient than subtracting weighted past 
 // output samples
}

void DampingFilter::calcLowShelfCoeffs()
{
 double G,                 // gain as factor
        G_B,               // gain at corner-frequency
        epsilon, beta,     // intermediate variables
        x, gainFactor;     // more intermediate variables

 double omegaAna,          // pre-warped analog radian corner freq.
        omegaDig;          // corner freq of the digital filter

 double pProto, zProto,    // prototype pole and zero
        pAna, zAna,        // denormalized analog pole and zero
        pDig, zDig;        // digital pole and zero

 // catch the special condition of zero dB-gain:
 if( fabs(lowGainFactor-1.0) < 0.000001 )   // equality check with margin
 {
  b0Ls = 1.0;
  b1Ls = 0.0;
  a1Ls = 0.0;

  // trigger re-calculation of the biquad-coefficients:
  calcBiquadCoeffs();

  return; // nothing more to do here, in this special case
 }

 // assign/calculate some intermediate variables:
 G   = lowGainFactor;
 G_B = lowCrossoverGainFactor;

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

 // trigger re-calculation of the biquad-coefficients:
 calcBiquadCoeffs();
}

void DampingFilter::calcHighShelfCoeffs()
{
 double G,                 // gain as factor
        G_B,               // gain at corner-frequency
        epsilon, beta,     // intermediate variables
        x, gainFactor;     // more intermediate variables

 double omegaAna,          // pre-warped analog radian corner freq.
        omegaDig;          // corner freq of the digital filter

 double pProto, zProto,    // prototype pole and zero
        pAna, zAna,        // denormalized analog pole and zero
        pDig, zDig;        // digital pole and zero

 // catch the special condition of zero dB-gain:
 if( fabs(highGainFactor-1.0) < 0.000001 )  // equality check with margin
 {
  b0Hs = 1.0;
  b1Hs = 0.0;
  a1Hs = 0.0;

  // trigger re-calculation of the biquad-coefficients:
  calcBiquadCoeffs();

  return; // nothing more to do here, in this special case
 }

 // assign/calculate some intermediate variables:
 G   = highGainFactor;
 G_B = G/highCrossoverGainFactor; 
  // we must take the reciprocal and scale by G, because the prototpye's 
  // frequency response will be inverted when transforming from low-shelf 
  // to high-shelf

 epsilon = sqrt( (G*G-G_B*G_B)/(G_B*G_B-1.0) );
 beta    = 1/epsilon;

 // make low-shelving prototype:
 zProto = -(G*beta);
 pProto = -beta;

 // convert the cutoff-frequency into a normalized cutoff-frequency in
 // radians/sec (a frequency between 0 and pi):
 omegaDig = (2.0*PI*highCrossoverFreq)/sampleRate;

 // prewarp the desired digital cutoff radian frequency to the desired analog
 // cutoff radian frequency:
 omegaAna = 2.0*sampleRate * tan(0.5*omegaDig);

 // transform the filter to the desired selectivity and cutoff-frequency via
 // an analog frequency transformation (LowShelf->High-Shelf):
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

void DampingFilter::resetBuffers()
{
 x1  = 0.0;
 x2  = 0.0;
 y1  = 0.0;
 y2  = 0.0;
 out = 0.0;
}
