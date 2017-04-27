#include "rosic_InfiniteImpulseResponseDesigner.h"
using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

InfiniteImpulseResponseDesigner::InfiniteImpulseResponseDesigner()
{
  sampleRate     = 44100.0;
  frequency      = 1000.0;
  setBandwidth(1.0); // sets up lowerFrequency and upperFrequency
  gain           = amp2dB(0.25);
  mode           = LOWPASS;
  prototypeOrder = 2;
}

InfiniteImpulseResponseDesigner::~InfiniteImpulseResponseDesigner()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void InfiniteImpulseResponseDesigner::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

void InfiniteImpulseResponseDesigner::setMode(int newMode)
{
  if( newMode >= BYPASS && newMode <= PEAK )
    mode = newMode;
  else
    DEBUG_BREAK;  // newMode must be one of the numbers defined in enum 'modes'
}

void InfiniteImpulseResponseDesigner::setPrototypeOrder(int newOrder)
{
  if( newOrder >= 1 )
    prototypeOrder = newOrder;
  else
    DEBUG_BREAK;  // newOrder must be at least 1
}

void InfiniteImpulseResponseDesigner::setApproximationMethod(int newApproximationMethod)
{
  prototypeDesigner.setApproximationMethod(newApproximationMethod);
}

void InfiniteImpulseResponseDesigner::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  calculateLowerAndUpperFrequency();
}

void InfiniteImpulseResponseDesigner::setBandwidth(double newBandwidth)
{
  bandwidth = newBandwidth;
  calculateLowerAndUpperFrequency();
}

void InfiniteImpulseResponseDesigner::setLowerFrequency(double newLowerFrequency)
{
  if( newLowerFrequency > 0.0 )
    lowerFrequency = newLowerFrequency;
  else
    DEBUG_BREAK;  // negative frequencies are not allowed
}

void InfiniteImpulseResponseDesigner::setUpperFrequency(double newUpperFrequency)
{
  if( newUpperFrequency > 0.0 )
    upperFrequency = newUpperFrequency;
  else
    DEBUG_BREAK;  // negative frequencies are not allowed
}

void InfiniteImpulseResponseDesigner::setGain(double newGain)
{
  gain = newGain;
}

void InfiniteImpulseResponseDesigner::setRipple(double newRipple)
{
  prototypeDesigner.setPassbandRipple(newRipple);
  prototypeDesigner.setPassbandGainRatio(1.0-0.01*newRipple);
  prototypeDesigner.setStopbandGainRatio(0.01*newRipple);
}

void InfiniteImpulseResponseDesigner::setStopbandRejection(double newStopbandRejection)
{
  prototypeDesigner.setStopbandRejection(newStopbandRejection);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

bool InfiniteImpulseResponseDesigner::doesModeDoubleTheOrder()
{
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    return true;
  else
    return false;
}

int InfiniteImpulseResponseDesigner::getFinalFilterOrder()
{
  if( doesModeDoubleTheOrder() == true )
    return 2*prototypeOrder;
  else
    return prototypeOrder;
}

int InfiniteImpulseResponseDesigner::getNumBiquadStages()
{
  int order = getFinalFilterOrder();
  if( isEven(order) )
    return order/2;
  else
    return (order+1)/2;
}

bool InfiniteImpulseResponseDesigner::hasCurrentModeBandwidthParameter()
{
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    return true;
  else
    return false;
}

bool InfiniteImpulseResponseDesigner::hasCurrentModeGainParameter()
{
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK )
    return true;
  else
    return false;
}

bool InfiniteImpulseResponseDesigner::hasCurrentModeRippleParameter()
{
  if( prototypeDesigner.hasCurrentMethodRippleParameter() )
  {
    if( mode != LOW_SHELV || mode != HIGH_SHELV || mode != PEAK )
      return true;
    else
      return false;
  }
  else
    return false;
}

bool InfiniteImpulseResponseDesigner::hasCurrentModeRejectionParameter()
{
  if( prototypeDesigner.hasCurrentMethodRejectionParameter() )
  {
    if( mode != LOW_SHELV || mode != HIGH_SHELV || mode != PEAK )
      return true;
    else
      return false;
  }
  else
    return false;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// coefficient retrieval:

void InfiniteImpulseResponseDesigner::getPolesAndZeros(Complex* poles, Complex* zeros)
{
  // calculate the required order and number of biquads for the filter:
  int finalOrder;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    finalOrder = 2*prototypeOrder;
  else
    finalOrder = prototypeOrder;

//  int numBiquads;
//  if( isEven(finalOrder) )
//    numBiquads = finalOrder/2;
//  else
//    numBiquads = (finalOrder+1)/2;

  // set up the type for the prototype (lowpass/low-shelv):
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK || mode == BYPASS )
    prototypeDesigner.setPrototypeMode(PrototypeDesigner::LOWSHELV_PROTOTYPE);
  else
    prototypeDesigner.setPrototypeMode(PrototypeDesigner::LOWPASS_PROTOTYPE);

  double  f1, f2, wd1, wd2, wa1, wa2;
  double  fs = sampleRate;
  int     k;

  // use sanity-checked local frequency variables here:
  f2 = upperFrequency;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    f1 = lowerFrequency;
    if( f2 > 0.999*0.5*fs )
      f2 = 0.999*0.5*fs;  // ensure upper frequency < sampleRate/2
    if( f1 > 0.999*f2  )
      f1 = 0.999*f2;      // ensure lower frequency < upper frequency
  }
  else
  {
    f1 = frequency;
    if( f1 > 0.999*0.5*fs )
      f1 = 0.999*0.5*fs;  // ensure frequency < sampleRate/2
  }

  // prewarp the frequencies to the desired frequencies required for the design of the
  // (unnormalized) analog prototype filter:
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    wd1 = 2.0*PI*f1/fs;        // normalized digital radian frequency 1
    wa1 = 2.0*fs*tan(0.5*wd1); // pre-warped analog radian frequency 1
    wd2 = 2.0*PI*f2/fs;        // normalized digital radian frequency 2
    wa2 = 2.0*fs*tan(0.5*wd2); // pre-warped analog radian frequency 2
  }
  else
  {
    wd1 = 2.0*PI*f1/fs;        // normalized digital radian frequency 1
    wa1 = 2.0*fs*tan(0.5*wd1); // pre-warped analog radian frequency 1
    wd2 = 0.0;                 // unused
    wa2 = 0.0;                 // unused
  }

  // set up the prototype-designer according to the specifications:
  prototypeDesigner.setOrder(prototypeOrder);
  prototypeDesigner.setGain(gain);
  if( mode == LOWPASS || mode == HIGHPASS || mode == BANDPASS || mode == BANDREJECT )
    prototypeDesigner.setReferenceGain(NEG_INF);
  else
    prototypeDesigner.setReferenceGain(0.0);

  // allocate temporary memory:
  Complex* protoPoles = new Complex[prototypeOrder];
  Complex* protoZeros = new Complex[prototypeOrder];

  // design the analog prototype filter:
  if( mode == HIGH_SHELV )
  {
    // we need to cope with exchange of roles of poles and zeros for high-shelving (inverse)
    // chebychevs because the low-shelv -> high-shelv frequency transform exchanges these roles
    // once again:
    if( prototypeDesigner.getApproximationMethod() == PrototypeDesigner::CHEBYCHEV )
    {
      prototypeDesigner.setApproximationMethod(PrototypeDesigner::INVERSE_CHEBYCHEV);
      prototypeDesigner.getPolesAndZeros(poles, zeros);
      prototypeDesigner.setApproximationMethod(PrototypeDesigner::CHEBYCHEV);
    }
    else if( prototypeDesigner.getApproximationMethod() == PrototypeDesigner::INVERSE_CHEBYCHEV )
    {
      prototypeDesigner.setApproximationMethod(PrototypeDesigner::CHEBYCHEV);
      prototypeDesigner.getPolesAndZeros(poles, zeros);
      prototypeDesigner.setApproximationMethod(PrototypeDesigner::INVERSE_CHEBYCHEV);
    }
    else
      prototypeDesigner.getPolesAndZeros(poles, zeros);
  }
  else
    prototypeDesigner.getPolesAndZeros(poles, zeros);

  // because the PrototypeDesigner returns only one representant for each pair of complex conjugate poles/zeros, we now create the full
  // set here:
  if( isOdd(prototypeOrder) )
  {
    // copy the real pole/zero to the end:
    protoPoles[prototypeOrder-1] = poles[prototypeOrder/2];
    protoZeros[prototypeOrder-1] = zeros[prototypeOrder/2];
  }
  // for each complex pole/zero create a pair of complex conjugates:
  for(k=0; k<prototypeOrder/2; k++)
  {
    protoPoles[2*k]   = poles[k];
    protoPoles[2*k+1] = poles[k].getConjugate();
    protoZeros[2*k]   = zeros[k];
    protoZeros[2*k+1] = zeros[k].getConjugate();
  }

  // s-plane frequency transformation:
  switch( mode )
  {
  case LOWPASS:    PoleZeroMapper::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case HIGHPASS:   PoleZeroMapper::sPlanePrototypeToHighpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case BANDPASS:   PoleZeroMapper::sPlanePrototypeToBandpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  case BANDREJECT: PoleZeroMapper::sPlanePrototypeToBandreject(protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  case LOW_SHELV:  PoleZeroMapper::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  case HIGH_SHELV:
    {
      // this is ugly - rewrite the prototype design code such that all approximation methods can be treated uniformly:
      if( prototypeDesigner.needsSpecialHighShelvTransform() )
        PoleZeroMapper::sPlanePrototypeToHighShelv( protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);
      else
        PoleZeroMapper::sPlanePrototypeToHighpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);
    }
    break;
  case PEAK:       PoleZeroMapper::sPlanePrototypeToBandpass(  protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1, wa2); break;
  default:         PoleZeroMapper::sPlanePrototypeToLowpass(   protoPoles, protoZeros, poles, zeros, prototypeOrder, wa1);      break;
  };

  // transform from the s-domain to the z-domain via the bilinear transform:
  double g; // not used actually
  PoleZeroMapper::bilinearAnalogToDigital(poles, finalOrder, zeros, finalOrder, fs, &g);

  // free dynamically allocated memory:
  delete[] protoPoles;
  delete[] protoZeros;
}

void InfiniteImpulseResponseDesigner::getBiquadCascadeCoefficients(double *b0, double *b1, double *b2, double *a1, double *a2)
{
  // calculate the required order and number of biquads for the filter:
  int finalOrder, numBiquads;

  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
    finalOrder = 2*prototypeOrder;
  else
    finalOrder = prototypeOrder;

  if( isEven(finalOrder) )
    numBiquads = finalOrder/2;
  else
    numBiquads = (finalOrder+1)/2;

  // set up the type for the prototype (lowpass/low-shelv):
  if( mode == LOW_SHELV || mode == HIGH_SHELV || mode == PEAK || mode == BYPASS )
  {
    prototypeDesigner.setPrototypeMode(PrototypeDesigner::LOWSHELV_PROTOTYPE);
    if( isCloseTo(gain, 0.0, 0.001) || mode == BYPASS ) // gains of zero yield a 'bypass' filter
    {
      for(int b = 0; b < numBiquads; b++)
      {
        b0[b] = 1.0;
        b1[b] = 0.0;
        b2[b] = 0.0;
        a1[b] = 0.0;
        a2[b] = 0.0;
        return;
      }
    }
  }
  else
    prototypeDesigner.setPrototypeMode(PrototypeDesigner::LOWPASS_PROTOTYPE);

  double  f1, f2, wd1, wd2, wa1, wa2;
  double  fs = sampleRate;

  // use sanity-checked local frequency variables here:
  f2 = upperFrequency;
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    f1 = lowerFrequency;
    if( f2 > 0.99*0.5*fs )
      f2 = 0.99*0.5*fs;  // ensure upper frequency < sampleRate/2
    if( f1 > 0.99*f2  )
      f1 = 0.99*f2;      // ensure lower frequency < upper frequency
  }
  else
  {
    f1 = frequency;
    if( f1 > 0.99*0.5*fs )
      f1 = 0.99*0.5*fs;  // ensure frequency < sampleRate/2
  }

  // prewarp the frequencies to the desired frequencies required for the design of the (unnormalized) analog prototype filter:
  if( mode == BANDPASS || mode == BANDREJECT || mode == PEAK )
  {
    wd1 = 2.0*PI*f1/fs;         // normalized digital radian frequency 1
    wa1 = 2.0*fs*tan(0.5*wd1);  // pre-warped analog radian frequency 1
    wd2 = 2.0*PI*f2/fs;         // normalized digital radian frequency 2
    wa2 = 2.0*fs*tan(0.5*wd2);  // pre-warped analog radian frequency 2
  }
  else
  {
    wd1 = 2.0*PI*f1/fs;         // normalized digital radian frequency 1
    wa1 = 2.0*fs*tan(0.5*wd1);  // pre-warped analog radian frequency 1
    wd2 = 0.0;                  // unused
    wa2 = 0.0;                  // unused
  }

  // allocate temporary memory:
  Complex* poles      = new Complex[finalOrder];
  Complex* zeros      = new Complex[finalOrder];

  getPolesAndZeros(poles, zeros);
  FilterCoefficientConverter::polesAndZerosToBiquadCascade(poles, zeros, finalOrder, b0, b1, b2, a1, a2);

  double deWarpedPassbandCenter = 2.0*atan(sqrt(tan(0.5*wd1)*tan(0.5*wd2)));
  normalizeGain(b0, b1, b2, a1, a2, deWarpedPassbandCenter, numBiquads);

  // free dynamically allocated memory:
  delete[] poles;
  delete[] zeros;
}

void InfiniteImpulseResponseDesigner::getDirectFormCoefficients(double *b, double *a)
{
  int numBiquads = getNumBiquadStages();
  double *b0 = new double[numBiquads];
  double *b1 = new double[numBiquads];
  double *b2 = new double[numBiquads];
  double *a1 = new double[numBiquads];
  double *a2 = new double[numBiquads];
  getBiquadCascadeCoefficients(b0, b1, b2, a1, a2);
  FilterCoefficientConverter::biquadCascadeToDirectForm(numBiquads, b0, b1, b2, a1, a2, b, a);
  delete[] b0;
  delete[] b1;
  delete[] b2;
  delete[] a1;
  delete[] a2;
}

void InfiniteImpulseResponseDesigner::calculateLowerAndUpperFrequency()
{
  lowerFrequency = frequency / pow(2.0, 0.5*bandwidth);
  upperFrequency = lowerFrequency * pow(2.0, bandwidth);
}

void InfiniteImpulseResponseDesigner::normalizeGain(double *b0, double *b1, double *b2,
                                               double *a1, double *a2, double wc, int numBiquads)
{
  double w;
  if( mode == LOWPASS || mode == BANDREJECT || mode == HIGH_SHELV || mode == PEAK )
    w = 0.0;           // normalize at DC
  else if( mode == HIGHPASS || mode == LOW_SHELV )
    w = PI;            // normalize at Nyquist frequency
  else if( mode == BANDPASS )
    w = wc;            // normalize at passband center frequency

  double normalizeFactor = 1.0;
  if( prototypeDesigner.getPrototypeMode() == PrototypeDesigner::LOWSHELV_PROTOTYPE )
  {
    if(  prototypeDesigner.getApproximationMethod() == PrototypeDesigner::INVERSE_CHEBYCHEV
      || prototypeDesigner.getApproximationMethod() == PrototypeDesigner::ELLIPTIC )
    {
      if( isEven(prototypeDesigner.getOrder()) )
      {
        double factor   = 1.0-prototypeDesigner.getPassbandGainRatio();
        double excessDb = -factor * gain;
        normalizeFactor = dB2amp(-excessDb/numBiquads);
      }
    }
  }
  else // prototype is lowpass
  {
    if(  prototypeDesigner.getApproximationMethod() == PrototypeDesigner::CHEBYCHEV
      || prototypeDesigner.getApproximationMethod() == PrototypeDesigner::ELLIPTIC )
    {
      if( isEven(prototypeDesigner.getOrder()) )
      {
        normalizeFactor = 1.0 / dB2amp(prototypeDesigner.getPassbandRipple());
        if( doesModeDoubleTheOrder() )
          normalizeFactor = pow(normalizeFactor, 1.0/prototypeDesigner.getOrder());
        else
          normalizeFactor = pow(normalizeFactor, 2.0/prototypeDesigner.getOrder());
      }
    }
  }

  FilterCoefficientConverter::normalizeBiquadStages(b0, b1, b2, a1, a2, w, numBiquads, normalizeFactor);
}
