//#include "rosic_FiniteImpulseResponseDesigner.h"
//using namespace rosic;

//-----------------------------------------------------------------------------------------------------------------------------------------
// construction/destruction:

FiniteImpulseResponseDesigner::FiniteImpulseResponseDesigner()
{
  sampleRate     = 44100.0;
  frequency      = 1000.0;
  setBandwidth(1.0); // sets up lowerFrequency and upperFrequency
  mode           = LOWPASS;
  windowType     = WindowDesigner::BLACKMAN;
}

FiniteImpulseResponseDesigner::~FiniteImpulseResponseDesigner()
{

}

//-----------------------------------------------------------------------------------------------------------------------------------------
// parameter settings:

void FiniteImpulseResponseDesigner::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;
}

void FiniteImpulseResponseDesigner::setMode(int newMode)
{ 
  mode = newMode;
}

void FiniteImpulseResponseDesigner::setFrequency(double newFrequency)
{
  frequency = newFrequency;
  calculateLowerAndUpperFrequency();
}

void FiniteImpulseResponseDesigner::setBandwidth(double newBandwidth)
{
  bandwidth = newBandwidth;
  calculateLowerAndUpperFrequency();
}

void FiniteImpulseResponseDesigner::setLowerFrequency(double newLowerFrequency)
{
  if( newLowerFrequency > 0.0 )
    lowerFrequency = newLowerFrequency;
  else
    DEBUG_BREAK;  // negative frequencies are not allowed
}

void FiniteImpulseResponseDesigner::setUpperFrequency(double newUpperFrequency)
{
  if( newUpperFrequency > 0.0 )
    upperFrequency = newUpperFrequency;
  else
    DEBUG_BREAK;  // negative frequencies are not allowed
}

void FiniteImpulseResponseDesigner::setWindowType(int newWindow)
{
  windowType = newWindow;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// inquiry:

bool FiniteImpulseResponseDesigner::hasCurrentModeBandwidthParameter()
{
  if( mode == BANDPASS || mode == BANDREJECT )
    return true;
  else
    return false;
}

void FiniteImpulseResponseDesigner::getImpulseResponse(double *impulseResponse, int length)
{
  switch( mode )
  {
  case LOWPASS:    
    getLowpassResponse(impulseResponse, length, 2*PI*frequency/sampleRate, windowType); 
    break;
  case HIGHPASS:   
    getHighpassResponse(impulseResponse, length, 2*PI*frequency/sampleRate, windowType);                                       
    break;
  case BANDREJECT: 
    getBandrejectResponse(impulseResponse, length, 2*PI*lowerFrequency/sampleRate, 2*PI*upperFrequency/sampleRate, windowType);  
    break;
  case BANDPASS:   
    getBandpassResponse(impulseResponse, length, 2*PI*lowerFrequency/sampleRate, 2*PI*upperFrequency/sampleRate, windowType);  
    break;
  case DIFFERENTIATOR:    
    getDifferentiatorResponse(impulseResponse, length, windowType);       
    break;
  case HILBERT:    
    getHilbertTransformerResponse(impulseResponse, length, windowType);       
    break;
  default: RAPT::rsArrayTools::fillWithZeros(impulseResponse, length);
  }
}

void FiniteImpulseResponseDesigner::getLowpassResponse(double *impulseResponse, int length, double omega, int windowType)
{
  rassert( length >= 1 && RAPT::rsIsOdd(length) ); // only odd lengths are supported

  // we use a similar naming-convention as in the eq. 16-4 in dspguide, page 290:
  double *h  = impulseResponse;
  int     M  = length-1; 
  double  M2 = (double) M / 2.0;
  for(int i=0; i<length; i++)
  {
    if( i == M2 )   // divide by zero - treat it separately 
     h[i] = omega;
    else
     h[i] = sin(omega*(i-M2)) / (i-M2);
  }

  WindowDesigner::applyWindow(impulseResponse, length, windowType);
  normalizeSumToUnity(impulseResponse, length); // normalizes to unit gain at DC
}

void FiniteImpulseResponseDesigner::getHighpassResponse(double *impulseResponse, int length, double omega, int windowType)
{
  getLowpassResponse(impulseResponse, length, omega, windowType);
  spectralInversion(impulseResponse, length);
}

void FiniteImpulseResponseDesigner::getBandpassResponse(double *impulseResponse, int length, double omegaLow, double omegaHigh, 
                                                        int windowType)
{
  getBandrejectResponse(impulseResponse, length, omegaLow, omegaHigh, windowType);
  spectralInversion(impulseResponse, length);
}

void FiniteImpulseResponseDesigner::getBandrejectResponse(double *impulseResponse, int length, double omegaLow, double omegaHigh, 
                                                          int windowType)
{
  // bandreject filters are obtained by adding a lowpass and a highpass output (which amounts to adding their impulse responses):
  getLowpassResponse(impulseResponse, length, omegaLow, windowType);

  double *tmpHighpass = new double[length];

  getHighpassResponse(tmpHighpass, length, omegaHigh, windowType);
  for(int i=0; i<length; i++)
    impulseResponse[i] += tmpHighpass[i];

  delete[] tmpHighpass;
}

void FiniteImpulseResponseDesigner::getDifferentiatorResponse(double *impulseResponse, int length, int windowType)
{
  double *h = impulseResponse;
  int    M  = length-1; 
  double M2 = (double) M / 2.0;
  for(int n=0; n<=M; n++)
  {
    double x = n-M2;
    if( x == 0.0 )
      h[n] = 0.0;
    else
      h[n] = (cos(PI*x) / x) - (sin(PI*x) / (PI*x*x)); 
      // \todo simplify this - the cosines and sines are actually just alternating -1s and +1s (i think so, at least).
  }
  WindowDesigner::applyWindow(impulseResponse, length, windowType);
}

void FiniteImpulseResponseDesigner::getHilbertTransformerResponse(double *impulseResponse, int length, int windowType)
{
  rassert( RAPT::rsIsOdd(length) ); // only odd lengths supported

  double *h   = impulseResponse;
  double gain = 2.0/PI; 
  int    M    = length-1; 
  int    M2   = M / 2;

  for(int i=M2; i<length; i++)
  {
    if( RAPT::rsIsEven(i-M2) )  
    {
      h[i]         = 0;
      h[M2-(i-M2)] = 0;
    }
    else 
    {
      h[i]         =  gain / (i-M2);
      h[M2-(i-M2)] = -gain / (i-M2);
    }
  }

  WindowDesigner::applyWindow(impulseResponse, length, windowType);
}

void FiniteImpulseResponseDesigner::spectralInversion(double *impulseResponse, int length)
{
  rassert( RAPT::rsIsOdd(length) ); // spectral inversion works only for odd lengths
  for(int i=0; i<length; i++)
    impulseResponse[i] = -impulseResponse[i];
  impulseResponse[(length-1)/2] += 1.0;
}

void FiniteImpulseResponseDesigner::spectralReversal(double *impulseResponse, int length)
{
  for(int i=0; i<length; i+=2)
    impulseResponse[i] = -impulseResponse[i];
}

void FiniteImpulseResponseDesigner::normalizeSumToUnity(double *impulseResponse, int length)
{
  double sum = 0.0;
  for(int i=0; i<length; i++)
    sum += impulseResponse[i];
  double normalizer = 1.0/sum;
  for(int i=0; i<length; i++)
    impulseResponse[i] *= normalizer;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// internal functions:

void FiniteImpulseResponseDesigner::calculateLowerAndUpperFrequency()
{
  lowerFrequency = frequency / pow(2.0, 0.5*bandwidth);  
  upperFrequency = lowerFrequency * pow(2.0, bandwidth); 
  // \todo move this function out of this class for sharing
}
