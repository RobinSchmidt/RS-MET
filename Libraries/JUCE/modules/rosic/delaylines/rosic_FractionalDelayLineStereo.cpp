//#include "rosic_FractionalDelayLineStereo.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FractionalDelayLineStereo::FractionalDelayLineStereo(int maximumDelayInSamples)
{
  // parameter initializition:
  tapIn               = 0;
  tapOut              = 0;
  delayInSeconds      = 0.001;
  sampleRate          = 44100.0; // default sample rate
  maxDelayInSamples   = maximumDelayInSamples;
  setDelayTime(delayInSeconds);

  // allocate memory for the actual delay-line:
  delayBufferL = new double[maxDelayInSamples+interpolatorMargin];
  delayBufferR = new double[maxDelayInSamples+interpolatorMargin];

  // initialize the content of the delay-lines with zeros:
  clearDelayBuffers();
}

FractionalDelayLineStereo::~FractionalDelayLineStereo()
{
  if( delayBufferL != NULL )
    delete[] delayBufferL;
  if( delayBufferR != NULL )
    delete[] delayBufferR;
}

//-------------------------------------------------------------------------------------------------
// parameter settings (set-functions):

void FractionalDelayLineStereo::setSampleRate(double newSampleRate)
{
  if(newSampleRate > 0.01)
    sampleRate = newSampleRate;

  setDelayTime(delayInSeconds);
}

// \todo inline this function
void FractionalDelayLineStereo::setDelayTime(double newDelayTime)
{
  if(newDelayTime > 0.0)
    delayInSeconds = newDelayTime;
  else
    delayInSeconds  = 0.0;

  delayInSamples = sampleRate*delayInSeconds;

  // cap the delay in samples to the maximum allowed delay-time:
  if( delayInSamples > (double) (maxDelayInSamples-1) )
    delayInSamples = (double) (maxDelayInSamples-1);
    ///< \todo optimization - avoid the two typecasts by means of maintaining a double-variable for
    // the maximum delay

  // depending on the chosen interpolation-method, we also need a minimum delay-time:
  if( delayInSamples < interpolatorMargin )
    delayInSamples = interpolatorMargin;
    ///< \todo make this minimum delay dependent on the interpolation-method

  // calculate the integer and fractional parts of the delay:
  double tmp                   = floor(delayInSamples);
  delayInSamplesFractionalPart = delayInSamples - tmp;
  frac                         = 1.0 - delayInSamplesFractionalPart;
  delayInSamplesIntegerPart    = (int) tmp;

  // adjust tapOut-pointer:
  tapOut = tapIn - delayInSamplesIntegerPart - 1;
  if( tapOut < 0 )
    tapOut += maxDelayInSamples;
}

void FractionalDelayLineStereo::setInterpolationMethod(int newMethod)
{
  interpolatorL.setInterpolationMethod(newMethod);
  interpolatorR.setInterpolationMethod(newMethod);
}

//-------------------------------------------------------------------------------------------------
// others:

void FractionalDelayLineStereo::clearDelayBuffers()
{
  for(int i=0; i<maxDelayInSamples+interpolatorMargin; i++)
  {
    delayBufferL[i] = 0.0;
    delayBufferR[i] = 0.0;
  }
  interpolatorL.reset();
  interpolatorR.reset();
}
