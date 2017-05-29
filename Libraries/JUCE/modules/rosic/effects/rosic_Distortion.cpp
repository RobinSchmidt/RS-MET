//#include "rosic_Distortion.h"
//using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction:

Distortion::Distortion()
{
  setSampleRate  (44100.0);    // sampleRate            = 44100 Hz by default


  setOverSampling(4);

  setDrive       (0.0);        // drive                 = 0 dB by default
  setLowThresh   (0.0);        // low  threshold        = 0 dB by default
  setHighThresh  (0.0);        // high threshold        = 0 dB by default
  setLowSat      (0.0);        // low  saturation level = 0 dB by default
  setHighSat     (0.0);        // high saturation level = 0 dB by default
  setLowClamp    (0.0);        // low  clamping level   = 0 dB by default
  setHighClamp   (0.0);        // high clamping level   = 0 dB by default
  setRectify     (0.0);        // no rectification by default
  lowSlopeInf  = true;       // threshold levels are equal to saturation levels
  highSlopeInf = true;
  mode         = 1;          // linear mapping betweem treshhold and clamp

  dcBlocker.setMode  (2);    // highpass mode
  dcBlocker.setCutoff(5.0);  // cutoff @ 5 Hz
}

Distortion::~Distortion()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void Distortion::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  dcBlocker.setSampleRate(sampleRate);
  //antiAliasFlt.setSampleRate(overSampling*sampleRate);
}

void Distortion::setOverSampling(int newOverSampling)
{
  if( newOverSampling > 0 )
    overSampling = newOverSampling;

  overSamplingRec = 1.0/overSampling;
  antiAliasFlt.setSubDivision(overSampling);
}

void Distortion::setMode(int newMode)
{
  if( newMode>=0 && newMode<=2)
    mode = newMode;
}

void Distortion::setDryWet(double newDryWet)
{
  if( newDryWet>=0 && newDryWet<=100 )
  {
    wetVolume = newDryWet/100; //DryWet is assumed to be in % wet
    dryVolume = 1 - wetVolume;
  }
}

void Distortion::setDrive(double newDrive)
{
  drive = DB2AMP(newDrive);
}

void Distortion::setDcOffset(double newDcOffset)
{
  dcOffset = newDcOffset;
}

void Distortion::setLowThresh(double newLowThresh)
{
  lowThresh = DB2AMP(newLowThresh);

  //calculate parameters for the linear mapping:
  if(lowSat==lowThresh)
    lowSlopeInf = true;
  else
  {
    lowSlopeInf = false;
    lowSlope    = (lowClamp-lowThresh)/(lowSat-lowThresh);
  }

  //calculate parameters for the tanh-mapping:
  lowAlpha = lowClamp - lowThresh;
  if(lowAlpha == 0)
    lowAlpha += DBL_MIN;
  lowBeta = 1/lowAlpha;
}

void Distortion::setHighThresh(double newHighThresh)
{
  highThresh = DB2AMP(newHighThresh);

  //calculate parameters for the linear mapping:
  if(highSat==highThresh)
    highSlopeInf = true;
  else
  {
    highSlopeInf = false;
    highSlope    = (highClamp-highThresh)/(highSat-highThresh);
  }

  //calculate parameters for the tanh-mapping:
  highAlpha = highClamp - highThresh;
  if(highAlpha == 0)
    highAlpha += DBL_MIN;
  highBeta = 1/highAlpha;
}

void Distortion::setLowSat(double newLowSat)
{
  lowSat = DB2AMP(newLowSat);

  if(lowSat==lowThresh)
    lowSlopeInf = true;
  else
  {
    lowSlopeInf = false;
    lowSlope    = (lowClamp-lowThresh)/(lowSat-lowThresh);
  }
}

void Distortion::setHighSat(double newHighSat)
{
  highSat = DB2AMP(newHighSat);

  if(highSat==highThresh)
    highSlopeInf = true;
  else
  {
    highSlopeInf = false;
    highSlope   = (highClamp-highThresh)/(highSat-highThresh);
  }
}

void Distortion::setLowClamp(double newLowClamp)
{
  lowClamp = DB2AMP(newLowClamp);

  //calculate parameters for the linear mapping:
  if(lowSat==lowThresh)
    lowSlopeInf = true;
  else
  {
    lowSlopeInf = false;
    lowSlope    = (lowClamp-lowThresh)/(lowSat-lowThresh);
  }

  //calculate parameters for the tanh-mapping:
  lowAlpha = lowClamp - lowThresh;
  if(lowAlpha == 0)
    lowAlpha += DBL_MIN;
  lowBeta = 1/lowAlpha;
}

void Distortion::setHighClamp(double newHighClamp)
{
  highClamp = DB2AMP(newHighClamp);

  //calculate parameters for the linear mapping:
  if(highSat==highThresh)
    highSlopeInf = true;
  else
  {
    highSlopeInf = false;
    highSlope   = (highClamp-highThresh)/(highSat-highThresh);
  }

  //calculate parameters for the tanh-mapping:
  highAlpha = highClamp - highThresh;
  if(highAlpha == 0)
    highAlpha += DBL_MIN;
  highBeta = 1/highAlpha;
}

void Distortion::setRectify(double newRectify)
{
  if( newRectify>=0 && newRectify<=1 )
    rectify = 2*newRectify;
}

