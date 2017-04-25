#include "DirectFormIir.h"

DirectFormIir::DirectFormIir()
{
 //init member variables:
 sampleRate  = 44100.0;
 order       = 4;

 filterDesigner.getDirectFormCoeffs(ffCoeffs, fbCoeffs);
}

DirectFormIir::~DirectFormIir()
{

}

//-----------------------------------------------------------------------------
//parameter settings:
void DirectFormIir::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);
}

void DirectFormIir::setSlope(int newSlope)
{
 if( newSlope>=1 && newSlope<=12 )
  slope = newSlope;

 filterDesigner.setSlope(slope);
 filterDesigner.getDirectFormCoeffs(ffCoeffs, fbCoeffs);

 if(mode>=3)
  order = 2*slope;
 else
  order = slope;

 //reset();
}

void DirectFormIir::setMode(int newMode)
{
 mode = newMode;

 filterDesigner.setMode(mode);
 filterDesigner.getDirectFormCoeffs(ffCoeffs, fbCoeffs);

 if(mode>=3)
  order = 2*slope;
 else
  order = slope;

 //reset();
}

void DirectFormIir::setApproximationMethod(int newMethod)
{
 method = newMethod;
 //...
}

void DirectFormIir::setFreq1(double newFreq1)
{
 if( (newFreq1>=0) && (newFreq1<=0.5*sampleRate) )
  freq1 = newFreq1;

 filterDesigner.setFreq1(freq1);
 filterDesigner.getDirectFormCoeffs(ffCoeffs, fbCoeffs);
}

void DirectFormIir::setFreq2(double newFreq2)
{
 if( (newFreq2>=0) && (newFreq2<=0.5*sampleRate) )
  freq2 = newFreq2;

 filterDesigner.setFreq2(freq2);
 filterDesigner.getDirectFormCoeffs(ffCoeffs, fbCoeffs);
}

//-----------------------------------------------------------------------------
//others:
void DirectFormIir::reset()
{
 for(long k=1; k<(maxOrder+1); k++)
 {
  xBuffer[k] = 0;
  yBuffer[k] = 0;
 }
}

void DirectFormIir::getMagnitudeResponse(float *frequencies, 
                                         float *magnitudes, 
                                         int    numBins, 
                                         bool   inDecibels)
{

}