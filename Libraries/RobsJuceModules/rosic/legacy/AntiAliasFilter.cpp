#include "AntiAliasFilter.h"

//----------------------------------------------------------------------------
//construction/destruction:
AntiAliasFilter::AntiAliasFilter()
{
 //init member variables:
 sampleRate = 44100.0;
 cutoff     = 18000.0;

 //set up the IIRDesigner:
 filtDesigner.setSampleRate(sampleRate);
 filtDesigner.setFreq1(cutoff);
 //fltDesigner.setMethod(1);               //chooses butterworth design
 filtDesigner.setMode(1);                 //chooses lowpass-mode
 filtDesigner.setSlope(4);                //chooses 4th order filter

 //let the fltDesigner calculate the filter-coefficients:
 filtDesigner.getDirectFormCoeffs(b, a);

 // reset the internal buffers
 reset();
}

AntiAliasFilter::~AntiAliasFilter()
{

}

//----------------------------------------------------------------------------
//parameter settings:
void AntiAliasFilter::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 //update sampleRate in the IirDesigner:
 filtDesigner.setSampleRate(sampleRate);

 //calculate new filter coefficients:
 filtDesigner.getDirectFormCoeffs(b, a);
}

void AntiAliasFilter::setCutoff(double newCutoff)
{
 if( newCutoff > 0 )
  cutoff = newCutoff;

 //update cutoff frequency in the fltDesigner:
 filtDesigner.setFreq1(cutoff);

 //calculate new filter coefficients:
 filtDesigner.getDirectFormCoeffs(b, a);
}

//----------------------------------------------------------------------------
//others:
void AntiAliasFilter::reset()
{
 int i;
 for(i=0; i<5; i++)
 {
  x[i] = 0.0;
  y[i] = 0.0;
 }
}