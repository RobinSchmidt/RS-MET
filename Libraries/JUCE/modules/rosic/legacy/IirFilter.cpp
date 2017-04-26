#include "IirFilter.h"

//----------------------------------------------------------------------------
//construction/destruction:

IirFilter::IirFilter()
{
	setSampleRate(44100.0); // samplerate = 44100 Hz by default
	setMode      (1);       // lowpass by default
	setSlope     (2);       // 2nd order filter
	setFreq1     (1000.0);  // (lower) cutoff frequency
	setFreq2     (2000.0);  // upper cutoff frequency
 setEqGain    (2.0);     // gain for shelving and peak modes
 setGlobalGain(1.0);     // global gain factor for the signal

	resetBuffers();         //reset the filters memory buffers (to 0)
}

IirFilter::~IirFilter()
{}

//----------------------------------------------------------------------------
//parameter settings:

void IirFilter::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 filter.setSampleRate(sampleRate);
 designer.setSampleRate(sampleRate);
 updateBiquadCoeffs();
}

void IirFilter::setMode(int newMode)
{
 mode = newMode;
 updateBiquadCoeffs();
}

void IirFilter::setSlope(int newSlope)
{
 slope = newSlope;
 updateBiquadCoeffs();
}

void IirFilter::setFreq1(double newFreq1)
{
 freq1 = newFreq1;
 updateBiquadCoeffs();
}

void IirFilter::setFreq2(double newFreq2)
{
 freq2 = newFreq2;
 updateBiquadCoeffs();
}

void IirFilter::setEqGain(double newEqGain)
{
 eqGain = newEqGain;
 updateBiquadCoeffs();
}

void IirFilter::setGlobalGain(double newGlobalGain)
{
 globalGain = newGlobalGain;
 updateBiquadCoeffs();
}

//----------------------------------------------------------------------------
//others:

void IirFilter::resetBiquadCoeffs()
{
 static intA b; // index for the biquad-stage 

 for(b=0; b<maxNumStages; b++)
 {
  b0[b] = 1.0;
  b1[b] = 0.0;
  b2[b] = 0.0;
  a0[b] = 1.0;
  a1[b] = 0.0;
  a2[b] = 0.0;
 }
}

void IirFilter::updateBiquadCoeffs()
{
 // reset the biquad-coefficients to neutral values:
 resetBiquadCoeffs();

 switch( mode )
 {
 case IirDesigner::BYPASS:
  {  
   designer.setMode(IirDesigner::BYPASS);
   filter.setNumStages(0);
  }
  break;

 case IirDesigner::LOWPASS:
  {
   designer.setMode(IirDesigner::LOWPASS);
   if( isEven(slope) )
    filter.setNumStages( slope/2 );
   else
    filter.setNumStages( (slope+1)/2 );
  }
  break;

 case IirDesigner::HIGHPASS:
  {
   designer.setMode(IirDesigner::HIGHPASS);
   if( isEven(slope) )
    filter.setNumStages( slope/2 );
   else
    filter.setNumStages( (slope+1)/2 );
  }
  break;

 case IirDesigner::BANDPASS:
  {
   designer.setMode(IirDesigner::BANDPASS);
   filter.setNumStages(slope);
  }
  break;

 case IirDesigner::BANDREJECT:
  {
   designer.setMode(IirDesigner::BANDREJECT);
   filter.setNumStages(slope);
  }
  break;

 case IirDesigner::LOW_SHELV:
  {
   designer.setMode(IirDesigner::LOW_SHELV);
   if( isEven(slope) )
    filter.setNumStages( slope/2 );
   else
    filter.setNumStages( (slope+1)/2 );
  }
  break;

 case IirDesigner::HIGH_SHELV:
  {
   designer.setMode(IirDesigner::HIGH_SHELV);
   if( isEven(slope) )
    filter.setNumStages( slope/2 );
   else
    filter.setNumStages( (slope+1)/2 );
  }
  break;

 case IirDesigner::PEAK:
  {
   designer.setMode(IirDesigner::PEAK);
   filter.setNumStages(slope);
  }
  break;

 default:
  {
   designer.setMode(IirDesigner::BYPASS);
   filter.setNumStages(0);
  }
 } // end of switch(mode)

 // set up the designer:
 designer.setSlope(slope);
 designer.setFreq1(freq1);
 designer.setFreq2(freq2);
 designer.setGain(eqGain);

 // let the designer calculate the filter-coefficients:
 designer.getBiquadCascadeCoeffs(b0, b1, b2, a1, a2);

 // pass the coefficients to the BiquadCascade-object:
 filter.setCoeffs(b0, b1, b2, a1, a2, globalGain);
}

void IirFilter::resetBuffers()
{
 filter.reset();
}



