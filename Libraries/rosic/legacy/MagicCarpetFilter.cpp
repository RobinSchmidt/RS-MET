#include "MagicCarpetFilter.h"

//----------------------------------------------------------------------------
//construction/destruction
MagicCarpetFilter::MagicCarpetFilter()
{
 mode          = LOWPASS_12;
 freq          = 1000.0;
 q             = 1.0;
 gain          = 0.0;
 A             = 1.0;
 twoStages     = true;
 sampleRate    = 44100.0;
 sampleRateRec = 1.0 / sampleRate;

	resetBuffers();          // reset the filters memory buffers (to 0)
	calcCoeffs();            // calculate coefficients from default specification
}

MagicCarpetFilter::~MagicCarpetFilter()
{}

//----------------------------------------------------------------------------
//parameter settings:
void MagicCarpetFilter::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.001 )
  sampleRate = newSampleRate;

 sampleRateRec = 1.0 / sampleRate;
 calcCoeffs();
}

void MagicCarpetFilter::setMode(int newMode)
{
 mode = newMode;
 calcCoeffs();
}

void MagicCarpetFilter::setTwoStages(bool newTwoStagesSwitch)
{
 twoStages = newTwoStagesSwitch;
}

//----------------------------------------------------------------------------
//others:
void MagicCarpetFilter::resetBuffers()
{
 x_s1_d2_L = 0.0;
 x_s1_d1_L = 0.0;
 y_s1_d2_L = 0.0;
 y_s1_d1_L = 0.0;

 x_s1_d2_R = 0.0;
 x_s1_d1_R = 0.0;
 y_s1_d2_R = 0.0;
 y_s1_d1_R = 0.0;

 //x_s2_d2_L = 0.0;
 //x_s2_d1_L = 0.0;
 y_s2_d2_L = 0.0;
 y_s2_d1_L = 0.0;

 //x_s2_d2_R = 0.0;
 //x_s2_d1_R = 0.0;
 y_s2_d2_R = 0.0;
 y_s2_d1_R = 0.0;
}



