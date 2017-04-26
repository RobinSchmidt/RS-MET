#include "rosic_CookbookFilter.h"
using namespace rosic;

//----------------------------------------------------------------------------
// construction/destruction
CookbookFilter::CookbookFilter()
{
 mode          = LOWPASS;
 freq          = 1000.0;
 oldFreq       = 1000.0;
 q             = 1.0;
 gain          = 0.0;
 A             = 1.0;
 numStages     = 2;
 sampleRate    = 44100.0;
 sampleRateRec = 1.0 / sampleRate;

	resetBuffers();          // reset the filters memory buffers (to 0)
	calcCoeffs();            // calculate coefficients from default specification
 convertDirectToLadder(); // init ladder-coeffs, too
}

CookbookFilter::~CookbookFilter()
{}

//----------------------------------------------------------------------------
// parameter settings:
void CookbookFilter::setSampleRate(double newSampleRate)
{
 if( newSampleRate > 0.0 )
  sampleRate = newSampleRate;

 sampleRateRec = 1.0 / sampleRate;
 calcCoeffs();
}

void CookbookFilter::setMode(int newMode)
{
 mode = newMode;
 calcCoeffs();
}

void CookbookFilter::setNumStages(int newNumStages)
{
 if( (newNumStages >=1 ) && (newNumStages <= maxNumStages) )
  numStages = newNumStages;
 resetBuffers();
}

//----------------------------------------------------------------------------
// others:
void CookbookFilter::resetBuffers()
{
 int i;
 for (i=0; i<maxNumStages; i++)
 {
  x2[i] = 0.0;
  x1[i] = 0.0;
  y2[i] = 0.0;
  y1[i] = 0.0;
  g2[i] = 0.0;
  g1[i] = 0.0;
 
  e2_0[i] = 0.0;
  e2_1[i] = 0.0;
  e1_0[i] = 0.0;
  e1_1[i] = 0.0;
  e0_0[i] = 0.0;
  e0_1[i] = 0.0;
  e0_2[i] = 0.0;
  e1t_0[i] = 0.0;
  e1t_1[i] = 0.0;
  e2t_0[i] = 0.0;
  e2t_1[i] = 0.0;
 }
}



