#include "rosic_Harmonics.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Harmonics::Harmonics() 
{
  driveFactor   = 1.0;
  dry           = 0.5;
  wet           = 0.5;
  chebychevMode = true;
  for(int i=0; i<numShapers; i++)
  {
    gainFactors[i] = 0.0;
    gainIsZero[i]  = true;
    subbandFiltersL[i].setSubDivision(i+2);
    subbandFiltersR[i].setSubDivision(i+2);
  }
}

Harmonics::~Harmonics()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void Harmonics::setSampleRate(double newSampleRate)
{
  inFilterL.setSampleRate(newSampleRate);
  inFilterR.setSampleRate(newSampleRate);
  outFilterL.setSampleRate(newSampleRate);
  outFilterR.setSampleRate(newSampleRate);
}

void Harmonics::setHarmonicGainFactor(int harmonicNumber, double newGain)
{
  if( harmonicNumber >= 2 && harmonicNumber < numShapers+2 )
  {
    gainFactors[harmonicNumber-2] = newGain;
    if( newGain == 0.f )
      gainIsZero[harmonicNumber-2] = true;
    else
      gainIsZero[harmonicNumber-2] = false;
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void Harmonics::reset()
{
  inFilterL.reset();
  inFilterR.reset();
  outFilterL.reset();
  outFilterR.reset();
  for(int i=0; i<numShapers; i++)
  {
    subbandFiltersL[i].reset();
    subbandFiltersR[i].reset();
  }
}

