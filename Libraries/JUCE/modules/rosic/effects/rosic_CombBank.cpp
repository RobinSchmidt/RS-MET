#include "rosic_CombBank.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

CombBank::CombBank()
{
  numActiveCombs = maxNumCombs;
  referencePitch = freqToPitch(110.0);
  detune         = 0.0;
  pan1           = 0.0;
  pan2           = 1.0;
  gain           = 1.0/sqrt((double)numActiveCombs);
  dry            = ONE_OVER_SQRT2;
  wet            = ONE_OVER_SQRT2;

  setStereoSpread(0.0);
  initPitchOffsets();
  reset();
  setReferenceFrequency(1000.0);
}

//-------------------------------------------------------------------------------------------------
// setup:

void CombBank::setSampleRate(double newSampleRate)
{
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setSampleRate(newSampleRate);
    combsR[i].setSampleRate(newSampleRate);
  }
}

void CombBank::setReferencePitch(double newPitch)
{
  referencePitch = newPitch;
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setFrequency( pitchToFreq(referencePitch + pitchOffsets[i] + 0.5*detune) );
    combsR[i].setFrequency( pitchToFreq(referencePitch + pitchOffsets[i] - 0.5*detune) );
  }
}

void CombBank::setDetune(double newDetune)
{
  detune = newDetune;
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setFrequency( pitchToFreq(referencePitch + pitchOffsets[i] + 0.5*detune) );
    combsR[i].setFrequency( pitchToFreq(referencePitch + pitchOffsets[i] - 0.5*detune) );
  }
}

void CombBank::initPitchOffsets()
{
  for(int i=0; i<maxNumCombs; i++)
    pitchOffsets[i] = (double) i;
}

void CombBank::setStereoSpread(double newStereoSpread)
{
  stereoSpread   = newStereoSpread;
  double panOdd  = 0.01*stereoSpread;
  double panEven = -panOdd;
  equalPowerGainFactors(panOdd,  &gainOddL,  &gainOddR,  -1.0, 1.0); 
  equalPowerGainFactors(panEven, &gainEvenL, &gainEvenR, -1.0, 1.0); 
}

void CombBank::setDecayTime(double newDecayTime)   
{ 
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setDecayTime(newDecayTime);
    combsR[i].setDecayTime(newDecayTime);
  }
}
    
void CombBank::setHighDecayScale(double newScale)
{ 
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setHighDecayScale(newScale); 
    combsR[i].setHighDecayScale(newScale); 
  }
}

void CombBank::setLowDecayScale(double newScale)
{ 
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setLowDecayScale(newScale); 
    combsR[i].setLowDecayScale(newScale); 
  }
}

void CombBank::setHighCrossoverFreq(double newFreq)    
{ 
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setHighCrossoverFreq(newFreq); 
    combsR[i].setHighCrossoverFreq(newFreq); 
  }
}

void CombBank::setLowCrossoverFreq(double newFreq)    
{ 
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setLowCrossoverFreq(newFreq); 
    combsR[i].setLowCrossoverFreq(newFreq); 
  }
}

void CombBank::setOddOnlyMode(bool shouldCreateOnlyOddHarmonics)
{
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setOddOnlyMode(shouldCreateOnlyOddHarmonics);
    combsR[i].setOddOnlyMode(shouldCreateOnlyOddHarmonics);
  }
}

void CombBank::setNegativePolarity(bool shouldBeNegative) 
{
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].setNegativePolarity(shouldBeNegative);
    combsR[i].setNegativePolarity(shouldBeNegative);
  }
}

//-------------------------------------------------------------------------------------------------
// others:

void CombBank::reset()
{
  for(int i=0; i<maxNumCombs; i++)
  {
    combsL[i].clearBuffer();
    combsR[i].clearBuffer();
  }
}

