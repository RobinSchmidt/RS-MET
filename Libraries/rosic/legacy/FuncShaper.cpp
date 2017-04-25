#include "FuncShaper.h"

//----------------------------------------------------------------------------
// construction/destruction:

FuncShaper::FuncShaper()
{
 // init member variables:
 sampleRate = 44100.0;
 a          = 0.0;
 b          = 0.0;
 c          = 0.0;
 d          = 0.0;
 drive      = 1.0;
 dcOffset   = 0.0;
 outVol     = 1.0;
 dryVol     = 0.0;
 wetVol     = 1.0;

 upsamplerL.setOversamplingFactor(4);
 upsamplerR.setOversamplingFactor(4);

 antiAliasFilterL.setSampleRate(4*sampleRate);
 antiAliasFilterR.setSampleRate(4*sampleRate);

 // init embedded objects:
 distortionCurve.setRange(-1.5, 1.5);
}

FuncShaper::~FuncShaper()
{

}

//----------------------------------------------------------------------------
// parameter settings:

void FuncShaper::setSampleRate(double newSampleRate)
{
 AudioModule::setSampleRate(newSampleRate);

 inputFilterL.setSampleRate(sampleRate);
 inputFilterR.setSampleRate(sampleRate);

 antiAliasFilterL.setSampleRate(4*sampleRate);
 antiAliasFilterR.setSampleRate(4*sampleRate);

 outputFilterL.setSampleRate(sampleRate);
 outputFilterR.setSampleRate(sampleRate);

 return;
}

bool FuncShaper::setFunctionString(char *newFunctionString)
{
 // set the new string in the TabulatedFunction object and return the
 // boolean result of this set-function to the calling function:
 return distortionCurve.setFunctionString(newFunctionString);
}

void FuncShaper::setA(double newA)
{
 a = newA;
 distortionCurve.assignVariable("a", a);
}

void FuncShaper::setB(double newB)
{
 b = newB;
 distortionCurve.assignVariable("b", b);
}

void FuncShaper::setC(double newC)
{
 c = newC;
 distortionCurve.assignVariable("c", c);
}

void FuncShaper::setD(double newD)
{
 d = newD;
 distortionCurve.assignVariable("d", d);
}

void FuncShaper::setDrive(double newDrive)
{
 drive = dB2amp(newDrive);
}

void FuncShaper::setInLpfCutoff(double newInLpfCutoff)
{
 inputFilterL.setLpfCutoff(newInLpfCutoff);
 inputFilterR.setLpfCutoff(newInLpfCutoff);
}

void FuncShaper::setInHpfCutoff(double newInHpfCutoff)
{
 inputFilterL.setHpfCutoff(newInHpfCutoff);
 inputFilterR.setHpfCutoff(newInHpfCutoff);
}

void FuncShaper::setDcOffset(double newDcOffset)
{
 dcOffset = newDcOffset;
}

void FuncShaper::setOutLpfCutoff(double newOutLpfCutoff)
{
 outputFilterL.setLpfCutoff(newOutLpfCutoff);
 outputFilterR.setLpfCutoff(newOutLpfCutoff);
}

void FuncShaper::setOutHpfCutoff(double newOutHpfCutoff)
{
 outputFilterL.setHpfCutoff(newOutHpfCutoff);
 outputFilterR.setHpfCutoff(newOutHpfCutoff);
}

void FuncShaper::setOutVol(double newOutVol)
{
 outVol = dB2amp(newOutVol);
}

void FuncShaper::setDryWet(double newDryWet)
{
 if( (newDryWet>=0.0) && (newDryWet<=100.0) )
 {
  wetVol = 0.01 * newDryWet;
  dryVol = 1 - wetVol;
 } 
}