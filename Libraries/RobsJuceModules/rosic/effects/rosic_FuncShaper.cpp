//#include "rosic_FuncShaper.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

FuncShaper::FuncShaper()
{
  // init member variables:
  sampleRate      = 44100.0;
  a               = 1.0;
  b               = 0.0;
  c               = 0.0;
  d               = 0.0;
  driveFactor     = 1.0;
  dcOffset        = 0.0;
  outVolFactor    = 1.0;
  dryVol          = 0.0;
  wetVol          = 1.0;
  minA            = 0.0;
  maxA            = 1.0;
  minB            = 0.0;
  maxB            = 1.0;
  minC            = 0.0;
  maxC            = 1.0;
  minD            = 0.0;
  maxD            = 1.0;

  drive           = 0.0;
  outVol          = 0.0;
  dryWet          = 100.0;

  inFilterActive  = false;
  outFilterActive = false;

  numFadeSamples  = (int) (sampleRate * 50.0/1000.0);
  fadeCountDown   = 0;
  deClickingFilter.setMode(rsOnePoleFilterStereo::LOWPASS_IIT);
  deClickingFilter.setSampleRate(sampleRate);
  deClickingFilter.setCutoff(20.0);

  // init embedded objects:
  oversampling    = 4;
  upsamplerL.setSubDivision(oversampling);
  upsamplerR.setSubDivision(oversampling);
  antiAliasFilterL.setSubDivision(oversampling);
  antiAliasFilterR.setSubDivision(oversampling);

  distortionCurve.setRange(-1.5, 1.5);

  setFunctionString("tanh(a*x);", false);
  setA(1.0, false);
  calculateTable();
}

FuncShaper::~FuncShaper()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void FuncShaper::setSampleRate(double newSampleRate)
{
  if( newSampleRate > 0.0 )
    sampleRate = newSampleRate;

  inputFilterL.setSampleRate(sampleRate);
  inputFilterR.setSampleRate(sampleRate);

  outputFilterL.setSampleRate(sampleRate);
  outputFilterR.setSampleRate(sampleRate);

  numFadeSamples  = (int) (sampleRate * 50.0/1000.0);
  deClickingFilter.setSampleRate(sampleRate);

  return;
}

bool FuncShaper::setFunctionString(const char *newFunctionString, bool reCalculateTable)
{
  // set the new string in the TabulatedFunction object and return the
  // boolean result of this set-function to the calling function:
  bool success = distortionCurve.setFunctionString(newFunctionString, reCalculateTable);
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
  return success;
}

void FuncShaper::setA(double newA, bool reCalculateTable)
{
  a = newA;
  distortionCurve.assignVariable("a", a, reCalculateTable);
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
}

/*
void FuncShaper::setMinA(double newMinA)
{
  if( newMinA <= (maxA-0.01) )
    minA = newMinA;
}

void FuncShaper::setMaxA(double newMaxA)
{
  if( newMaxA >= (minA+0.01) )
    maxA = newMaxA;
}
*/

void FuncShaper::setB(double newB, bool reCalculateTable)
{
  b = newB;
  distortionCurve.assignVariable("b", b, reCalculateTable);
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
}

/*
void FuncShaper::setMinB(double newMinB)
{
  if( newMinB <= (maxB-0.01) )
    minB = newMinB;
}

void FuncShaper::setMaxB(double newMaxB)
{
  if( newMaxB >= (minB+0.01) )
    maxB = newMaxB;
}
*/

void FuncShaper::setC(double newC, bool reCalculateTable)
{
  c = newC;
  distortionCurve.assignVariable("c", c, reCalculateTable);
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
}

/*
void FuncShaper::setMinC(double newMinC)
{
  if( newMinC <= (maxC-0.01) )
    minC = newMinC;
}

void FuncShaper::setMaxC(double newMaxC)
{
  if( newMaxC >= (minC+0.01) )
    maxC = newMaxC;
}
*/

void FuncShaper::setD(double newD, bool reCalculateTable)
{
  d = newD;
  distortionCurve.assignVariable("d", d, reCalculateTable);
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
}

/*
void FuncShaper::setMinD(double newMinD)
{
  if( newMinD <= (maxD-0.01) )
    minD = newMinD;
}

void FuncShaper::setMaxD(double newMaxD)
{
  if( newMaxD >= (minD+0.01) )
    maxD = newMaxD;
}
*/

void FuncShaper::useInputFilter(bool shouldBeUsed)
{
  inFilterActive = shouldBeUsed;
  fadeCountDown  = numFadeSamples;
}

void FuncShaper::setInLowpassCutoff(double newInLowpassCutoff)
{
  inputFilterL.setLowpassCutoff(newInLowpassCutoff);
  inputFilterR.setLowpassCutoff(newInLowpassCutoff);
}

void FuncShaper::setInHighpassCutoff(double newInHighpassCutoff)
{
  inputFilterL.setHighpassCutoff(newInHighpassCutoff);
  inputFilterR.setHighpassCutoff(newInHighpassCutoff);
}

void FuncShaper::setDrive(double newDrive)
{
  drive         = newDrive;
  driveFactor   = RAPT::rsDbToAmp(drive);
  //fadeCountDown = numFadeSamples;
}

void FuncShaper::setDcOffset(double newDcOffset)
{
  dcOffset      = newDcOffset;
  //fadeCountDown = numFadeSamples;
}

void FuncShaper::setOversampling(int newOversamplingFactor)
{
  oversampling = max(1, newOversamplingFactor);
  if( oversampling > 1 )
  {
    upsamplerL.setSubDivision(oversampling);
    upsamplerR.setSubDivision(oversampling);
    antiAliasFilterL.setSubDivision(oversampling);
    antiAliasFilterR.setSubDivision(oversampling);
  }
}

void FuncShaper::useOutputFilter(bool shouldBeUsed)
{
  outFilterActive = shouldBeUsed;
  fadeCountDown   = numFadeSamples;
}

void FuncShaper::setOutLowpassCutoff(double newOutLowpassCutoff)
{
  outputFilterL.setLowpassCutoff(newOutLowpassCutoff);
  outputFilterR.setLowpassCutoff(newOutLowpassCutoff);
}

void FuncShaper::setOutHighpassCutoff(double newOutHighpassCutoff)
{
  outputFilterL.setHighpassCutoff(newOutHighpassCutoff);
  outputFilterR.setHighpassCutoff(newOutHighpassCutoff);
}

void FuncShaper::setOutVol(double newOutVol)
{
  outVol        = newOutVol;
  outVolFactor  = RAPT::rsDbToAmp(outVol);
  //fadeCountDown = numFadeSamples;
}

void FuncShaper::setDryWet(double newDryWet)
{
  if( (newDryWet >= 0.0) && (newDryWet <= 100.0) )
  {
    dryWet        = newDryWet;
    wetVol        = 0.01 * dryWet;
    dryVol        = 1.0 - wetVol;
    //fadeCountDown = numFadeSamples;
  } 
}

//-------------------------------------------------------------------------------------------------
// others:

void FuncShaper::calculateTable()
{
  distortionCurve.calculateTable();
  distortionCurve.clipTableValues(-4.0, 4.0);
  fadeCountDown = numFadeSamples;
}



/*=================================================================================================

ToDo:

-Let the user set up the range of the lookup-table. Maybe the resolution and interpolation method 
 also.
-Make the usage of the lookup table optional.
-Rationale: some distortion functions like (foldovers, i.e. sin, etc) may not play well with 
 tabulation...although, the different table range could be simulated by pre- and post scaling of 
 the signal, so maybe it's not really necessary introduce extra parameter for that

 See: 
 https://www.youtube.com/watch?v=oIChUOV_0w4  Ivan Cohen - Fifty shades of distortion (ADC'17)

*/
