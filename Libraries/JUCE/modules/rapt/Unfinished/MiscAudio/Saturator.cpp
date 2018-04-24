namespace RSLib
{

//-------------------------------------------------------------------------------------------------
// class rsHalfWaveSaturator:

// Construction/Destruction:

rsHalfWaveSaturator::rsHalfWaveSaturator()
{
  setThreshold(0.0); 
  setSaturationFunction(&rsTanh);
}

// Setup:

void rsHalfWaveSaturator::setThreshold(double newThreshold)
{
  rsAssert(newThreshold >= 0.0);
  rsAssert(newThreshold <= 1.0);
  thresh   = newThreshold;
  outScale = 1-thresh;
  if(outScale >= RS_EPS(double))
    inScale = 1/outScale;
  else
    inScale = 1.0;
}

void rsHalfWaveSaturator::setSaturationFunction(double (*newFunction)(double))
{
  saturate = newFunction;
}

//-------------------------------------------------------------------------------------------------
// class rsSaturator:

// Setup:

void rsSaturator::setLowerThreshold(double newThreshold)
{
  loSaturator.setThreshold(newThreshold);
}

void rsSaturator::setUpperThreshold(double newThreshold)
{
  upSaturator.setThreshold(newThreshold);
}

void rsSaturator::setThresholds(double newThreshold)
{
  setUpperThreshold(newThreshold);
  setLowerThreshold(newThreshold);
}

void rsSaturator::setLowerSaturationFunction(double (*newFunction)(double))
{
  loSaturator.setSaturationFunction(newFunction);
}

void rsSaturator::setUpperSaturationFunction(double (*newFunction)(double))
{
  upSaturator.setSaturationFunction(newFunction);
}

void rsSaturator::setSaturationFunctions(double (*newFunction)(double))
{
  setLowerSaturationFunction(newFunction);
  setUpperSaturationFunction(newFunction);
}

//-------------------------------------------------------------------------------------------------
// class rsSidechainSaturator:

void rsSidechainSaturator::setMode(int newMode)
{
  mode = newMode;
}

}