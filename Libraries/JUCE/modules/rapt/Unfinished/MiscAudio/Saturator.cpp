// class rsHalfWaveSaturator:
// Construction/Destruction:

template<class TSig, class TPar>
rsHalfWaveSaturator<TSig, TPar>::rsHalfWaveSaturator()
{
  setThreshold(0.0); 
  setSaturationFunction(&rsTanh);
}

// Setup:

template<class TSig, class TPar>
void rsHalfWaveSaturator<TSig, TPar>::setThreshold(TPar newThreshold)
{
  rsAssert(newThreshold >= 0.0);
  rsAssert(newThreshold <= 1.0);
  thresh   = newThreshold;
  outScale = 1-thresh;
  if(outScale >= RS_EPS(TPar))
    inScale = 1/outScale;
  else
    inScale = 1.0;
}

template<class TSig, class TPar>
void rsHalfWaveSaturator<TSig, TPar>::setSaturationFunction(TSig (*newFunction)(TSig))
{
  saturate = newFunction;
}

//-------------------------------------------------------------------------------------------------
// class rsSaturator:

// Setup:

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setLowerThreshold(TPar newThreshold)
{
  loSaturator.setThreshold(newThreshold);
}

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setUpperThreshold(TPar newThreshold)
{
  upSaturator.setThreshold(newThreshold);
}

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setThresholds(TPar newThreshold)
{
  setUpperThreshold(newThreshold);
  setLowerThreshold(newThreshold);
}

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setLowerSaturationFunction(TSig (*newFunction)(TSig))
{
  loSaturator.setSaturationFunction(newFunction);
}

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setUpperSaturationFunction(TSig (*newFunction)(TSig))
{
  upSaturator.setSaturationFunction(newFunction);
}

template<class TSig, class TPar>
void rsSaturator<TSig, TPar>::setSaturationFunctions(TSig (*newFunction)(TSig))
{
  setLowerSaturationFunction(newFunction);
  setUpperSaturationFunction(newFunction);
}

//-------------------------------------------------------------------------------------------------
// class rsSidechainSaturator:

template<class TSig, class TPar>
void rsSidechainSaturator<TSig, TPar>::setMode(int newMode)
{
  mode = newMode;
}
