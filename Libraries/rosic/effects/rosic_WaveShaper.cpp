#include "rosic_WaveShaper.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

WaveShaper::WaveShaper()
{
  transferFunctionIndex =   1;     // tanh by default
  oversampling          =   4;
  driveFactor           =   1.0;
  dcOffset              =   0.0;
  outVolFactor          =   1.0;
  amount                = 100.0;
  slope                 =   1.0;
  intercept             =   0.0;
  calculateQuinticCoefficients();  // to initialize the a-coefficients
  upsamplerL.setSubDivision(oversampling);
  upsamplerR.setSubDivision(oversampling);
  antiAliasFilterL.setSubDivision(oversampling);
  antiAliasFilterR.setSubDivision(oversampling);
  reset();
}

WaveShaper::~WaveShaper()
{

}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void WaveShaper::setTransferFunction(int newTransferFunction)
{
  if( newTransferFunction >= 0 )
    transferFunctionIndex = newTransferFunction;
}

void WaveShaper::setOversampling(int newOversamplingFactor)
{
  oversampling = rmax(1, newOversamplingFactor);
  if( oversampling > 1 )
  {
    upsamplerL.setSubDivision(oversampling);
    upsamplerR.setSubDivision(oversampling);
    antiAliasFilterL.setSubDivision(oversampling);
    antiAliasFilterR.setSubDivision(oversampling);
  }
}

//-------------------------------------------------------------------------------------------------
// others:
    
void WaveShaper::calculateQuinticCoefficients()
{
  a0 = intercept;
  a1 = slope; 
  a2 = -2*intercept;
  a3 = -0.5*(4*a1-5); 
  a4 = intercept;
  a5 = 0.5*(2*a1-3);
}

void WaveShaper::reset()
{
  upsamplerL.reset();
  upsamplerR.reset();
  antiAliasFilterL.reset();
  antiAliasFilterR.reset();
}




