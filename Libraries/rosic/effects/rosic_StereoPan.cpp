#include "rosic_StereoPan.h"
using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

StereoPan::StereoPan()
{
  g                 = 1.0;
  p                 = 0.0;
  panLaw            = SINCOS;
  normalizeAtCenter = true;
  calculateGainFactors();
}

StereoPan::~StereoPan()
{

}

//-------------------------------------------------------------------------------------------------
// internal functions:

void StereoPan::calculateGainFactors()
{
  switch( panLaw )
  {
    // pan-laws without cross-mixing:
  case LINEAR:                   
    panLawLinear(p, gLL, gRL, gLR, gRR, normalizeAtCenter);                         break;
  case SINCOS:                   
    panLawSinCos(p, gLL, gRL, gLR, gRR, normalizeAtCenter);                         break;
  case SQUARE_ROOT:              
    panLawSquareRoot(p, gLL, gRL, gLR, gRR, normalizeAtCenter);                     break;
  case LINEAR_CLIPPED:           
    panLawLinearClipped(p, gLL, gRL, gLR, gRR, normalizeAtCenter);                  break;
  case LINEAR_SQUARE_NORMALIZED: 
    panLawLinearSquareNormalized(p, gLL, gRL, gLR, gRR, normalizeAtCenter);         break;
  case LINEAR_CLIPPED_SQUARE_NORMALIZED: 
    panLawLinearClippedSquareNormalized(p, gLL, gRL, gLR, gRR, normalizeAtCenter);  break;

    // pan-laws with cross-mixing:
  case LINEAR_CROSSMIX_SQUARE_NORMALIZED: 
    panLawLinearCrossMixSquareNormalized(p, gLL, gRL, gLR, gRR, normalizeAtCenter); break;
  case SINCOS_CROSSMIX_SQUARE_NORMALIZED: 
    panLawSinCosCrossMixSquareNormalized(p, gLL, gRL, gLR, gRR, normalizeAtCenter); break;
  case SQRT_CROSSMIX_SQUARE_NORMALIZED: 
    panLawSqrtCrossMixSquareNormalized(p, gLL, gRL, gLR, gRR, normalizeAtCenter);   break;

  default: panLawSinCos(p, gLL, gRL, gLR, gRR, false);
  }

  gLL *= g; gRL *= g; gLR *= g; gRR *= g;
}

