#ifndef RAPT_FILTERS_H_INCLUDED
#define RAPT_FILTERS_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Filters or Filter

#include "Basic/MovingWindowFilters.h"
#include "Basic/OnePoleFilter.h"
#include "Basic/SmoothingFilter.h"   // maybe move to a Tool
//#include "Basic/OnePoleOneZero.h"
//#include "Basic/TwoPoleOneZero.h"  // decaying sine, etc.
//#include "Basic/Biquad.h"
//#include "Basic/BiquadChain.h"
//#include "Basic/CookBookFilter.h"
//...

#include "Scientific/PrototypeDesigner.h"               // unit cutoff analog prototypes
#include "Scientific/PoleZeroMapper.h"                  // LP->LP, LP->HP, bilinear S->Z, etc.
#include "Scientific/FilterCoefficientConverter.h"      // pole-zero-to-biquad, etc
#include "Scientific/InfiniteImpulseResponseDesigner.h" // coeffs for Butterworth, elliptic,..
#include "Scientific/FilterAnalyzer.h"                  // transfer-function, frequency response,..
#include "Scientific/BiquadCascade.h"                   // rename to BiquadChain
#include "Scientific/EngineersFilter.h" 
#include "Scientific/LinkwitzRileyCrossOver.h" 
#include "Scientific/CrossOver4Way.h" 
#include "Scientific/DirectFormFilter.h"
#include "Scientific/EllipticSubBandFilter.h"
#include "Scientific/EllipticSubBandFilterDirectForm.h"
#include "Scientific/QuadratureNetwork.h"
#include "Scientific/FilterSpecifications.h"
#include "Scientific/QuantileFilter.h"

//#include "Convolution/WindowedSinc.h" 

#include "Musical/LadderFilter.h"
#include "Musical/PhasorFilter.h"
#include "Musical/StateVariableFilter.h"
//include "Musical/AttackDecaySineFilter.h" // maybe get rid of writing "Filter" all the time
//#include "Basic/Equalizer.h"

}

#endif