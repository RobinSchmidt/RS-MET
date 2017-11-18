#ifndef RAPT_FILTERS_H_INCLUDED
#define RAPT_FILTERS_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Filters or Filter

//#include "Basic/OnePoleFilter.h"
#include "Basic/SmoothingFilter.h"
//#include "Basic/OnePoleOneZero.h"
//#include "Basic/TwoPoleOneZero.h"  // decaying sine, etc.
//#include "Basic/Biquad.h"
//#include "Basic/BiquadChain.h"
//...

#include "Musical/LadderFilter.h"
#include "Musical/PhasorFilter.h"
#include "Musical/StateVariableFilter.h"
//include "Musical/AttackDecaySineFilter.h" // maybe get rid of writing "Filter" all the time

#include "Scientific/PrototypeDesigner.h"               // unit cutoff analog prototypes
#include "Scientific/PoleZeroMapper.h"                  // LP->LP, LP->HP, bilinear S->Z, etc.
#include "Scientific/FilterCoefficientConverter.h"      // pole-zero-to-biquad, etc
#include "Scientific/InfiniteImpulseResponseDesigner.h" // coeffs for Butterworth, elliptic,..
#include "Scientific/FilterAnalyzer.h"                  // transfer-function, frequency response,..
#include "Scientific/BiquadCascade.h" // rename to BiquadChain


//#include "Convolution/WindowedSinc.h" 

}

#endif