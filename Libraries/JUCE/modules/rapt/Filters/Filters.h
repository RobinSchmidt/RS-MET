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

//files to come:
//#include "Scientific/Prototypes.h"
//#include "Scientific/PoleZeroMapping.h" // LP->LP, LP->HP, S->Z bilinear, etc.
//#include "Scientific/CoefficientConversion.h" // pole-zero-to-biquad, 
//#include "Scientific/WindowedSinc.h"
//...

}

#endif