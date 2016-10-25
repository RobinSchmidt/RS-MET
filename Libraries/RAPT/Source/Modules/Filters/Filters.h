#ifndef RAPT_FILTERS_H_INCLUDED
#define RAPT_FILTERS_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Filters or Filter

//files to come:
//#include "Basics/OnePole.h"
//#include "Basics/OnePoleOneZero.h"
//#include "Basics/TwoPoleOneZero.h"  // decaying sine, etc.
//#include "Basics/Biquad.h"
//#include "Basics/BiquadChain.h"
//...

#include "Musical/LadderFilter.h"
//include "Musical/StateVariableFilter.h"
//include "Musical/AttackDecaySineFilter.h" // maybe get rid of writing "Filter" all the time

//files to come:
//#include "Scientific/Prototypes.h"
//#include "Scientific/PoleZeroMapping.h" // LP->LP, LP->HP, S->Z bilinear, etc.
//#include "Scientific/CoefficientConversion.h" // pole-zero-to-biquad, 
//#include "Scientific/WindowedSinc.h"
//...

}

#endif