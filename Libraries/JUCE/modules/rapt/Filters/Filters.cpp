#include "Filters.h"

namespace RAPT
{

//#include "Basic/OnePoleFilter.cpp"
#include "Basic/SmoothingFilter.cpp"

#include "Musical/LadderFilter.cpp"
#include "Musical/PhasorFilter.cpp"
#include "Musical/StateVariableFilter.cpp"

#include "Scientific/PrototypeDesigner.cpp"
#include "Scientific/PoleZeroMapper.cpp" 
#include "Scientific/FilterCoefficientConverter.cpp"
#include "Scientific/InfiniteImpulseResponseDesigner.cpp"
#include "Scientific/FilterAnalyzer.cpp"
#include "Scientific/BiquadCascade.cpp"
#include "Scientific/EngineersFilter.cpp"
#include "Scientific/LinkwitzRileyCrossOver.cpp"
#include "Scientific/CrossOver4Way.cpp"
#include "Scientific/DirectFormFilter.cpp" 

}