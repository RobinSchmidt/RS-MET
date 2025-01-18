#include "Filters.h"

namespace RAPT
{

#include "General/MovingWindowFilters.cpp"
#include "General/OnePoleFilter.cpp"
#include "General/SmoothingFilter.cpp"
#include "General/Interpolator.cpp" 
#include "General/SparseFilter.cpp" 

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
#include "Scientific/EllipticSubBandFilter.cpp"
#include "Scientific/EllipticSubBandFilterDirectForm.cpp"
#include "Scientific/QuadratureNetwork.cpp"
#include "Scientific/FilterSpecifications.cpp"
#include "Scientific/QuantileFilter.cpp"
#include "Scientific/WindowedFilterDesigner.cpp"
#include "Scientific/Convolver.cpp"
#include "Scientific/HilbertFilter.cpp"

#include "Musical/LadderFilter.cpp"
#include "Musical/PhasorFilter.cpp"
#include "Musical/StateVariableFilterOld.cpp"
#include "Musical/StateVariableFilter.cpp"
#include "Musical/DelayLine.cpp"
#include "Musical/Allpasses.cpp"

}