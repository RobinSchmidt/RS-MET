/** We create simplified typenames for some explicit instantiations of the RAPT template classes. 
We use the following convention: we append one or more letters to template class name that 
indicate the template parameters. For example, PhaseScopeBuffer<float, float, double> becomes: 
PhaseScopeBufferFFD. For each typedef here, we also need the corresponding explicit 
instantiation in RaptInstantiations.cpp, otherwise linker errors will occur when the defined types 
are actually used somewhere. On the other hand, we don't necessarily need a typedef for every 
explicit instantiation in RaptInstatiations.cpp, but then we have to always write out the full 
template: we would have to write PhaseScopeBuffer<float, float, double> everywhere, if the 
PhaseScopeBufferFFD typedef wouldn't exist. */

#ifndef RAPT_INSTANTIATIONS_H
#define RAPT_INSTANTIATIONS_H

#include "../../../../../Libraries/JUCE/modules/rapt/rapt.h"
using namespace RAPT;

// Math:
typedef RAPT::rsConicSection<float> rsConicSectionF;
typedef RAPT::rsEllipse<float> rsEllipseF;
typedef RAPT::rsSinCosTable<float> rsSinCosTableF;

// Filters:
typedef RAPT::LadderFilter<float, float> LadderFilterFF;
typedef RAPT::PhasorFilter<float, float> PhasorFilterFF;
typedef RAPT::PhasorStateMapper<float> PhasorStateMapperF;
typedef RAPT::StateVariableFilter<float, float> StateVariableFilterFF; 

// Generators:
typedef RAPT::rsRayBouncer<float> rsRayBouncerF;

// Modulation:
typedef RAPT::rsBreakpointModulator<float> rsBreakpointModulatorF;

// Visualization:
typedef RAPT::Image<float> ImageF;
typedef RAPT::AlphaMask<float> AlphaMaskF;
typedef RAPT::ImagePainter<float, float, float> ImagePainterFFF;
typedef RAPT::ImageDrawer<float, float, float> ImageDrawerFFF;
typedef RAPT::LineDrawer<float, float, float> LineDrawerFFF;
typedef RAPT::PhaseScopeBuffer<float, float, double> PhaseScopeBufferFFD;

#endif