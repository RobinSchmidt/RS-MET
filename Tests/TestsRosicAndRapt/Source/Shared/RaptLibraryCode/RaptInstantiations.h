/** We create simplified typenames for some explicit instantiations of the RAPT template classes. 
We use the following convention: we append one or more letters to template class name that 
indicate the template parameters. For example, rsPhaseScopeBuffer<float, float, double> becomes: 
PhaseScopeBufferFFD. For each typedef here, we also need the corresponding explicit 
instantiation in RaptInstantiations.cpp, otherwise linker errors will occur when the defined types 
are actually used somewhere. On the other hand, we don't necessarily need a typedef for every 
explicit instantiation in RaptInstatiations.cpp, but then we have to always write out the full 
template: we would have to write rsPhaseScopeBuffer<float, float, double> everywhere, if the 
PhaseScopeBufferFFD typedef wouldn't exist. */

#ifndef RAPT_INSTANTIATIONS_H
#define RAPT_INSTANTIATIONS_H

#include "../../../../../Libraries/JUCE/modules/rapt/rapt.h"
using namespace RAPT;

// Math:
typedef RAPT::rsConicSection<float> rsConicSectionF;
typedef RAPT::rsEllipse<float> rsEllipseF;
typedef RAPT::rsSinCosTable<float> rsSinCosTableF;
typedef RAPT::rsSinCosTable<double> rsSinCosTableD;

typedef RAPT::rsVector3D<float> rsVector3DF;



// Filters:
typedef RAPT::rsSmoothingFilter<float, float> rsSmoothingFilterFF;
typedef RAPT::rsLadderFilter<float, float> rsLadderFilterFF;
typedef RAPT::rsPhasorFilter<float, float> rsPhasorFilterFF;
typedef RAPT::rsPhasorStateMapper<float> rsPhasorStateMapperF;
typedef RAPT::rsStateVariableFilter<float, float> rsStateVariableFilterFF; 

// Physics:
typedef RAPT::rsParticleSystem<float> rsParticleSystemF;

// Generators:
typedef RAPT::rsRayBouncer<float> rsRayBouncerF;
typedef RAPT::rsRayBouncer1D<float> rsRayBouncer1DF;

// Modulation:
typedef RAPT::rsBreakpointModulator<float> rsBreakpointModulatorF;

// Visualization:
typedef RAPT::rsImage<float> rsImageF;
typedef RAPT::rsAlphaMask<float> rsAlphaMaskF;
typedef RAPT::rsImagePainter<float, float, float> rsImagePainterFFF;
typedef RAPT::rsImageDrawer<float, float, float> rsImageDrawerFFF;
typedef RAPT::rsLineDrawer<float, float, float> rsLineDrawerFFF;
typedef RAPT::rsPhaseScopeBuffer<float, float, double> rsPhaseScopeBufferFFD;

#endif