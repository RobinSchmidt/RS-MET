/** This cpp file includes all templated code of the rapt library. It is included, for example, by 
rosic/basics/TemplateInstantiations.cpp where the templates get instantiated and compiled as part 
of the rosic module. The rapt module itself has very little code to compile - those bits that 
actually are compiled into rapt.obj are in rapt.cpp. You should make sure, that any template 
instantiation gets compiled onlyinto one compilation unit, otherwise you may get "symbol already 
defined" linker errors. */

// todo: maybe move this file inot a subdirectory - otherwise the projucer compiles it - apparently
// it compiles all top-level .cpp files for a module?

#include "rapt.h"

#include "Basics/Basics.cpp"
#include "Data/Data.cpp"
#include "Math/Math.cpp"
#include "AudioBasics/AudioBasics.cpp" 
#include "Filters/Filters.cpp"
#include "Analysis/Analysis.cpp"
#include "Physics/Physics.cpp"
//#include "Circuits/Circuits.cpp"
//#include "Spectral/Spectral.cpp"
#include "Visualization/Visualization.cpp"
#include "Generators/Generators.cpp"
#include "Modulators/Modulators.cpp"
//#include "Effects/Effects.cpp"
//#include "Framework/Framework.cpp".
//#include "Music/Music.cpp"
//#include "Instruments/Instruments.cpp"
#include "Unfinished/Unfinished.cpp" 
#include "Spectral/Spectral.cpp"


//void fixLinkerWarningDummyFunction() {} // dummy function to fix linker warning (duh!)
                                          // nope - it actually adds a 2nd warning
/*

// We request some explicit instantiations here - later, when we add modules to the jura framework
// which use these classes, they may be deleted. At the moment, they are needed for Elan's
// Chaosfly but are nowhere instantiatied within jura. It's not a very elegant solution, but it's
// supposed to be temporary anyway:

//template double RAPT::rsLinToLin(double x, double inMin, double inMax, double outMin, double outMax);

template class RAPT::rsMatrixView<double>;
template class RAPT::rsMatrixOld<double>;
template class RAPT::rsNodeBasedFunction<double>;
template class RAPT::rsParametricBellFunction<double>;
template class RAPT::rsPositiveBellFunctions<double>;
template class RAPT::rsNormalizedSigmoids<double>;
template class RAPT::rsSinCosTable<double>;
template class RAPT::rsScaledAndShiftedSigmoid<double>;
template class RAPT::rsEllipse<double>;
template class RAPT::rsRotationXYZ<double>;
template class RAPT::rsCoordinateMapper<double>;
template class RAPT::rsCoordinateMapper2D<double>;


template class RAPT::rsStateVariableFilter<double, double>;

template class RAPT::rsAlphaMask<float>;
template class RAPT::rsImagePainter<float, float, float>;

// for PhaseScope:
template class RAPT::rsAlphaMask<double>;
template class RAPT::rsScopeScreenScanner<float>;  // do we need this?
template class RAPT::rsScopeScreenScanner<double>;
template class RAPT::rsPhaseScopeBuffer<double, float, double>;
template class RAPT::rsPhaseScopeBuffer2<double, float, double>;

// needed for the release build of Chaosfly on Linux - without them, apparently the compiler
// generates the classes only partially - some member functions are missing probably because they
// not called from anywhere inside jura or rosic:
template double RAPT::rsAbs(double x);
template class RAPT::rsBreakpointModulator<double>;
template class RAPT::rsSmoothingFilter<double, double>;
template class RAPT::rsLadderFilter<double, double>;
template class RAPT::rsLadderFilter<rsFloat64x2, double>;
//template class RAPT::rsLadderFilter<std::complex<double>, double>; // needed for TestPluginJUCE
template class RAPT::rsPhasorFilter<double, double>;
template class RAPT::rsPhasorStateMapper<double>;
// todo: get rid of directly using rapt classes in jura and/or products - create instantiations for
// double in rosic and use these instantiations only

template class RAPT::rsRayBouncer<double>;
template class RAPT::rsRayBouncerDriver<double>;
template class RAPT::rsLissajousOscillator3D<double>;

template class RAPT::rsMultiBandSplitter<double, double>;

// hmm...it seems, we need all these explicit instantiations anyway - maybe clean up the build
// system...rename the files to raptJuceModule.h/cpp
*/

