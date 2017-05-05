#ifdef JURA_PROCESSORS_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_processors.h"

// obsolete - moved to rosic module:
//#include "../../RAPT/Source/Modules/RAPT.cpp"
//// i think, in order to safely avoid "multiple definition" linker errors, this here needs to be the
//// one and only place where RAPT.cpp gets included
//
//// We request some explicit instantiations here - later, when we add modules to the jura framework
//// which use these classes, they may be deleted. At the moment, they are needed for Elan's
//// Chaosfly but are nowhere instantiatied within jura. It's not a very elegant solution, but it's
//// supposed to be temporary anyway:
//
//template class RAPT::rsParametricBellFunction<double>;
//template class RAPT::rsPositiveBellFunctions<double>;
//template class RAPT::NormalizedSigmoids<double>;
//template class RAPT::ScaledAndShiftedSigmoid<double>;
//
//template class RAPT::StateVariableFilter<double, double>;
//
//template class RAPT::AlphaMask<float>;
//template class RAPT::ImagePainter<float, float, float>;
//
//// for Elan's PrettyScope (may be irrelevant now):
//template class RAPT::AlphaMask<double>;
//template class RAPT::PhaseScopeBuffer2<double, float, double>;
//
//
//// needed for the release build of ChaosFly on Linux - withou them, apparently the compiler
//// generates the classes only partially - some member functions are missing probably they called
//// from nowhere inside JURA:
//template double RAPT::rsAbs(double x);
//template class RAPT::rsBreakpointModulator<double>;
//template class RAPT::LadderFilter<double, double>;
//// ..i really should copy over the rosic code and only use actual rosic classes in products (not
//// class templates) - these classes can themselves, one-by-one, be made instantiations of
//// templates -  i propagate up the code intio RAPT - which them becomes
//// Rob's Audio Processing Templates


namespace jura
{

#include "analyzers/jura_PhaseScope.cpp"
//#include "analyzers/jura_PhaseScopeMultiColor.cpp"

#include "filters/jura_LadderFilter.cpp"
#include "filters/jura_PhasorFilter.cpp"
#include "filters/jura_EngineersFilter.cpp"

#include "modulators/jura_BreakpointModulatorAudioModule.cpp"
#include "modulators/jura_ModulatorCurveEditor.cpp"
#include "modulators/jura_BreakpointModulatorEditor.cpp"
#include "modulators/jura_ModulatorCurveEditorMulti.cpp"
#include "modulators/jura_BreakpointModulatorEditorMulti.cpp"

#include "effects/jura_Enveloper.cpp"
#include "effects/jura_FuncShaper.cpp"
#include "effects/jura_AlgoVerb.cpp"

#include "instruments/jura_AcidSequencer.cpp"
#include "instruments/jura_AciDevil.cpp"

#include "misc/jura_AudioModuleChain.cpp"

}
