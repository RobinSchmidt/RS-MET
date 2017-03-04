#ifdef JURA_PROCESSORS_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_processors.h"

#include "../../RAPT/Source/Modules/RAPT.cpp"
// i think, in order to safely avoid "multiple definition" linker errors, this here needs to be the 
// one and only place where RAPT.cpp gets included

// We request some explicit instantiations here - later, when we add modules to the jura framework 
// which use these classes, they may be deleted. At the moment, they are needed for Elan's 
// Chaosfly but are nowhere instantiatied within jura. It's not a very elegant solution, but it's 
// supposed to be temporary anyway:

template RAPT::rsParametricBellFunction<double>;
template RAPT::rsPositiveBellFunctions<double>;
template RAPT::NormalizedSigmoids<double>;
template RAPT::ScaledAndShiftedSigmoid<double>;

template RAPT::StateVariableFilter<double, double>;

template RAPT::AlphaMask<float>;
template RAPT::ImagePainter<float, float, float>;

// for Elan's PrettyScope:
template RAPT::AlphaMask<double>;
template RAPT::PhaseScopeBuffer2<double, float, double>;


namespace jura
{

#include "baseclasses/jura_AudioModule.cpp"
#include "baseclasses/jura_AudioPlugin.cpp"

#include "analyzers/jura_PhaseScope.cpp"
#include "analyzers/jura_PhaseScopeMultiColor.cpp"

#include "filters/jura_LadderFilter.cpp"
#include "filters/jura_PhasorFilter.cpp"

#include "modulators/jura_BreakpointModulatorAudioModule.cpp"
#include "modulators/jura_ModulatorCurveEditor.cpp"
#include "modulators/jura_BreakpointModulatorEditor.cpp"
#include "modulators/jura_ModulatorCurveEditorMulti.cpp"
#include "modulators/jura_BreakpointModulatorEditorMulti.cpp"

#include "effects/jura_Enveloper.cpp"

#include "misc/jura_AudioModuleChain.cpp"

}
