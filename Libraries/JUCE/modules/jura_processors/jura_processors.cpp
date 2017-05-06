#ifdef JURA_PROCESSORS_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_processors.h"

namespace jura
{

#include "analyzers/jura_PhaseScope.cpp"
//#include "analyzers/jura_PhaseScopeMultiColor.cpp"
#include "analyzers/jura_MultiAnalyzer.cpp"

#include "filters/jura_LadderFilter.cpp"
#include "filters/jura_PhasorFilter.cpp"
#include "filters/jura_EngineersFilter.cpp"
#include "filters/jura_CrossOver.cpp"
#include "filters/jura_Equalizer.cpp"
#include "filters/jura_MultiModeFilter.cpp"

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
#include "misc/jura_ChannelMatrix2x2.cpp"
#include "misc/jura_DspWorkbench.cpp"

}
