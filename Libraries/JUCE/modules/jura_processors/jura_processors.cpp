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
//using namespace RAPT;
// i think, in order to safely avoid "multiple definition" linker errors, this here needs to be the 
// one and only place where RAPT.cpp gets included

namespace jura
{

#include "baseclasses/jura_AudioModule.cpp"
#include "baseclasses/jura_AudioPlugin.cpp"

#include "filters/jura_LadderFilter.cpp"

#include "modulators/jura_BreakpointModulatorAudioModule.cpp"
#include "modulators/jura_ModulatorCurveEditor.cpp"
#include "modulators/jura_BreakpointModulatorEditor.cpp"
#include "modulators/jura_ModulatorCurveEditorMulti.cpp"
#include "modulators/jura_BreakpointModulatorEditorMulti.cpp"

#include "effects/jura_Enveloper.cpp"

}
