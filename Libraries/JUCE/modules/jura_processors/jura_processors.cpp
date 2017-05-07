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

#include "tools/jura_TuningFileManager.cpp"

#include "custom_widgets/jura_CustomComboBoxes.cpp"
#include "custom_widgets/jura_CustomSliders.cpp"
#include "custom_widgets/jura_EffectSelectionPopup.cpp"

#include "basics/jura_WaveTable.cpp"

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
#include "modulators/jura_LowFrequencyOscillator.cpp"

#include "generators/jura_OscillatorStereo.cpp"
#include "generators/jura_FourOscSection.cpp"
#include "generators/jura_SamplePlayer.cpp"
#include "generators/jura_VectorMixer.cpp"         // maybe move to basics
#include "generators/jura_VectorSamplePlayer.cpp"

#include "effects/jura_Enveloper.cpp"
#include "effects/jura_FuncShaper.cpp"
#include "effects/jura_AlgoVerb.cpp"
#include "effects/jura_EchoLab.cpp"
#include "effects/jura_PitchShifter.cpp"
#include "effects/jura_CombStereoizer.cpp"
#include "effects/jura_StereoPanPlotEditor.cpp"

#include "instruments/jura_AcidSequencer.cpp"
#include "instruments/jura_AciDevil.cpp"
#include "instruments/jura_PolyphonicInstrument.cpp"
#include "instruments/jura_Straightliner.cpp"
#include "instruments/jura_SimpleSampler.cpp"
#include "instruments/jura_KeyShot.cpp"

#include "misc/jura_AudioModuleChain.cpp"
#include "misc/jura_ChannelMatrix2x2.cpp"
#include "misc/jura_DspWorkbench.cpp"

}
