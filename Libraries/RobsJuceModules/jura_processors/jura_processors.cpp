#ifdef JURA_PROCESSORS_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_processors.h"

using namespace rosic;

namespace jura
{

#include "tools/jura_ClassConversions.cpp"
#include "tools/jura_TuningFileManager.cpp"

#include "custom_widgets/jura_CustomComboBoxes.cpp"
#include "custom_widgets/jura_CustomSliders.cpp"
#include "custom_widgets/jura_EffectSelectionPopup.cpp"
#include "custom_widgets/jura_MeteringDisplay.cpp"

#include "basics/jura_RoutingMatrix.cpp"
#include "basics/jura_WaveTable.cpp"

#include "analyzers/jura_MidiMessageFilter.cpp"  // maybe move to basics
#include "analyzers/jura_MidiMonitor.cpp"
#include "analyzers/jura_TrackMeter.cpp"
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
#include "modulators/jura_VariousModulators.cpp"

#include "generators/jura_VariousGenerators.cpp"
#include "generators/jura_Oscillator3D.cpp"
#include "generators/jura_OscillatorStereo.cpp"
#include "generators/jura_FourOscSection.cpp"
#include "generators/jura_SamplePlayer.cpp"
#include "generators/jura_VectorMixer.cpp"         // maybe move to basics
#include "generators/jura_VectorSamplePlayer.cpp"
#include "generators/jura_RayBouncer.cpp"
#include "generators/jura_Snowflake.cpp"

#include "effects/jura_Enveloper.cpp"
#include "effects/jura_FuncShaper.cpp"
#include "effects/jura_NodeShaper.cpp"
#include "effects/jura_AlgoVerb.cpp"
#include "effects/jura_EchoLab.cpp"
#include "effects/jura_PitchShifter.cpp"
#include "effects/jura_CombStereoizer.cpp"
#include "effects/jura_StereoPanPlotEditor.cpp"
#include "effects/jura_VariousModules.cpp"
#include "effects/jura_Quadrifex.cpp"
#include "effects/jura_Moduluxury.cpp"
#include "effects/jura_StereoDelay.cpp"
#include "effects/jura_MultiBandEffect.cpp"

#include "instruments/jura_AcidSequencer.cpp"
#include "instruments/jura_AciDevil.cpp"
#include "instruments/jura_PolyphonicInstrument.cpp"
#include "instruments/jura_Straightliner.cpp"
#include "instruments/jura_SimpleSampler.cpp"
#include "instruments/jura_KeyShot.cpp"
#include "instruments/jura_MagicCarpet.cpp"
#include "instruments/jura_Quadrigen.cpp" 
#include "instruments/jura_Quadriga.cpp"
#include "instruments/jura_Workhorse.cpp"
#ifdef _MSC_VER
#include "instruments/jura_Liberty.cpp" // Liberty not yet available on gcc 
#include "instruments/jura_LibertyModules.cpp"
#endif

#include "unfinished/jura_ModalSynth.cpp"

#include "misc/jura_ChannelMatrix2x2.cpp"
#include "misc/jura_DspWorkbench.cpp"
#include "misc/jura_DebugAudioModule.cpp"
#include "misc/jura_AudioModuleSlotArray.cpp"
#include "misc/jura_ToolChain.cpp"

//#include "unfinished/jura_QuadSource.cpp"
//#include "unfinished/jura_DualFilter.cpp"
//#include "unfinished/jura_NewSynth.cpp"

#include "unfinished/jura_OscArrays.cpp"

}
