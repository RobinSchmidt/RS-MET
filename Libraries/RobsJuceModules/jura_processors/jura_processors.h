/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               jura_processors
  vendor:           RS-MET
  version:          0.0.1
  name:             JUCE wrappers for RAPT and rosic
  description:      JUCE wrappers and GUIs for RAPT's and rosic's DSP algorithms
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     juce_core, juce_audio_basics, juce_graphics, juce_gui_basics,
                    juce_audio_formats, juce_audio_processors, jura_framework,
					          rapt, rosic, romos
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


#ifndef JURA_PROCESSORS_H_INCLUDED
#define JURA_PROCESSORS_H_INCLUDED

//#include <juce_core/juce_core.h>
//#include <juce_audio_basics/juce_audio_basics.h>
//#include <juce_audio_formats/juce_audio_formats.h>
//#include <juce_audio_processors/juce_audio_processors.h>
//#include <juce_graphics/juce_graphics.h>
//#include <juce_gui_basics/juce_gui_basics.h>
#include <jura_framework/jura_framework.h>

//#include <rapt/rapt.h>

#ifdef _MSC_VER
#include <romos/romos.h>
#else
#include <rosic/rosic.h>  // included now by romos.h, but not for gcc bcs romos doesn't compile there
#endif


using namespace juce;
// do we actually need all these includes? - most of them are included by jura_framework already

//#include "../../RAPT/Source/Modules/RAPT.h" // get rid of that
//using namespace RAPT;

// disable warnings related to "inherits ... via dominance", todo: try to get rid of virtual
// inheritance and reactivate the warning
#if defined _MSC_VER
#pragma warning (disable : 4250)
#endif


namespace jura
{

#include "tools/jura_ClassConversions.h"
#include "tools/jura_TuningFileManager.h"

#include "custom_widgets/jura_CustomComboBoxes.h"
#include "custom_widgets/jura_CustomSliders.h"
#include "custom_widgets/jura_EffectSelectionPopup.h"
#include "custom_widgets/jura_MeteringDisplay.h"

#include "basics/jura_RoutingMatrix.h"
#include "basics/jura_WaveTable.h"

#include "analyzers/jura_MidiMessageFilter.h"  // maybe move to basics
#include "analyzers/jura_MidiMonitor.h"
#include "analyzers/jura_TrackMeter.h"
#include "analyzers/jura_PhaseScope.h"
//#include "analyzers/jura_PhaseScopeMultiColor.h"
#include "analyzers/jura_MultiAnalyzer.h"

#include "filters/jura_LadderFilter.h"
#include "filters/jura_PhasorFilter.h"
#include "filters/jura_EngineersFilter.h"
#include "filters/jura_CrossOver.h"
#include "filters/jura_Equalizer.h"
#include "filters/jura_MultiModeFilter.h"

#include "modulators/jura_BreakpointModulatorAudioModule.h"
#include "modulators/jura_ModulatorCurveEditor.h"
#include "modulators/jura_BreakpointModulatorEditor.h"
#include "modulators/jura_ModulatorCurveEditorMulti.h"
#include "modulators/jura_BreakpointModulatorEditorMulti.h"
#include "modulators/jura_LowFrequencyOscillator.h"
#include "modulators/jura_VariousModulators.h"

#include "generators/jura_VariousGenerators.h"
#include "generators/jura_Oscillator3D.h"
#include "generators/jura_OscillatorStereo.h"
#include "generators/jura_FourOscSection.h"
#include "generators/jura_SamplePlayer.h"
#include "generators/jura_VectorMixer.h"  // maybe move to basics
#include "generators/jura_VectorSamplePlayer.h"
#include "generators/jura_RayBouncer.h"
#include "generators/jura_Snowflake.h"
 
#include "effects/jura_Enveloper.h"
#include "effects/jura_FuncShaper.h"
#include "effects/jura_NodeShaper.h"
#include "effects/jura_AlgoVerb.h"
#include "effects/jura_EchoLab.h"
#include "effects/jura_PitchShifter.h"
#include "effects/jura_CombStereoizer.h"
#include "effects/jura_StereoPanPlotEditor.h"
#include "effects/jura_VariousModules.h"
#include "effects/jura_Quadrifex.h"
#include "effects/jura_Moduluxury.h"
#include "effects/jura_StereoDelay.h"
#include "effects/jura_MultiBandEffect.h"

#include "instruments/jura_AcidSequencer.h"  // maybe move to a "sequencers" folder someday
#include "instruments/jura_AciDevil.h"
#include "instruments/jura_PolyphonicInstrument.h"
#include "instruments/jura_Straightliner.h"
#include "instruments/jura_SimpleSampler.h"
#include "instruments/jura_KeyShot.h"
#include "instruments/jura_MagicCarpet.h"
#include "instruments/jura_Quadrigen.h"  // move to generators, requires to split relevant parts out of jura_VariousModules
#include "instruments/jura_Quadriga.h"
#include "instruments/jura_Workhorse.h"
#ifdef _MSC_VER
#include "instruments/jura_Liberty.h"    // romos currently doesn't compile on gcc
#include "instruments/jura_LibertyModules.h"
#endif

#include "unfinished/jura_ModalSynth.h"

#include "misc/jura_ChannelMatrix2x2.h"
#include "misc/jura_DspWorkbench.h"
#include "misc/jura_DebugAudioModule.h"
#include "misc/jura_AudioModuleSlotArray.h"
#include "misc/jura_ToolChain.h"

//#include "unfinished/jura_QuadSource.h"
//#include "unfinished/jura_DualFilter.h"
//#include "unfinished/jura_NewSynth.h"

#include "unfinished/jura_OscArrays.h"

}

#endif   // JURA_PROCESSORS_H_INCLUDED
