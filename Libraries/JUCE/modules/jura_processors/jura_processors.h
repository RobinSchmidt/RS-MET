/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               jura_processors
  vendor:           RS-MET
  version:          0.0.1
  name:             JUCE wrappers for RAPT 
  description:      JUCE wrappers and GUIs for RAPT's audio signal processing algorithms
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     juce_core, juce_audio_basics, juce_graphics, juce_gui_basics, 
                    juce_audio_formats, juce_audio_processors, jura_framework, rosic
  OSXFrameworks:    
  iOSFrameworks:    

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


#ifndef JURA_PROCESSORS_H_INCLUDED
#define JURA_PROCESSORS_H_INCLUDED

#include <juce_core/juce_core.h> 
#include <juce_audio_basics/juce_audio_basics.h> 
#include <juce_audio_formats/juce_audio_formats.h> 
#include <juce_audio_processors/juce_audio_processors.h> 
#include <juce_graphics/juce_graphics.h> 
#include <juce_gui_basics/juce_gui_basics.h> 
#include <jura_framework/jura_framework.h> 
#include <rosic/rosic.h>
using namespace juce;
// do we actually need all these includes? - most of them are included by jura_framework already

//#include "../../RAPT/Source/Modules/RAPT.h" // get rid of that
//using namespace RAPT;

// disable warnings related to "inherits ... via dominance", todo: try to get rid of virtual 
// inheritance and reactivate the warning
#pragma warning (disable : 4250)  
                                  

namespace jura
{

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

#include "effects/jura_Enveloper.h"
#include "effects/jura_FuncShaper.h"
#include "effects/jura_AlgoVerb.h"

#include "instruments/jura_AcidSequencer.h"  // maybe move to a "sequencers" folder someday
#include "instruments/jura_AciDevil.h"

#include "misc/jura_AudioModuleChain.h"
#include "misc/jura_ChannelMatrix2x2.h"
#include "misc/jura_DspWorkbench.h"

}

#endif   // JURA_PROCESSORS_H_INCLUDED
