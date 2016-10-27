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
                    juce_audio_formats, juce_audio_processors, jura_framework
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
using namespace juce;
// do we actually need all these includes? - most of them are included by jura_framework already

#include "../../RAPT/Source/Modules/RAPT.h"
using namespace RAPT;

namespace jura
{

#include "baseclasses/jura_AudioModule.h"
#include "baseclasses/jura_AudioPlugin.h"

#include "filters/jura_LadderFilter.h"

#include "modulators/jura_BreakpointModulatorAudioModule.h"
//#include "modulators/jura_BreakpointModulatorEditor.h"

}

#endif   // JURA_PROCESSORS_H_INCLUDED
