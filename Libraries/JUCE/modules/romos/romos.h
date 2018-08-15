/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               romos
  vendor:           RS-MET
  version:          0.0.1
  name:             Rob's Modular Synthesizer
  description:      Modular software synthesizer
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     rosic
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


/** This is the main include file for the RoMoS (Rob's Modular Synthesizer) library. */

#ifndef ROMOS_H_INCLUDED
#define ROMOS_H_INCLUDED

#ifdef _MSC_VER  // it currently doesn't compile on gcc

#include <rosic/rosic.h> // for dsp algorithms


#include "Algorithms/romos_FilterDesign.h"
#include "Algorithms/romos_Interpolation.h"
#include "Framework/romos_AudioConnection.h"
#include "Framework/romos_ModuleTypeRegistry.h"
#include "Framework/romos_ProcessingStatus.h"
#include "Framework/romos_Module.h"
#include "Framework/romos_ModuleAtomic.h"
#include "Framework/romos_ModuleContainer.h"

#include "Framework/romos_ModuleFactory.h"
#include "Framework/romos_NoteEvent.h"
#include "Framework/romos_TopLevelModule.h"
#include "Framework/romos_VoiceAllocator.h"
#include "Framework/romos_WorkArea.h"

#include "Modules/romos_ArithmeticModules.h"
#include "Modules/romos_DelayModules.h"
#include "Modules/romos_FilterModules.h"
#include "Modules/romos_FunctionModules.h"
#include "Modules/romos_InfrastructuralModules.h"
#include "Modules/romos_ModulationModules.h"
#include "Modules/romos_SoundGeneratorModules.h"


#include "Framework/romos_Liberty.h"


#include "TestSuite/romos_ModuleBuildCodeGenerator.h"
#include "TestSuite/romos_TestModuleBuilder.h" // not needed anymore when new testsuite is complete
#include "TestSuite/romos_UnitTestRunner.h"
#include "TestSuite/romos_PerformanceTestRunner.h"
#include "TestSuite/romos_InteractiveTestRunner.h"
// maybe more?

#include "Framework/romos_Liberty.h"

#endif   // #ifdef _MSC_VER
#endif   // #ifndef ROMOS_H_INCLUDED







