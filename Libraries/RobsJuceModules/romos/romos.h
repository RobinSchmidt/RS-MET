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

//#ifdef _MSC_VER  // it currently doesn't compile on gcc

#include <rosic/rosic.h> // for dsp algorithms - maybe include it in the cpp file?
//#include <map>           // for std::map, used to save/restore module states

//#define RS_BUILD_OLD_MODULE_FACTORY    // for transition to new factory implementation

namespace romos
{
#include "Algorithms/romos_FilterDesign.h"
#include "Algorithms/romos_Interpolation.h"

#include "Framework/romos_AudioConnection.h"
#include "Framework/romos_ModuleFactory.h"
#include "Framework/romos_VoiceAllocator.h"
#include "Framework/romos_ProcessingStatus.h"
#include "Framework/romos_Module.h"
#include "Framework/romos_AtomicModule.h"
#include "Framework/romos_ContainerModule.h"
#include "Framework/romos_NoteEvent.h"
#include "Framework/romos_TopLevelModule.h"
#include "Framework/romos_WorkArea.h"

#include "Modules/romos_ModuleDefinitionMacros.h"
#include "Modules/romos_ArithmeticModules.h"
#include "Modules/romos_DelayModules.h"
#include "Modules/romos_FilterModules.h"
#include "Modules/romos_FunctionModules.h"
#include "Modules/romos_FormulaModules.h"
#include "Modules/romos_InfrastructuralModules.h"
#include "Modules/romos_ModulationModules.h"
#include "Modules/romos_SoundGeneratorModules.h"

#include "Framework/romos_Liberty.h"  // should not be in the framework folder - maybe top-level

}



// maybe these should go into a separate include file:
#include "TestSuite/romos_TestEventGenerator.h"
#include "TestSuite/romos_UnitTest.h"
#include "TestSuite/romos_ProcessingTest.h"
#include "TestSuite/romos_PerformanceTest.h"
#include "TestSuite/romos_GenerateDesiredOutput.h"
#include "TestSuite/romos_UnitTest.h"
#include "TestSuite/romos_ModularSystemTest.h"
#include "TestSuite/romos_ModuleBuildCodeGenerator.h"
#include "TestSuite/romos_TestModuleBuilder.h" // not needed anymore when new testsuite is complete
#include "TestSuite/MiscTests.h"
#include "TestSuite/romos_UnitTestRunner.h"
//#include "TestSuite/romos_ConcretePerformanceTests.h"
#include "TestSuite/romos_PerformanceTestRunner.h"
#include "TestSuite/romos_InteractiveTestRunner.h"
// maybe more?


#ifdef _MSC_VER

namespace romos
{
struct PostExitMemLeakChecker
{
  ~PostExitMemLeakChecker()
  {
    if(_CrtDumpMemoryLeaks() == 1)
    {
      std::cout << "\n\n!!! Memory leaks detected (post exit of main) !!! \n";
      getchar();
    }
  }
};
extern PostExitMemLeakChecker postExitMemLeakChecker;
}
// there is a way to prevent the debugger to trigger false memory leaks, explained here:
// http://www.cplusplus.com/forum/general/9614/ it says:
// So the way to get around the problem is this:
// And then instantiate MemCheck as the first variable in the block/function so that it gets
// destroyed last, after everything else.

// maybe this check should be moved up to rosic or even rapt so that it catches memleaks that
// occur there also

//#include "Framework/romos_Liberty.h"

#endif   // #ifdef _MSC_VER

#endif   // #ifndef ROMOS_H_INCLUDED







