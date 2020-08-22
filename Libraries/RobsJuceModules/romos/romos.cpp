#ifdef ROMOS_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using.
*/
#error "Incorrect use of JUCE cpp file"
#endif


//#ifdef _MSC_VER

#include "romos.h"

using namespace rosic;  // get rid

//romos::PostExitMemLeakChecker romos::postExitMemLeakChecker;

// just include all files in alphabetical order - it's impractical to try to include them in
// dependency order due to a messy dependency network

namespace romos
{
#include "Algorithms/romos_FilterDesign.cpp"
#include "Algorithms/romos_Interpolation.cpp"

#include "Framework/romos_AudioConnection.cpp"
#include "Framework/romos_Module.cpp"
#include "Framework/romos_AtomicModule.cpp"
#include "Framework/romos_ContainerModule.cpp"
#include "Framework/romos_ModuleFactory.cpp"
#include "Framework/romos_NoteEvent.cpp"
#include "Framework/romos_ProcessingStatus.cpp"
#include "Framework/romos_TopLevelModule.cpp"
#include "Framework/romos_VoiceAllocator.cpp"
#include "Framework/romos_WorkArea.cpp"

#include "Modules/romos_ArithmeticModules.cpp"
#include "Modules/romos_DelayModules.cpp"
#include "Modules/romos_FilterModules.cpp"
#include "Modules/romos_FunctionModules.cpp"
#include "Modules/romos_FormulaModules.cpp"
#include "Modules/romos_InfrastructuralModules.cpp"
#include "Modules/romos_ModulationModules.cpp"
#include "Modules/romos_SoundGeneratorModules.cpp"

#include "Framework/romos_Liberty.cpp"
}

using namespace romos; 
// get rid - it's still needed in the code below

// todo: wrap the includes below alos into a namespace and get rid of the namespace declarations 
// in all the files:

#include "TestSuite/AutomaticTests.cpp"
#include "TestSuite/InteractiveTests.cpp"
#include "TestSuite/romos_ConcreteModularSystemTests.cpp"
#include "TestSuite/romos_ConcretePerformanceTests.cpp"
#include "TestSuite/romos_ConcreteProcessingTests.cpp"
#include "TestSuite/romos_ContainerManipulationTests.cpp"
#include "TestSuite/romos_GenerateDesiredOutput.cpp"
#include "TestSuite/romos_GlobalFrameworkTests.cpp"
#include "TestSuite/romos_InteractivePlottingTests.cpp"
#include "TestSuite/romos_InteractiveTestRunner.cpp"
#include "TestSuite/romos_ModularSystemTest.cpp"
#include "TestSuite/romos_ModuleBuildCodeGenerator.cpp"
#include "TestSuite/romos_PerformanceTest.cpp"
#include "TestSuite/romos_PerformanceTestRunner.cpp"
#include "TestSuite/romos_ProcessingTest.cpp"
#include "TestSuite/romos_TestEventGenerator.cpp"
#include "TestSuite/romos_TestModuleBuilder.cpp"
#include "TestSuite/romos_UnitTest.cpp"
#include "TestSuite/MiscTests.cpp"
#include "TestSuite/romos_UnitTestRunner.cpp"
#include "TestSuite/TestFilter.cpp"    // these two do not yet compile - maybe they use old
#include "TestSuite/TestHelpers.cpp"   // functions that have been removed
#include "TestSuite/TestsMain.cpp"
