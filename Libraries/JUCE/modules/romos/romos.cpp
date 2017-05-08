#ifdef ROMOS_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using.
*/
#error "Incorrect use of JUCE cpp file"
#endif


#ifdef _MSC_VER

#include "romos.h"

//namespace romos { // it'S still all wrapped into a namespace inside the included files

// just include all files in alphabetical order - it's impractical to try to include them in
// dependency order due to a messy dependency network

#include "Algorithms/romos_FilterDesign.cpp"
#include "Algorithms/romos_Interpolation.cpp"

#include "Framework/romos_AudioConnection.cpp"
#include "Framework/romos_Liberty.cpp"
#include "Framework/romos_Module.cpp"
#include "Framework/romos_ModuleAtomic.cpp"
#include "Framework/romos_ModuleContainer.cpp"
#include "Framework/romos_ModuleFactory.cpp"
#include "Framework/romos_ModuleTypeRegistry.cpp"
#include "Framework/romos_NoteEvent.cpp"
#include "Framework/romos_ProcessingStatus.cpp"
#include "Framework/romos_TopLevelModule.cpp"
#include "Framework/romos_VoiceAllocator.cpp"
#include "Framework/romos_WorkArea.cpp"

#include "Modules/romos_ArithmeticModules.cpp"
#include "Modules/romos_DelayModules.cpp"
#include "Modules/romos_FilterModules.cpp"
#include "Modules/romos_FunctionModules.cpp"
#include "Modules/romos_InfrastructuralModules.cpp"
#include "Modules/romos_ModulationModules.cpp"
#include "Modules/romos_SoundGeneratorModules.cpp"

#include "TestSuite/romos_TestModuleBuilder.cpp" // seems to be needed by ModuleFactory - get rid of this dependency

#endif






//#include "Framework/romos_ModuleTypeRegistry.cpp"
//#include "Framework/romos_ProcessingStatus.cpp"

//// if everything else is commented, up to here, we get no linker errors when building
//
//#include "Framework/romos_AudioConnection.cpp"
//#include "Framework/romos_VoiceAllocator.cpp"
//#include "Framework/romos_WorkArea.cpp"
//#include "Framework/romos_Module.cpp"
//#include "Framework/romos_ModuleContainer.cpp"
//#include "Framework/romos_ModuleAtomic.cpp"
//
//#include "Modules/romos_ArithmeticModules.cpp"
//#include "Modules/romos_DelayModules.cpp"
//#include "Modules/romos_FilterModules.cpp"
//#include "Modules/romos_FunctionModules.cpp"
//#include "Modules/romos_InfrastructuralModules.cpp"
//#include "Modules/romos_ModulationModules.cpp"
//#include "Modules/romos_SoundGeneratorModules.cpp"
//
//
//
//#include "TestSuite/romos_TestModuleBuilder.cpp"  // seems to need all the atomic modules
//#include "Framework/romos_ModuleFactory.cpp"      // need TestModuleBuilder



//#include "Framework/romos_TopLevelModule.cpp"


//}
