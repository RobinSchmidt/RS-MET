#ifdef ROMOS_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using.
*/
#error "Incorrect use of JUCE cpp file"
#endif

#include "romos.h"

//namespace romos { // it'S still all wrapped into a namespace inside the included files

#include "Framework/romos_ModuleTypeRegistry.cpp"
#include "Framework/romos_ProcessingStatus.cpp"

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