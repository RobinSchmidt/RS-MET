#ifdef ROMOS_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using.
*/
#error "Incorrect use of JUCE cpp file"
#endif

#include "romos.h"

namespace romos
{

#include "Framework/romos_ModuleTypeRegistry.cpp"


//#include "Framework/romos_AudioConnection.cpp"

//#include "Framework/romos_Module.cpp"

//#include "Framework/romos_ModuleAtomic.cpp"

//#include "Framework/romos_ModuleFactory.cpp"
//#include "Framework/romos_ModuleContainer.cpp" 
//#include "Framework/romos_TopLevelModule.cpp"


}