#ifdef JURA_BASICS_H_INCLUDED
 /* When you add this cpp file to your project, you mustn't include it in a file where you've
    already included any other headers - just put it inside a file on its own, possibly with your config
    flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
    header files that the compiler may be using.
 */
 #error "Incorrect use of JUCE cpp file"
#endif

#include "jura_basics.h"

//#include "../../RAPT/Code/Library/RAPT.h"
#include "../../RAPT/Code/Library/RAPT.cpp"
//using namespace RAPT;
// i think, in order to safely avoid "multiple definition" linker errors, this here needs to be the 
// one and only place where RAPT.cpp gets included

namespace jura
{

#include "control/jura_Parameter.cpp"
#include "control/jura_AutomatableParameter.cpp"
#include "control/jura_AutomatableModule.cpp"
// maybe this control stuff can be done in RAPT to make it available to non-JUCE based client code, 
// too



}
