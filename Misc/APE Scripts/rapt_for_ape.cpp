#ifndef RAPT_FOR_APE_CPP_INCLUDED
#define RAPT_FOR_APE_CPP_INCLUDED

/* 

Include file for RAPT for use in APE (Audio Programming Environment). In APE scripts, it has to be 
included as relative path from "[PathToAPE]/Audio Programming Environment/includes". Example:

  #include "../../../../RS-MET/Misc/APE Scripts/rapt_for_ape.cpp"

We make it a .cpp file because it also inlcudes the relevant cpp files from rapt, such that the
relevant code actually gets compiled - otherwise, we would get linker errors.  

*/

// headers:
#include "../../Libraries/RobsJuceModules/rapt/Basics/Basics.h"
#include "../../Libraries/RobsJuceModules/rapt/Data/Data.h"
#include "../../Libraries/RobsJuceModules/rapt/Math/Math.h"
#include "../../Libraries/RobsJuceModules/rapt/AudioBasics/AudioBasics.h"
//#include "../../Libraries/RobsJuceModules/rapt/Music/Music.h"
#include "../../Libraries/RobsJuceModules/rapt/Filters/Filters.h"
#include "../../Libraries/RobsJuceModules/rapt/Analysis/Analysis.h"
#include "../../Libraries/RobsJuceModules/rapt/Physics/Physics.h"
//#include "../../Libraries/RobsJuceModules/rapt/Circuits/Circuits.h"
//#include "../../Libraries/RobsJuceModules/rapt/Spectral/Spectral.h"
#include "../../Libraries/RobsJuceModules/rapt/Visualization/Visualization.h"
#include "../../Libraries/RobsJuceModules/rapt/Generators/Generators.h"
#include "../../Libraries/RobsJuceModules/rapt/Modulators/Modulators.h"
//#include "../../Libraries/RobsJuceModules/rapt/Effects/Effects.h"
//#include "../../Libraries/RobsJuceModules/rapt/Framework/Framework.h"
//#include "../../Libraries/RobsJuceModules/rapt/Instruments/Instruments.h"
#include "../../Libraries/RobsJuceModules/rapt/Unfinished/Unfinished.h"
#include "../../Libraries/RobsJuceModules/rapt/Spectral/Spectral.h"

// implementations:
#include "../../Libraries/RobsJuceModules/rapt/Basics/Basics.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Data/Data.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Math/Math.cpp"
#include "../../Libraries/RobsJuceModules/rapt/AudioBasics/AudioBasics.cpp" 
#include "../../Libraries/RobsJuceModules/rapt/Filters/Filters.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Analysis/Analysis.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Physics/Physics.cpp"
//#include "../../Libraries/RobsJuceModules/rapt/Circuits/Circuits.cpp"
//#include "../../Libraries/RobsJuceModules/rapt/Spectral/Spectral.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Visualization/Visualization.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Generators/Generators.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Modulators/Modulators.cpp"
//#include "../../Libraries/RobsJuceModules/rapt/Effects/Effects.cpp"
//#include "../../Libraries/RobsJuceModules/rapt/Framework/Framework.cpp".
//#include "../../Libraries/RobsJuceModules/rapt/Music/Music.cpp"
//#include "../../Libraries/RobsJuceModules/rapt/Instruments/Instruments.cpp"
#include "../../Libraries/RobsJuceModules/rapt/Unfinished/Unfinished.cpp" 
#include "../../Libraries/RobsJuceModules/rapt/Spectral/Spectral.cpp"

#endif
