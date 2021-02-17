When a script needs access to some DSP class from rapt or rosic, it must 
include the corresponding .cpp file like:

#include "../../../../RS-MET/Libraries/RobsJuceModules/rapt/rapt_for_ape.cpp"

The path must be relative from "Audio Programming Environment/includes", 
whereever that folder happens to be installed. It must be the .cpp and not a 
.h or .hpp because the .cpp files contain the actual implementations of some
member functions of the class templates

ToDo:
-Prepare the library such that we can include also rosic and rs_testing into 
 APE scripts. In particular, we want acces to the prototypes in the rs_testing
 module.