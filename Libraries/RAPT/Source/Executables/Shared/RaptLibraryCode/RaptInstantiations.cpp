#include "../../../../Source/Modules/RAPT.cpp"

/** Many RAPT functions and classes are written as C++ templates in order to be type independent. 
This implies that the compiler does not actually compile any code unless the template gets 
instantiated from somewhere inside the client code. That means, if you include the RAPT.h file from 
somewhere in your client code, your code may compile correctly, but you'll get linker errors of the 
type "unresolved external symbol". So, in order to make your code also link correctly, you'll have 
to include the RAPT.cpp somewhere - but if you include it in multiple places (because you need RAPT
classes in many of your source files), you may get linker errors of the type "multiple definition" 
or "symbol already defined" or something like that. It means that a particular template was 
instantiated in multiple compilation units and now the linker does not know which one to link to.
The only resolution i could think of was to consolidate all template instantiations into one 
single .cpp file in the client code. That file must include RAPT.cpp (and it must be the only file
that does so), and you must request the compiler explicitly here to instantiate the templates that
your client code will need. This way (1) instantiations are available to the linker in your project
and (2) there will be only one such instantiation of a particular template for a particular type. 
If somewhere you'll get an "undefined external symbol" error, add the respective request for 
explicit template instantiation here. 

For the demos here, we instantiate all classes and functions for "float" (because "double" would 
not show all the trunction-warnings that we get for "float").

\todo: provide a predefined template instantiation file within the RAPT library folder, that
instantiates all templates for double. That file can be used by client code by default but client
code may also define its own instantiation file. */

// Data:
template void RAPT::ArrayTools::rsFillWithRangeLinear(float* x, int N, float min, float max);
template void RAPT::ArrayTools::rsFillWithRandomValues(float* x, int N, double min, double max, int seed);

// Math:
template RAPT::NormalizedSigmoids<float>;
template RAPT::ParametricSigmoid<float>;
template RAPT::ScaledAndShiftedSigmoid<float>;
template RAPT::rsPositiveBellFunctions<float>;    // get rid of rs-prefixes
template RAPT::rsParametricBellFunction<float>;

// Filters:
template RAPT::LadderFilter<float, float>;
template RAPT::PhasorFilter<float, float>;
template RAPT::PhasorStateMapper<float>;
template RAPT::StateVariableFilter<float, float>; 

// Visualization
template RAPT::Image<float>;
template RAPT::AlphaMask<float>;
template RAPT::ImagePainter<float, float, float>;
template RAPT::PhaseScopeBuffer<float, float, double>;

// Modulators:
//template RAPT::rsBreakpointModulator<double>; // will be needed, when the class is templatized
