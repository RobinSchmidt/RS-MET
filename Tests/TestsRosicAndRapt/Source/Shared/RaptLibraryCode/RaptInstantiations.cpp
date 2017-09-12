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

#include "../../../../../Libraries/JUCE/modules/rapt/rapt.cpp"

// Basics:


// Data:
template void RAPT::ArrayTools::rsFillWithRangeLinear(float* x, int N, float min, float max);
template void RAPT::ArrayTools::rsFillWithRandomValues(float* x, int N, double min, double max, int seed);
template void RAPT::ArrayTools::rsFillWithRandomValues(double* x, int N, double min, double max, int seed);

// Math:
template int RAPT::rsLimit(int x, int min, int max);
template RAPT::rsConicSection<float>;
template RAPT::rsEllipse<float>;
template RAPT::rsNormalizedSigmoids<float>;
template RAPT::rsParametricSigmoid<float>;
template RAPT::rsScaledAndShiftedSigmoid<float>;
template RAPT::rsPositiveBellFunctions<float>;    // get rid of rs-prefixes
template RAPT::rsParametricBellFunction<float>;
template RAPT::rsSinCosTable<float>;
template RAPT::rsSinCosTable<double>;

template void RAPT::rsStatistics::linearRegression(int N, float* x, float* y, float& a, float& b);
template float RAPT::rsStatistics::proportionalRegression(int N, float* x, float* y);

// Filters:
template RAPT::rsSmoothingFilter<float, float>;
template RAPT::LadderFilter<float, float>;
template RAPT::PhasorFilter<float, float>;
template RAPT::PhasorStateMapper<float>;
template RAPT::StateVariableFilter<float, float>; 

// Generators:
template RAPT::rsRayBouncer<float>;

// Modulation:
template RAPT::rsBreakpointModulator<float>;

// Visualization:
template RAPT::Image<float>;
template RAPT::AlphaMask<float>;
template RAPT::ImagePainter<float, float, float>;
template RAPT::LineDrawer<float, float, float>;
template RAPT::PhaseScopeBuffer<float, float, double>;

// Modulators:
//template RAPT::rsBreakpointModulator<double>; // will be needed, when the class is templatized
