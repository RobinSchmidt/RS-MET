#ifdef RAPT_H_INCLUDED
/* When you add this cpp file to your project, you mustn't include it in a file where you've
already included any other headers - just put it inside a file on its own, possibly with your config
flags preceding it, but don't include anything else. That also includes avoiding any automatic prefix
header files that the compiler may be using. */
#error "Incorrect use of JUCE cpp file"
#endif


// This file is where all the non-templated code of the rapt-library goes

#include "Basics/Plotting.cpp"
// only the plotting code gets compiled into rapt.obj - everything else is template code and must
// be included and instantiated elsewhere, for example in rosic
