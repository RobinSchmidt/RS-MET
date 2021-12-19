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



/*

ToDo:
-Make sure that all classes where this is possible are trivial. Don't define empty destructors
 because this will render types nontrivial, see here (at around 20 min):
   https://www.youtube.com/watch?v=ZDKfKotvEZ4
 -> figure out what else needs to be done to make classes trivial. Maybe we should avoid 
 constructors, too? ...yep - and avoid assigment operators, too - see here:
   https://docs.microsoft.com/en-us/cpp/cpp/trivial-standard-layout-and-pod-types?view=msvc-170
   https://www.geeksforgeeks.org/trivial-classes-c/
   https://mariusbancila.ro/blog/2020/08/10/no-more-plain-old-data/
 ensure the triviality in unit-tests via this:
   https://en.cppreference.com/w/cpp/types/is_trivial

*/