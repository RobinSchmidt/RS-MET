#ifndef PROTOTYPES_FOR_APE_H_INCLUDED
#define PROTOTYPES_FOR_APE_H_INCLUDED

/** This is the header for some prototype DSP classes that need to be accessible from the 
TestsRosicAndRapt project and also from APE scripts. This is necessarry because APE 
currently can't use rosic and therefore also not the rs_testing module due to lack of 
support for some standard C/C++ headers that are used in rosic (intrin.h, etc.). So, currently we 
are limited to use only RAPT classes in APE but not rosic and also not the Prototypes defined in
rs_testing.  */



#endif