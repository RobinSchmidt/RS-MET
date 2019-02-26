Some of the BLAS and LAPACK routines originally needed code from the files in
the LibF2C and XBlas directories. The relevant code has been copied into 
LibF2C.hpp/cpp and XBlas.hpp/cpp and adapted (i.e. templatized, where 
necessary), so these subdirectories are just for reference and to grab new code
from, when the files must be updated to include new functionality.

To build it, just add the files LaPack.hpp/cpp to your project. If you add the
files Blas.cpp, XBlas.cpp etc. also to your project, make sure to exclude them 
from build. They are included from LaPack.cpp in a sort of unity-build way.