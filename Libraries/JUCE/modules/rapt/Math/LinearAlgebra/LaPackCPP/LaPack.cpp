#pragma once


#include "LaPack.hpp"

//#include <stdarg.h>
//#include <cmath> 
//#include <limits>      // uses in lamch to inquire numeric parameters
#include <algorithm>   // for min/max - remove - re-implement

#include "LibF2C.cpp"
#include "Blas.cpp"
#include "XBlas.cpp"
#include "LapackBanded.cpp"

#include "TypedFunctions.cpp"
#include "TemplateInstantiations.cpp"  
// remove - it should be compiled by itself...but then *it* must include the implementation files
// above - not good
