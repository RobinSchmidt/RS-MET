#ifndef RAPT_BASICS_H_INCLUDED
#define RAPT_BASICS_H_INCLUDED

namespace RAPT
{

#include "Constants.h"
#include "TypeDefinitions.h"
#include "MacroDefinitions.h"
//#include "MathBasics.h"  // add from RSLib
#include "SortAndSearch.h"

// move to some other file:
template<class T>
inline void rsSwap(T& x, T& y)
{
  T t = x;
  x = y;
  y = t;
}

}

#endif