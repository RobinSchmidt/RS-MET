#ifndef RAPT_MATH_H_INCLUDED
#define RAPT_MATH_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Math

#include "Functions/BasicFunctions.h"
#include "Functions/IntegerFunctions.h"
#include "Functions/RealFunctions.h"
#include "Functions/AudioFunctions.h"     // maybe move into a Music or Audio submodule
#include "Functions/ComplexFunctions.h"
#include "Functions/FunctionIterators.h"
#include "Functions/MoebiusTransform.h"
#include "Functions/BellFunctions.h"
#include "Functions/FunctionObjects.h"
#include "Functions/Sigmoids.h"
#include "Functions/SinCosTable.h"

#include "Types/Matrix.h"
#include "Types/Polynomial.h"
#include "Types/Vector.h"

#include "Geometry/Line2D.h"
#include "Geometry/ConicSection.h"
#include "Geometry/Ellipse.h"

#include "Misc/Statistics.h"

}

#endif