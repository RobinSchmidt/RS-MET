#ifndef RAPT_MATH_H_INCLUDED
#define RAPT_MATH_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Math

#include "Misc/LinearAlgebra.h"
#include "Misc/Statistics.h"
#include "Misc/CurveFitting.h"

#include "Functions/BasicFunctions.h"
#include "Functions/IntegerFunctions.h"
#include "Functions/InterpolatingFunction.h"
#include "Functions/NodeBasedFunction.h"
#include "Functions/RealFunctions.h"
#include "Functions/ComplexFunctions.h"
#include "Functions/FunctionIterators.h"
#include "Functions/Mappers.h"
#include "Functions/MoebiusTransform.h"
#include "Functions/BellFunctions.h"
#include "Functions/FunctionObjects.h" // obsolete thx to std::function?
#include "Functions/Sigmoids.h"
#include "Functions/SinCosTable.h"

#include "Types/Matrix.h"
#include "Types/Polynomial.h"
#include "Types/Vector.h"

#include "Geometry/Line2D.h"
#include "Geometry/ConicSection.h"
#include "Geometry/Ellipse.h"
#include "Geometry/GeometricTransformations.h"

#include "Numerics/RootFinder.h"
#include "Numerics/NumericCalculus.h"
#include "Numerics/Interpolation.h"
// todo: Optimizer, CurveFitter, Interpolator, Differentiator, InitialValueSolver, 

}

#endif