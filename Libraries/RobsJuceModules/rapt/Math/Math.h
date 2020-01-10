#ifndef RAPT_MATH_H_INCLUDED
#define RAPT_MATH_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Math

//#include "LinearAlgebra/LaPackCPP/LaPack.hpp"
#include "LinearAlgebra/BandDiagonalSolver.hpp"
#include "Types/Vector.h"  // move to LinearAlgebra
#include "Types/Matrix.h"
#include "LinearAlgebra/LinearAlgebra.h"

#include "Misc/FourierTransformer.h"
#include "Misc/Statistics.h"
#include "Misc/CurveFitting.h"
#include "Misc/RatioGenerator.h"

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

#include "Types/Polynomial.h"  // maybe move to Functions

#include "Geometry/Line2D.h"
#include "Geometry/ConicSection.h"
#include "Geometry/Ellipse.h"
#include "Geometry/GeometricTransformations.h"

#include "Numerics/RootFinder.h"
#include "Numerics/NumericCalculus.h"
#include "Numerics/Interpolation.h"
#include "Numerics/Optimization.h"
// todo: Optimizer, CurveFitter, Interpolator, Differentiator, InitialValueSolver, 

#include "Functions/WindowFunctions.h"   // may use FFT stuff for Dolph/Chebychev window later

#include "Functions/FunctionOperators.h" // may use stuff from NumericCalculus later - move it down, then
                                         // perhaps this should not yet be in the library


}

#endif