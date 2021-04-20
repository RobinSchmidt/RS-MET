#ifndef RAPT_MATH_H_INCLUDED
#define RAPT_MATH_H_INCLUDED

namespace RAPT
{
  
// make a nested namespace Math

// Why does it compile with Fraction.h being included before IntegerFunctions.h? rsFraction needs
// rsGcd...Try it in gcc as well...yeah...in gcc, it fails

#include "Functions/IntegerFunctions.h"

//#include "LinearAlgebra/LaPackCPP/LaPack.hpp"
#include "LinearAlgebra/BandDiagonalSolver.hpp"
#include "LinearAlgebra/LinearAlgebra.h"

#include "Types/Vector.h"  // move to LinearAlgebra
#include "Types/Matrix.h"
#include "Types/Polynomial.h"
#include "Types/RationalFunction.h"
#include "Types/Fraction.h"

#include "LinearAlgebra/LinearAlgebraNew.h"
// the required order of the includes is a bit messy - maybe matrices, vectors, tensors, 
// polynomials and rational functions should be in an "Algebra" folder




#include "Misc/FourierTransformer.h"
#include "Misc/Statistics.h"
#include "Misc/CurveFitting.h"
#include "Misc/RatioGenerator.h"

#include "Functions/BasicMathFunctions.h"
//#include "Functions/IntegerFunctions.h"  // moved up bcs rsFraction needs rsGcd
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

#include "Geometry/Line2D.h"
#include "Geometry/ConicSection.h"
#include "Geometry/Ellipse.h"
#include "Geometry/GeometricTransformations.h"

#include "Numerics/RootFinder.h"
#include "Numerics/NumericCalculus.h"
#include "Numerics/Interpolation.h"
#include "Numerics/Optimization.h"
// todo: Optimizer, CurveFitter, Interpolator, Differentiator, InitialValueSolver, 

#include "Functions/WindowFunctions.h"

#include "Functions/FunctionOperators.h" 
// may use stuff from NumericCalculus later - move it down, then perhaps this should not yet be in
// the library but rather in the prototypes section - it's not yet used anywhere anyway and perhaps
// never will


}

#endif