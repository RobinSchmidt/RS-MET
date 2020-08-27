#include "Math.h"

namespace RAPT
{

//#if !defined LAPACK_COMPILED
//#include "LinearAlgebra/LaPackCPP/LaPack.cpp"  // multiple definition errors (*)
//#define LAPACK_COMPILED
//#endif
// ...nope - that conditional compilation also doesn't work

#include "LinearAlgebra/BandDiagonalSolver.cpp"
#include "LinearAlgebra/LinearAlgebra.cpp"
#include "LinearAlgebra/LinearAlgebraNew.cpp"

#include "Misc/FourierTransformer.cpp"
#include "Misc/Statistics.cpp"
#include "Misc/CurveFitting.cpp"
#include "Misc/RatioGenerator.cpp"

#include "Functions/BasicMathFunctions.cpp"
#include "Functions/IntegerFunctions.cpp"
#include "Functions/InterpolatingFunction.cpp"
#include "Functions/NodeBasedFunction.cpp"
#include "Functions/RealFunctions.cpp"
#include "Functions/ComplexFunctions.cpp"
#include "Functions/FunctionIterators.cpp"
#include "Functions/Mappers.cpp"
#include "Functions/BellFunctions.cpp"
#include "Functions/Sigmoids.cpp"
#include "Functions/SinCosTable.cpp"

#include "Types/Matrix.cpp"
#include "Types/Polynomial.cpp"
#include "Types/Vector.cpp"

#include "Geometry/Line2D.cpp"
#include "Geometry/ConicSection.cpp"
#include "Geometry/Ellipse.cpp"
#include "Geometry/GeometricTransformations.cpp"

#include "Numerics/RootFinder.cpp"
#include "Numerics/NumericCalculus.cpp"
#include "Numerics/Interpolation.cpp"
#include "Numerics/Optimization.cpp"

#include "Functions/WindowFunctions.cpp"
#include "Functions/FunctionOperators.cpp" 

}

// how to fix the LaPack linking problems
// for the multiple definition errors, :
// -if the function is small: inline it
// -if it's larger: maybe (artificially) templatize it, even if that doesn't make much sense
// or: collect the non-templatized functions in some other cpp file, such that it can be included 
// by rosic and not by rapt (so the functions only go into the rosic compilation unit
// ...hmm...somehow the problem seems to arise only because of the conflict between Shared.obj and
// rosic.obj
// -maybe don't include LaPack.cpp here in rapt/Math.cpp - instead include it in rosic.cpp (via
//  TemplateInstatiations.cpp)
// ...nope - that doesn't work either
// i need a way of being able to have some non-templated functions compiled into the rapt.obj
// compilation unit without having them compiled again into rosic.obj and/or Shared.obj
// the problem is that rosic.cpp includes rapt.cpp because it must isntantiate the rapt templates
// ...maybe use conditional compilation in rapt.cpp