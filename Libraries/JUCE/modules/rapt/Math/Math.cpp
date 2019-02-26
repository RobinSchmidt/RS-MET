#include "Math.h"

namespace RAPT
{

//#include "LinearAlgebra/LaPackCPP/LaPack.cpp"  // multiple definition errors
#include "LinearAlgebra/BandDiagonalSolver.cpp"
#include "LinearAlgebra/LinearAlgebra.cpp"

#include "Misc/Statistics.cpp"
#include "Misc/CurveFitting.cpp"

#include "Functions/BasicFunctions.cpp"
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

}