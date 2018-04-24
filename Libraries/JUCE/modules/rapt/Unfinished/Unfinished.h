#ifndef RAPT_UNFINISHED_H_INCLUDED
#define RAPT_UNFINISHED_H_INCLUDED

/* Here are classes and functions that are not yet finished and not yet ready for use in production
code. It's a sort of construction yard for new code. */

namespace RAPT
{

#include "Data/Flags.h"
#include "Data/MultiArray.h"

#include "Math/Vector.h"
#include "Math/Matrix.h"
#include "Math/FunctionObjects.h"
#include "Math/GradientBasedMinimizer.h"
#include "Math/MultiLayerPerceptron.h"
#include "Math/DifferentialEquationSystem.h"
#include "Math/FourierTransformer.h"
#include "Math/Interpolation.h"
#include "Math/NumberTheory.h"
#include "Math/NumericCalculus.h"
#include "Math/Transforms.h"
#include "Math/Statistics.h"  // merge with the Statistics.h aleady present in the not "Unfinished" Math folder
#include "Math/GeometricFunctions.h"
#include "Math/Point2D.h"
#include "Math/Polygon2D.h"
#include "Math/Rectangle2D.h"
#include "Math/Triangle2D.h"
#include "Math/AffineTransform2D.h"
#include "Math/ModularInteger.h"

// still missing math files from RSLib:  BigInt/BigFloat (should go into rosic)

#include "Audio/BandSplitter.h" // should go to Filters



}

#endif