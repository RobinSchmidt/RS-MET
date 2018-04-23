#ifndef RAPT_UNFINISHED_H_INCLUDED
#define RAPT_UNFINISHED_H_INCLUDED

/* Here are classes and functions that are not yet finished and not yet ready for use in production
code. It's a sort of construction yard for new code. */

namespace RAPT
{

#include "Flags.h" // should go to Data

#include "BandSplitter.h" // should go to Filters

// should go to Math:
#include "Vector.h"
#include "Matrix.h"
#include "FunctionObjects.h"
#include "GradientBasedMinimizer.h"
#include "MultiLayerPerceptron.h"
#include "DifferentialEquationSystem.h"
#include "FourierTransformer.h"
#include "Interpolation.h"
#include "NumberTheory.h"
#include "NumericCalculus.h"
#include "Transforms.h"
#include "Statistics.h"  // merge with the Statistics.h aleady present in the Math folder
#include "GeometricFunctions.h"
#include "Point2D.h"
#include "Polygon2D.h"
#include "Rectangle2D.h"
#include "Triangle2D.h"
#include "AffineTransform2D.h"

// for these, updates/adaptions for RAPT are under construction:
#include "ModularInteger.h"
#include "MultiArray.h"

// still missing math files from RSLib:  BigInt/BigFloat (should go into rosic)


}

#endif