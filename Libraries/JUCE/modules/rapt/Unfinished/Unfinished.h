#ifndef RAPT_UNFINISHED_H_INCLUDED
#define RAPT_UNFINISHED_H_INCLUDED

/* Here are classes and functions that are not yet finished and not yet ready for use in production
code. It's a sort of construction yard for new code. */

namespace RAPT
{

#include "BandSplitter.h"

#include "Vector.h"
#include "Matrix.h"
#include "FunctionObjects.h"
#include "GradientBasedMinimizer.h"
#include "MultiLayerPerceptron.h"

// for these, updates/adaptions for RAPT are under construction:
#include "DifferentialEquationSystem.h"
#include "FourierTransformer.h"
#include "Interpolation.h"
#include "NumberTheory.h"
#include "NumericCalculus.h"
#include "Transforms.h"
//#include "Statistics.h"  // merge with the Statistics.h aleady present in the Math folder


}

#endif