/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               rosic
  vendor:           RS-MET
  version:          0.0.1
  name:             Rob's Signal Processing Classes
  description:      Library of audio DSP algorithms
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/

#ifndef ROSIC_H_INCLUDED
#define ROSIC_H_INCLUDED

/** ToDo: re-order the includes in the order in which they depend on each other */

//-------------------------------------------------------------------------------------------------
// new, centralized header includes (when finsihed, delete old ones):

#include <malloc.h>  // for alloca - try to get rid..
//#include <math.h>

// basics:
#include "basics/GlobalDefinitions.h"
#include "basics/GlobalFunctions.h"
#include "basics/rosic_Constants.h"
#include "basics/rosic_ChannelMatrix2x2.h"
#include "basics/rosic_HelperFunctions.h"
#include "basics/rosic_NumberManipulations.h"
#include "basics/rosic_FunctionTemplates.h"
// we need to intersperse some includes from other directories before we can finish the includes 
// from basics (todo: fix the dependency structure/layering):
#include "infrastructure/rosic_MutexLock.h"
#include "math/rosic_CephesDeclarations.h"
#include "math/rosic_RealFunctionEvaluationAlgorithms.h"
#include "math/rosic_IntegerFunctions.h"
#include "math/rosic_ElementaryFunctionsReal.h"
#include "math/rosic_SpecialFunctionsReal.h"
#include "math/rosic_Complex.h"
#include "math/rosic_ComplexFunctions.h"
#include "scripting/rosic_ExpressionEvaluatorFunctions.h"
#include "scripting/rosic_ExpressionEvaluatorComplexFunctions.h"
#include "scripting/rosic_ExpressionEvaluator.h"
// ..now we can finish the includes from "basics":
#include "basics/rosic_Interpolator.h"                // needs ElementaryFunctionsReal
#include "basics/rosic_SampleBuffer.h"                // needs MutexLock
#include "basics/rosic_SamplePlaybackParameters.h"    // needs ElementaryFunctionsReal
#include "basics/rosic_TabulatedFunction.h"           // needs ExpressionEvaluator, Mutexlock
#include "basics/rosic_WarpedAllpassInterpolator.h"   // needs ElementaryFunctionsReal
#include "basics/rosic_WindowDesigner.h"              // needs SpecialFunctionsReal

// datastructures:
#include "datastructures/rosic_Array.h"
#include "datastructures/rosic_String.h"
#include "datastructures/rosic_KeyValueMap.h"
#include "datastructures/rosic_ExtensionsForSTL.h"

// math:
#include "math/rosic_Interpolation.h"
#include "math/rosic_LinearAlgebra.h"
#include "math/rosic_Matrix.h"
#include "math/rosic_Vector.h"
#include "math/rosic_MatrixVectorFunctions.h"
#include "math/rosic_PolynomialAlgorithms.h"          // needs Array
#include "math/rosic_PrimeNumbers.h"
//#include "math/rosic_PrimeArray.h"  
#include "math/rosic_Transformations.h"

// good, until here


//-------------------------------------------------------------------------------------------------
// old, de-centralized header includes:

#include "analysis/rosic_Analysis.h"
//#include "basics/rosic_Basics.h"
//#include "datastructures/rosic_DataStructures.h"
#include "delaylines/rosic_DelayLines.h"
#include "dynamics/rosic_Dynamics.h"
#include "effects/rosic_Effects.h"
#include "filters/rosic_Filters.h"
#include "generators/rosic_Generators.h"
#include "infrastructure/rosic_Infrastructure.h"
#include "instruments/rosic_Instruments.h"
//#include "math/rosic_Math.h"
#include "modulators/rosic_Modulators.h"
#include "neural/rosic_Neural.h"
#include "numerical/rosic_Numerical.h"
#include "others/rosic_Others.h"
#include "rendering/rosic_Rendering.h"
#include "scripting/rosic_Scripting.h"
//#include "plugins/rosic_PlugIns.h"
#include "transforms/rosic_Transforms.h"

#endif // #ifndef ROSIC_H_INCLUDED







