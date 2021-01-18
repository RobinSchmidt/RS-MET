/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               rs_testing
  vendor:           RS-MET
  version:          0.0.1
  name:             Test Tools
  description:      Test tools for rapt and rosic modules
  website:          http://www.rs-met.com
  license:          Custom

  dependencies:     rapt, rosic
  OSXFrameworks:
  iOSFrameworks:

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/

#ifndef RS_TESTING_H_INCLUDED
#define RS_TESTING_H_INCLUDED


#include <rapt/rapt.h>
#include <rosic/rosic.h>



using namespace RAPT;  // get rid of this

namespace RAPT
{
#include "Legacy/FunctionObjects.h"
#include "Legacy/GradientBasedMinimizer.h"
#include "Legacy/MultiLayerPerceptron.h"
}

#include "Prototypes/OscDrivers.h"
#include "Prototypes/PartialDifferentialEquations.h"

#include "Prototypes/BivariatePolynomial.h"
#include "Prototypes/TrivariatePolynomial.h"
#include "Prototypes/PiecewisePolynomial.h"

#include "Prototypes/Prototypes.h"



#include "TestTools/DSPPlotters.h"
#include "TestTools/Plotting.h"

#include "TestTools/Utilities/TestInputCreation.h"
#include "TestTools/Utilities/FileWriting.h"
#include "TestTools/Utilities/PerformanceTestTools.h"
#include "TestTools/Utilities/TestUtilities.h"



#include "Experiments/MiscExperiments.h"

#include "Misc/TestClasses.h"


#endif // #ifndef RS_TESTING_H_INCLUDED
