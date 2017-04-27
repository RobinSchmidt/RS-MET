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




/** ToDo: re-order the includes in the order in whcih they depend on each other */

#ifndef ROSIC_H_INCLUDED
#define ROSIC_H_INCLUDED

#include <malloc.h>  // for alloca - try to get rid..
//#include <math.h>

#include "analysis/rosic_Analysis.h"
#include "basics/rosic_Basics.h"
#include "datastructures/rosic_DataStructures.h"
#include "delaylines/rosic_DelayLines.h"
#include "dynamics/rosic_Dynamics.h"
#include "effects/rosic_Effects.h"
#include "filters/rosic_Filters.h"
#include "generators/rosic_Generators.h"
#include "infrastructure/rosic_Infrastructure.h"
#include "instruments/rosic_Instruments.h"
#include "math/rosic_Math.h"
#include "modulators/rosic_Modulators.h"
#include "neural/rosic_Neural.h"
#include "numerical/rosic_Numerical.h"
#include "others/rosic_Others.h"
#include "rendering/rosic_Rendering.h"
#include "scripting/rosic_Scripting.h"
//#include "plugins/rosic_PlugIns.h"
#include "transforms/rosic_Transforms.h"

#endif // #ifndef ROSIC_H_INCLUDED







