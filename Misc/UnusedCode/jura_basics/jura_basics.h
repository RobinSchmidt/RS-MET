/*******************************************************************************
 The block below describes the properties of this module, and is read by
 the Projucer to automatically generate project code that uses it.
 For details about the syntax and how to create or use a module, see the
 JUCE Module Format.txt file.


 BEGIN_JUCE_MODULE_DECLARATION

  ID:               jura_basics
  vendor:           RS-MET
  version:          0.0.1
  name:             basic tools for JUCE/RAPT glue classes
  description:      basic tools for JUCE/RAPT glue classes
  website:          http://www.rs-met.com
  license:          GPL/Commercial

  dependencies:     juce_core
  OSXFrameworks:    
  iOSFrameworks:    

 END_JUCE_MODULE_DECLARATION

*******************************************************************************/


#ifndef JURA_BASICS_H_INCLUDED
#define JURA_BASICS_H_INCLUDED

#include <juce_core/juce_core.h> 
using namespace juce;

#include "../../RAPT/Code/Library/RAPT.h"
using namespace RAPT;


namespace jura  // can we use our own namespace "jura" here or do we have to stick to juce? 
{

#include "control/jura_ScopedPointerLock.h"
#include "control/jura_Parameter.h"
#include "control/jura_AutomatableParameter.h"
#include "control/jura_AutomatableModule.h"

}

#endif   // JURA_BASICS_H_INCLUDED
