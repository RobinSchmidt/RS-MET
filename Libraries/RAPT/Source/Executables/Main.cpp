#include "Demos/ArrayDemos.h"
#include "Demos/MathDemos.h"
#include "Demos/FilterDemos.h"
#include "Demos/ModulatorDemos.h"

int main(int argc, char** argv)
{
  // In this main function, you can select, which demo, test or whatever will be run by 
  // uncommenting the respective function call (and maybe commenting the one which was previously
  // active).

  //-----------------------------------------------------------------------------------------------
  // Demos:

  // Array Demos:
  //convolutionDemo();

  // Math Demos:
  parametricBell();
  sigmoids();

  // Filter Demos:
  ladderImpulseResponse();
  svfImpulseResponse();
  // add ladderMagnitudeResponse... - maybe create a template that plots 
  // impulse/step/magnitude/phase - responses in a single window

  // Filter Demos:
  breakpointModulatorDefault();


  // ToDo: check for memory leaks here
}