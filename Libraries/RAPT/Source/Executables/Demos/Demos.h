#ifndef DEMOS_INCLUDED
#define DEMOS_INCLUDED

#include "ArrayDemos.h"
#include "MathDemos.h"
#include "FilterDemos.h"
#include "ModulatorDemos.h"

/** Runs all the demos one at a time. Each will typically create some plot, show it in GNUPlot 
(if installed) and when you close the window, it runs the next and gives the next plot. */
void runDemos();
// todo maybe make it hierarchical: runnAllDemos calls runMathDemos, runFilterDemos, etc


#endif