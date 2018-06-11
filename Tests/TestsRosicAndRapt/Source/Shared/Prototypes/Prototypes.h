#ifndef RS_PROTOTYPES_H
#define RS_PROTOTYPES_H

//#include "rapt/rapt.h"
#include "rosic/rosic.h"

// new implementation of classic IIR filter design:
#include "ClassicFilterDesign/PoleZeroPrototype.h"
#include "ClassicFilterDesign/PoleZeroMapper.h"
#include "ClassicFilterDesign/PoleZeroDesignerAnalog.h"
#include "ClassicFilterDesign/PoleZeroDesignerDigital.h"

#include "Projection3Dto2D.h"


/** This file contains prototypical implementations of algorithms. These prototypes are not meant 
to be used for production code but are useful for a more readable proof-of-concept (because of lack 
of optimizations), for tweaking an algorithm's internal parameters which might not be even exposed 
in the production-code versions, and to create reference output for the unit-tests for production 
code. */

/** Prototype for rsResampler::signalValueViaSincAt(). It provides as additional parameters for 
tweaking: 
-pointer to a window-function
-parameter for the window (if applicable)
-switch for normalizing the output by the sum of the tap weights 
*/
double signalValueViaSincAt(double *x, int N, double t, double sincLength, double stretch,
  //FunctionPointer3DoublesToDouble windowFunction = rsExactBlackmanWindow, 
  double (*windowFunction)(double,double,double) = RAPT::rsExactBlackmanWindow,
  double windowParameter = 0.0, bool normalizeDC = true);

/** Generates polynomial coefficients of the polynomial used in Halpern filters. It's the T^2(w) 
polynomial in Eq. 8.18 in Paarmann: Design and Analysis of Analog Filters. */
void halpernT2(double *c, int N);

/** Generates polynomial coefficients of the polynomial used in Papoulis filters. It's the L^2(w) 
polynomial in Eq. 8.14 in Paarmann: Design and Analysis of Analog Filters */
void papoulisL2(double *c, int N);



#endif
