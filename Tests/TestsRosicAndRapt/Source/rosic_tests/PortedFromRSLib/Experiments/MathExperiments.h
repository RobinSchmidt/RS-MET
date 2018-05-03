#ifndef RS_MATHEXPERIMENTS_H
#define RS_MATHEXPERIMENTS_H

//#include "../ExperimentsUtilities.h" // old
#include "../../../Shared/Shared.h"  // new


void binomialDistribution();
void sineParameters();
void bandLimitedStep();
void cubicInterpolationNonEquidistant();   // move to unit tests
void hyperbolicFunctions(); 
void splineInterpolationNonEquidistant();
void rationalInterpolation();
void splineInterpolationAreaNormalized();
void numericDerivative();
void shiftPolynomial();
//void stretchPolynomial();
void monotonicPolynomials();
void parametricBell();         // move to MathExperiments
void partialFractionExpansion();
void partialFractionExpansion2();
void partialFractionExpansionQuadratic();
void dampedSineEnergy();
void sineIntegral();
void logarithmQuotient();
void stirlingNumbers();
void bernoulliNumbers();
void sequenceSquareRoot();
void conicSystem();
void logisticMapNoise();

//void findPatternForVectorMultiply();
void nestedComplexOrder2();
void nestedComplexOrder3();


// these are some ideas, there's not much code yet:
void bigFloatErrors();

void primeRecursion();
void primeSieveSchmidt1();
void primeSieveSchmidt2();
void primeSieveAtkin();
void primeSieve();
void primeDistribution();
void numberTheoreticTransform();



#endif 
