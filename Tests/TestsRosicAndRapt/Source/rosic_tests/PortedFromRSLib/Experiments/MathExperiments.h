#ifndef RS_MATHEXPERIMENTS_H
#define RS_MATHEXPERIMENTS_H

// maybe move to RaptTests.h and get rid of this file

void determinant();
void characteristicPolynomial();
bool testSubSpaces();
bool testSigularValueDecomp();
void linearIndependence();
void eigenstuff();
void iterativeLinearSolvers();
void linearSolverPrecision();

void bandMatrix();
void pentaDiagnonalMatrix();
void pentaDiagnonalMatrix2();

void minSqrdDifsForFixSums();
void minSqrdCurvForFixSums();

void binomialDistribution();
void sineParameters();
void bandLimitedStep();
void cubicInterpolationNonEquidistant();   // move to unit tests
void chebychevInterpolant();
void cubicSplines();
void hyperbolicFunctions(); 
void splineInterpolationNonEquidistant();
void rationalInterpolation();
void splineInterpolationAreaNormalized();
void numericDifferentiation();
void numericIntegration();
void nonUniformArrayDiffAndInt();
void uniformArrayDiffAndInt();
void iteratedNumDiff();
//void vertexMeshGradient();
//void vertexMeshHessian();
//void vertexMeshLaplacian();

void convolvePolynomials();
void convolvePiecewise();
void shiftPolynomial();
//void stretchPolynomial();
void monotonicPolynomials();
void mixedPolynomialRoots();
void parametricBell();
void partialFractionExpansion();
void partialFractionExpansion2();
void partialFractionExpansion3();
//void partialFractionExpansion4();

void partialFractionExpansionQuadratic();
void dampedSine();
void sineIntegral();
void logarithmQuotient();
void stirlingNumbers();
void bernoulliNumbers();
void sequenceSquareRoot();
void conicSystem();
void logisticMapNoise();
void variousFunctions(); 
void functionOperators();

// these are some ideas, there's not much code yet:
void bigFloatErrors();

void primeRecursion();
void primeSieveSchmidt1();
void primeSieveSchmidt2();
void primeSieveAtkin();
void primeSieve();
void primeDistribution();
void numberTheoreticTransform();
void numberTheoreticTrafoModuli();

void powerIterator();
void gaussianIterator();
void expPolyIterator();
void reciprocalIterator();

//-------------------------------------------------------------------------------------------------
// Functions for dealing with meshes - their implementations are in a a separate .cpp file 
// MeshExperiments.cpp because the amount of code is quite large

void derivativeFormulas();
void vertexMeshGradient();
void vertexMeshHessian();
void vertexMeshLaplacian();



#endif 
