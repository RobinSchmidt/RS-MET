/** Many RAPT functions and classes are written as C++ templates in order to be type independent.
This implies that the compiler does not actually compile any code unless the template gets
instantiated from somewhere inside the client code. That means, if you include the RAPT.h file from
somewhere in your client code, your code may compile correctly, but you'll get linker errors of the
type "unresolved external symbol". So, in order to make your code also link correctly, you'll have
to include the RAPT.cpp somewhere - but if you include it in multiple places (because you need RAPT
classes in many of your source files), you may get linker errors of the type "multiple definition"
or "symbol already defined" or something like that. It means that a particular template was
instantiated in multiple compilation units and now the linker does not know which one to link to.
The only resolution i could think of was to consolidate all template instantiations into one
single .cpp file in the client code. That file must include RAPT.cpp (and it must be the only file
that does so), and you must request the compiler explicitly here to instantiate the templates that
your client code will need. This way (1) instantiations are available to the linker in your project
and (2) there will be only one such instantiation of a particular template for a particular type.
If somewhere you'll get an "undefined external symbol" error, add the respective request for
explicit template instantiation here.

For the demos here, we instantiate all classes and functions for "float" (because "double" would
not show all the trunction-warnings that we get for "float").

\todo: provide a predefined template instantiation file within the RAPT library folder, that
instantiates all templates for double. That file can be used by client code by default but client
code may also define its own instantiation file. */

//#include "../../../../../Libraries/RobsJuceModules/rapt/rapt_templates.cpp"
//#include "../../../../../Libraries/JUCE/modules/rosic/rosic.cpp"

#include "../rapt/rapt_templates.cpp"
using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

// we need to make sure to provide only instantiations that are not already there in rosic.cpp
// todo: move instantiations to rosic, where appropriate

typedef std::complex<double> cmplxD;

// Basics:


template bool RAPT::rsGreater(const int&, const int&);
template bool RAPT::rsGreater(const double&, const double&);

template bool RAPT::rsLess(const double&, const double&);

template class rsMovingMaximumFilter<int>;
template class rsMovingMaximumFilter<double>;
template class rsQuantileFilter<double>;

template bool RAPT::rsLess(const int& left, const int& right);

template void RAPT::rsHeapSort(int *buffer, int length,
  bool (*lessThen)(const int& left, const int& right));

template void RAPT::rsHeapSort(rsRange<double>* buffer, int length,
  bool (*lessThen)(const rsRange<double>& left, const rsRange<double>& right));

//template bool RAPT::rsIsSortedAscending(int *buffer, int length);
template bool rsArrayTools::isSortedAscending(const int *buffer, int length);


template std::vector<int> RAPT::rsFindAllOccurencesOf(char* buffer, int bufferLength,
  char* pattern, int patternLength);
template int RAPT::rsFindHighestPeakIndex(double *buffer, int length); // move to rsArrayTools
// maybe move these into rsArrayTools

//-------------------------------------------------------------------------------------------------
// Data:

// rsArrayTools<char>
template void rsArrayTools::fillWithValue(char* buffer, int length, char value);

// rsArrayTools<int>
template int rsArrayTools::findSplitIndex(const int* A, int N, int key);
template int rsArrayTools::copyIfMatching(const int *, int *, int, const int *, int);
template int rsArrayTools::copyIfNotMatching(const int *, int *, int, const int *, int);
template void rsArrayTools::copySection(const int *source, int sourceLength, int *destination, int copyStart, int copyLength);
template void rsArrayTools::cumulativeSum(const int *x, int *y, int length, int order);
template void rsArrayTools::fillWithRangeLinear(int* x, int N, int min, int max);
template void rsArrayTools::fillWithRandomValues(int* x, int N, double min, double max, int seed);
template void rsArrayTools::allocateSquareArray2D(int**& theArray, int size);
template void rsArrayTools::deAllocateSquareArray2D(int**& theArray, int size);
template void rsArrayTools::rightShift(int *buffer, int length, int numPlaces);

// rsArrayTools<rsUint32>
template bool rsArrayTools::contains(const rsUint32 *buffer, int length, rsUint32 elementToFind);
template void rsArrayTools::copy(const rsUint32 *src, rsUint32 *dst, int N);
template void rsArrayTools::fillWithRangeLinear(rsUint32* x, int N, rsUint32 min, rsUint32 max);
template int rsArrayTools::firstIndexWithNonZeroValue(const rsUint32 *a, int N);
template void rsArrayTools::fillWithZeros(rsUint32 *buffer, int length);
template rsUint32 rsArrayTools::maxValue(const rsUint32 *x, int length);

// rsArrayTools<float>
template void rsArrayTools::fillWithRange(float* x, int N, float min, float max, float shape);
template void rsArrayTools::fillWithRandomValues(float* x, int N, double min, double max, int seed); // ?why the mix of float and double?
template void rsArrayTools::fillWithRandomIntegers(double* x, int N, int min, int max, int seed);
template float rsArrayTools::generalizedMean(const float *x, int length, float p);
template float rsArrayTools::maxAbs(const float *x, int length);
template float rsArrayTools::minValue(const float *x, int length);
template float rsArrayTools::maxValue(const float *x, int length);
template float rsArrayTools::maxDeviation(const float *buffer1, const float *buffer2, int length);
//template void rsArrayTools::negate(    const float *x, float *y, int N);
template void rsArrayTools::negateEven(const float *x, float *y, int N);
template void rsArrayTools::negateOdd( const float *x, float *y, int N);
template int rsArrayTools::maxDeviationIndex(const float *buffer1, const float *buffer2, int length);
template void rsArrayTools::reverse(const float* x, float* y, int length);
//template void rsArrayTools::reverse(float* x, int length);

// rsArrayTools<double>
template void rsArrayTools::applyFunction(const double *x, double *y, int N, double (*f) (double));
template void rsArrayTools::cumulativeSum(const double *x, double *y, int length, int order);
template void rsArrayTools::deConvolve(const double *y, int yLength, const double *h, int hLength, double *x);
template void rsArrayTools::deInterleave(double*, int, int);
template void rsArrayTools::difference(double *buffer, int length, int order, bool periodic);
template void rsArrayTools::divide(const double *buffer1, const double *buffer2, double *result, int length);
//template void rsArrayTools::fillWithRangeLinear(double* x, int N, double min, double max);
//template void rsArrayTools::fillWithRangeExponential(double* x, int N, double min, double max);
template void rsArrayTools::fillWithRange(double* x, int N, double min, double max, double shape);
template void rsArrayTools::fillWithRandomValues(double* x, int N, double min, double max, int seed);
template void rsArrayTools::filter(const double *x, int xLength, double *y, int yLength,
  const double *b, int bOrder, const double *a, int aOrder);
template void rsArrayTools::filterBiDirectional(const double *x, int xLength, double *y, int yLength,
  const double *b, int bOrder, const double *a, int aOrder, int numRingOutSamples);
template double rsArrayTools::generalizedMean(const double *x, int length, double p);
template double rsArrayTools::maxAbs(const double *x, int length);
template int rsArrayTools::maxAbsIndex(const double* const buffer, int length);
template double rsArrayTools::maxValue(const double *x, int length);
template double rsArrayTools::mean(const double *x, int length);
template double rsArrayTools::median(const double *x, int length);
template double rsArrayTools::maxDeviation(const double *buffer1, const double *buffer2, int length);
template double rsArrayTools::minValue(const double *x, int length);
template void rsArrayTools::movingAverage3pt(const double* x, int N, double* y, bool endsFixed);
template void rsArrayTools::movingMedian3pt(const double* x, int N, double* y);
template void rsArrayTools::negate(const double *x, double *y, int N);
template void rsArrayTools::normalize(double *buffer, int length, double maximum, bool subtractMean);
template void rsArrayTools::normalizeMean(double *x, int N, double newMean);
template void rsArrayTools::reverse(const double* x, double* y, int length);
template void rsArrayTools::rightShift(double *buffer, int length, int numPlaces);
template void rsArrayTools::sequenceSqrt(const double *y, int yLength, double *x);
template void rsArrayTools::shift(double *buffer, int length, int numPlaces);
template void rsArrayTools::transposeSquareArray(double **in, double **out, int size);
template void rsArrayTools::unwrap(double* a, int N, double p);

// rsArrayTools<complex<double>>
template void rsArrayTools::convert(const std::complex<double> *source, std::complex<double> *destination, int length);
template void rsArrayTools::fillWithRandomValues(std::complex<double>* x, int N, double min, double max, int seed);
//template int rsArrayTools::maxAbsIndex(const std::complex<double>* const buffer, int length);

// rsArrayTools<rsRange>
template int rsArrayTools::maxIndex(const rsRange<double>*, int length);

template double rsArrayTools::maxAbs(const std::complex<double> *buffer, int length);


// rsArrayTools<rosic::rsFloat32x4>
//template void rsArrayTools::fillWithRandomValues(rosic::rsFloat32x4* x, int N, double min, double max, int seed);

template class RAPT::rsBinaryHeap<int>;
template class RAPT::rsDoubleHeap<int>;


template void rsMatrixTools::initMatrix(double** A, int N, int M, double value);
template void rsMatrixTools::copyMatrix(double** source, double **destination, int N, int M);
template bool rsMatrixTools::areMatricesApproximatelyEqual(double **A, double **B, int N, int M, double tol);

template void rsMatrixTools::matrixMultiply(double **A, double **B, double **C, int N, int M, int P);
template void rsMatrixTools::matrixMultiplyFirstTransposed(double **A, double **B, double **C, int N, int M, int P);
template void rsMatrixTools::matrixMultiplySecondTransposed(double **A, double **B, double **C, int N, int M, int P);
template void rsMatrixTools::matrixMultiplyBothTransposed(double **A, double **B, double **C, int N, int M, int P);
template void rsMatrixTools::matrixInPlaceMultiply(double **A, double **B, int N, int M);


//-------------------------------------------------------------------------------------------------
// Math:

//template RAPT::rsLinearAlgebra<float>; // doens't work bcs the template parameters are decalred in the member functions
template void rsLinearAlgebra::rsSolveLinearSystem2x2(const double A[2][2], double x[2], const double y[2]);
template void rsLinearAlgebra::rsSolveLinearSystem3x3(const double A[3][3], double x[3], const double y[3]);

template double rsLinearAlgebra::eigenvalue2x2_1(double, double, double, double);
template double rsLinearAlgebra::eigenvalue2x2_2(double, double, double, double);
template void rsLinearAlgebra::eigenvector2x2_1(double, double, double, double, double*, double*, bool);
template void rsLinearAlgebra::eigenvector2x2_2(double, double, double, double, double*, double*, bool);

template bool rsLinearAlgebra::rsSolveLinearSystem(double **A, double *x, const double *b, int N);
template bool rsLinearAlgebra::rsInvertMatrix(double **A, int N);
template bool rsLinearAlgebra::rsSolveTridiagonalSystem(double *lowerDiagonal, double *mainDiagonal,
  double *upperDiagonal, double *rightHandSide, double *solution, int N);
template bool rsLinearAlgebra::rsChangeOfBasisColumnWise(double **A, double **B, double *va,
  double *vb, int N);
template bool rsLinearAlgebra::rsChangeOfBasisRowWise(double **A, double **B, double *va,
  double *vb, int N);
template bool rsLinearAlgebra::rsChangeOfBasisMatrixColumnWise(double **A, double **B, double **C, int N);
template bool rsLinearAlgebra::rsChangeOfBasisMatrixRowWise(   double **A, double **B, double **C, int N);
template bool rsLinearAlgebra::rsSolveLinearSystem(cmplxD **A, cmplxD *x, const cmplxD *b, int N);

template std::vector<double> rsLinearAlgebraNew::solveOld(rsMatrix<double> A, std::vector<double> b);

template class RAPT::rsMatrixOld<double>;  // try to get rid


template class RAPT::rsPolynomial<float>;
template class RAPT::rsPolynomial<double>;
template class RAPT::rsPolynomial<std::complex<double>>;
//template class RAPT::rsPolynomial<std::complex<float>>;  // template doesn't compile with float
//template  class RAPT::rsPolynomial<int>;                 // template doesn't compile with int
// todo: instantiate rsPolynomial also for float, int and rsFraction<int>, maybe also for 
// rsMatrix<float>, etc.

template void RAPT::rsPolynomial<double>::divideByMonomialInPlace(double*, int, double, double*);
  // needs separate instantiation because function itself has a (second) template parameter

//template void RAPT::rsPolynomial<std::complex<double>>::subtract(
//  const std::complex<double>* p, int pN, const std::complex<double>* q, int qN,
//  std::complex<double>* r);
// duplicate explicit instantiation in gcc


template void  RAPT::rsPolynomial<float>::rootsQuadraticComplex(
  const std::complex<float>& a0, const std::complex<float>& a1, const std::complex<float>& a2,
  std::complex<float>* x1, std::complex<float>* x2);
template float RAPT::rsPolynomial<float>::cubicDiscriminant(
  const float& a0, const float& a1, const float& a2, const float& a3);

template void RAPT::rsPolynomial<float>::rootsCubicComplex(
  std::complex<float> a0, std::complex<float> a1,
  std::complex<float> a2, std::complex<float> a3,
  std::complex<float>* r1, std::complex<float>* r2, std::complex<float>* r3);

template void RAPT::rsPolynomial<std::complex<double>>::roots(
  const std::complex<double>* a, int degree, std::complex<double>* roots);

template std::vector<std::complex<double>> RAPT::rsRationalFunction<double>::partialFractions(
  const std::vector<std::complex<double>>& numerator,
  const std::vector<std::complex<double>>& denominator,
  const std::vector<std::complex<double>>& poles);

template std::vector<std::complex<double>> RAPT::rsRationalFunction<double>::partialFractions(
  const std::vector<std::complex<double>>& numerator,
  const std::vector<std::complex<double>>& denominator,
  const std::vector<std::complex<double>>& poles,
  const std::vector<int>& multiplicities);

template class RAPT::rsRationalFunction<std::complex<double>>;


template class RAPT::rsMatrixView<double>;
template class RAPT::rsMatrix<double>;

template std::vector<double> RAPT::rsLinearAlgebraNew::solve(
  const rsMatrixView<double>& A, const std::vector<double>& B);



template std::vector<std::complex<double>> RAPT::rsLinearAlgebraNew::solve(
  const rsMatrixView<std::complex<double>>& A,
  const std::vector<std::complex<double>>& B);

// doesn't compile because the > comparison in the pivoting doesn't work with complex numbers
// possible solution: implement a > operator for std::complex numbers - compare real parts first,
// then imag
// or replace
//  if(rsAbs(A(j, i)) > maxAbs)
// with
//  if(rsAbs(A(j, i)) > rsAbs(maxAbs))
// or use comparison function rsGreater with an explicit specialization for complex<double>
// or use rsGreaterAbs

// figure out!

// these instantiations need some more operations defined on rsRationalFunction
//template std::vector<RAPT::rsRationalFunction<double>> RAPT::rsLinearAlgebraNew::solveLinearSystem(
//  rsMatrixView<RAPT::rsRationalFunction<double>>& A, std::vector<RAPT::rsRationalFunction<double>>& B);

/*
template void RAPT::rsLinearAlgebraNew::makeTriangularNoPivot(
  rsMatrixView<RAPT::rsRationalFunction<double>>& A,
  rsMatrixView<RAPT::rsRationalFunction<double>>& B);
  */

template RAPT::rsMatrix<double> RAPT::rsLinearAlgebraNew::inverse(
  const RAPT::rsMatrixView<double>& A);

template int RAPT::rsLinearAlgebraNew::makeDiagonal(
  RAPT::rsMatrixView<double>& A, RAPT::rsMatrixView<double>& B);


// should go into an "Interpolation" class...or maybe CurveFitter
//template void RAPT::fitCubicWithDerivative(double x1, double x2, double y1, double y2, double yd1,
//  double yd2, double *a3, double *a2, double *a1, double *a0);
template void RAPT::fitCubicWithDerivativeFixedX(double y0, double y1, double yd0, double yd1,
  double *a3, double *a2, double *a1, double *a0);


template class RAPT::UnivariateScalarFunction<double>;
template class RAPT::UnivariateScalarFunctionViaPointer<double>;
template class RAPT::QuadraticTestErrorFunction<double>;

template class RAPT::GradientBasedMinimizer<double>;

template class RAPT::MultiLayerPerceptron<double>;
template class RAPT::MultiLayerPerceptronErrorFunction<double>;
template class RAPT::MultiLayerPerceptronTrainer<double>;

template double RAPT::rsInterpolateWrapped(double, double, double, double, double);

template bool RAPT::rsCurveFitter::fitExponentialSum(double* y, int numValues, double* A, double* a,
  int numExponentials);

template rsUint32 RAPT::rsBinomialCoefficient(rsUint32 n, rsUint32 k);
template rsUint32 RAPT::rsBinomialCoefficientUpTo20(rsUint32 n, rsUint32 k);
template int RAPT::rsMultinomialCoefficient(int* n, int k);
template long long RAPT::rsGcd(long long x, long long y);  // maybe use rsInt64
//template int RAPT::rsGcd(int, int);
template int RAPT::rsLcm(int, int);
template rsUint32 RAPT::rsLcm(rsUint32, rsUint32);
template rsUint32 RAPT::rsMultinomialCoefficient(rsUint32* n, rsUint32 k);
template rsUint32 RAPT::rsMultinomialCoefficientUpTo12(rsUint32* n, rsUint32 k);
template int RAPT::rsLeviCivita(int indices[], int N);

template void RAPT::rsStirlingNumbersFirstKind(int **s, int nMax);
template void RAPT::rsFillPrimeTable(int*, rsUint32, rsUint32);
template void RAPT::rsEGCD(int x, int y, int& a, int& b, int& g);
template int RAPT::rsModularInverse(const int& x, const int& m);
template int RAPT::rsPrimeModularInverse(const int& x, const int& p);
template int RAPT::rsPrimeModularInverse2(const int& x, const int& p);
template int RAPT::rsChineseRemainderTheorem(int* remainders, int* moduli, rsUint32 count);
//template rsUint32 RAPT::rsChineseRemainderTheorem(rsUint32* remainders, rsUint32* moduli, rsUint32 count);
template void RAPT::rsFillPrimeTable(rsUint32 *primes, rsUint32 numPrimes, rsUint32 bufferSize);

//template void RAPT::rsNextPascalTriangleLine(const double* x, double* y, int N);
template void RAPT::rsPascalTriangleLine(double* y, int N);
template void RAPT::rsPascalTriangleLine(float*  y, int N);
template void RAPT::rsPascalTriangleLine(std::complex<double>* y, int N);


template void RAPT::smbFft(float *fftBuffer, long fftFrameSize, long sign);
template void RAPT::rsDFT(std::complex<double> *buffer, int N);
template void RAPT::rsFFT(std::complex<double> *buffer, int N);
template void RAPT::rsIFFT(std::complex<double> *buffer, int N);
template void RAPT::rsMagnitudeAndPhase(double *signal, int N, double *magnitudes, double *phases);
template void RAPT::rsCrossCorrelation(float x[], float y[], int N, float r[], bool removeBias);
template void RAPT::rsAutoCorrelationFFT(float x[], int N, float r[]);


template class RAPT::rsFourierTransformerRadix2<float>;
template class RAPT::rsFourierTransformerBluestein<float>;

template class RAPT::rsFourierTransformerRadix2<double>;
template class RAPT::rsFourierTransformerBluestein<double>;



template class RAPT::rsSineIterator<double>;
template class RAPT::rsComplexExponentialIterator<double>;

template double RAPT::rsNormalizedSinc(double x);
template double RAPT::rsSineIntegral(double x);

template int RAPT::rsClip(int x, int min, int max);
template int RAPT::rsFactorial(int n);

//template double RAPT::rsWrapToInterval(double x, double min, double max); // rename to rsWrap

template void RAPT::rsBinomialDistribution(double*, int, double);


template class RAPT::rsConicSection<float>;
template class RAPT::rsRotationXY<double>;
template class RAPT::rsGeometricTransforms<double>;

template class RAPT::rsCoordinateMapper<float>;
template class RAPT::rsCoordinateMapper<double>;
template class RAPT::rsCoordinateMapper2D<float>;
template class RAPT::rsCoordinateMapper2D<double>;

template class RAPT::rsEllipse<float>;
template class RAPT::rsLine2D<float>;
template class RAPT::rsNormalizedSigmoids<float>;
template class RAPT::rsParametricSigmoid<float>;
template class RAPT::rsScaledAndShiftedSigmoid<float>;
template class RAPT::rsPositiveBellFunctions<float>;
template class RAPT::rsParametricBellFunction<float>;
template class RAPT::rsRootFinder<float>;
template class RAPT::rsInterpolatingFunction<float, double>;
template class RAPT::rsInterpolatingFunction<double, double>;
template class RAPT::rsNodeBasedFunction<float>;
template class RAPT::rsSinCosTable<float>;
template class RAPT::rsSinCosTable<double>;
template class RAPT::rsVector2D<float>;
template class RAPT::rsVector3D<float>;

//template class RAPT::rsFourierTransformerRadix2<double>;  // duplicate instantiation

template void RAPT::rsStatistics::linearRegression(int N, const float* x, const float* y, float& a, float& b);
template float RAPT::rsStatistics::proportionalRegression(int N, const float* x, const float* y);
template void RAPT::rsRemoveCorrelationBias(double x[], int N, double r[]);

template rsPolynomial<double> RAPT::rsCurveFitter::fitPolynomial(
  double* x, double* y, int numDataPoints, int degree);

template class RAPT::rsMultiArrayOld<float>;

// numeric calculus
template void RAPT::rsNumericDifferentiator<double>::stencilCoeffs(
  const double* x, int N, int d, double* c);

template void RAPT::rsNumericDifferentiator<double>::hessian(
  const std::function<double(double*)>& f, double* x, int N, double* H, const double* h);

template void RAPT::rsNumericDifferentiator<double>::derivative(
  const double *x, const double *y, double *yd, int N, bool extrapolateEnds);

template void RAPT::rsNumericDifferentiator<float>::gradient2D(
  const rsGraph<rsVector2D<float>, float>& mesh, const std::vector<float>& u, 
  std::vector<float>& u_x, std::vector<float>& u_y);

template void RAPT::rsNumericDifferentiator<double>::gradient2D(
  const rsGraph<rsVector2D<double>, double>& mesh, const std::vector<double>& u, 
  std::vector<double>& u_x, std::vector<double>& u_y);

template void RAPT::rsNumericDifferentiator<double>::laplacian2D(
  const rsGraph<rsVector2D<double>, double>& mesh, const double* u, double* L);

template void RAPT::rsNumericDifferentiator<double>::laplacian2D_2(
  const rsGraph<rsVector2D<double>, double>& mesh, const std::vector<double>& u, 
  std::vector<double>& L);


template double RAPT::rsBandwidthConverter::bandedgesToCenterFrequency(double fl, double fu);
template double RAPT::rsBandwidthConverter::bandedgesToAbsoluteBandwidth(double fl, double fu);
template double RAPT::rsBandwidthConverter::absoluteBandwidthToQ(double fl, double fu);

// Filters-Basic:
template class RAPT::rsOnePoleFilter<float, float>;

// Filters-Musical:
template class RAPT::rsSmoothingFilter<float, float>;
template class RAPT::rsLadderFilter<float, float>;
template class RAPT::rsLadderFilter<double, double>;
//template class RAPT::rsLadderFilter<rsFloat64x2, double>;
template class RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2>;
template class RAPT::rsPhasorFilter<float, float>;
template class RAPT::rsPhasorStateMapper<float>;
template class RAPT::rsStateVariableFilter<float, float>;

// Filters-Scientific:
template class RAPT::rsPrototypeDesigner<float>;
//template class RAPT::rsPrototypeDesigner<double>;
template class RAPT::rsPoleZeroMapper<float>;
template class RAPT::rsFilterCoefficientConverter<float>;
template class RAPT::rsInfiniteImpulseResponseDesigner<float>;
template class RAPT::rsFilterAnalyzer<float>;
template class RAPT::rsBiquadCascade<float, float>;
template class RAPT::rsEngineersFilter<float, float>;
template class RAPT::rsLinkwitzRileyCrossOver<float, float>;
template class RAPT::rsCrossOver4Way<float, float>;
template class RAPT::rsDirectFormFilter<float, float>;
template class RAPT::rsEllipticSubBandFilter<float, float>;
template class RAPT::rsEllipticSubBandFilterDirectForm<float, float>;
template class RAPT::rsQuadratureNetwork<float, float>;
template void RAPT::rsBiDirectionalFilter::applyLowpass(
  const double *x, double *y, int N, double fc, double fs, int numPasses, double gc);

//template class RAPT::rsQuantileFilter<double>;



// Physics:
template class RAPT::rsParticleSystem<float>;

// Generators:
template class RAPT::rsTriSawOscillator<float>;
template class RAPT::rsBouncillator<float>;
template class RAPT::rsRayBouncer<float>;
template class RAPT::rsNoiseGenerator<float>;

// Modulation:
template class RAPT::rsBreakpointModulator<float>;

// Analysis:
template class RAPT::rsAutoCorrelationPitchDetector<double>;
template class RAPT::rsEnvelopeFollower2<double>;

// Visualization:
template class RAPT::rsImage<float>;
template class RAPT::rsImage<rsPixelRGB>;
template class RAPT::rsAlphaMask<float>;
template class RAPT::rsImagePainter<float, float, float>;
template class RAPT::rsImageDrawer<float, float, float>;
template class RAPT::rsLineDrawer<float, float, float>;
template class RAPT::rsImagePlotter<float, double>;
template class RAPT::rsImageProcessor<float>;
template class RAPT::rsImageContourPlotter<float, float>;
template class RAPT::rsPhaseScopeBuffer<float, float, double>;

// Modulators:
//template class RAPT::rsBreakpointModulator<double>; // will be needed, when the class is templatized

// Offline:
template class RAPT::rsResampler<double, double>;

// Unfinished:
template class RAPT::rsTwoBandSplitter<float, float>;
template class RAPT::rsMultiBandSplitter<float, float>;
template class RAPT::rsExponentialEnvelopeMatcher<double>;
template class RAPT::rsPeakPicker<double>;
template class RAPT::rsEnvelopeExtractor<double>;


// misc audio functions
template void RAPT::rsFadeOut(double* buffer, int start, int end);
template void RAPT::rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p);
template double RAPT::rsSineFrequency(double y0, double y1, double y2, double smalll);
template double RAPT::rsSineFrequencyAtCore(const double *x, int N, int n0, double smalll);
template double RAPT::rsSineFrequencyAt(const double *x, int N, int n0, bool refine);
template double RAPT::rsSinePhaseAt(double *x, int N, int n0, double w);
template double RAPT::rsSinePhaseAtViaZeros(double *x, int N, int n0, int precision);
template double RAPT::rsGetShiftForBestMatch(double *x1, double *x2, int N, bool deBias);

template double RAPT::rsCentroid(double *x, int N);
template double RAPT::rsCentroidOfEnergy(double *x, int N);
template double RAPT::rsWindowFunction::raisedCosine(double x, double length, double p);
template double RAPT::rsWindowFunction::exactBlackman(double x, double length, double dummy);
template double RAPT::rsWindowFunction::windowedSinc(double x, double length, double stretch);


template void RAPT::rsSineQuadraturePart(double *x, double *y, int N, double f, double fs,
  bool backward);
template void RAPT::rsSineEnvelopeViaQuadrature(double *x, double *y, int N, double f, double fs,
  double smooth);
template void RAPT::rsSineEnvelopeViaAmpFormula(double *x, double *y, int N, double f, double fs,
  double smooth);
template void RAPT::rsEnvelopedSine(double *y, int N, double f, double fs, double p, double *a);
template void RAPT::rsEnvelopedPhaseCatchSine(double *y, int N, double f, double fs, double p0,
  double pk, int k, double *a, int sweepDirection);
template void RAPT::rsRecreateSine(double *x, double *y, int N, double fx, double fy, double fs,
  double p0, double smooth);
template void RAPT::rsRecreateSineWithPhaseCatch(double *x, double *y, int N, double fx, double fy,
  double fs, double p0, double pk, int k, double smooth, int sweepDirection);
template double RAPT::rsSineShiftAmount(double *x, int N, int n0, double p0, double w);
template double RAPT::rsSineShiftAmount(double *x, int N, int n0, double p0);
template double RAPT::rsEnvelopeMatchOffset(const double* x, int Nx, const double* y, int Ny, int D);
template std::vector<double> RAPT::rsExpDecayTail(const RAPT::rsSinusoidalPartial<double>& partial,
  int spliceIndex, double sampleRate);

template std::vector<double> RAPT::rsFlattenPitch(const double *x, int N, double fs, double ft);

// move to rsFilterAnalyzer:
template double RAPT::analogBiquadMagnitudeSquaredAt(double B0, double B1, double B2, double A0,
  double A1, double A2, double w);


template void RAPT::rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(const double *x, double *y,
  int N, double fc, double bw, double fs, int numPasses, double gc);
template void RAPT::rsBiDirectionalFilter::applyButterworthBandpassBwInHz(const double *x, double *y,
  int N, double fc, double bw, double fs, int order, int numPasses, double gc);


//template class RAPT::rsSinusoidalSynthesizer<double>;
//template class RAPT::rsHarmonicAnalyzer<double>;

//=================================================================================================
// Instantiations of (prototype) classes that are not (yet) in the RAPT namespace

template class rsBivariatePolynomial<double>;
template class rsBivariatePolynomial<std::complex<double>>;

// it's really annoying that we have to instantiate these member functions separately - maybe move 
// their code to the header file - it is short anyway:
template rsPolynomial<double> rsBivariatePolynomial<double>::integralX(
  rsPolynomial<double> a, rsPolynomial<double> b) const;
template rsPolynomial<double> rsBivariatePolynomial<double>::integralY(
  rsPolynomial<double> a, rsPolynomial<double> b) const;
template rsPolynomial<double> rsBivariatePolynomial<double>::integralY(
  double a, rsPolynomial<double> b) const;
template rsPolynomial<double> rsBivariatePolynomial<double>::integralY(
  rsPolynomial<double> a, double b) const;
template rsPolynomial<double> rsBivariatePolynomial<double>::integralX(
  rsPolynomial<double> a, double b) const;
template rsPolynomial<double> rsBivariatePolynomial<double>::integralX(
  double a, rsPolynomial<double> b) const;

template class rsTrivariatePolynomial<double>;
template class rsPiecewisePolynomial<double>;