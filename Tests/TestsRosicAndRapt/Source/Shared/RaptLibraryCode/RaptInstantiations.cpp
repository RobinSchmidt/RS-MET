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

#include "../../../../../Libraries/JUCE/modules/rapt/rapt.cpp"
//#include "../../../../../Libraries/JUCE/modules/rosic/rosic.cpp"
using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

// we need to make sure to provide only instantiations that are not already there in rosic.cpp
// todo: move instantiations to rosic, where appropriate

typedef std::complex<double> cmplxD;

// Basics:


template bool RAPT::defaultLess(const int& left, const int& right);
template void RAPT::rsHeapSort(int *buffer, int length, 
  bool (*lessThen)(const int& left, const int& right));
template bool RAPT::rsIsSortedAscending(int *buffer, int length);
template std::vector<int> RAPT::rsFindAllOccurencesOf(char* buffer, int bufferLength, 
  char* pattern, int patternLength);
template int RAPT::rsFindHighestPeakIndex(double *buffer, int length); // move to rsArray
// maybe move these into rsArray

// Data:

template void rsArray::fillWithRangeLinear(double* x, int N, double min, double max);
template void rsArray::fillWithRangeExponential(double* x, int N, double min, double max);
template void rsArray::fillWithRandomValues(double* x, int N, double min, double max, int seed);
template void rsArray::deConvolve(double *y, int yLength, double *h, int hLength, double *x);
template void rsArray::sequenceSqrt(double *y, int yLength, double *x);
template void rsArray::transposeSquareArray(double **in, double **out, int size);
template void rsArray::cumulativeSum(double *x, double *y, int length, int order);
template void rsArray::difference(double *buffer, int length, int order, bool periodic);
template double rsArray::minValue(double *x, int length);
template double rsArray::maxValue(double *x, int length);
template double rsArray::maxAbs(double *x, int length);
template void rsArray::normalize(double *buffer, int length, double maximum, bool subtractMean);
template double rsArray::mean(double *x, int length);
template double rsArray::maxDeviation(double *buffer1, double *buffer2, int length);
template void rsArray::applyFunction(double *x, double *y, int N, double (*f) (double));
template void rsArray::negate(double *source, double *destination, int length);
template void rsArray::filter(double *x, int xLength, double *y, int yLength, 
  double *b, int bOrder, double *a, int aOrder);
template void rsArray::filterBiDirectional(double *x, int xLength, double *y, int yLength, 
  double *b, int bOrder, double *a, int aOrder, int numRingOutSamples);
template void rsArray::fillWithRangeLinear(float* x, int N, float min, float max);
template void rsArray::fillWithRandomValues(float* x, int N, double min, double max, int seed); // ?
template void rsArray::divide(double *buffer1, double *buffer2, double *result, int length);
template int rsArray::maxAbsIndex(const double* const buffer, int length);
template void rsArray::reverse(double* x, double* y, int length);
template void rsArray::rightShift(double *buffer, int length, int numPlaces);
template void rsArray::shift(double *buffer, int length, int numPlaces);
template void rsArray::unwrap(double* a, int N, double p);

template float rsArray::maxDeviation(float *buffer1, float *buffer2, int length);

template void rsArray::fillWithRangeLinear(int* x, int N, int min, int max);
template void rsArray::fillWithRandomValues(int* x, int N, double min, double max, int seed);
template void rsArray::allocateSquareArray2D(int**& theArray, int size);
template void rsArray::deAllocateSquareArray2D(int**& theArray, int size);
template void rsArray::rightShift(int *buffer, int length, int numPlaces);
template void rsArray::copySection(int *source, int sourceLength, int *destination, int copyStart, int copyLength);
template void rsArray::cumulativeSum(int *x, int *y, int length, int order);

template rsUint32 rsArray::maxValue(rsUint32 *x, int length);
template void rsArray::fillWithRangeLinear(rsUint32* x, int N, rsUint32 min, rsUint32 max);
template void rsArray::copyBuffer(const rsUint32 *src, rsUint32 *dst, int N);
template int rsArray::firstIndexWithNonZeroValue(rsUint32 *a, int N);
template bool rsArray::contains(rsUint32 *buffer, int length, rsUint32 elementToFind);
template void rsArray::fillWithZeros(rsUint32 *buffer, int length);

template void rsArray::convertBuffer(std::complex<double> *source, std::complex<double> *destination, int length);
template void rsArray::fillWithRandomValues(std::complex<double>* x, int N, double min, double max, int seed);






template void MatrixTools::rsInitMatrix(double** A, int N, int M, double value);
template void MatrixTools::rsCopyMatrix(double** source, double **destination, int N, int M);
template bool MatrixTools::rsAreMatricesApproximatelyEqual(double **A, double **B, int N, int M, double tol);

template void MatrixTools::rsMatrixMultiply(double **A, double **B, double **C, int N, int M, int P);
template void MatrixTools::rsMatrixMultiplyFirstTransposed(double **A, double **B, double **C, int N, int M, int P);
template void MatrixTools::rsMatrixMultiplySecondTransposed(double **A, double **B, double **C, int N, int M, int P);
template void MatrixTools::rsMatrixMultiplyBothTransposed(double **A, double **B, double **C, int N, int M, int P);
template void MatrixTools::rsMatrixInPlaceMultiply(double **A, double **B, int N, int M);



// Math:
//template RAPT::rsLinearAlgebra<float>; // doens't work bcs the template parameters are decalred in the member functions
template void rsLinearAlgebra::rsSolveLinearSystem2x2(const double A[2][2], double x[2], const double y[2]);
template void rsLinearAlgebra::rsSolveLinearSystem3x3(const double A[3][3], double x[3], const double y[3]);
template bool rsLinearAlgebra::rsSolveLinearSystem(double **A, double *x, double *b, int N);
template bool rsLinearAlgebra::rsInvertMatrix(double **A, int N);
template bool rsLinearAlgebra::rsSolveTridiagonalSystem(double *lowerDiagonal, double *mainDiagonal, 
  double *upperDiagonal, double *rightHandSide, double *solution, int N);
template bool rsLinearAlgebra::rsChangeOfBasisColumnWise(double **A, double **B, double *va, 
  double *vb, int N);
template bool rsLinearAlgebra::rsChangeOfBasisRowWise(double **A, double **B, double *va, 
  double *vb, int N);
template bool rsLinearAlgebra::rsChangeOfBasisMatrixColumnWise(double **A, double **B, double **C, int N);
template bool rsLinearAlgebra::rsChangeOfBasisMatrixRowWise(   double **A, double **B, double **C, int N);
template bool rsLinearAlgebra::rsSolveLinearSystem(cmplxD **A, cmplxD *x, cmplxD *b, int N);


template RAPT::rsMatrix<double>;

template RAPT::rsPolynomial<float>;
template RAPT::rsPolynomial<double>;
//template RAPT::rsPolynomial<int>; // template doesn't compile with int
template void RAPT::rsPolynomial<double>::dividePolynomialByMonomialInPlace(double*, int, double, double*);
  // needs separate instantiation because function itself has a (second) template parameter

// should go into an "Interpolation" class...or maybe CurveFitter
//template void RAPT::fitCubicWithDerivative(double x1, double x2, double y1, double y2, double yd1, 
//  double yd2, double *a3, double *a2, double *a1, double *a0);
template void RAPT::fitCubicWithDerivativeFixedX(double y0, double y1, double yd0, double yd1, 
  double *a3, double *a2, double *a1, double *a0);


template RAPT::UnivariateScalarFunction<double>;
template RAPT::UnivariateScalarFunctionViaPointer<double>;
template RAPT::QuadraticTestErrorFunction<double>;

template RAPT::GradientBasedMinimizer<double>;

template RAPT::MultiLayerPerceptron<double>;
template RAPT::MultiLayerPerceptronErrorFunction<double>;
template RAPT::MultiLayerPerceptronTrainer<double>;

template bool RAPT::rsFitSumOfExponentials(double* y, int numValues, double* A, double* a, 
  int numExponentials);

template rsUint32 RAPT::rsBinomialCoefficient(rsUint32 n, rsUint32 k);
template rsUint32 RAPT::rsBinomialCoefficientUpTo20(rsUint32 n, rsUint32 k);
template int RAPT::rsMultinomialCoefficient(int* n, int k);
template rsUint32 RAPT::rsMultinomialCoefficient(rsUint32* n, rsUint32 k);
template rsUint32 RAPT::rsMultinomialCoefficientUpTo12(rsUint32* n, rsUint32 k);
template int RAPT::rsLeviCivita(int indices[], int N);

template void RAPT::rsStirlingNumbersFirstKind(int **s, int nMax);


template void RAPT::rsFillPrimeTable(int*, rsUint32, rsUint32);
template void RAPT::rsFindPrimesUpTo(std::vector<rsUint32>& primes, rsUint32 upperLimit);
template void RAPT::rsFindPrimesUpTo(std::vector<rsUint64>& primes, rsUint64 upperLimit);
template void RAPT::rsPrimeFactors(rsUint32 x, std::vector<rsUint32>& factors, 
  std::vector<rsUint32>& exponents, std::vector<rsUint32>* primeTable);

template void RAPT::rsEGCD(int x, int y, int& a, int& b, int& g);
template int RAPT::rsModularInverse(const int& x, const int& m);
template int RAPT::rsPrimeModularInverse(const int& x, const int& p);
template int RAPT::rsPrimeModularInverse2(const int& x, const int& p);
template int RAPT::rsChineseRemainderTheorem(int* remainders, int* moduli, rsUint32 count);
//template rsUint32 RAPT::rsChineseRemainderTheorem(rsUint32* remainders, rsUint32* moduli, rsUint32 count);
template void RAPT::rsFillPrimeTable(rsUint32 *primes, rsUint32 numPrimes, rsUint32 bufferSize);

template void RAPT::smbFft(float *fftBuffer, long fftFrameSize, long sign);
template void RAPT::rsDFT(std::complex<double> *buffer, int N);
template void RAPT::rsFFT(std::complex<double> *buffer, int N);
template void RAPT::rsIFFT(std::complex<double> *buffer, int N);
template void RAPT::rsMagnitudeAndPhase(double *signal, int N, double *magnitudes, double *phases);
template void RAPT::rsCrossCorrelation(float x[], float y[], int N, float r[], bool removeBias);
template void RAPT::rsAutoCorrelationFFT(float x[], int N, float r[]);




template RAPT::rsSineIterator<double>;
template RAPT::rsComplexExponentialIterator<double>;

template double RAPT::rsNormalizedSinc(double x);
template double RAPT::rsSineIntegral(double x);

template int RAPT::rsClip(int x, int min, int max);
template int RAPT::rsFactorial(int n);

//template double RAPT::rsWrapToInterval(double x, double min, double max); // rename to rsWrap

template void RAPT::rsBinomialDistribution(double*, int, double);


template RAPT::rsConicSection<float>;
template RAPT::rsRotationXY<double>;

template RAPT::rsCoordinateMapper<float>;
template RAPT::rsCoordinateMapper<double>;
template RAPT::rsCoordinateMapper2D<float>;
template RAPT::rsCoordinateMapper2D<double>;

template RAPT::rsEllipse<float>;
template RAPT::rsLine2D<float>;
template RAPT::rsNormalizedSigmoids<float>;
template RAPT::rsParametricSigmoid<float>;
template RAPT::rsScaledAndShiftedSigmoid<float>;
template RAPT::rsPositiveBellFunctions<float>;
template RAPT::rsParametricBellFunction<float>;
template RAPT::rsRootFinder<float>;
template RAPT::rsNodeBasedFunction<float>;
template RAPT::rsSinCosTable<float>;
template RAPT::rsSinCosTable<double>;
template RAPT::rsVector2D<float>;
template RAPT::rsVector3D<float>;

template RAPT::rsFourierTransformerRadix2<double>;

template void RAPT::rsStatistics::linearRegression(int N, float* x, float* y, float& a, float& b);
template float RAPT::rsStatistics::proportionalRegression(int N, float* x, float* y);
template void RAPT::rsRemoveCorrelationBias(double x[], int N, double r[]);


template RAPT::rsMultiArray<float>; 

template double RAPT::rsBandwidthConverter::bandedgesToCenterFrequency(double fl, double fu);
template double RAPT::rsBandwidthConverter::bandedgesToAbsoluteBandwidth(double fl, double fu);
template double RAPT::rsBandwidthConverter::absoluteBandwidthToQ(double fl, double fu);

// Filters-Musical:
template RAPT::rsSmoothingFilter<float, float>;
template RAPT::rsLadderFilter<float, float>;
template RAPT::rsLadderFilter<double, double>;
//template RAPT::rsLadderFilter<rsFloat64x2, double>;
//template RAPT::rsLadderFilter<rsFloat64x2, rsFloat64x2>;
template RAPT::rsPhasorFilter<float, float>;
template RAPT::rsPhasorStateMapper<float>;
template RAPT::rsStateVariableFilter<float, float>; 

// Filters-Scientific:
template RAPT::rsPrototypeDesigner<float>;
//template RAPT::rsPrototypeDesigner<double>;
template RAPT::rsPoleZeroMapper<float>;
template RAPT::rsFilterCoefficientConverter<float>;
template RAPT::rsInfiniteImpulseResponseDesigner<float>;
template RAPT::rsFilterAnalyzer<float>;
template RAPT::rsBiquadCascade<float, float>;
template RAPT::rsEngineersFilter<float, float>;
template RAPT::rsLinkwitzRileyCrossOver<float, float>;
template RAPT::rsCrossOver4Way<float, float>;
template RAPT::rsDirectFormFilter<float, float>;
template RAPT::rsEllipticSubBandFilter<float, float>;
template RAPT::rsEllipticSubBandFilterDirectForm<float, float>;
template RAPT::rsQuadratureNetwork<float, float>;

// Physics:
template RAPT::rsParticleSystem<float>;

// Generators:
template RAPT::rsBouncillator<float>;
template RAPT::rsRayBouncer<float>;

// Modulation:
template RAPT::rsBreakpointModulator<float>;

// Analysis:
template RAPT::rsAutoCorrelationPitchDetector<double>;

// Visualization:
template RAPT::rsImage<float>;
template RAPT::rsAlphaMask<float>;
template RAPT::rsImagePainter<float, float, float>;
template RAPT::rsLineDrawer<float, float, float>;
template RAPT::rsPhaseScopeBuffer<float, float, double>;

// Modulators:
//template RAPT::rsBreakpointModulator<double>; // will be needed, when the class is templatized

// Offline:
template RAPT::rsResampler<double, double>;

// Unfinished:
template RAPT::rsTwoBandSplitter<float, float>;
template RAPT::rsMultiBandSplitter<float, float>;

// misc audio functions
template void RAPT::rsFadeOut(double* buffer, int start, int end);
template void RAPT::rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p);
template double RAPT::rsSineFrequency(double y0, double y1, double y2, double smalll);
template double RAPT::rsSineFrequencyAtCore(double *x, int N, int n0, double smalll);
template double RAPT::rsSineFrequencyAt(double *x, int N, int n0, bool refine);
template double RAPT::rsSinePhaseAt(double *x, int N, int n0, double w);
template double RAPT::rsSinePhaseAtViaZeros(double *x, int N, int n0, int precision);
template double RAPT::rsGetShiftForBestMatch(double *x1, double *x2, int N, bool deBias);

template double RAPT::rsCentroid(double *x, int N);
template double RAPT::rsCentroidOfEnergy(double *x, int N);
template double RAPT::rsRaisedCosineWindow(double x, double length, double p);
template double RAPT::rsExactBlackmanWindow(double x, double length, double dummy);
template double RAPT::rsWindowedSinc(double x, double length, double stretch);


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

// move to rsFilterAnalyzer:
template double RAPT::analogBiquadMagnitudeSquaredAt(double B0, double B1, double B2, double A0, 
  double A1, double A2, double w);


template void RAPT::rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(double *x, double *y, 
  int N, double fc, double bw, double fs, int numPasses, double gc);
template void RAPT::rsBiDirectionalFilter::applyButterworthBandpassBwInHz(double *x, double *y, 
  int N, double fc, double bw, double fs, int order, int numPasses, double gc);