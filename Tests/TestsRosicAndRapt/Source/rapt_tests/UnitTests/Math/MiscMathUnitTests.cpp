#include "MiscMathUnitTests.h"

bool testMiscMath()
{
  std::string dummy;
  bool testResult = true;

  testResult &= testExponentialCurveFitting(  dummy);
  testResult &= testRootFinding(              dummy);
  testResult &= testGradientBasedOptimization(dummy);
  testResult &= testMinSqrDifFixSum(          dummy);
  testResult &= testPhaseUnwrapStuff(         dummy);

  //testResult &= testMultiLayerPerceptronOld(  dummy); // produces verbose output
  //testResult &= testMultiLayerPerceptron(     dummy); // maybe move to experiments

  return testResult;
}

bool testExponentialCurveFitting(std::string &reportString)
{
  std::string testName = "ExponentialCurveFitting";
  bool testResult = true;

  // example from Hamming's "Numerical Methods ..." page 620
  static const int N = 4;
  static const int k = N/2;
  //double x[N] = { 0.0,  1.0,  2.0,  3.0};
  double y[N] = {32.0, 20.0, 14.0, 11.0};

  // find exponents and weights:
  double a[k];   // exponents
  double A[k];   // weights
  bool success = rsCurveFitter::fitExponentialSum(y, N, A, a, k);
  rsAssert(success);

  // now we should have: y[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n) - verify this:
  double yc[N];
  int n;
  for(n = 0; n < N; n++)
    yc[n] = A[0]*exp(a[0]*n) + A[1]*exp(a[1]*n);
  for(n = 0; n < N; n++)
    testResult &= (yc[n] == y[n]);

  return testResult;
}

bool testRootFinding(std::string &reportString)
{
  std::string testName = "RootFinding";
  bool testResult = true;

  double r;
  double tol = std::numeric_limits<double>::epsilon();

  UnivariateScalarFunctionViaPointer<double> sine(&sin, &cos);

  // Newton iteration:
  r = sine.findRootViaNewtonNonRobust(3.0);
  testResult &= rsIsCloseTo(r, PI, tol);
  r = sine.findRootViaNewtonNonRobust(6.0);
  testResult &= rsIsCloseTo(r, 2*PI, tol);

  // Chebychev method (generalization of Newton iteration using also 2nd derivative):
  // ...stopping criterion still not working - so it's commented:
  //r = sine.findRootViaChebychevNonRobust(3.0);
  //testResult &= rsIsCloseTo(r, PI, tol);
  //r = sine.findRootViaChebychevNonRobust(6.0);
  //testResult &= rsIsCloseTo(r, 2*PI, tol);

  // Ridders' method does not yet converge here, so it's commented
  //r = sine.findRootViaRidders(6.2, 6.3);
  //r = sine.findRootViaRidders(3.1, 3.2);
  //r = sine.findRootViaRidders(2.5, 3.5);

  return testResult;
}

bool testGradientBasedOptimization(std::string &reportString)
{
  std::string testName = "GradientBasedOptimization";
  bool testResult = true;

  // set up the minimizer:
  double tol = 0.00001; // tolerance - later, pass this to the minimizer
  GradientBasedMinimizer<double> minimizer;
  minimizer.setBetaFormula(GradientBasedMinimizer<double>::POLAK_RIBIERE);
  //minimizer.setBetaFormula(GradientBasedMinimizer::FLETCHER_REEVES);
  //minimizer.setBetaFormula(GradientBasedMinimizer::HESTENES_STIEFEL);
  //minimizer.setPrintInfo(true);

  // create the error function object:
  QuadraticTestErrorFunction<double> error;
  double aMin[2]  = { 2, -2};  // x = [ 2, -2] is the desired minimum
  double aInit[2] = {-2, -2};  // x = [-2, -2] is the intial guess
  rsVectorDbl xMin( 2, aMin);
  rsVectorDbl xInit(2, aInit);
  rsVectorDbl xFinal;

  // do the minimization using different algorithms:
  minimizer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  minimizer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  xFinal = minimizer.minimizeFunction(&error, xInit);
  testResult &= (xFinal-xMin).getEuclideanNorm() < tol;

  return testResult;
}

bool testMinSqrDifFixSum(std::string &reportString)
{
  // code moved to experiments - todo: implement actual unit-tests
  std::string testName = "MinSqrDifFixSum";
  bool testResult = true;




  return testResult;
}

bool testPhaseUnwrapStuff(std::string &reportString)  // rename to testUnwrapping
{
  bool r = true;  // test result

  // we consider the range -3..+7 with wrap-around - all numbers ineach column from an equivalence
  // class and numbers in the leftmost column are identified with numbers in the rightmost column 
  // (for example -PI and PI would represent the same angle):
  // -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13   k = -2
  // -13 -12 -11 -10 -09 -08 -07 -06 -05 -04 -03   k = -1
  // -03 -02 -01 +00 +01 +02 +03 +04 +05 +06 +07   k =  0 (base range)
  // +07 +08 +09 +10 +11 +12 +13 +14 +15 +16 +17   k = +1
  // +17 +18 +19 +20 +21 +22 +23 +24 +25 +26 +27   k = +2
   

  double rangeMin = -3.0;
  double rangeMax = +7.0;

  double val;



  val = rsFindCosistentUnwrappedValue(-10.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 9.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 8.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 7.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 6.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 5.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 4.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 3.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 2.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;
  val = rsFindCosistentUnwrappedValue(- 1.0, 5.0, rangeMin, rangeMax); r &= val == -5.0;

  val = rsFindCosistentUnwrappedValue(  0.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  1.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  2.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  3.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  4.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  5.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  6.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  7.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  8.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue(  9.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;
  val = rsFindCosistentUnwrappedValue( 10.0, 5.0, rangeMin, rangeMax); r &= val ==  5.0;

  val = rsFindCosistentUnwrappedValue(11.0, 5.0, rangeMin, rangeMax); r &= val == 15.0;
  // hangs


  // tests various functions that have to do with phase-unwrapping

  // todo: test wrapped interpolation, rsWarpTointerval, unwrap, ...


  return r;
}


bool testMultiLayerPerceptronOld(std::string &reportString)
{
  std::string testName = "MultiLayerPerceptron";
  bool testResult = true;

  // this should actually be in the "Experiments" suite

  //int hiddenNeuronsInLayers[2] = {2, 4};
  //MultiLayerPerceptron mlp(3, 2, 2, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[3] = {2, 4, 3};
  //MultiLayerPerceptron mlp(1, 1, 3, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[3] = {4, 3, 5};
  //MultiLayerPerceptron mlp(1, 1, 3, hiddenNeuronsInLayers);

  int hiddenNeuronsInLayers[1] = {10};
  MultiLayerPerceptron<double> mlp(1, 1, 1, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[2] = {3, 5};
  //MultiLayerPerceptron mlp(1, 1, 2, hiddenNeuronsInLayers);

  //int hiddenNeuronsInLayers[2] = {2,4};
  //MultiLayerPerceptron mlp(1, 1, 2, hiddenNeuronsInLayers);

  // create the error-function object and pass some training data to it:
  MultiLayerPerceptronErrorFunction<double> mlpError(&mlp);
  static const int N = 200;
  rsVectorDbl x[N];
  rsVectorDbl y[N];
  for(int n=0; n<N; n++)
  {
    x[n].setDimensionality(1);
    y[n].setDimensionality(1);
  }
  x[0].v[0] = rsRandomUniform(-3.0, 3.0, 7);
  for(int n=0; n<N; n++)
  {
    x[n].v[0] = rsRandomUniform(-3.0, 3.0);
    //x[n].v[0] = (double) n / (double) (N-1);
    //x[n].v[0] = linToLin(x[n].v[0], 0.0, 1.0, -3.0, 3.0);
    y[n].v[0] = x[n].v[0] * x[n].v[0];            // y = x^2
    //y[n].v[0] = fabs(x[n].v[0]);                    // y = |x|
    //y[n].v[0] = sign(x[n].v[0]);                    // y = sign(x)
    //y[n].v[0] = sin(x[n].v[0] * x[n].v[0]);            // y = sin(x^2)
  }
  mlpError.setTrainingData(x, y, N);

  // create the minimizer and minimize the training error:
  GradientBasedMinimizer<double> mlpTrainer;
  mlpTrainer.setPrintInfo(true);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::POLAK_RIBIERE);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::FLETCHER_REEVES);
  //mlpTrainer.setBetaFormula(GradientBasedMinimizer<double>::HESTENES_STIEFEL);
  rsVectorDbl w = mlpTrainer.minimizeFunction(&mlpError, mlp.getWeightsAsVector());
  mlp.setWeightVector(w);

  int dummy = 0;


  //double *test = new double[0];
  //test[0] = 1.0;

  //double x[3]  = {3,2,1};
  //double yt[2] = {1,-1};
  //double y[2];
  //mlp.computeNetworkOutput(x, y);
  //mlp.printWeights();
  //mlp.printActivations();

  //mlpError.computePatternGradientByWeightPerturbation(Vector(2, yt));
  //mlpError.printPatternGradient();

  //mlpError.computePatternGradient(Vector(2, yt));
  //mlpError.printPatternGradient();



  //Vector wv = mlp.getWeightsAsVector();
  //wv.print();
  //mlp.setWeightVector(wv);
  //mlp.printWeights();
  //Vector gv = mlpTrainer.getGradientVector();
  //gv.print();


  double xTest[N];
  double yTest[N];
  for(int n=0; n<N; n++)
  {
    xTest[n] = (double) n / (double) (N-1);
    xTest[n] = rsLinToLin(xTest[n], 0.0, 1.0, -6.0, 6.0);
    mlp.computeNetworkOutput(&xTest[n], &yTest[n]);
    rsVectorDbl yVec = mlp.getOutput();
    int dummy = 0;
  }
  //Plotter::plotData(N, xTest, yTest);
    
  printf("%s", "\n Press Key \n");
  getchar();

  return testResult;
}


rsVectorDbl mlpTestFunction(rsVectorDbl x, double noiseAmplitude)
{
  rsAssert(x.dim == 3);
  rsVectorDbl y(2);

  // compute output values according to deterministic function:
  y[0]  =  2*x[0] + 3*x[1] - 1*x[2] - 2*x[0]*x[1] + 1*x[0]*x[2] + 3*x[1]*x[2];
  y[1]  = -3*x[0] + 1*x[1] + 2*x[2] + 3*x[0]*x[1] - 2*x[0]*x[2] + 1*x[1]*x[2];

  // add some noise:
  y[0] += noiseAmplitude * rsRandomUniform(-1.0, 1.0);
  y[1] += noiseAmplitude * rsRandomUniform(-1.0, 1.0);

  return y;
}
rsVectorDbl randomVector(int numDimensions, double min, double max)
{
  rsVectorDbl x(numDimensions);
  for(int i = 0; i < numDimensions; i++)
    x[i] = rsRandomUniform(min, max);
  return x;
}
bool testMultiLayerPerceptron(std::string &reportString)
{
  std::string testName = "MultiLayerPerceptron";
  bool testResult = true;

  // we try to approximate the function f: R^3 -> R^2, defined by:
  // y1 =  2*x1 + 3*x2 - 1*x3 - 2*x1*x2 + 1*x1*x3 + 3*x2*x3
  // y2 = -3*x1 + 1*x2 + 2*x3 + 3*x1*x2 - 2*x1*x3 + 1*x2*x3
  // by using noisy datapoints

  // ...that particular function seems to be quite difficult to approximate by an MLP
  // maybe use a function that is easier to represent - possibly by using outputs of an MLP
  // itself (with randomized weights) - see if the weights can be retrieved by training (or an 
  // equivalent weight-vector - there are weight-space symmetries)

  // create test data:
  int n;
  static const int numPatterns    = 200;
  double noiseAmplitude = 0.2;
  rsVectorDbl x[numPatterns];
  rsVectorDbl y[numPatterns];
  for(n = 0; n < numPatterns; n++)
  {
    x[n].setDimensionality(3);
    y[n].setDimensionality(2);
  }
  rsRandomUniform(-1.0, 1.0, 7); // init PRNG
  for(n = 0; n < numPatterns; n++)
  {
    x[n] = randomVector(3, -2.0, 2.0);
    y[n] = mlpTestFunction(x[n], noiseAmplitude);
  }

  // create an MLP with 3 hidden layers with 4, 3 and 5 neurons respectively:
  //int hiddenNeuronsInLayers[3] = {4, 3, 5};
  //MultiLayerPerceptron mlp(3, 2, 3, hiddenNeuronsInLayers);

  int hiddenNeuronsInLayers[1] = {10};
  MultiLayerPerceptron<double> mlp(3, 2, 1, hiddenNeuronsInLayers);

  // create the error-function object and pass the training data to it:
  MultiLayerPerceptronErrorFunction<double> mlpError(&mlp);
  mlpError.setTrainingData(x, y, numPatterns);

  // create the minimizer for the error-function find the weights that minimize the training error:
  GradientBasedMinimizer<double> mlpTrainer;
  mlpTrainer.setPrintInfo(true);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::GRADIENT_DESCENT);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::BOLD_DRIVER_WITH_MOMENTUM);
  //mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::CONJUGATE_GRADIENT);
  mlpTrainer.setAlgorithm(GradientBasedMinimizer<double>::SCALED_CONJUGATE_GRADIENT);
  rsVectorDbl w = mlpTrainer.minimizeFunction(&mlpError, mlp.getWeightsAsVector());
  mlp.setWeightVector(w);  // set up the network with the optimal weight-vector

  return testResult;
}