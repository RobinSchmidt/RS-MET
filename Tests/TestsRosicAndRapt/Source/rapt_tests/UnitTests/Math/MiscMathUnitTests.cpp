#include "MiscMathUnitTests.h"

bool testMiscMath()
{
  std::string dummy;
  bool testResult = true;

  testResult &= testExponentialCurveFitting(  dummy);
  testResult &= testRootFinding(              dummy);
  testResult &= testGradientBasedOptimization(dummy);
  testResult &= testMinSqrDifFixSum(          dummy);

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

// for production code later use plain c-arrays..
// N: length of v, v: computed output values, s: desired sums, w: weights for the squared 
// differences (optional, if not given, all weights are unity)
//template<class T>
//void rsMinSqrDifFixSum(int N, T* v, const T* s, const T* w = nullptr)
//{
//
//}


bool testMinSqrDifFixSum(std::string &reportString)
{
  // Test minimization of the sum of the squares of the differences between adjacent array elements
  // when their sums are given fixed values.

  // not yet complete - maybe move to experiments

  std::string testName = "MinSqrDifFixSum";
  bool testResult = true;

  std::vector<double> v, s, w; // values, sums and weights



  s = { 12, 24, 36 };  
  v = rsMinSqrDifFixSum(s);       // v = 4, 8, 16, 20, d = 4,4,4

  s = { 12, 24, 36, 48 };  
  v = rsMinSqrDifFixSum(s);       // v = 3, 9, 15, 21, 27, d = 6,6,6,6

  w = { 2, 3, 4, 5 };
  v = rsMinSqrDifFixSum(s, w); 

  s = { 12, 24, 36, 48, 60 };  
  w = {  2,  3,  4,  5,  6 };
  v = rsMinSqrDifFixSum(s, w);   
  // 3.6, 8.4, 15.6, 20.4, 27.6, 32.4 -> d = 4.8, 7.2, 4.8, 7.2, 4.8 -> ?
  // the result seems wrong - the increasing weights should make the differences between adjacent
  // elements go down for later values...instead, they alternate - wtf?


  s = { 20 };        
  v = rsMinSqrDifFixSum(s);       // 10, 10
  s = { 20, 30 };    
  v = rsMinSqrDifFixSum(s);       // 7.5, 12.5, 17.5
  s = { 20, 30, 40 };  
  v = rsMinSqrDifFixSum(s);       // 8.33, 11.66, 18.33, 21.66
  s = { 20, 30, 40, 50 }; 
  v = rsMinSqrDifFixSum(s);       // 7.5, 12.5, 17.5, 22.5, 27.5
  s = { 20, 30, 40, 50, 60 }; 
  v = rsMinSqrDifFixSum(s);       // 8, 12, 18, 22, 28, 32
  s = { 20, 30, 40, 50, 60, 70 }; 
  v = rsMinSqrDifFixSum(s);       // 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5
  s = { 20, 30, 40, 50, 60, 70, 80 }; 
  v = rsMinSqrDifFixSum(s);       // 7.86, 12.14, 17.86, 22.14, 27.86, 32.14, 37.86, 42.14
  return testResult;



  // maybe when the length of the sum vector s is odd (i.e. length of v is even) we need to do 
  // something special? the computed solution looks suboptimal - calculate N=4 case by hand - see,
  // if the matrix looks different - use v1,v2,v3,v4,w1,w2,w3
  // it seems, the longer the vector, the less strong the suboptimality...maybe it sort of averages
  // out for longer vectors...figure out

  // no - for s = { 20, 30, 40 } the result v = { 8.33, 11.66, 18.33, 21.66 } gives indeed a lower
  // sum-of-squared-differences than for example v = { 7.5, 12.5, 17.5, 22.5 } - so the math works
  // out correctly. maybe my minimization criterion is not well suited for even N - maybe i should
  // give the two outermost squared differences less weight (maybe 1/2)

  // probably the best and most flexible way to deal with this is to give the user the option to 
  // apply a weight to each squared difference - we will need that anyway for the freq-estimation.
  // is there a meaningful way to apply weights to the constraints as well? or rather some sort of
  // slack-variables maybe instead of requiring the hard constraint 
  // v1 + v2 = s1 -> v1 + v2 - s1 = 0, require k * (v1 + v2 - s1)^2 = min ...but i think, for the 
  // frequency-estimation problem, hard constraints are actually appropriate - but maybe for other
  // problems...wait...that makes the whole thing nonlinear - so screw it

  // todo: make production version of the function that includes weights for the squared 
  // differences, maybe make a version that uses the simple penta-solver and one that uses 
  // rsbandDiagonalSolver  - maybe call the simple one "Fast"

  // to figure out, if the results are really correct, verify a couple of cases for small N 
  // symbolically in sage

  // maybe instead of minimizing the squared differences, it may also make sense to minimize the 
  // squared (numeric) 2nd derivative or a mix of the two or the squared difference of the squared
  // difference

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