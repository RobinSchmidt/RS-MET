using namespace RAPT;
using namespace std;

void ellipseLineIntersections()
{
  // create and set up ellipse:
  //rsEllipseF ellipse;
  rsEllipse<float> ellipse;
  ellipse.setParameters(1.5f, 3.f, float(PI/4), 0.5f, 1.f); 

  // create line parameters:
  float x, y, dx, dy;
  x  = 0.5f;
  y  = 0.5f;
  dx = 0.5f;
  dy = 0.2f;

  // create data for drawing the ellipse:
  static const int Ne = 500;
  float xe[Ne], ye[Ne];
  float err;
  for(int n = 0; n < Ne; n++)
  {
    float phi = float(2*PI*n) / (Ne-1);
    ellipse.getPointOnEllipse(phi, &xe[n], &ye[n]);
    err = ellipse.evaluate(xe[n], ye[n]);
  }

  // create data for drawing the line:
  float tMin = -3.5;
  float tMax = +3.0;
  float xl[2], yl[2];
  xl[0] = x + tMin*dx;
  xl[1] = x + tMax*dx;
  yl[0] = y + tMin*dy;
  yl[1] = y + tMax*dy;

  // find intersection points between line and ellipse:
  float ti1, ti2, xi1, yi1, xi2, yi2;
  ellipse.lineIntersectionParameter(x, dx, y, dy, &ti1, &ti2);
  xi1 = x + ti1*dx;
  yi1 = y + ti1*dy;
  xi2 = x + ti2*dx;
  yi2 = y + ti2*dy;

  // find tangent line to intersection point:
  float A, B, C, a, b;
  float xt1[2], yt1[2], xt2[2], yt2[2];
  ellipse.getTangentCoeffs(xi1, yi1, &A, &B, &C); // implicit  A*x + B*y + C = 0
  a = -A/B;     // explicit y = a*x + b
  b = -C/B;
  xt1[0] = -1.2f;   // use points where x=0, x=2 to draw the tangent:
  yt1[0] = a*xt1[0] + b;   
  xt1[1] = -0.5f;
  yt1[1] = a*xt1[1] + b; 

  // same for the 2nd intersection:
  ellipse.getTangentCoeffs(xi2, yi2, &A, &B, &C);
  a = -A/B;
  b = -C/B;
  xt2[0] = 0;
  yt2[0] = a*xt2[0] + b;   
  xt2[1] = 2;
  yt2[1] = a*xt2[1] + b; 

  GNUPlotter plt;
  plt.setRange(-1.5, +2.5, -1, +3);
  plt.setPixelSize(600, 600);
  plt.addCommand("set size square");   // set aspect ratio to 1:1 ..encapsulate in GNUPlotter
  plt.addDataArrays(Ne, xe,  ye);   // ellipse
  plt.addDataArrays(2,  xl,  yl);   // line
  plt.addDataArrays(2,  xt1, yt1);  // tangent at 1st intersection
  plt.addDataArrays(2,  xt2, yt2);  // tangent at 2nd intersection
  plt.plot();
}

void finiteDifferenceStencilCoeffs()
{
  // Computation of coefficients for arbitrary finite difference stencils, see:
  // http://web.media.mit.edu/~crtaylor/calculator.html

  typedef std::vector<double> Vec;
  Vec s = {-2, -1, 0, 1, 2};         // normalized stencil offsets
  int d = 4;                         // order of derivative to be approximated

  // establish matrix:
  int N = (int) s.size();    // stencil length
  double** A;                // matrix data
  rsMatrixTools::allocateMatrix(A, N, N);
  for(int i = 0; i < N; i++)
    for(int j = 0; j < N; j++)
      A[i][j] = pow(s[j], i);
  //rsMatrix<double> A_dbg(N, N, A);  // for debug

  // establish right-hand-side vector:
  Vec rhs(N);
  rsFill(rhs, 0.0);
  rhs[d] = rsFactorial(d);

  // compute coeffs by solving the linear system:
  Vec c(N);
  rsLinearAlgebra::rsSolveLinearSystem(A, &c[0], &rhs[0], N);
  // In practice, the resulting coefficients have to be divided by h^d where h is the step-size and
  // d is the order of the derivative to be approximated. The stencil values in s are actually 
  // multipliers for some basic step-size h, i.e. a stencil -2,-1,0,1,2 means that we use values
  // f(x-2h),f(x-h),f(x),f(x+h),f(x+2h) to approximate the d-th derivative of f(x) at x=0

  // todo: 
  // -move to library
  // -add unit test
  // -try the example with the 4-th order derivative and 5-point stencil that is presented
  //  on the website
  // -try different examples and compare results with results from the website - use also
  //  asymmetrical and/or non-equidistant stencils
  // -if it all works, maybe implement it also in sage to get rid of roundoff errors
  // -maybe we should round the final coeffs? are they supposed to be integer? ...maybe only, if
  //  the stencil offsets are all integers?

  // stencil = -2,-1,0,1,2, d = 4:

  // stencil = -2,-1,0,3,4, d = 3:
  // f_xxx = (-6*f[i-2]+15*f[i-1]-10*f[i+0]+1*f[i+3]+0*f[i+4])/(10*1.0*h**3)

  // -2,-1,0,0.5,2, d=3
  // f_xxx = (-27*f[i-2]+40*f[i-1]+90*f[i+0]-128*f[i+0.5]+25*f[i+2])/(60*1.0*h**3)

  // -2,-1,0,0.5:
  // f_xx = (-1*f[i-2]+10*f[i-1]-25*f[i+0]+16*f[i+0.5])/(5*1.0*h**2)
  // f_xxx = (-6*f[i-2]+20*f[i-1]-30*f[i+0]+16*f[i+0.5])/(5*1.0*h**3)

  rsMatrixTools::deallocateMatrix(A, N, N);
}

void interpolatingFunction()
{
  typedef RAPT::rsInterpolatingFunction<float, double> IF;

  // create data to interpolate:
  int N = 5;
  float  x[5] = { 2, 4, 5, 7, 8 };
  double y[5] = { 1, 3, 1, 2, 3 };

  // create and set up interpolating function object:
  IF intFunc;
  //intFunc.setMode(IF::LINEAR);
  intFunc.setMode(IF::CUBIC_HERMITE);
  intFunc.setPreMap( &log);
  intFunc.setPostMap(&exp);
  intFunc.setPreMap( nullptr);
  intFunc.setPostMap(nullptr);


  // do extra/interpolation:
  static const int M = 500; // number of interpolated values
  float  xi[M];  
  double yi[M];
  float  xiMin = 0;
  float  xiMax = 10;
  RAPT::rsArrayTools::fillWithRangeLinear(xi, M, xiMin, xiMax);
  intFunc.interpolate(x, y, N, xi, yi, M);

  // convert xi to double for plotter and plot:
  double xid[M];
  RAPT::rsArrayTools::convert(xi, xid, M);
  GNUPlotter plt;
  plt.addDataArrays(M, xid, yi);
  plt.setRange(xiMin, xiMax, 0.0, 4.0);
  plt.plot();
  // todo: plot the original data as points
  // todo: plot all the different possible interpolants in one plot for comparison
}

void linearRegression()
{
  static const int N = 500;  // number of samples
  float minDist = 0.1f;       // minimum distance between successive samples
  float maxDist = 1.0f;      // maximum ...
  float a       = 0.8f;
  float b       = 20.0f;
  float noise   = 10.0f;
  int   seed    = 0;

  // create input data:
  float x[N], y[N];
  float dx;      // delta between samples
  x[0] = 10.f;
  y[0] = a*x[0] + b + noise * (float)round(rsRandomUniform(-1, +1, seed)); 
  int n;
  for(n = 1; n < N; n++){
    dx = (float)rsRandomUniform(minDist, maxDist);
    x[n] = x[n-1] + dx;
    y[n] = a*x[n] + b + noise * (float)rsRandomUniform(-1, +1); 
  }

  // retrieve a,b from the data:
  float ar, br = 0;
  rsStatistics::linearRegression(N, x, y, ar, br);
  //ar = Statistics::proportionalRegression(N, x, y);

  // create lines for correct and retrieved line parameters:
  float yc[N], yr[N];
  for(n = 0; n < N; n++){
    yc[n] = a *x[n] + b;
    yr[n] = ar*x[n] + br; }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0], &yc[0], &yr[0]);
  plt.plot();
}

void multipleRegression()
{
  // under construction - we use this example:
  // https://www.youtube.com/watch?v=vvay3t19HU8&list=PLb0zKSynM2PCmp5J5LWM3PcZXBaCoQkXj&index=40
  // The price of pizza is assumed to have a linear relationship with the diameter and the number
  // of extra ingredients, so the model is y = b0 + b1*x1 + b2*x2 + noise.

  static const int N = 7;                                   // number of datapoints (pizza orders)
  double x0[N] = {   1,   1,   1,    1,    1,    1,    1 }; // 1st regressor: constant term 
  double x1[N] = {  24,  25,  28,   30,   40,   36,   32 }; // 2nd regressor: pizza diameters
  double x2[N] = {   0,   1,   0,    1,    2,    1,    3 }; // 3rd regressor: number of extras
  double y[N]  = { 5.5, 8.9, 7.4, 10.9, 19.4, 13.8, 14.0 }; // regressand:    pizza prices

  // create input data matrix:
  double* X[3];
  X[0] = x0;
  X[1] = x1;
  X[2] = x2;

  // Compute right hand side for linear equation system (X * X^T) * b = X * Y:
  double Y[3];

  // on wikipedia, the formula reads (X^T * X) * b = X^T * Y because there, the data arrays with 
  // the independent variables are stored in the columns of the X-matrix whereas we store them in 
  // the rows

  //RAPT::rsMatrixTools::transposedMatrixVectorMultiply(X, y, Y, N, 3);
  // something is wrong - maybe we don't need the transposed version bcs our X matrix is defined 
  // differently? - yes - we have the regressors in the rows, wikipedia has them in the columns

  //RAPT::rsMatrixTools::matrixVectorMultiply(X, y, Y, N, 3); // access violation
  RAPT::rsMatrixTools::matrixVectorMultiply(X, y, Y, 3, N);   // works but shouldn't according to doc

  // establish matrix on left-hand-side:
  double flatXX[9];     // flat storage array
  double* XX[3];        // array of pointers to rows
  XX[0] = &flatXX[0];
  XX[1] = &flatXX[3];
  XX[2] = &flatXX[6];
  RAPT::rsMatrixTools::matrixMultiplySecondTransposed(X, X, XX, 3, 7, 3); // is this correct?

  // estimate parameter vector by solving the linear system:
  double b[3];  
  RAPT::rsLinearAlgebra::rsSolveLinearSystemInPlace(XX, b, Y, 3); // get rid of rs prefix
  // verify, if the result is correct!

  // Use the model to predict the pizza prices from their diameters and extras:
  double yp[N], ea[N], er[N];
  for(int n = 0; n < N; n++) {
    yp[n] = b[0]*x0[n] + b[1]*x1[n] + b[2]*x2[n]; // predicted price
    ea[n] = y[n]  - yp[n];                        // absolute error
    er[n] = ea[n] / y[n];                         // relative error
  }

  // plot - draw true datapoints and prdicted ones:
  GNUPlotter plt;
  plt.addDataArrays(N, x1, x2, y);  // the triplets should be interpreted as 3D points
  plt.addDataArrays(N, x1, x2, yp);
  plt.setGraphStyles("points", "points");  // maybe select style - or use better default
  plt.plot3D(); 
  // that doesn't look too bad actually :-) ...but verify and clean up
  // todo: plot the plane (semitransparently?) that is defined by our model
}

// https://en.wikipedia.org/wiki/Linear_regression#Least-squares_estimation_and_related_techniques
// https://en.wikipedia.org/wiki/Weighted_least_squares

// ToDo: generalize and factor out the code for multiple regression and move it to rapt 
// (Numerics/DataFitting.h/cpp), then 
// implement polynomial regression on top of that (use 1,x,x^2,x^3,... as regressors for a given
// x-array) - this may e generalized to use a set of arbitrary functions of x - the model is a 
// weighted sum of these functions of x - for example, we could model a signal with sines/cosines
// of given frequencies - maybe create experiments for that (polynomialRegression, 
// sinusoidalRegression, exponentialRegression) - problem: the frequnecies and decays must be known
// in advance - can these be estimated too? ...for exponential decays, there's code somewhere but 
// what about sine frequencies? maybe the exponential can be made complex? -> figure out

template<class T>
RAPT::rsPolynomial<T> fitPolynomial(int numDataPoints, T* x, T* y, int degree)
{
  typedef RAPT::rsArrayTools AT;
  //typedef RAPT::rsMatrixTools MT;


  // create MxN data matrix X:
  int M = degree+1;       // # rows
  int N = numDataPoints;  // # cols
  rsMatrix<T> X(M, N);
  AT::fillWithValue(X.getRowPointer(0), N, 1.0);   // 1st row is all ones
  for(int i = 1; i < M; i++)                       // i-th row is (i-1)th row times x
    AT::multiply(X.getRowPointer(i-1), x, X.getRowPointer(i), N);

  //plotMatrixRows(X); // ok - looks good

  // compute M-vector Y

  // compute MxM matrix X * X^T that appears in the system (X * X^T) * b = X * Y:

  rsMatrix<T> XX = X * X.getTranspose();
  // optimize: the product of X * X^T can be allocated and filled without explicitly
  // creating the transposed matrix as temporary object


  // compute target vector Y = X*y

  // Compute right hand side for linear equation system (X * X^T) * b = X * Y:


  RAPT::rsPolynomial<T> p(degree);
  //T* coeffs = p.coeffs(); // private - but we are a friend - friend stuff doesn't compile
  T* coeffs = p.getCoeffPointer();  // temporary solution - getCoeffPointer should not exist!




  /*
  T** X;
  MT::allocateMatrix(X, M, N);
  AT::fillWithValue(X[0], N, 1.0);      // 1st row is all ones
  for(int i = 1; i < M; i++)
    AT::multiply(X[i-1], x, X[i], N);   // i-th row is (i-1)th row times x
    */





  //...

  //MT::deallocateMatrix(X, M, N);
  return p;
}
// todo: use the new matrix stuff - we should adapt the Gaussian elimination algorithm so it may
// work with a matrix given in flat storage format...or maybe it should get a pointer to a 
// MatrixView

void polynomialRegression()
{
  // not yet finished

  // We create data from a polynomial function with added noise and try to estimate the polynomial
  // coefficients from the data.

  // User parameters:
  double noise         = 1.0;   // amount of noise
  int    numDataPoints = 100;   // number of data points
  int    modelDegree   = 5;     // degree of our model polynomial
  double xMin          = -1.0;
  double xMax          = +2.0;


  typedef std::vector<double> Vec;
  typedef RAPT::rsPolynomial<double> Poly;
  typedef RAPT::rsArrayTools AT;

  // The polynomial to generate our data is 1 - 2x + 3x^2 - 4x^3 + 5x^4:
  Poly p(Vec({ 5,-4,3,-2,1 }));  // todo: have a constructor that takes an inititalizer list

  // Generate data:
  int N = numDataPoints;
  Vec x(N), yc(N), yn(N);  // x-values, clean and noisy y-values
  AT::fillWithRangeLinear(&x[0], N, xMin, xMax);
  RAPT::rsNoiseGenerator<double> prng;
  for(int n = 0; n < N; n++) {
    yc[n] = p(x[n]);
    yn[n] = yc[n] + noise * prng.getSample();
  }

  // estimate polynomial - make a convenience function that takes the data as input and returns
  // an rsPolynomial object:
  Poly qc = fitPolynomial(N, &x[0], &yc[0], modelDegree); // fit to the clean data



  // plot:
  rsPlotVectorsXY(x, yc, yn);
  //GNUPlotter plt;
}

// todo: make a similar function for sinusoidalRegression


/*
Idea: we model the signal x(t) by a polynomial such that:
x(t) = a0 + a1*t + a2*t^2 + a3*t^3 + ... + aN*t^N
and we have signal values x(0), x(-1), x(-2), x(-3),... available. We want to find the 
polynomial coeffs and then evaluate x(1) to predict/extrapolate the signal value at t=1 
(conveniently, at t=1, the evaluation of the polynomial boils down to a sum over the coeffs)
We get a linear system of equations:
x( 0) = a0
x(-1) = a0 -   a1 +    a2 -    a3 +     a4 -      a5 + ...
x -2) = a0 - 2*a1 +  4*a2 -  8*a3 +  16*a4 -   32*a5 + ... 
x(-3) = a0 - 3*a1 +  9*a2 - 27*a3 +  81*a4 -  243*a5 + ...
x(-4) = a0 - 4*a1 + 16*a2 - 64*a3 + 256*a4 - 1024*a5 + ...
...
hmm...the matrix entries grow exponentially - that's probably not good for numerical 
conditioning - maybe the prediction can be made without explictly solving the system by 
somehow constructing the Netwon polynomial?....figure that out

but maybe a rational function is better than a polynomial?
        a0 + a1*t + a2*t^2 + a3*t^3 + ... + aN*t^N
x(t) = --------------------------------------------
        b0 + b1*t + b2*t^2 + b3*t^3 + ... + bM*t^M
this also leads to a linear system of equations:
x( 0) *  b0                             = a0
x(-1) * (b0 -   b1 +   b2 -   b3 + ...) = a0 -   a1 +    a2 -    a3 + ...
x(-2) * (b0 - 2*b1 + 4*b2 - 8*b3 + ...) = a0 - 2*a1 +  4*a2 -  8*a3 + ...
...
the advantage of the rational function might be that it doesn't necessarily grow so much, in 
fact, if N==M, it approaches aN/bM as t -> inf and if M > N, it approaches zero
maybe try a biquadratic -> 6 coeffs a0,a1,a2,b0,b1,b2 but actually we have only 5 degrees of 
freedom due to being able to normalize the polynomials (make them monic and drag a "gain" 
factor before the fraction or something) so we would use 5 samples
*/

void polynomialPrediction()
{

}


/*
        a0 + a1*t + a2*t^2
x(t) = --------------------
        b0 + b1*t + b2*t^2

x( 0) *  b0                 = a0
x(-1) * (b0 -   b1 +    b2) = a0 -   a1 +    a2
x(-2) * (b0 - 2*b1 +  4*b2) = a0 - 2*a1 +  4*a2
x(-3) * (b0 - 3*b1 +  9*b2) = a0 - 3*a1 +  9*a2
x(-4) * (b0 - 4*b1 + 16*b2) = a0 - 4*a1 + 16*a2

let's normalize, such that b2 = 1 (make denominator monic, overall scale factor is a2), we get 5
equations for 5 unknowns and solve them with sage:

var("x0 x1 x2 x3 x4 a0 a1 a2 b0 b1 b2")
e1 = x0 *  b0              == a0
e2 = x1 * (b0 -   b1 +  1) == a0 -   a1 +    a2
e3 = x2 * (b0 - 2*b1 +  4) == a0 - 2*a1 +  4*a2
e4 = x3 * (b0 - 3*b1 +  9) == a0 - 3*a1 +  9*a2
e5 = x4 * (b0 - 4*b1 + 16) == a0 - 4*a1 + 16*a2
solve([e1,e2,e3,e4,e5],[a0,a1,a2,b0,b1])

which gives the result:

a0 = 12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)*x0
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
a1 = (x2*x3*x4 + (x1*(7*x2 - 36*x3 + 30*x4) + x2*(45*x3 - 64*x4) + 18*x3*x4)*x0 - (x2*(16*x3 - 27*x4) + 12*x3*x4)*x1)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
a2 = (x2*x3*x4 + (x1*(x2 - 6*x3 + 6*x4) + x2*(9*x3 - 16*x4) + 6*x3*x4)*x0 - (x2*(4*x3 - 9*x4) + 6*x3*x4)*x1)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
b0 = 12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4), 
b1 = (x0*(x1 - 12*x2 + 27*x3 - 16*x4) + x1*(18*x2 - 64*x3 + 45*x4) + 6*x2*(5*x3 - 6*x4) + 7*x3*x4)
     /(x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4)

*/
double biquadraticPrediction(double x0, double x1, double x2, double x3, double x4)
{
  double s, a0, a1, a2, b0, b1;
  s  = 1.0 / (x0*(x1 - 6*x2 + 9*x3 - 4*x4) + x1*(6*x2 - 16*x3 + 9*x4) + 6*x2*(x3 - x4) + x3*x4);
  a0 = s * (12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4)*x0);
  a1 = s * ((x2*x3*x4 + (x1*(7*x2 - 36*x3 + 30*x4) + x2*(45*x3 - 64*x4) + 18*x3*x4)*x0 
            - (x2*(16*x3 - 27*x4) + 12*x3*x4)*x1));
  a2 = s * ((x2*x3*x4 + (x1*(x2 - 6*x3 + 6*x4) + x2*(9*x3 - 16*x4) + 6*x3*x4)*x0 
            - (x2*(4*x3 - 9*x4) + 6*x3*x4)*x1));
  b0 = s * (12*(x1*(x2 - 4*x3 + 3*x4) + x2*(3*x3 - 4*x4) + x3*x4));
  b1 = s * ((x0*(x1 - 12*x2 + 27*x3 - 16*x4) + x1*(18*x2 - 64*x3 + 45*x4) 
            + 6*x2*(5*x3 - 6*x4) + 7*x3*x4));
  // maybe optimize ..the product x3*x4 appears often

  return (a0 + a1 + a2) / (b0 + b1 + 1.0); // evaluate rational function at x=1
}
//....but: what if the function has a pole at or near x=1? or what if s = inf ...maybe rational 
// extrapolation is not such a good idea after all...hmmm


// maybe move to RAPT into the Statistics section
double variance(double *x, int N)
{
  double mx  = RAPT::rsArrayTools::mean(x, N); // mean of x
  double sum = 0;
  for(int n = 0; n < N; n++) {
    double d = x[n] - mx;
    sum += d*d;
  }
  return sum / (N-1);
  // todo: make an optimized version that takes the mean as argument
}
double standardDeviation(double *x, int N)
{
  return sqrt(variance(x, N));
}
double covariance(double *x, double *y, int N)
{
  double mx = RAPT::rsArrayTools::mean(x, N); // mean of x
  double my = RAPT::rsArrayTools::mean(y, N); // mean of y
  double sum = 0;
  for(int n = 0; n < N; n++)
    sum += (x[n]-mx) * (y[n]-my);
  return sum / (N-1);
}
double correlation(double *x, double *y, int N)
{
  double vxy = covariance(x, y, N);
  double vx  = variance(x, N);
  double vy  = variance(y, N);
  return vxy / sqrt(vx*vy);
}
double conditionalProbability(double* a, double* b, int N)
{
  // conditional probability of "a given b"
  // a and b must be arrays of boolean values (0 or 1) but as double data type
  int na = 0, nb = 0;
  for(int n = 0; n < N; n++) {
    if(b[n] == 1) {
      nb++;
      if(a[n] == 1)
        na++;
    }
  }
  return (double) na / (double) nb;

  // can this somehow be generalized...like sum(prod(a, b)) / sum(b) ...sum over the elementwise 
  // product of the realizations of a,b divided by sum of realizations of b - should allow
  // realizations between 0..1
}
// or maybe it's called conditional relative frequency?
// https://mathbitsnotebook.com/Algebra1/StatisticsReg/ST2TwoWayTable.html

double jointProbability(double* a, double* b, int N)
{
  int sum = 0;
  for(int n = 0; n < N; n++)
    if(a[n] == 1 && b[n] == 1)
      sum += 1;
  return sum / (double) N;
  // maybe generalize by using a sum-of-products (divided by N)
}
// maybe rename to andProbability and write also an orProbability


void probabilityLogic()
{
  // consider an example: produce random number x in 0..1, 

  // let: A:  x < 0.5, B: 0.3 < x < 0.6,
  // -> P(A)=0.5, P(B)=0.3, P(A|B)=2/3, P(B|A)=2/5=0.4 -> what's the correlation C(A,B)?


  // setup:
  int N = 100000;   // number of realizations to produce
  double aL = 0.0;  // lower limit for x, such that event A is considered to have occured
  double aU = 0.5;  // upper limit ....
  double bL = 0.3;  // ...same for event B
  double bU = 0.6;


  typedef std::vector<double> Vec;
  Vec A(N), B(N);   // realizations of events A and B (true/false, represented as 0/1)
  Vec x(N);
  int n;
  RAPT::rsNoiseGenerator<double> prng;
  prng.setRange(0, 1);
  for(n = 0; n < N; n++)
  {
    double xn = prng.getSample();
    x[n] = xn;

    if(xn > aL && xn < aU)  
      A[n] = 1;
    else         
      A[n] = 0;

    if(xn > bL && xn < bU)
      B[n] = 1;
    else
      B[n] = 0;
  }

  // compute relative frequencies of events A and B (should approximate their probabilities):
  //double fA = RAPT::rsArrayTools::sum(&A[0], N) / N;
  //double fB = RAPT::rsArrayTools::sum(&B[0], N) / N;
  // are actually the mean values

  // compute sample mean values for event A and B:
  double mA = RAPT::rsArrayTools::mean(&A[0], N);
  double mB = RAPT::rsArrayTools::mean(&B[0], N);

  // compute sample variances:
  double vA = variance(&A[0], N);
  double vB = variance(&B[0], N);

  // compute sample covariance and correlation:
  double cov = covariance( &A[0], &B[0], N);
  double cor = correlation(&A[0], &B[0], N);

  // compute empirical probabilities (by relative frequencies):
  double pa  = RAPT::rsArrayTools::sum(&A[0], N) / N;        // P(A), empirical prob of event A
  double pb  = RAPT::rsArrayTools::sum(&B[0], N) / N;        // P(B), empricial prob of event B
  double cab = conditionalProbability(&A[0], &B[0], N); // P(A|B), empirical prob of A given B
  double cba = conditionalProbability(&B[0], &A[0], N); // P(B|A), empirical prob of B given A
  double jab = jointProbability(      &A[0], &B[0], N); // P(A,B), empirical prob of A and B

  // compute joint probability by formulas via conditional probability:
  //double jab1 = cab * pb;   // P(A,B) = P(A|B) * P(B)
  //double jab2 = cba * pa;   // P(A,B) = P(B|A) * P(A)
  // ok, this works - both are equal to jab

  // soo, let's say P(A) and P(B) are known and let's assume that the correlation (or maybe 
  // covariance) between a and b is also known - can we compute P(A,B)?
  // how does the correlations relate to the conditional probabilities P(A|B) or P(B|A)?

  // https://math.stackexchange.com/questions/1751950/from-correlation-coefficient-to-conditional-probability
  // says:
  // Corr(A,B) = (P(A,B) - P(A)*P(B)) / sqrt( P(A)*(1-P(A)) * P(B)*(1-P(B)) )

  // let's try it:
  double cor2 = (jab-pa*pb) / sqrt(pa*(1-pa)*pb*(1-pb));
  // yes, looks good

  // from this, we may build a general continuous "and" formula that incorporates correlation
  // ...and from that, we can also create an "or" formula

  int dummy = 0;
}

double productLog(const double z) 
{
  // Evaluates the product-log function W(z) defined as the inverse of the function f(W) = W * e^W.
  // This function is also known as Lambert-W or Omega function.

  // adapted from http://keithbriggs.info/software/LambertW.c
  // here's more: http://keithbriggs.info/software.html

  if(0.0 == z) 
    return 0.0;

  double eps = 4.0e-16;
  double em1 = 0.3678794411714423215955237701614608; // 1/Euler
  double p, e, t, w;

  if(z < -em1+1e-4) // series near -em1 in sqrt(q)
  { 
    double q = z+em1, r = sqrt(q), q2 = q*q, q3 = q2*q;
    return 
      -1.0
      +2.331643981597124203363536062168*r
      -1.812187885639363490240191647568*q
      +1.936631114492359755363277457668*r*q
      -2.353551201881614516821543561516*q2
      +3.066858901050631912893148922704*r*q2
      -4.175335600258177138854984177460*q3
      +5.858023729874774148815053846119*r*q3
      -8.401032217523977370984161688514*q3*q;  // error approx 1e-16
  }

  // initial approximation for Halley iteration:
  if(z < 1.0)     // series near 0
  { 
    p = sqrt(2.0*(2.7182818284590452353602874713526625*z+1.0));                 // euler-number
    w = -1.0+p*(1.0+p*(-0.333333333333333333333+p*0.152777777777777777777777)); // -1/3, 11/72
  } 
  else 
    w = log(z);   // asymptotic

  if(z > 3.0) 
    w -= log(w);  // useful?

  // Halley iteration:
  for(int i = 0; i < 10; i++) 
  { 
    e  = exp(w); 
    t  = w*e-z;
    p  = w+1.0;
    t /= e*p-0.5*(p+1.0)*t/p; 
    w -= t;
    if(fabs(t) < eps*(1.0+fabs(w)))  // rel-abs error
      return w;
  }

  //rsAssertFalse(); // no convergence
  return 0.0;
}
void productLogPlot()
{
  int N = 1000;        // number of values
  float xMin = -3.0;
  float xMax = +15.0;

  // evaluate function:
  vector<float> x(N), y(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
    y[n] = (float) productLog(x[n]);

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}






// try something based on recursive subdivision of an interval at a ratio r: Take the interval
// 0..1, split it into to sections of length r and 1-r, apply the same subdivision to both
// sub-intervals (starting with the larger one), maybe alternate which side get r and 1-r
// ...but what if the number of desired splits is not a power of 2?...always subdivide the
// currently largest subinterval...and the decision which side gets r and which gets 1-r is based
// on which of the neighbouring sections is larger
// example: r=0.8
// 0...1 -> 0....0.8..1 ...that's good for the supersaw because when changing the number of
// saws (intervals/splits), the old splits stay where they are and it is just a new split/saw
// introduced - that means we can just fade it in when implementing a continuous density
// parameter (we may have to do something special for the transition between 1 and 2, though - like
// moving the first saw in frequency while fading in the second - but as soon as the 2 outer freqs
// are established at their places at the interval limits, more saws are just faded in when 
// increasing density continuously)
std::vector<double> intervalSplitting(int numSegments, double midpoint)
{
  int N = numSegments;
  double r = midpoint;
  std::vector<double> a(N+1);
  a[0] = 0.0;
  a[N] = 1.0;
  int n = 2;    // number of points
  int kl = 0;   // lower index
  int ku = N;   // upper index
  while(n < N)
  {
    int km = (int)round(r*(ku-kl));
    a[km] = a[kl] + r * (a[ku] - a[kl]);

    // update kl, ku
    // ...

    n++;
  }
  return a;
}


bool rangeStartLess(const rsRange<double>& r1, const rsRange<double>& r2)
{
  return r1.getMin() < r2.getMin();
}
void splitRange(const rsRange<double>& r, rsRange<double>& rl, rsRange<double>& ru, double ratio)
{
  rl = rsRange<double>(r.getMin(),  r.getMin() + ratio*r.getSize()); // lower part of range r
  ru = rsRange<double>(rl.getMax(), r.getMax());                     // upper part of range r
}
std::vector<double> intervalSplittingProto(int numSegments, double midpoint, int splitStrategy = 1)
{
  // prototype implementation - may not be optimal in terms of efficiency, but is algorithmically
  // easier to understand

  int N = numSegments;
  double r = midpoint;
  typedef rsRange<double> Range;
  std::vector<Range> s(N);  // array/set of segments/intervals/ranges (initially empty)
  s[0] = Range(0.0, 1.0);   // seed set-of-ranges is { [0,1) }
  int n = 1;                // current number of intervals
  while(n < N) 
  {
    int k = rsArrayTools::maxIndex(&s[0], n);    // index of largest range in the current set

    // split s[k] into two ranges:
    Range rl, ru;
    if(splitStrategy == 0) splitRange(s[k], rl, ru, r);  // always use r
    else if(splitStrategy == 1)                          // alternating, ...
      if(rsIsOdd(n))  splitRange(s[k], rl, ru, r);       // ...odd n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...even n uses 1-r
    else if(splitStrategy == 2)                          // alternating, ...
      if(rsIsEven(n)) splitRange(s[k], rl, ru, r);       // ...even n uses r
      else            splitRange(s[k], rl, ru, 1-r);     // ...odd n uses 1-r
    // is it somehow possible to continuously morph between the odd and even version of this
    // strategy? ...maybe just "crossfade" the result split-point arrays? that would add another
    // potentially interesting dimension for tweaking ...maybe even vector-crossfade between
    // skew-right/skew-left/alternate-odd/alternate-even?

    // the lower part of s[k] replaces sk, the upper part gets appended to the array:
    s[k] = rl;
    s[n] = ru;
    n++;
  }
  // the cost of this algorithm is O(N^2) because in each iteration, we search for the maximum 
  // which is itself an O(N) operation - can we avoid this search by somehow keeping the array of
  // ranges sorted (by size), thereby turning it into an O(N) algorithm?

  // sort the array ranges by their start-point:
  rsHeapSort(&s[0], N, &rangeStartLess);
  // this is an O(N*log(N)) operation - so if we can turn the above into an O(N) operation, the 
  // overall comlexity of range splitting would still be O(N*log(N)) - certainly much better than
  // O(N^2) - but to achieve O(N), we would have to avoid the final sorting too - maybe by always
  // keeping a version sorted by size and another version sorted by start around?

  // turn the ranges into a vector of split-points:
  std::vector<double> a(N+1); // split-points
  for(n = 0; n < N; n++)
    a[n] = s[n].getMin();
  a[n] = 1.0;
  return a;
}

void ratioGenerator()
{
  int numRatios = 6;     // number of ratios (i.e. "density")
  int numParams = 200;   // number of sampe values for the parameter of the ratio algo
  int numPrimes = numRatios+1;

  std::vector<rsUint32> primes(numPrimes);
  rsFillPrimeTable(&primes[0], numPrimes);

  typedef rsRatioGenerator<double>::RatioKind RK;
  rsRatioGenerator<double> ratGen;
  ratGen.setPrimeTable(&primes);
  ratGen.setRatioKind(RK::linToExp);
  ratGen.setParameter1(0.5);

  std::vector<double> ratios(numRatios);
  ratGen.fillRatioTable(&ratios[0], numRatios);

  std::vector<double> p(numParams);  // rename to params
  std::vector<double> r(numRatios);  // rename to tmp
  rsArrayTools::fillWithRangeLinear(&p[0], numParams, 0.0, 1.0);
  double** y; // rename to r(atios)
  rsMatrixTools::allocateMatrix(y, numRatios, numParams);

  for(int i = 0; i < numParams; i++) {
    ratGen.setParameter1(p[i]);
    ratGen.fillRatioTable(&r[0], numRatios);
    rsArrayTools::transformRange(&r[0], &r[0], numRatios, 1., 2.);  // ratios in 1...2
    for(int j = 0; j < numRatios; j++)
      y[j][i] = r[j];
  }



  GNUPlotter plt;
  plt.addDataArrays(numParams, &p[0], numRatios, y);
  plt.setRange(-0.1, 1.1, 0.9, 2.1);
  plt.plot();


  rsMatrixTools::deallocateMatrix(y, numRatios, numParams);

  // primePowerDiff has graphs that cross each other - that means that there are values for the
  // parameter for which not all ratios are distinct - but we want them to be distict

  // make an LinToExp mapping where p=0 means linear spacing, p=1 means exponential spacing 
  // (lin-spacing of the logs) and in between, there's some intermediate setting (crossfade)

  // f = i / (numRatios-1); // fraction 
  // linVal = 1 + i / (numRatios-1);
  // expVal = exp(linVal) / exp(2);
  // mixVal = (1-p)*linVal + p*expVal;




  // try the self-similar interval splitting algorithm:
  //std::vector<double> a = intervalSplittingProto(100, 1.0/GOLDEN_RATIO, 1);
  //std::vector<double> a = intervalSplittingProto(100, 0.9, 2);
  std::vector<double> a = intervalSplittingProto(16, 0.75, 1);
  double mean = rsMean(a);;
  //rsPlotMarkers(&a[0], (int) a.size());
  // on the gui for the supersaw audio module, we should visualize the spread in a similar way
  // maybe with vertical lines instead of + markers - the visualization may or may not use a 
  // normalized interval - if it doesn't, it will also visualize the amount of spread ..maybe
  // also show the sawtooth amplitudes by the height of the lines - especially useful when we later
  // introduce amp envelopes (bell-shaped curves, but mayby also inverse bell curves emphasizing 
  // the outer freqs - how about: 0: flat, 1: bell, -1: inverted bell)

  // todo: plot the ratios as function of the parameter, let the parameter go from 0 to 1
}

void ratiosLargeLcm()
{
  // We want to produce rational numbers r between 1 and 2 that have the property that when their
  // decimal expansion (or any other expansion) is truncated, the truncated numbers have a large
  // lowest common multiple

  // We want to see all LCMs of a all pairs in a range of numbers in order to pick out those, for
  // which we have a resonably large LCM for any pair of numbers between some nMin and nMax. When 
  // nMin is fixed, we actually only need to consider numbers up to nMax = 2*nMin-1 - but we plot 
  // more...
  //  ...maybe it makes sense to choose prime numbers from the upper half
  // ..to make them irrational, square them, add 1 and take the square root

  unsigned int nMin = 100;
  unsigned int nMax = 2*nMin-1;
  unsigned int N = nMax - nMin;

  // fill LCM matrix:
  std::vector<double> axes(N);
  RAPT::rsArrayTools::fillWithRangeLinear(&axes[0], N, double(nMin), double(nMax));
  double** lcmMatrix;
  RAPT::rsMatrixTools::allocateMatrix(lcmMatrix, N, N);
  for(unsigned int i = 0; i < N; i++)
    for(unsigned int j = 0; j < N; j++)
      lcmMatrix[i][j] = (double) RAPT::rsLcm(nMin+i, nMin+j);


  GNUPlotter plt;
  plt.addDataMatrix(N, N, &axes[0], &axes[0], lcmMatrix);
  plt.setPixelSize(800, 800);
  plt.addCommand("set size square");

  // factor out into fucntion plotHeatMap
  plt.addGraph("i 0 nonuniform matrix w image notitle");
  //p.addCommand("set palette color");                  // this is used by default
  //p.addCommand("set palette color negative");         // reversed colors
  //p.addCommand("set palette gray negative");          // maximum is black
  //p.addCommand("set palette gray");                   // maximum is white
  plt.addCommand("set palette rgbformulae 30,31,32");     // colors printable as grayscale
  plt.plot();
  //plt.plot3D();

  RAPT::rsMatrixTools::deallocateMatrix(lcmMatrix, N, N);
}

void ratiosMetallic()
{
  // The metallic ratio rn for integer n is given by (n + sqrt(n^2+4))/2 - so the ratio of two
  // metallic ratios for integers m, n is: rn/rm = (n + sqrt(n^2+4)) / (m + sqrt(m^2+4))
  // ...maybe obtain the continued fraction expansion of these numerically and see, if we get small
  //  coefficients


  int dummy = 0;
}


void sinCosTable()
{
  // A test for the rsSinCosTable class.

  rsSinCosTable<float> table(8); // parameter is the table size 

  // create data:
  int N = 2000;  // number of values to plot
  float xMin = -15.0;
  //float xMin =   0.0;
  float xMax = +15.0;
  vector<float> x(N), ySin(N), yCos(N), ySinTbl(N), yCosTbl(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    ySin[n] = sin(x[n]);
    yCos[n] = cos(x[n]);
    //table.getValuesRounded(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesTruncated(x[n], &ySinTbl[n], &yCosTbl[n]);
    //table.getValuesLinear(x[n], &ySinTbl[n], &yCosTbl[n]);
    table.getValuesCubic(x[n], &ySinTbl[n], &yCosTbl[n]);
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &ySin[0], &yCos[0], &ySinTbl[0], &yCosTbl[0]);
  plt.plot();
}

void expBipolar()
{
  // find coeffs for a parametrized exponential function y = f(x) = a * e^(b*x) + c such that 
  // f(0) = y0, f(1) = y1, f'(0) = s. We have f'(x) = a*b*e^(b*x). The equation system is:
  // y0 = a + c         
  // y1 = a * e^b + c
  // s  = a * b
  // solving 1 and 3 for c and b and plugging into 2 gives:
  // d = y1 - y0 == a e^(s/a) - a
  // which can be given to wolfram alpha as:
  // Solve[d == a e^(s/a) - a, a ]

  double y0 = 0.2;
  double y1 = 8.0;
  double s  = 1.0;

  // compute coeffs:
  //double dy  = y0 - y1;
  //double tmp = (exp(s/dy)*s)/dy;
  //tmp = productLog(tmp);
  //double a = s*dy / (-dy * tmp + s); // nope - formula must be wrong 

  double d = y1 - y0;
  double w = productLog(- (exp(-s/d)*s) / d);
  double a = -d*s / (d*w+s);
  double b = s  / a;
  double c = y0 - a;

  // compute values at the endpoints for test:
  double f0  = a * exp(b * 0) + c; // f0 is wrong
  double f1  = a * exp(b * 1) + c;
  double fp0 = a * b;

  int dummy = 0;
}

void expGaussBell()
{
  // This is still wrong. The idea is to find a,b,c,d parameters for the function
  // f(x) = a * exp(b*x) + c * exp(d*x^2)
  // and impose the constraints f(0) = 1, f'(0) = s0, f(1) = 0, f'(1) = s1.
  // we have:
  // f (x) = a * exp(b*x) +   c  *  exp(d*x^2)
  // f'(x) = a*b*exp(b*x) + 2*c*d*x*exp(d*x^2)
  // so:
  // (1) f (0) = 1  = a + c
  // (2) f'(0) = s0 = a * b
  // (3) f (1) = 0  = a*exp(b) + c*exp(d)
  // (4) f'(1) = s1 = a*b*exp(b) + 2*c*d*exp(d)
  // solve(1=a+c,s_0=a*b,0=a*exp(b)+c*exp(d),s_1=a*b*exp(b)+2*c*d*exp(d);a,b,c,d) -> nope
  // we still need to solve the system - maybe an expression that contains only one unknown 
  // - say a - can be derived and we can st up a Newton iteration to the 1D root finding problem
  // we have: (1) c = (1-a), (2) b = (s0/a), 
  // (3) d = ln(-a*exp(b)/c) = ln(-a*exp(s0/a)/(1-a))
  // (4) s1 = a*b*exp(b) + 2*c*d*exp(d)
  // s1 = s0*exp(s0/a) + 2*(1-a)*d*exp(d)
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a))*exp(ln(-a*exp(s0/a)/(1-a))) - s1
  //  0 = s0*exp(s0/a) + 2*(1-a)*ln(-a*exp(s0/a)/(1-a)) * -a*exp(s0/a)/(1-a) - s1
  // ...more math to do....


  // settings:
  int N = 300;        // number of values
  float xMin = 0.0;
  float xMax = 1.0;
  double s0   = -0.1f; // slope at 0
  double s1   = -0.2f; // slope at 1

                       // compute coefficients:
  double a, b, c, d;
  a = -s0 / productLog(0.5*s0*s1);
  b = -s0 / a;
  c = 1 - a;
  d = b; // wrong!

         // evaluate function:
  vector<float> x(N), y(N);
  rsArrayTools::fillWithRangeLinear(&x[0], N, xMin, xMax);
  for(int n = 0; n < N; n++)
  {
    double t = x[n]; // temporary
    t = a*exp(b*t) + c*exp(d*t);
    y[n] = (float)t;
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &x[0], &y[0]);
  plt.plot();
}


void spreadAmounts(double amount, double spread, double* amount1, double *amount2)
{
  double absAmt = amount;
  double target = 1;
  if(amount < 0) {
    absAmt = -amount;
    target = -target;
  }
  *amount1 = RAPT::rsLinToLin(absAmt, 0.0, 1.0,  spread, target);
  *amount2 = RAPT::rsLinToLin(absAmt, 0.0, 1.0, -spread, target);
}

double spreadAmounts1(double x, double y)
{ 
  double out1, out2; spreadAmounts(x, y, &out1, &out2); 
  return out1;
}
double spreadAmounts2(double x, double y)
{
  double out1, out2; spreadAmounts(x, y, &out1, &out2);
  return out2;
}
void twoParamRemap()
{
  // ...figure out, if for every output pair (z1,z2) there is an input pair (x,y) that leads to
  // the desired output pair...maybe try to find the inverse function 
  // we have (z1,z2) = f(x,y) -> find (x,y) = g(z1,z2) where g is the inverse of z

  // maybe plot contour lines on both functions z1(x,y), z2(x,y) and see if every possible pair
  // of such contour lines has at least one intersection - done - yes - that seems to work

  int N = 21;
  GNUPlotter plt;
  plt.addDataBivariateFunction(N, -1.0, +1.0, N, -1.0, +1.0, &spreadAmounts1); 
  plt.addDataBivariateFunction(N, -1.0, +1.0, N, -1.0, +1.0, &spreadAmounts2);
  plt.addCommand("set hidden3d");
  //plt.addCommand("set pm3d");
  plt.addCommand("set palette gray");
  plt.addCommand("set lmargin 0");         // margin between plot and left border
  plt.addCommand("set tmargin 0");         // margin between plot and top border
  plt.addCommand("set contour");       // draws contours on the base plane
  //plt.addCommand("set contour surface"); // draws contours on the surface
  //plt.addCommand("set contour both"); // draws contours on the surface and base plane
  plt.addCommand("set cntrparam levels incr -1,0.2,1");

  //plt.addCommand("set palette rgbformulae 8, 9, 7");
  //plt.addCommand("set style fill solid 1.0 noborder");
  //plt.addCommand("set pm3d depthorder noborder");
  //plt.addCommand("set pm3d lighting specular 0.6");

  //plt.addCommand("set view 0, 0");
  plt.addCommand("set view 60, 45");
  plt.plot3D();



  //http://gnuplot.sourceforge.net/demo/contours.html
}

// fun stuff:

//=================================================================================================

rsGroupString add(rsGroupString A, rsGroupString B)
{
  return A + B;
  // associative, not distributive, also, it would be weird to have addition and multipication 
  // doing the same thing (although, it's not explicitly forbidden, i think)
}
rsGroupString mul1(rsGroupString A, rsGroupString B)
{
  return (-A) + (-B); // reverse and add
  // neither associative nor distributive :-(
}
rsGroupString mul2(rsGroupString A, rsGroupString B)
{
  return -((-A) + (-B));  // reverse, add, reverse result
  // maybe associative, not distributive
}
// other ideas
// -ab*cde = (a+c)(a+d)(a+e)(b+c)(b+d)(b+e) where + is modular addition of the char-numbers, so we 
//  have 6=2*3 letters and each is given by modular addition
// -compute the (modular) sum of the letters of the first input and mod-add them to each letter in
//  the 2nd input
// -mod-add element--wise...but what if inputs have different length?
// -some sort of "convolution"

rsGroupString mul3(rsGroupString A, rsGroupString B)
{
  rsGroupString C = A;
  for(int i = 0; i < B.length(); i++) {
    if(C.length() > 0 && C[C.length()-1] != B[i])
      C.append(B[i]);
    else if(C.length() > 0)     // avoid popping on empty vector
      C.removeLast();
  }
  return C;
}
// is the same as the current "+" operator of the class - but is probably better to use as mul
// maybe we should give the operations neutral names here, such as op1, op2, etc. or name them
// what they do - concatDel, modAdd


rsGroupString add2(rsGroupString A, rsGroupString B, unsigned int m) // m: modulus
{
  // do modular addition of of the characters where the shorter string is appropriately zero-padded
  // (zero standing for a "blank" character). if the result contains trailing 0s, these will be 
  // removed. the inverse element is obtained by taking charcter-wise m-x where x is the input char
  // and m is the modulus
  // could this be distributive with this concatenation and delete pairs thing that we use as the
  // other addition? if so, we could switch to using concatenation as multiplication

  //std::vector<unsigned int> s1 = A.get(), s1 = B.get();

  int LA = A.length();
  int LB = B.length();
  //int m  = 5;              // modulus

  int i;                   // loop index 
  rsGroupString C;
  if(LA >= LB) {           // A is longer than B or has same length
    C.resize(LA);
    for(i = 0; i < LB; i++)   C[i] = (A[i] + B[i]) % m;
    for(i = LB; i < LA; i++)  C[i] =  A[i];
  }
  else {                   // A is shorter than B
    C.resize(LB);
    for(i = 0; i < LA; i++)   C[i] = (A[i] + B[i]) % m;
    for(i = LA; i < LB; i++)  C[i] =  B[i];
  }
  while(C.last() == 0)  C.removeLast();  // remove trailing zeros
  return C;
}
rsGroupString add2(rsGroupString A, rsGroupString B) // to make compatible with isAsso..., etc.
{
  return add2(A, B, 7);  // uses fixed modulus
}


// tests associativity of the given operation for the given triple of arguments
bool isAssociative(rsGroupString (*op) (rsGroupString A, rsGroupString B),
  rsGroupString A, rsGroupString B, rsGroupString C)
{
  rsGroupString s1 = op(op(A, B), C); // (A+B)+C
  rsGroupString s2 = op(A, op(B, C)); // A+(B+C)
  return s1 == s2; 
}

// tests, if the given "mul" operation is distributive over the given "add" operation for the given
// triple of input arguments
bool isDistributive(
  rsGroupString (*add) (rsGroupString A, rsGroupString B),
  rsGroupString (*mul) (rsGroupString A, rsGroupString B),
  rsGroupString A, rsGroupString B, rsGroupString C)
{
  rsGroupString s1 = mul(A,add(B,C));         // A*(B+C)
  rsGroupString s2 = add(mul(A,B), mul(A,C)); // A*B + A*C
  return s1 == s2; 
}

bool testStringMultiplication(rsGroupString (*mul) (rsGroupString A, rsGroupString B))
{
  //bool r = true;
  typedef rsGroupString2 GS;

  // test associativity:
  GS abc("abc"), cde("cde"), efg("efg");
  GS A = abc, B = cde, C = efg;

  std::string t1, t2;

  // do this in a loop with various A, B, C (maybe random or by systematically checking all 
  // possible strings up to a given length) - it should be a function that automatically
  // checks a lot of strings and also takes the operations as inputs (as pointers), so we can try
  // various things


  bool asso = true;
  GS s1 = mul(mul(A, B), C); // (A*B)*C
  GS s2 = mul(A, mul(B, C)); // A*(B*C)
  t1 = s1.toString();        // for inspection/debugging
  t2 = s2.toString();
  asso &= s1 == s2; 

  // distributivity:
  bool dist = true;
  s1 = mul(A,(B+C));        // A*(B+C)
  s2 = mul(A,B) + mul(A,C); // A*B + A*C
  t1 = s1.toString();       // for inspection/debugging
  t2 = s2.toString();
  dist &= s1 == s2;
  // should also use a loop over many strings

  // maybe try all strings of length up to 4 from the alphabet a,b,c,d - each of the 3 inputs
  // A,B,C should take on all possible values - so we need a doubly nested loop


  return asso && dist;
}

bool groupString()
{
  bool r = true;  // test result - todo: turn into unit test
  typedef rsGroupString2 GS;
  GS abc("abc"), cde("cde"), efg("efg");
  GS abde = abc + cde;    r &= abde == "abde"; 
  GS edba = -abde;        r &= edba == "edba";
  GS empty = abde - abde; r &= empty == "";

  // test associativity of addition:
  GS s1 = (abc + cde) + efg;
  GS s2 = abc + (cde  + efg);
  r &= s1 == s2;  // of course, this is only an example - do more random tests

  // test neutral element of addition:
  s1 = abc + GS(""); r &= s1 == "abc";
  s1 = GS("") + abc; r &= s1 == "abc";


  // test multiplication (with the various candidate rules):
  //r &= testStringMultiplication(&add); // asso, not distri
  //r &= testStringMultiplication(&mul1); // not asso, not distri - bad!
  r &= testStringMultiplication(&mul2);


  typedef std::vector<unsigned int> Vec;

  unsigned int m = 7; // the modulus
  rsGroupString s2314 = (Vec({2, 3, 1, 4}));
  rsGroupString s546  = (Vec({5, 4, 6}));
  rsGroupString s3613 = (Vec({3, 6, 1, 3}));
  rsGroupString t1, t2, t3;     // temporaries
  t1 = add2(s2314, s546,  m);   // 0004
  t2 = add2(s546, s2314,  m);   // 0004 -> commutative
  t1 = add2(add2(s2314, s546, m), s3613, m);  // 361
  t2 = add2(s2314, add2(s546, s3613, m), m);  // 361 -> associative
  // OK - this looks good - addition is commutative and associative - at least for the tried 
  // examples - todo: try more examples - if it works out, try inverse elements and then 
  // distributivity with the concat-delete operation ...hmm - this addition actually allows
  // pairs of equal characters - which is a good thing

  bool asso = isAssociative(&add2, s2314, s3613, s546);
  bool distri = isDistributive(&add2, &mul3, s2314, s3613, s546);



  // test multiplication: define various candidate multiplication functions
  // rsGroupString mul1(rsGroupString A, rsGroupString B), mul2, etc. and make a function 
  // testMultiplication passing a function pointer to one of these candidates. Inside this 
  // function, test associativity, distributivity

  // -maybe first try to get a ring with a couple of candidate multiplication rules and then try to 
  //  find among the candiates a rule that turns the ring into a field
  // -if no suitable multiplication rule between strings can be found, mybe try a multiplication 
  //  rule between strings and single characters - some sort of "scalar-multiplication" just like 
  //  in a vector space - for example modular-addition of the scalar value to all chars in the 
  //  string ('a' is the the neutral element and the is no null element)

  return r;
}

//=================================================================================================

void primeAlternatingSums()
{
  // Take the array of prime numbers with alternating signs
  // 2 -3 5 -7 11 -13 17 -19 ...
  // and take running sums of various order ...just for fun to see what happens
  // ...do they also all have alternating signs? what happens, if we add them to one of the
  // previous arrays...just mess around a little - may interesting patterns emerge...
  // -what changes, if we give the even-indexed primes a negative sign?
  // -what happens, if we flip every k-th sign instead of every 2nd?
  // -what happens, if we filter the arrays with integer-coefficient IIR filters? this is a 
  //  generalization of (iterated) running sums (i think)
  // -what about nonlinear and/or time-variant recursive filters with coeffs depending on index?
  // -what about trying to use non/linear prediction methods on the time series?
  //  -maybe the equations (correlation coeffs, etc) should use rational numbers?
  //  -test on number sequences with known structure first (like fibonacci numbers), see, if they
  //   pick it up correctly
  // -what if we sign invert not every second number but every second pair of numbers? ..and then 
  //  combine the resulting sequences with those obtained from inverting every other number ...then
  //  generalize: invert every triple, quadruple, etc...

  int N = 200; // number of primes

  //typedef std::vector<int> IntVec;
  typedef RAPT::rsArrayTools AR;

  // create array of primes:
  //std::vector<int> primes(N);
  std::vector<int> primes(N);
  RAPT::rsFillPrimeTable(&primes[0], N);  // make convenience function getPrimes

  // create arrays of primes with alternating signs:
  int n;
  std::vector<int> ev, od; // ev: have even-numberd signs flipped, od: odd numbered signs flipped
  ev = primes; od = primes;  
  for(n = 0; n < N; n += 2) ev[n] = -ev[n];
  for(n = 1; n < N; n += 2) od[n] = -od[n];

  // 1st order running sums:
  std::vector<int> ev1(N), od1(N);
  AR::cumulativeSum(&ev[0], &ev1[0], N);
  AR::cumulativeSum(&od[0], &od1[0], N);

  // 2nd order running sums:
  std::vector<int> ev2(N), od2(N);
  AR::cumulativeSum(&ev1[0], &ev2[0], N);
  AR::cumulativeSum(&od1[0], &od2[0], N);

  // 3rd order running sums:
  std::vector<int> ev3(N), od3(N);
  AR::cumulativeSum(&ev2[0], &ev3[0], N);
  AR::cumulativeSum(&od2[0], &od3[0], N);

  // ...at this point, we really should use a loop



  // plot:
  GNUPlotter plt;
  plt.setPixelSize(1000, 500);
  //plt.addDataArrays(N, &primes[0]);
  //plt.addDataArrays(N, &ev[0]);
  plt.addDataArrays(N, &od[0]);
  //plt.addDataArrays(N, &ev1[0]);
  plt.addDataArrays(N, &od1[0]);
  //plt.addDataArrays(N, &ev2[0]);
  plt.addDataArrays(N, &od2[0]);
  //plt.addDataArrays(N, &ev3[0]);
  plt.addDataArrays(N, &od3[0]);
  plt.plot();
  // todo: maybe plot with stems or at least, get rid of the linear interpolation (draw steps)

  // the 3rd order sum is very smooth ... subtract even and odd 3rd order sums...

  // ..but maybe it's actually faster to noodle around with that stuff in python or sage
}

// naive implementation - can probably be optimized, if we know the prime factorization of n
// see: https://www.quora.com/What-is-the-fastest-way-to-find-the-divisors-of-a-number#
void rsFindNonTrivialDivisors(rsUint32 n, std::vector<rsUint32>& d)
{
  d.clear();
  if(n == 0) return; // or should we say that 0 is divisible by any number?
  for(rsUint32 i = 2; i <= rsIntSqrt(n); i++)
    if(n % i == 0) {
      d.push_back(i);
      rsUint32 j = n/i;
      if(j != i)
        d.push_back(j);
      // todo: optimize by using a divmod operation: divmod(n, j, i)
    }
  //RAPT::rsHeapSort(&d[0], (int)d.size());
  std::sort(d.begin(), d.end());
  //d.push_back(n);
}
// we do not add the trivial divisors 1 and n to the array...or should we?
// the definition of divisors would include them: https://en.wikipedia.org/wiki/Divisor
// ...unless we qualify them as nontrivial

// computes the number of divisors of a number from the exponents of its prime-factorization
// which must be passed as argument.
rsUint32 rsNumDivisors(std::vector<rsUint32>& exponents)
{
  rsUint32 nd = 1;
  for(size_t i = 0; i < exponents.size(); i++)
    nd *= exponents[i] + 1;
  return nd;
}

class rsNumberDivisibilityInfo
{
public:
  rsNumberDivisibilityInfo(rsUint32 number)
  {
    this->number = number;
    RAPT::rsPrimeFactors(number, factors, exponents);
    rsFindNonTrivialDivisors(number, divisors);
  }

  rsUint32 number;
  std::vector<rsUint32> factors;
  std::vector<rsUint32> exponents;
  std::vector<rsUint32> divisors;

protected:

};

void divisibility()
{
  rsUint32 max = 5040;
  std::vector<rsNumberDivisibilityInfo> numInfos;
  numInfos.reserve(max+1);
  for(rsUint32 i = 0; i <= max; i++)
    numInfos.push_back(rsNumberDivisibilityInfo(i));

  // todo: find highly composite and largely composite numbers...they can be useful for GUI sizes

  // plot the number of non-trivial divisors as function of the number itself:
  std::vector<int> numDivisors(numInfos.size());
  for(size_t i = 0; i < numInfos.size(); i++)
    numDivisors[i] = (int) numInfos[i].divisors.size();

  GNUPlotter plt;
  plt.addDataArrays((int)numDivisors.size(), &numDivisors[0]);
  plt.plot();
  // todo: mark highly composite numbers with a big mark and largely composite numbers with a 
  // smaller mark
}


// computes the "arithmetic derivative" of given natural number (todo: generalize to integers and
// rationals)
rsUint32 numDeriv(rsUint32 n)  // number derivative
{
  if(n == 0 || n == 1)
    return 0;

  std::vector<rsUint32> p, e; // prime factors and exponents
  RAPT::rsPrimeFactors(n, p, e);

  double s = 0; // sum of fractions - use a rational number class later
  for(int i = 0; i < p.size(); i++)
    s += (double)e[i] / (double)p[i];

  double ns = n*s;
  rsUint32 d = (rsUint32) round(ns); // the derivative

  return d;
}
void arithmeticDerivative()
{
  int N = 100;
  std::vector<rsUint32> x(N), d(N); // numbers and their derivatives
  for(int i = 0; i < N; i++) {
    x[i] = i;
    d[i] = numDeriv(i);
  }
  int dummy = 0;
}
/*
https://oeis.org/A003415/list
[0,0,1,1,4,1,5,1,12,6,7,1,16,1,9,8,32,1,21,1,24,
10,13,1,44,10,15,27,32,1,31,1,80,14,19,12,60,1,21,
16,68,1,41,1,48,39,25,1,112,14,45,20,56,1,81,16,
92,22,31,1,92,1,33,51,192,18,61,1,72,26,59,1,156,
1,39,55,80,18,71]
*/
// https://en.wikipedia.org/wiki/Arithmetic_derivative
// https://web.archive.org/web/20050426071741/http://web.mit.edu/lwest/www/intmain.pdf
// http://oeis.org/wiki/Arithmetic_derivative

