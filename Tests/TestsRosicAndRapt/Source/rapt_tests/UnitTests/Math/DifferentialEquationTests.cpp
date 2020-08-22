
// create some arrays for the time-axis and the 3 time series created by the system:
namespace rsTestODE
{
static const int maxNumDimensions = 3;
static const int numValues         = 1000;
double x[numValues], y[maxNumDimensions][numValues], error[maxNumDimensions][numValues];
}

class rsLorentzSystem : public rsDifferentialEquationSystemDbl
{

public:

  rsLorentzSystem()
  {
    sigma = 10.0;
    rho   = 28.0;
    beta  = 8.0/3.0;

    y.setDimensionality(3);
    y.initWithZeros();
    y.v[0] = 0.5;

    x = 0.0;
  }

  virtual rsVector<double> f( const double &x, const rsVector<double> &y)
  {
    rsVector<double> v(3);

    v[0] = sigma*(y[1]-y[0]);        // dx/dt = sigma*(y-x)
    v[1] = y[0]*(rho-y[2]) - y[1];   // dy/dt = x*(rho-z) - y;
    v[2] = y[0]*y[1] - beta*y[2];    // dz/dt = x*y - beta*z

    return v;
  }

  double sigma, rho, beta;

};


class rsBesselSystem : public rsDifferentialEquationSystemDbl
{

public:

  rsBesselSystem()
  {
    n = 0.0;

    y.setDimensionality(2);
    y.initWithZeros();
    y.v[0] = 1.0;

    x = 0.0;
  }

  virtual rsVector<double> f( const double &x, const rsVector<double> &y)
  {
    rsVector<double> v(2);

    v[0] = y[1];
    if( fabs(x) < EPS )
      v[1] = -y[0];  // avoid division by zero (use limit)
    else
      v[1] = -(x*y[1] + (x*x-n*n)*y[0]) / (x*x);

    return v;
  }

  double n;

};

// function: y = e^x -> ODE: y' = y   rename to ExponentialSystem
class rsTestSystem1 : public rsDifferentialEquationSystemDbl
{

public:

  rsTestSystem1()
  {
    y.setDimensionality(1);
    y.v[0] = 1.0;
    x = 0.0;
  }

  virtual rsVector<double> f( const double &x, const rsVector<double> &y)
  {
    rsVector<double> v(1);
    v[0] = y[0];
    return v;
  }

};


// 1-dimenstional "system" for test purposes with y' = n*((y-c)/x) - (y-c)/tau. The analytic
// solution is y(x) = x^2 * e^(-x/tau) + c.  -> check this - it doesn't seem to work:

class rsTestSystem : public rsDifferentialEquationSystemDbl
{

public:

  rsTestSystem()
  {
    n   = 2.0;
    tau = 1.0;
    c   = 1.0;

    y.setDimensionality(1);
    y.v[0] = 1.0;

    x = 0.0;
  }

  void setParameters(double newN, double newTau, double newC) // remove
  {
    n   = newN;
    tau = newTau;
    c   = newC;
  }

  virtual rsVector<double> f( const double &x, const rsVector<double> &y)
  {
    rsVector<double> v(1);

    if( fabs(x) < EPS )
      v[0] = -(y[0]-c)/tau;  // avoid division by zero (use limit)
    else
      v[0] = n*((y[0]-c)/x) - (y[0]-c)/tau;

    return v;
  }

  double n, tau, c;

};















void retrieveVariables(rsDifferentialEquationSystemDbl &theSystem, int n)
{
  rsTestODE::x[n] = theSystem.getX();
  for(int d = 0; d < theSystem.getNumDimensions(); d++)
    rsTestODE::y[d][n] = theSystem.getElementOfY(d);
}

enum integrationMethods
{
  Euler = 0,
  Midpoint,
  Heun2,
  Heun3,
  RungeKutta4,
  CashKarp5,
};

void iterateState(rsDifferentialEquationSystemDbl &theSystem, int integrationMethod, double h)
{
  switch( integrationMethod )
  {
  case Euler:       theSystem.stepEuler(h);                     break;
  case Midpoint:    theSystem.stepMidpoint(h);                  break;
  case Heun2:       theSystem.stepHeun2(h);                     break;
  case RungeKutta4: theSystem.stepRungeKutta4(h);               break;
  case CashKarp5:   theSystem.stepCashKarpWithErrorEstimate(h); break;
  }


  //
  //theSystem.stepMidpoint(h);
  //theSystem.stepRungeKutta(h);
  // \todo include a switch for the method
}

void runDifferentialEquationSystem(rsDifferentialEquationSystemDbl &theSystem,
                                   int integrationMethod, double h)
{
  for(int n = 0; n < rsTestODE::numValues; n++)
  {
    retrieveVariables(theSystem, n);
    iterateState(theSystem, integrationMethod, h);
  }
}

void testLorentzSystem()
{
  // create and initialize the system:
  rsLorentzSystem ls;
  ls.setElementOfY(0, 1.0);  // initial x-coordinate
  ls.setElementOfY(1, 2.0);  // initial y-coordinate
  ls.setElementOfY(2, 3.0);  // initial z-coordinate
  ls.setX(0.0);              // intial time

    // note that the notation used here is inconsistent with the notation used inside
    // rsDifferentialEquationSystem: our state-vector is here is defined as v = (x,y,z) and the
    // independent variable here is t (time) - there, the independent variable is called x,
    // and the whole state-vector is called y

  runDifferentialEquationSystem(ls, RungeKutta4, 0.01);

  // maybe use the 1 dimensional system y' = f(x,y) = n*(y/x) - y/tau with the solution
  // y(x) = x^n * e^(-x/tau) for automatic correctness tests - compare the numerical solution
  // against the true solution
}

void testBesselSystem()
{
  rsBesselSystem bs;
  runDifferentialEquationSystem(bs, RungeKutta4, 0.005);
}

void testTestSystem1()
{
  rsTestSystem1 ts;
  double h = 0.01;

  //runDifferentialEquationSystem(ts, Euler, h);
  //runDifferentialEquationSystem(ts, Midpoint, h);
  //runDifferentialEquationSystem(ts, Heun2, h);
  //runDifferentialEquationSystem(ts, Heun3, h);
  //runDifferentialEquationSystem(ts, RungeKutta4, h);
  runDifferentialEquationSystem(ts, CashKarp5, h);



  // compare this to the analytical solution:
  for(int i = 0; i < rsTestODE::numValues; i++)
  {
    double x = i*h;
    rsTestODE::error[0][i] = exp(x) - rsTestODE::y[0][i];
  }


  //double maxError = maxAbs(error[0], numValues);
  // if numValues = 1000, the measured maximum (accumulated) errors are:
  // Euler:       1055.6595528701000
  // Midpoint:    3.6034945509309182
  // Heun2:       3.6034945509309182    // hmm - is this method equivalent to Midpoint?
  // RungeKutta4: 1.8003898730967194e-005
  // CashKarp5:   3.0340743251144886e-009

  // \todo - auto-check that error is below some threshold.

  //int dummy = 0;




  /*
  rsTestSystem ts;

  double h = 0.01;

  double n   = 2.0;
  double tau = 1.0;
  double c   = 1.0;  // initial value -> integration constant
  ts.setParameters(n, tau, c);

  runDifferentialEquationSystem(ts, h);


  // compare this to the analytical solution:
  for(int i = 0; i < numValues; i++)
  {
    double x = i*h;
    yT[0][i] = pow(x, n) * exp(-x/tau) + c;
  }

  int dummy = 0;
  */
}


bool testDifferentialEquationSystem()
{
  std::string testName = "DifferentialEquationSystem";
  bool testResult = true;

  //testLorentzSystem();
  //rsNormalize(y[0], numValues, 1.0);
  //writeToMonoWaveFile("LorentzTest.wav", y[0], numValues, 44100, 16);

  /*
  testBesselSystem();
  RSCore::normalize(y[0], numValues, 1.0);
  RSCore::writeToMonoWaveFile("D:\\TmpData\\BesselTest.wav", y[0], numValues, 44100, 16);
  */

  testTestSystem1();

  //appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}
