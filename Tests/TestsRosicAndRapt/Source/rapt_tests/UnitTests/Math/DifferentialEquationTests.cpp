
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


// 1-dimensional "system" for test purposes with y' = n*((y-c)/x) - (y-c)/tau. The analytic
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

bool testNewOdeSolver()
{
  bool ok = true;

  int  N = 2000;                         // Number of samples to produce

  // Generate a reference signal with our old, known to work, Lorenz-system implementation:
  rosic::LorentzSystem lorentzSystem;    
  lorentzSystem.setPseudoFrequency(500); // Determines step-size (together with sample rate) 
  lorentzSystem.setState(1, 1, 1);       // Set up initial state, (0,0,0) doesn't work - is fixed point
  std::vector<double> x(N), y(N), z(N);  // Signal arrays
  for(int n = 0; n < N; n++)
  {
    lorentzSystem.getState(&x[n], &y[n], &z[n]);
    lorentzSystem.iterateState();
  }
  //rsPlotVectors(x, y, z);  // OK - looks good.

  // Retrieve system and solver parameters from reference implementation:
  double h     = lorentzSystem.getStepSize();
  double sigma = lorentzSystem.getSigma();
  double rho   = lorentzSystem.getRho();
  double beta  = lorentzSystem.getBeta();


  // Create some temporary workspace variables that are used by the low-level API of the solver:
  double p[3];             // Vector in phase space
  p[0] = p[1] = p[2] = 1;  // We use again (1,1,1) as initial state  --is this redundant?


  double wrk[6];           // Workspace for the solver (for storing computed derivatives)
  // ToDo: document the required size. For Euler-steps, it's equal to N, for midpoint steps it's
  // 2*N, etc.


  // Now let's try to reproduce it with the old general ODE solver. The old implementation required 
  // client cod to create a subsclass of rsDifferentialEquationSystem and there implement the 
  // derivative computation method in a "template method" design pattern. The class rsLorentzSystem 
  // is such a subclass that implements the Lorenz equations in the overriden method:

  rsLorentzSystem lsOld;
  std::vector<double> x2(N), y2(N), z2(N);
  rsVector<double> state(3);
  state[0] = 1;
  state[1] = 1;
  state[2] = 1;
  lsOld.setX(0.0);
  lsOld.setY(state);
  for(int n = 0; n < N; n++)
  {
    x2[n] = lsOld.getElementOfY(0);
    y2[n] = lsOld.getElementOfY(1);
    z2[n] = lsOld.getElementOfY(2);
    lsOld.stepEuler(h);
  }
  //rsPlotVectors(x2, y2, z2); 
  ok &= x2 == x;
  ok &= y2 == y;
  ok &= z2 == z;


  // Now let's try to reproduce that using the new general ODE solver. This one has a different 
  // API. Instead of subclassing a solver class, the user must pass it a std::function object for
  // the derivative computations. The low-level API takes a reference to such a funtion as argument
  // to the stepper methods (along with pointers to a state-vector and a workspace).

  // Define the std::function object that computes the derivatives (i.e. the phase-space velocity)
  // from a given position y in phase space. The phase space point where we currently are in passed
  // in y, the derivative (velocity) should be stored in dy. That's hoW the API of the solver 
  // works.
  using Func = std::function<void(const double* y, double* dy)>;
  Func f = [&](const double* y, double* dy)
  {
    dy[0] = sigma*(y[1]-y[0]);                    // dx/dt = sigma*(y-x)
    dy[1] = y[0]*(rho-y[2]) - y[1];               // dy/dt = x*(rho-z) - y
    dy[2] = y[0]*y[1] - beta*y[2];                // dz/dt = x*y - beta*z;
  };

  // Solve the Lorenz system with the genral solver. The reference implementation uses the simple
  // forward Euler method, so we use that here, too.
  using ODES2 = rsInitialValueSolver2<double>;    // Why the 2? I don't see any with a 1.
  p[0] = p[1] = p[2] = 1;                         // Reset state
  std::vector<double> x3(N), y3(N), z3(N);
  for(int n = 0; n < N; n++)
  {
    // Retrieve current state and write into output signals:
    x3[n] = p[0];
    y3[n] = p[1];
    z3[n] = p[2];

    // Iterate state in phase space:
    ODES2::stepForwardEuler(f, 3, p, h, wrk);
  }
  //rsPlotVectors(x2, y2, z2);  // OK - looks good.

  // Check that the solver produced the same result as the reference implementation:
  ok &= x3 == x;
  ok &= y3 == y;
  ok &= z3 == z;


  // OK - we have compared the forward Euler results of all 3 implementations. Now we want to check 
  // the higher order solver methods like Runge-Kutta, etc....

  // First, the midpoint method:

  // Old solver:
  state[0] = 1;
  state[1] = 1;
  state[2] = 1;
  lsOld.setX(0.0);
  lsOld.setY(state);
  for(int n = 0; n < N; n++)
  {
    x2[n] = lsOld.getElementOfY(0);
    y2[n] = lsOld.getElementOfY(1);
    z2[n] = lsOld.getElementOfY(2);
    lsOld.stepMidpoint(h);
  }
  //rsPlotVectors(x2, y2, z2); 

  // New solver:
  p[0] = p[1] = p[2] = 1;              // Reset state
  for(int n = 0; n < N; n++)
  {
    // Retrieve current state and write into output signals:
    x3[n] = p[0];
    y3[n] = p[1];
    z3[n] = p[2];

    // Iterate state in phase space:
    //ODES2::stepRungeKutta4(f, 3, p, h, wrk);
    ODES2::stepMidpoint(f, 3, p, h, wrk);
  }
  //rsPlotVectors(x3, y3, z3); 

  ok &= x3 == x2;
  ok &= y3 == y2;
  ok &= z3 == z2;



  // Now, the 4th order Runge-Kutta method:

  //h = 0.1;  // For debugging

  // Old solver:
  state[0] = 1;
  state[1] = 1;
  state[2] = 1;
  lsOld.setX(0.0);
  lsOld.setY(state);
  for(int n = 0; n < N; n++)
  {
    x2[n] = lsOld.getElementOfY(0);
    y2[n] = lsOld.getElementOfY(1);
    z2[n] = lsOld.getElementOfY(2);
    lsOld.stepRungeKutta4(h);
  }
  //rsPlotVectors(x2, y2, z2); 

  //rsPlotVectors(x2, x3);   
  // Compare RK4 to midpoint. They stay pretty close for 1200 samples and then they go off into
  // different directions.

  // New solver:

  // New solver:
  p[0] = p[1] = p[2] = 1;              // Reset state
  for(int n = 0; n < N; n++)
  {
    // Retrieve current state and write into output signals:
    x3[n] = p[0];
    y3[n] = p[1];
    z3[n] = p[2];

    // Iterate state in phase space:
    ODES2::stepRungeKutta4(f, 3, p, h, wrk);
  }
  //rsPlotVectors(x3, y3, z3); 

  rsPlotVectors(x2, x3);
  rsPlotVectors(x2-x3,y2-y3,z2-z3); 
  // Wrong!


  ok &= x3 == x2;
  ok &= y3 == y2;
  ok &= z3 == z2;




  // Compare results of both implementattions:
  //rsPlotVectors(x2, x3); // Yes - we see a match!


  //rsPlotVectors(x, x3);   // Compare result of x-coordinate Euler and midpoint method


  // The look very different! I guess that shouldn't be surprising. Unfortunately, we have no 
  // reference signal for the midpoint method solution.
  // Oh - but we could use rsLorentzSystem as reference. It uses the old solver code



  rsAssert(ok);
  return ok;

  // ToDo:
  //
  // - Implement and test a more convenient high-level API for the new solver. It should be used 
  //   like solver.setDerivativeFunction(f); solver.setMethod(R); solver.doStep(); 
  //   solver.setState(); solver.getState(); and it should manage its workspace memory internally.
  //
  // - Maybe do some tests with a simple system that has an analytic solution such that we can 
  //   compare the results of different solvers with the analytic solution. Maybe a damped 
  //   sinuosoid would be a good example system. The ODE is given by: ...
}


bool testDifferentialEquationSystem()
{
  std::string testName = "DifferentialEquationSystem";
  bool ok = true;

  ok &= testNewOdeSolver();


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
  return ok;
}
