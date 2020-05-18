#ifndef RAPT_GRADIENTBASEDMINIMIZER_H
#define RAPT_GRADIENTBASEDMINIMIZER_H

/** This class implements gradient based minimization algorithms. It uses a function-object of 
class MultivariateErrorFunction to represent the function to be minimized. If you want to maximize 
some function instead, you may simply define the error function as the negative of the function to
be maximized. It implements a couple of different gradient based optimization methods, ranging
from the crude gradient descent with fixed stepsize to the sophisticated scaled conjugate
gradient algorithm (the latter of which is often a good choice which is why this is the default
algorithm).

References:
-(1) Christopher M. Bishop: Neural Networks for Pattern Recognition
-(2) Jonathan Richard Shewchuck: An Introduction to the Conjugate Gradient Method Without the
 Agonizing Pain

\todo: manage maximization internally by having a 'sign' member - if +1, we minimize, if -1, we
maximize ....or probably better vice versa -> rename class to GradientBasedOptimizer
maybe: incorporate entirely different optimization strategies (possibly not gradient based)
like simulated annealing, differential evolution, etc. -> rename class to Optimizer

\todo: design change instead of switching between algorithms inside this class, make subclasses
for the different algorithms, baseclass: Optimizer, subclasses OptimizerGradientDescent -
use template method pattern, let the user also take single optimization steps himself instead
of just optimizing in one go

\in many of the algorithms, we use ad-hoc values for some values like an (initial) step-size,
stopping criterion and so on - make these values available as user parameters and/or use some
more scientifically derived default values

\todo: maybe use std:vector<T> for a parameter vector and std::function for the error function,
if that makes sense ..do these changes after unit tests are in place...make it possible to
re-interpret a std::vector as a matrix and as a set of matrices (have a class rsMatrixSet, use
it to represent neural network synaptic layers). */

template<class T>
class GradientBasedMinimizer
{

public:

  enum algorithms
  {
    GRADIENT_DESCENT,
    BOLD_DRIVER_WITH_MOMENTUM,
    CONJUGATE_GRADIENT,
    SCALED_CONJUGATE_GRADIENT
  };

  enum betaFormulas
  {
    POLAK_RIBIERE,
    FLETCHER_REEVES,
    HESTENES_STIEFEL,
  };

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  GradientBasedMinimizer();

  /** Destructor. */
  ~GradientBasedMinimizer();

  //---------------------------------------------------------------------------------------------
  // \name Setup
  // todo: remove the virtual declarations

  /**  If you want to see the progress of the algorithm at the standard output, set this to
  true. */
  virtual void setPrintInfo(bool shouldPrint) { printInfo = shouldPrint; }

  /** Sets the threshold for the norm of the gradient (divided by its dimensionality) - if it
  falls below this value in the optimization algorithm, the algorithm will be considered to have
  converged and return. */
  virtual void setConvergenceThreshold(T newThreshold)
  {
    convergenceThreshold = newThreshold;
  }

  /** Sets the maximum number of steps that the algorithm will take in order not run indefinitely
  in cases where it doesn't converge. */
  virtual void setMaxNumSteps(int newMaximum) { maxNumSteps = newMaximum; }

  /** Selects one of the optimization algorithms. @see: algorithms */
  virtual void setAlgorithm(int newAlgorithm) { algorithm = newAlgorithm; }

  /** Selects the formula by which the 'beta' constant (the weight for the previous direction)
  in the (scaled) conjugate gradient algorithm is calculated. There are 3 different forms for
  this constant which are all equivalent for a quadratic function but will give rise to different
  convergence behavior when conjugate gradient is applied to non-quadratic functions. These forms
  are known as Fletcher/Reeves, Hestenes/Stiefel and Polak/Ribiere forms respectively.
  Polak/Ribiere is the default choice as the literature ascribes the best convergence properties
  to this form. You may, however, experiment with the other forms - perhaps, the optimal form
  depends on the problem at hand. @see: betaFormulas */
  virtual void setBetaFormula(int newFormula) { betaFormula = newFormula; }

  /** Sets the momentum constant for algorithms that use momentum. The value should be between
  0.0 (inclusive) and 1.0 (exclusive). */
  virtual void setMomentum(T newMomentum)
  {
    rsAssert(newMomentum >= 0.0 && newMomentum < 1.0, "Momentum out of range");
    momentum = rsClip(newMomentum, 0.0, 1.0);
  }

  void setStepSize(T newStepSize)
  {
    stepsize = newStepSize;
  }

  //---------------------------------------------------------------------------------------------
  // inquiry - not really useful - get rid of them:

  /** Returns the threshold for convergence. @see: setConvergenceThreshold */
  virtual T getConvergenceThreshold() const { return convergenceThreshold; }

  /** Sets the maximum number of steps. @see: setMaxNumSteps */
  virtual int getMaxNumSteps() const { return maxNumSteps; }

  /** Returns the index of the selected algorithm. @see: setAlgorithm */
  virtual int getAlgorithm() const { return algorithm; }

  /** Returns the index of the formula for 'beta'. @see: setBetaFormula */
  virtual int getBetaFormula() const { return betaFormula; }

  /** Returns the momentum constant. @see: setMomentum */
  virtual T getMomentum() const { return momentum; }

  //---------------------------------------------------------------------------------------------
  // optimization:

  /** Minimizes the function 'functionToMinimize' starting at the 'initialGuess' and returns the
  result. */
  virtual rsVectorDbl minimizeFunction(MultivariateErrorFunction<T> *functionToMinimize,
    rsVectorDbl initialGuess);

  //=============================================================================================

protected:

  /** Minimizes the function via a simple gradient descent procedure with a fixed stepsize. It is
  not really recommended to use this procedure for production code as it has rather slow
  convergence properities and may sometimes even diverge. It is mainly there for experimentation
  and benchmarking other algorithms against it during development. The algorithm requires
  1 gradient evaluation per step. */
  virtual void minimizeViaGradientDescent();

  /** Minimizes the function via gradient descent with varying stepsize and momentum. The
  momentum constant can be set via setMomentum and the adaption of the stepsize is done with two
  fixed factors 'rho' (for error decrease) and 'sigma' (for error increase). The algorithm
  requires 1 gradient evaluation and one function evaluation per step.*/
  virtual void minimizeViaBoldDriverWithMomentum();

  /** Minimizes the function via the conjugate gradient algorithm. The algorithm is implemented
  in its pure form, that is, with a computation of the optimal stepsize from the (approximate)
  Hessian matrix (as opposed to invoking a line search). This implies that it will converge only
  for problems where the Hessian is positive definite, so it's actually not recommended to use
  this algorithm for production code but rather for experimentation. The algorithm requires
  2 gradient evaluations per step. */
  virtual void minimizeViaConjugateGradient();

  /** Minimizes the function via the scaled conjugate gradient algorithm - this is often the best
  choice, so it's the default algorithm. The algorithm requires 2 gradient evaluations and 2
  function evaluations per step. */
  virtual void minimizeViaScaledConjugateGradient();

  // some internal function to print out some information during development:
  virtual void printStartInfo();
  virtual void printProgressInfo();
  virtual void printEndInfo();

  // pointer to the objective function to be minimized:
  MultivariateErrorFunction<T> *functionToMinimize;

  int algorithm;               // choice of the optimization algorithm
  int maxNumSteps;             // maximum number of steps
  int step;                    // current step
  int betaFormula;             // formula for beta in the conjugate gradient algoritm
  bool converged;              // convergence flag
  bool printInfo;              // flag to indicate to print information to the standard output
  T e;                         // error function value
  T stepsize;                  // stepsize 
  T momentum;                  // momentum constant for algorithms that use it
  T convergenceThreshold;      // threshold (for the normalized gradient norm) for considering 
                               // the minimization as converged
  rsVectorDbl p;               // parameter vector
  rsVectorDbl g;               // gradient vector

};


#endif 
