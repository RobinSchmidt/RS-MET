//#include "rosic_GradientBasedMinimizer.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------    
// construction/destruction:

GradientBasedMinimizer::GradientBasedMinimizer()
{
  functionToMinimize   = NULL;
  algorithm            = SCALED_CONJUGATE_GRADIENT;
  printInfo            = false;   
  convergenceThreshold = 0.00001;  
  momentum             = 0.9;
  maxNumSteps          = 10000;
  step                 = 0;
  converged            = false;
  e                    = 0.0;
  stepsize             = 0.1;
  betaFormula          = POLAK_RIBIERE;
}

GradientBasedMinimizer::~GradientBasedMinimizer()
{

}

//-------------------------------------------------------------------------------------------------    
// optimization:

rosic::Vector GradientBasedMinimizer::minimizeFunction(
  MultivariateErrorFunction *functionToMinimize, rosic::Vector initialGuess)
{
  if( functionToMinimize == NULL )
    return initialGuess;

  this->functionToMinimize = functionToMinimize;
  p = initialGuess;

  switch( algorithm )
  {
  case GRADIENT_DESCENT:               minimizeViaGradientDescent();               break;
  case BOLD_DRIVER_WITH_MOMENTUM:      minimizeViaBoldDriverWithMomentum();        break;
  case CONJUGATE_GRADIENT:             minimizeViaConjugateGradient();             break;
  case SCALED_CONJUGATE_GRADIENT:      minimizeViaScaledConjugateGradient();       break;
  }

  return p;
}

void GradientBasedMinimizer::minimizeViaGradientDescent()
{
  stepsize     = 0.1;                                  
  g            = functionToMinimize->getGradient(p);   // gradient
  converged    = false;
  step         = 0;
  printStartInfo();
  while( converged == false && step <= maxNumSteps )
  {
    printProgressInfo();

    // do the update step:
    p -= stepsize * g;

    // compute new gradient:
    g = functionToMinimize->getGradient(p);

    // done - next step:
    step++;
    if( g.getEuclideanNorm()/g.dim < convergenceThreshold )
      converged = true;
  }
  printEndInfo();
}

void GradientBasedMinimizer::minimizeViaBoldDriverWithMomentum()
{
  stepsize     = 0.2;                                  // initial stepsize
  double rho   = 1.01;                                 // growth factor for the stepsize
  double sigma = 0.95;                                 // decay factor for the stepsize
  double mu    = momentum;                             // momentum constant
  g            = functionToMinimize->getGradient(p);   // gradient
  Vector d     = -g;                                   // current direction
  d.initWithZeros();
  e            = functionToMinimize->getValue(p);
  double eOld  = e;
  converged    = false;
  step         = 0;
  Vector pTmp;
  double eTmp;
  printStartInfo();
  while( converged == false && step <= maxNumSteps )
  {
    printProgressInfo();

    // calculate new direction:
    d  = (1-mu)*(-g) + mu*d;    
    //d  = -g/g.getEuclideanNorm() + mu*d;

    // do the update step, but remember the old parameter vector in order to restore it in case of
    // error increase:
    pTmp  = p;
    p    += stepsize * d;

    // check, whether the error has actually decreased due to the update step - if yes, increase 
    // the learning rate for the next step, if not, undo the update step, decrease learning rate 
    // and reset the direction to the  negative gradient:
    eTmp = e;
    e    = functionToMinimize->getValue(p);  
    if( e < eOld )
      stepsize *= rho;                             // increase learning rate
    else
    {                      
      p         = pTmp;                            // undo step
      e         = eTmp;
      stepsize *= sigma;                           // decrease stepsize
      d         = -g;                              // reset direction to the negative gradient
    }
    stepsize  = clip(stepsize, 0.000001, 100.0);
    eOld = e;

    // evaluate gradient for next iteration:
    g  = functionToMinimize->getGradient(p);

    // done - next step:
    step++;
    if( g.getEuclideanNorm()/g.dim < convergenceThreshold )
      converged = true;
  }
  printEndInfo();
}

void GradientBasedMinimizer::minimizeViaConjugateGradient()
{
  Vector gOld;              // gradient from previous iteration
  Vector gTmp;              // for temporary storage of a gradient vector
  Vector d;                 // current direction
  Vector dtH;               // approximation of d^T * H (with Hessian matrix H)
  double eps = 0.0001;      // epsilon for approximating d^T * H
  double eps2;              // epsilon scaled by norm of direction vector
  double beta;              // weight for the old direction
  g = functionToMinimize->getGradient(p);
  d = -g;
  converged    = false;
  step         = 0;
  printStartInfo();
  while( converged == false && step <= maxNumSteps )
  {
    // approximate d^T * H via a finite difference:
    eps2 = eps / d.getEuclideanNorm();
    gTmp = functionToMinimize->getGradient(p + eps2*d);
    dtH  = (gTmp-g) / eps2;

    // compute optimal stepsize (this formula assumes the Hessian to be positive definite - for 
    // more robust optimization, a line search would have to be used instead):
    stepsize = - (d*g) / (dtH*d);

    printProgressInfo();

    // do the update step:
    p += stepsize * d;

    // compute new gradient, while remembering the old one for computation of beta:
    gOld = g;
    g    = functionToMinimize->getGradient(p);

    // compute the weight for the old direction:
    if( betaFormula == FLETCHER_REEVES )
      beta = (g*g)        / (gOld*gOld);   
    else if( betaFormula == HESTENES_STIEFEL )
      beta = (g*(g-gOld)) / (d*(g-gOld));   
    else
      beta = g*(g-gOld)   / (gOld*gOld);        // Polak/Ribiere 
    beta = rmax(beta, 0.0);

    // compute new search direction:
    d = -g + beta*d;

    // done - next step:
    step++;
    if( g.getEuclideanNorm()/g.dim < convergenceThreshold )
      converged = true;
  }
  printEndInfo();
}

void GradientBasedMinimizer::minimizeViaScaledConjugateGradient()
{
  Vector gOld;              // gradient from previous iteration
  Vector gTmp;              // for temporary storage of a gradient vector
  Vector d;                 // current direction
  Vector dtH;               // approximation of d^T * H (with Hessian matrix H)
  double eNew;              // error at the new p where we (possibly) go to
  double eps = 0.0001;      // epsilon for approximating d^T * H
  double eps2;              // epsilon scaled by norm of direction vector
  double beta;              // weight for the old direction
  double lambda = 0.1;      // scale factor for the unit matrix
  double delta;             // denominator in equation for alpha - should be > 0
  double Delta;             // comparison parameter between predicted and actual error decrease
  double norm = 0;          // Euclidean norm of current direction vector
  bool   success   = true;  // flag to indicate a successful step - if false in some iteration, we 
                            // re-use the gradient and error value from the previous iteration
  g = functionToMinimize->getGradient(p);
  d = -g;
  converged    = false;
  step         = 0;
  printStartInfo();
  while( converged == false && step <= maxNumSteps )
  {
    // approximate d^T * H via a finite difference:
    if( success == true )
    {
      norm = d.getEuclideanNorm();
      eps2 = eps / norm;
      gTmp = functionToMinimize->getGradient(p + eps2*d); 
      dtH  = (gTmp-g) / eps2;
    }
    else
    {
      // the old values from the previous iteration are still valid
    }

    // compute the denominator 'delta' for the formula for alpha:
    delta = dtH*d + lambda * norm*norm;

    // check if this denominator is positive (corresponding to a positive definite modified 
    // Hessian), if not, increase lambda and re-evaluate delta:
    if( delta < 0.0 )
    {
      lambda = 2.0 * ( lambda - delta / (norm*norm) );
      delta  = dtH*d + lambda * norm*norm;
    }

    // compute the stepsize:
    stepsize = - (d*g) / delta;

    printProgressInfo();

    // compute the comparison parameter which measures the quality of the quadratic approximation:
    if( success == true )
      e = functionToMinimize->getValue(p);
    else
    {
      // the old value from the previous iteration is still valid
    }     
    eNew  = functionToMinimize->getValue(p + stepsize*d);
    Delta = -2.0*(e-eNew) / (stepsize*(d*g));

    // according to the value of this comparison parameter, different actions take place...
    if( Delta > 0.0 )
    {
      // Delta > 0 indicates a decrease of the error function for the proposed step, so we 
      // actually do the step:
      p += stepsize*d;

      // Delta also mesures the quality of the quadratic approximation - if it is large 
      // (Delta > 0.75), then the quadratic approximation is good and we may decrease lambda, if,
      // on the other hand, Delta is small (Delta < 0.25), the quadractic approximation is bad and
      // we should increase lambda:
      if( Delta > 0.75 )
        lambda *= 0.5;
      else if( Delta < 0.25 )
        lambda *= 4.0;
      lambda  = clip(lambda, DBL_MIN, DBL_MAX);

      success = true;
    }
    else
    {
      // Delta < 0 indicates an increase of the error function for the proposed step (or no change 
      // at all for Delta == 0 (?)), in this case, we increase lambda in the same way as above, but 
      // we don't actually perform the step - instead we jump out of the current iteration, such 
      // that the next iteration will attempt the step with the new value value of lambda:
      lambda *= 4.0;
      lambda  = clip(lambda, DBL_MIN, DBL_MAX);
      success = false;
      step++;
      if( printInfo == true )
        printf("%s %.4f %s",  "increase of error predicted, increasing lambda to: ", lambda, "\n");
      continue;
    }

    // O.K. - we have done a successful update step, so we now choose a new direction according to 
    // the rules of the standard conjugate gradient algorithm (see comments there for details):
    gOld = g;
    g    = functionToMinimize->getGradient(p);
    if( betaFormula == FLETCHER_REEVES )
      beta = (g*g)        / (gOld*gOld);   
    else if( betaFormula == HESTENES_STIEFEL )
      beta = (g*(g-gOld)) / (d*(g-gOld));   
    else
      beta = g*(g-gOld)   / (gOld*gOld);        // Polak/Ribiere 
    beta = rmax(beta, 0.0);
    d    = -g + beta*d;

    // done - next step:
    step++;
    if( g.getEuclideanNorm()/g.dim < convergenceThreshold )
      converged = true;
  }
  printEndInfo();
}

//-------------------------------------------------------------------------------------------------
// printing:

void GradientBasedMinimizer::printStartInfo()
{
  if( printInfo == true )
  {
    switch( algorithm )
    {
    case GRADIENT_DESCENT:
      printf("%s", "starting gradient descent algorithm with momentum with initial parameters: \n");
      break;
    case BOLD_DRIVER_WITH_MOMENTUM:
      printf("%s", "starting bold driver algorithm with momentum with initial parameters: \n");
      break;
    case CONJUGATE_GRADIENT:
      printf("%s", "starting conjugate gradient algorithm with initial parameters: \n");
      break;
    case SCALED_CONJUGATE_GRADIENT:
      printf("%s", "starting scaled conjugate gradient algorithm with initial parameters: \n");
      break;
    }
    p.print();
  }
}

void GradientBasedMinimizer::printProgressInfo()
{
  if( printInfo == true )
  {
    e = functionToMinimize->getValue(p);
    printf("%s %d %s",  "step: ", step, " ");

    //if( algorithm == BOLD_DRIVER_WITH_MOMENTUM )
      
    printf("%s %4f %s", "stepsize: ",   stepsize, " ");

    e = functionToMinimize->getValue(p);
    printf("%s %4f %s", "gradient norm: ", g.getEuclideanNorm(), " ");
    printf("%s %4f %s", "error: ",         e,                    " \n");
  }
}

void GradientBasedMinimizer::printEndInfo()
{
  if( printInfo == true )
  {
    //char *algoString;
    std::string algoString; // todo: don't mix cout/printf
    switch( algorithm )
    {
    case GRADIENT_DESCENT:           algoString = "gradient descent";            break;
    case BOLD_DRIVER_WITH_MOMENTUM:  algoString = "bold driver with momentum";   break;
    case CONJUGATE_GRADIENT:         algoString = "conjugate gradient";          break;
    case SCALED_CONJUGATE_GRADIENT:  algoString = "scaled conjugate gradient";   break;
    }
    if( converged == true )
    {
      e = functionToMinimize->getValue(p);
      printf("%s %d %s",  "step: ", step, "   ");  
      printf("%s %4f %s", "error: ",   e, " \n");
      cout << algoString;
      printf("%s %d %s", " algorithm converged at step: ", step, "\n");
      printf("%s", "optimized parameters: \n");
      p.print();
    }
    else
    {
      cout << algoString;
      printf("%s %d %s", " algorithm did not converge and aborted at step: ", step, "\n");
    }
  }
}