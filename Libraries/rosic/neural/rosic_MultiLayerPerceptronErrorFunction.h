#ifndef rosic_MultiLayerPerceptronErrorFunction_h
#define rosic_MultiLayerPerceptronErrorFunction_h

// rosic-indcludes:
#include "rosic_MultiLayerPerceptron.h"
#include "../numerical/rosic_FunctionObjects.h"

namespace rosic
{

  /**

  This class implements a function object that represents the error-function of a multilayer 
  perceptron. It implements algorithms for calculating the derivative of the error-function
  with respect to the synaptic weights (most notably, the celebrated backpropagation algorithm)
  and provides functions to retrieve this gradient and manipulate the weights as monolithic 
  vectors so as to facilitate the use of generic gradient based optimization algorithms, which 
  is also the reason why this class is derived from MultivariateErrorFunction.

  References: 
   -(1) Christopher M. Bishop: Neural Networks for Pattern Recognition

   \todo: define an activation function of the form: f(z) = x / (1+abs(x)) with x being defined
   as: x = k*z^3 + (1-k)*z and optimize k for each neuron. this activation function allows for 
   concave and convex regions and has therefore potentially greater representational flexibility
   ....maybe - but this will possibly make it necessary to keep track of the inputs to the 
   activation before the function is applied unless we find a formula to express the derivative
   in terms of the function value

  */

  class MultiLayerPerceptronErrorFunction : public MultivariateErrorFunction
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor.  */
    MultiLayerPerceptronErrorFunction(MultiLayerPerceptron *mlpToTrain);  

    /** Destructor. */
    ~MultiLayerPerceptronErrorFunction();

    //---------------------------------------------------------------------------------------------
    // overrides:

    /** Returns the training error. */
    virtual double getValue(Vector p);

    /** Returns the gradient of the training error with respect to the network weights. */
    virtual Vector getGradient(Vector p);

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the input data vectors and corresponding target vectors of the training data 
    set. To save memory, the arrays are not copied into local variables but referenced directly 
    - so keep them valid as long as you use this object. */
    void setTrainingData(Vector *inputs, Vector *targets, int numPatterns);

    //---------------------------------------------------------------------------------------------
    // inquiry:

 
    //---------------------------------------------------------------------------------------------
    // weight initialization:

    /** Prints all the partial derivative matrices for one pattern to the standard output - 
    mainly for debugging. */
    void printPatternGradient();

    //=============================================================================================

    // to be made protected later:

    /** Computes the gradient of the per-pattern error function with respect to the weights for one
    particular target vector by means of the backpropagation algorithm and stores the results in 
    the member dwn. The target vector has to be passed here (which has to have a dimensionality 
    equal to the number of output nodes), the corresponding input pattern and network output are 
    taken to be the ones that are still present in the members x, y inside the 'mlp' member - thus, 
    a call to mlp->forwardPropagate with the corresponding input pattern should take place 
    before calling this function. */
    void computePatternGradient(const Vector& yTarget);

    /** Naive gradient computation by perturbing each synaptic weight in turn - for check and 
    debug pruposes only. */
    void computePatternGradientByWeightPerturbation(const Vector& yTarget);

    /** Returns all relevant per-pattern partial derivatives of the weights (excluding the dummy 
    weights to the bias nodes) in a single vector suitable for convenient optimization. */
    Vector getPatternGradient();

  protected:

    /** Returns the derivative of the activation function, given the function value (which 
    represents the activation itself). */
    INLINE double activationDerivative(double activation);

    /** Computes the deltas for all neurons, given some target vector. The corresponding input 
    pattern and network output are taken to be the ones that are still present in the members 
    x, y inside the 'mlp' member - thus, a call to mlp->computeNetworkOutput with the corresponding 
    input pattern should take place before calling this function. The results will be stored in the 
    member deltas. */
    void computeDeltas(const Vector& yTarget);

    /** Computes the error function for some desired output vector - it assumes that the network 
    output is available in mlp->y, so you may want to call mlp->forwardPropagate before. */
    double getPatternError(const Vector& yTarget);

    /** Computes the error function of the network with respect to the training data and returns 
    the result. */
    double getTrainingError();

    MultiLayerPerceptron *mlp;      // pointer to the MLP to be trained

    Matrix *dwn;      // per-pattern derivatives of the error with respect to synaptic weights
    Matrix *dwN;      // avaraged derivatives of the error with respect to synaptic weights
    Vector *deltas;   // the delta values for each layer of neurons
    Vector *xTrain;   // input vectors for training
    Vector *yTrain;   // target vectors for training
    int    nTrain;    // number of input/target pairs for training

  };

  double MultiLayerPerceptronErrorFunction::activationDerivative(double z)
  {
    double gp;
    switch( mlp->activationFunctionIndex )
    {
    case MultiLayerPerceptron::LOGISTIC:        gp = z * (1.0-z);         break; // g' = g * (1-g)
    case MultiLayerPerceptron::TANH:            gp = 1.0 - z*z;           break; // g' = 1 - g^2
    case MultiLayerPerceptron::LINEAR_RATIONAL: gp = 1.0-fabs(z); gp*=gp; break; // g' = (1-|g|)^2
    }
    return gp;
  }

} // end namespace rosic

#endif // rosic_MultiLayerPerceptronErrorFunction_h
