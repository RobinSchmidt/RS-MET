#ifndef rosic_MultiLayerPerceptronTrainer_h
#define rosic_MultiLayerPerceptronTrainer_h

//// rosic-indcludes:
//#include "rosic_MultiLayerPerceptron.h"

namespace rosic
{

  /**

  This class implements the training algorithms for the MultiLayerPerceptron class. It maintains a
  pointer to the MLP that is to be trained which must be passed to the constructor.

  References: 
   -(1) Christopher M. Bishop: Neural Networks for Pattern Recognition

  */

  class MultiLayerPerceptronTrainer  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor.  */
    MultiLayerPerceptronTrainer(MultiLayerPerceptron *mlpToTrain);  

    /** Destructor. */
    ~MultiLayerPerceptronTrainer();

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

    /** Initializes all synaptic weights with zeros. */
    void initializeWeightsToZeros();

    /** Initializes all synaptic weights with pseudo-random numbers between min and max. */
    void initializeWeightsRandomly(double min = -1.0, double max = +1.0, int seed = 0);

    /** Prints all the partial derivative matrices for one pattern to the standard output - 
    mainly for debugging. */
    void printPatternGradient();

    //=============================================================================================

    /** Computes the gradient of the per-pattern error function with respect to the weights for one
    particular target vector by means of the backpropagation algorithm and stores the results in 
    the member dwn. The target vector has to be passed here (which has to have a dimensionality 
    equal to the number of output nodes), the corresponding input pattern and network output are 
    taken to be the ones that are still present in the members x, y inside the 'mlp' member - thus, 
    a call to mlp->computeNetworkOutput with the corresponding input pattern should take place 
    before calling this function. */
    //void computePatternGradient(double *yTarget);
    void computePatternGradient(const Vector& yTarget);

    /** Naive gradient computation by perturbing each synaptic weight in turn - for check and 
    debug pruposes only. */
    //void computePatternGradientByWeightPerturbation(double *yTarget);
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
    input pattern should take place before calling this function. The results will preliminarily be
    stored in the member db, as this has the right dimensionality and is not needed otherwise 
    during this step. */
    //void computeDeltas(double *yTarget);
    void computeDeltas(const Vector& yTarget);

    /** Computes the error function for some desired output vector - it assumes that the network 
    output is available in mlp->y, so you may want to call mlp->forwardPropagate before. */
    //double computePatternError(double *yTarget);
    double computePatternError(const Vector& yTarget);

    // computeError

    //ErrorFunction errorFunction;
    MultiLayerPerceptron *mlp;      // pointer to the MLP to be trained

    Matrix *dwn;      // per-pattern derivatives of the error with respect to synaptic weights
    Matrix *dwN;      // summed derivatives of the error with respect to synaptic weights
    Vector *deltas;   // the delta values for each layer of neurons

    Vector *xTrain;   // input vectors for training
    Vector *yTrain;   // target vectors for training
    int    nTrain;    // number of input/target pairs for training

  };

  double MultiLayerPerceptronTrainer::activationDerivative(double z)
  {
    double gp = 0;
    switch( mlp->activationFunctionIndex )
    {
    case MultiLayerPerceptron::LOGISTIC:
      {
        gp  = z * (1.0-z);     // g' = g * (1-g)
      }
      break;
    case MultiLayerPerceptron::TANH:
      {
        gp  = 1.0 - z*z;       // g' = 1 - g^2
      }
      break;
    case MultiLayerPerceptron::LINEAR_RATIONAL:
      {
        gp  = 1.0 - fabs(z);   // g' = (1 - |g|)^2
        gp *= gp;
      }
      break;
    }
    return gp;
  }

} // end namespace rosic

#endif // rosic_MultiLayerPerceptronTrainer_h
