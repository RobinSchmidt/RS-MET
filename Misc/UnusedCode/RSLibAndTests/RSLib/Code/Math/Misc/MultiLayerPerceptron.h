#ifndef RS_MULTILAYERPERCEPTRON_H
#define RS_MULTILAYERPERCEPTRON_H

namespace RSLib
{

  /**

  This class implements a multilayer perceptron to be used for nonlinear function approximation. 
  This class here only implements the neural network 'at work', that is, the forward propagation
  of activations to produce an output vector from an input vector. The network learning/training 
  should be done by means of the MultipLayerPerceptronErrorFunction class (which establishes the 
  error function itself and realizes the gradient computation via backpropagation) and the 
  GradientBasedMinimizer class (which realizes the required optimization algorithms in a generic 
  way).

  \todo move this class to somewhere else

  */

  class RSLib_API MultiLayerPerceptron  
  {

  public:

    /** The available activation functions. */
    enum activationFunctions
    {
      LOGISTIC,
      TANH,
      LINEAR_RATIONAL,
      CUBIC_RATIONAL,

      NUM_ACTIVATION_FUNCTIONS
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass the desired number input- and output neurons, the number of 
    hidden layers and the number of neurons (excluding bias neurons) for each hidden layer as array 
    'numNeuronsInHiddenLayers' which should be of length 'numHiddenLayers'. */
    MultiLayerPerceptron(int numInputs, int numOutputs, int numHiddenLayers, 
      int *numNeuronsInHiddenLayers);  

    /** Destructor. */
    ~MultiLayerPerceptron();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Selects one of the activation functions @see: activationFunctions. */
    void setActivationFunction(int newActivationFunction);

    /** Sets up all the weight-matrices from a vector representation of the relevant weights as 
    returned by getWeightsAsVector. */
    void setWeightVector(const rsVectorDbl& newWeightVector);

    /** Saves the current state of the network (number of layers, neurons per layer, weights, 
    activation-function, ...) into a file. */
    void saveStateToFile(char *fileName);

    /** Loads a state of the network (number of layers, neurons per layer, weights, 
    activation-function, ...) from a file. */
    void loadStateFromFile(char *fileName);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of weight layers in the network. */
    int getNumWeightLayers() const { return numWeightLayers; }

    /** Returns the total number of relevant network weights (excluding the dummy weights). */
    int getNumWeights() const { return numWeights; }

    /** Returns all relevant weights (excluding the dummy weights to the bias nodes) in a single 
    vector suitable for convenient optimization. */
    rsVectorDbl getWeightsAsVector();

    //---------------------------------------------------------------------------------------------
    // computation:

    /** Sets the input vector. After setting the input vector and calling forwardPropagate, you may
    retrieve the corresponding output vector via getOutput. */
    void setInput(const rsVectorDbl& xIn);

    /** Computes the ouput activations 'y' from the input activations 'x'. */
    void forwardPropagate();

    /** Returns the output vector that has most recently be computed. */
    rsVectorDbl getOutput();

    /** Computes a network output vector y from an input vector x - the lengths of the two arrays
    must match the number of inputs and outputs that have been passed to the constructor. This 
    function is mainly made for simplified use, when no learning is taking place - it will also
    avoid some overhead computations and data copying that would be necessary in learning mode.
    For learning, the three functions setInput, forwardPropagate and getOutput should be used 
    instead. */
    void computeNetworkOutput(double *xIn, double *yOut);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Initializes all synaptic weights with zeros. */
    void initializeWeightsToZeros();

    /** Initializes all synaptic weights with pseudo-random numbers between min and max. */
    void initializeWeightsRandomly(double min = -1.0, double max = +1.0, int seed = 0);

    /** Prints all the synaptic weight-matrices to the standard output - mainly for debugging. */
    void printWeights();

    /** Prints the activations of all neurons (including input and output neurons) to the standard 
    output - mainly for debugging. */
    void printActivations();

    //=============================================================================================

  protected:

    /** Applies the selected activation function to some value x and returns the result. */
    inline double activationFunction(double x) const;

    /** Iterates through all the weights and applies the passed function to them. */
    void forEachWeight( double (*f) (double) );

    /** Computes the number of relevant weights (excluding the dummy weights to the bias nodes) and
    stores it in the member 'numWeights' - this is the dimensionality of the vector representation 
    of the network weights. */
    void computeNumberOfWeights();

    /** Amalgamates all relevant weights (excluding the dummy weights to the bias nodes) into a 
    single vector suitable for convenient optimization. The weight matrices should be passed as an 
    array of matrices with a number of elements equal to the member 'numWeightLayers' - not 
    referencing our member 'w' directly will facilitate reuse of the function for vectorizing the 
    matrices of partial derivatives in the optimization process. To avoid unnecessary data moving, 
    the function is realized as void and will assign the passed vector pointer. */
    void vectorizeWeightMatrices(rsMatrixDbl *weightMatrices, rsVectorDbl *weightVector);

    int activationFunctionIndex;
    int numInputs, numOutputs, numWeightLayers, numWeights;

    rsVectorDbl *z;   // array of vectors representing the neurons (including bias nodes for next layer)
    rsMatrixDbl *w;   // array of matrices representing synaptic weights
    rsVectorDbl  x;   // vector of network inputs (excluding bias)
    rsVectorDbl  y;   // vector of network outputs (exluding bias)

    friend class MultiLayerPerceptronErrorFunction;
    friend class MultiLayerPerceptronTrainer;

  };

  inline double MultiLayerPerceptron::activationFunction(double x) const
  {
    switch( activationFunctionIndex )
    {
    case LOGISTIC:        return 1.0 / (1.0 + exp(-x));
    case TANH:            return tanh(x);
    case LINEAR_RATIONAL: return x / (1.0 + fabs(x));
    default:              return 0.0;
    }
  }

  //===============================================================================================

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

  class RSLib_API MultiLayerPerceptronErrorFunction : public MultivariateErrorFunction
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
    virtual double getValue(rsVectorDbl p);

    /** Returns the gradient of the training error with respect to the network weights. */
    virtual rsVectorDbl getGradient(rsVectorDbl p);

    //---------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the input data vectors and corresponding target vectors of the training data 
    set. To save memory, the arrays are not copied into local variables but referenced directly 
    - so keep them valid as long as you use this object. */
    void setTrainingData(rsVectorDbl *inputs, rsVectorDbl *targets, int numPatterns);

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
    void computePatternGradient(const rsVectorDbl& yTarget);

    /** Naive gradient computation by perturbing each synaptic weight in turn - for check and 
    debug pruposes only. */
    void computePatternGradientByWeightPerturbation(const rsVectorDbl& yTarget);

    /** Returns all relevant per-pattern partial derivatives of the weights (excluding the dummy 
    weights to the bias nodes) in a single vector suitable for convenient optimization. */
    rsVectorDbl getPatternGradient();

  protected:

    /** Returns the derivative of the activation function, given the function value (which 
    represents the activation itself). */
    inline double activationDerivative(double activation);

    /** Computes the deltas for all neurons, given some target vector. The corresponding input 
    pattern and network output are taken to be the ones that are still present in the members 
    x, y inside the 'mlp' member - thus, a call to mlp->computeNetworkOutput with the corresponding 
    input pattern should take place before calling this function. The results will be stored in the 
    member deltas. */
    void computeDeltas(const rsVectorDbl& yTarget);

    /** Computes the error function for some desired output vector - it assumes that the network 
    output is available in mlp->y, so you may want to call mlp->forwardPropagate before. */
    double getPatternError(const rsVectorDbl& yTarget);

    /** Computes the error function of the network with respect to the training data and returns 
    the result. */
    double getTrainingError();

    MultiLayerPerceptron *mlp;      // pointer to the MLP to be trained

    rsMatrixDbl *dwn;      // per-pattern derivatives of the error with respect to synaptic weights
    rsMatrixDbl *dwN;      // avaraged derivatives of the error with respect to synaptic weights
    rsVectorDbl *deltas;   // the delta values for each layer of neurons
    rsVectorDbl *xTrain;   // input vectors for training
    rsVectorDbl *yTrain;   // target vectors for training
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

  //===============================================================================================

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
    void setTrainingData(rsVectorDbl *inputs, rsVectorDbl *targets, int numPatterns);



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
    void computePatternGradient(const rsVectorDbl& yTarget);

    /** Naive gradient computation by perturbing each synaptic weight in turn - for check and 
    debug pruposes only. */
    //void computePatternGradientByWeightPerturbation(double *yTarget);
    void computePatternGradientByWeightPerturbation(const rsVectorDbl& yTarget);

    /** Returns all relevant per-pattern partial derivatives of the weights (excluding the dummy 
    weights to the bias nodes) in a single vector suitable for convenient optimization. */
    rsVectorDbl getPatternGradient();

  protected:

    /** Returns the derivative of the activation function, given the function value (which 
    represents the activation itself). */
    inline double activationDerivative(double activation);

    /** Computes the deltas for all neurons, given some target vector. The corresponding input 
    pattern and network output are taken to be the ones that are still present in the members 
    x, y inside the 'mlp' member - thus, a call to mlp->computeNetworkOutput with the corresponding 
    input pattern should take place before calling this function. The results will preliminarily be
    stored in the member db, as this has the right dimensionality and is not needed otherwise 
    during this step. */
    //void computeDeltas(double *yTarget);
    void computeDeltas(const rsVectorDbl& yTarget);

    /** Computes the error function for some desired output vector - it assumes that the network 
    output is available in mlp->y, so you may want to call mlp->forwardPropagate before. */
    //double computePatternError(double *yTarget);
    double computePatternError(const rsVectorDbl& yTarget);

    // computeError

    //ErrorFunction errorFunction;
    MultiLayerPerceptron *mlp;      // pointer to the MLP to be trained

    rsMatrixDbl *dwn;      // per-pattern derivatives of the error with respect to synaptic weights
    rsMatrixDbl *dwN;      // summed derivatives of the error with respect to synaptic weights
    rsVectorDbl *deltas;   // the delta values for each layer of neurons

    rsVectorDbl *xTrain;   // input vectors for training
    rsVectorDbl *yTrain;   // target vectors for training
    int    nTrain;    // number of input/target pairs for training

  };

  double MultiLayerPerceptronTrainer::activationDerivative(double z)
  {
    double gp;
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

}

#endif
