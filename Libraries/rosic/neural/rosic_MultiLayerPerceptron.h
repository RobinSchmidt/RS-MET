#ifndef rosic_MultiLayerPerceptron_h
#define rosic_MultiLayerPerceptron_h

// rosic-indcludes:
#include "../math/rosic_MatrixVectorFunctions.h"

namespace rosic
{

  /**

  This class implements a multilayer perceptron to be used for nonlinear function approximation. 
  This class here only implements the neural network 'at work', that is, the forward propagation
  of activations to produce an output vector from an input vector. The network learning/training 
  should be done by means of the MultipLayerPerceptronErrorFunction class (which establishes the 
  error function itself and realizes the gradient computation via backpropagation) and the 
  GradientBasedMinimizer class (which realizes the required optimization algorithms in a generic 
  way).

  */

  class MultiLayerPerceptron  
  {

    friend class MultiLayerPerceptronErrorFunction;

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
    void setWeightVector(const Vector& newWeightVector);

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
    Vector getWeightsAsVector();

    //---------------------------------------------------------------------------------------------
    // computation:

    /** Sets the input vector. After setting the input vector and calling forwardPropagate, you may
    retrieve the corresponding output vector via getOutput. */
    void setInput(const Vector& xIn);

    /** Computes the ouput activations 'y' from the input activations 'x'. */
    void forwardPropagate();

    /** Returns the output vector that has most recently be computed. */
    Vector getOutput();

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
    INLINE double activationFunction(double x) const;

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
    void vectorizeWeightMatrices(Matrix *weightMatrices, Vector *weightVector);

    int activationFunctionIndex;
    int numInputs, numOutputs, numWeightLayers, numWeights;

    Vector *z;   // array of vectors representing the neurons (including bias nodes for next layer)
    Matrix *w;   // array of matrices representing synaptic weights
    Vector  x;   // vector of network inputs (excluding bias)
    Vector  y;   // vector of network outputs (exluding bias)

  };

  INLINE double MultiLayerPerceptron::activationFunction(double x) const
  {
    switch( activationFunctionIndex )
    {
    case LOGISTIC:        return 1.0 / (1.0 + exp(-x));
    case TANH:            return tanh(x);
    case LINEAR_RATIONAL: return x / (1.0 + fabs(x));
    default:              return 0.0;
    }
  }

} // end namespace rosic

#endif // rosic_MultiLayerPerceptron_h
