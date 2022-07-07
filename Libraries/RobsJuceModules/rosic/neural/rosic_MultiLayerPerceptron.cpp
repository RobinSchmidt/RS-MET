//#include "rosic_MultiLayerPerceptron.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiLayerPerceptron::MultiLayerPerceptron(int numInputs, int numOutputs, int numHiddenLayers,
                                           int *numNeuronsInHiddenLayers)
{
  activationFunctionIndex = LINEAR_RATIONAL;
  //activationFunctionIndex = LOGISTIC;
  //activationFunctionIndex = TANH;
  this->numInputs         = numInputs;
  this->numOutputs        = numOutputs;
  this->numWeightLayers   = numHiddenLayers+1;

  x.setDimensionality(numInputs);
  y.setDimensionality(numOutputs);

  // create the neuron-layers - the '+1' is for the bias nodes:
  z    = new Vector[numWeightLayers+1];
  z[0] = Vector(numInputs+1);
  for(int i=1; i<numWeightLayers; i++)
  {
    z[i] = Vector(numNeuronsInHiddenLayers[i-1]+1);
    //Vector* zDbg = &(z[i]);  // debug
  }
  z[numWeightLayers] = Vector(numOutputs+1); // output layer has a dummy-bias, too
  for(int i=0; i<=numWeightLayers; i++)
    z[i].v[0] = 1.0;  // init bias nodes

  // create the weight layers:
  w = new Matrix[numWeightLayers];
  for(int i=0; i<numWeightLayers; i++)
  {
    int numColumns = z[i].dim;     // numColumns == number of input nodes for this weight matrix
    int numRows    = z[i+1].dim;   // numRows == number of output nodes for this weight matrix
    w[i] = Matrix(numRows, numColumns);
  }
  // the first rows of the weight matrices are actually irrelevant as they go into the bias nodes
  // of the subsequent layer - we put up with this slight overhead to simplify the implementation

  computeNumberOfWeights();
  initializeWeightsRandomly(-1.0, +1.0, 1);
}

MultiLayerPerceptron::~MultiLayerPerceptron()
{
  delete[] z;
  delete[] w;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void MultiLayerPerceptron::setActivationFunction(int newActivationFunction)
{
  if( newActivationFunction >= 0 && newActivationFunction < NUM_ACTIVATION_FUNCTIONS )
    activationFunctionIndex = newActivationFunction;
}

void MultiLayerPerceptron::setWeightVector(const rosic::Vector &newWeightVector)
{
  int i = 0;
  for(int layer=0; layer < numWeightLayers; layer++)
  {
    Matrix *wm  = &(w[layer]);

    // assign the irrelevant weights to zeros:
    for(int c=0; c<wm->numColumns; c++)
      wm->m[0][c] = 0.0;

    // retrieve the relevant weights from the vector:
    for(int r=1; r<wm->numRows; r++)
    {
      for(int c=0; c<wm->numColumns; c++)
      {
        wm->m[r][c] = newWeightVector.v[i];
        i++;
      }
    }
  }
}

void MultiLayerPerceptron::saveStateToFile(char* /*fileName*/)
{
  //....
}

void MultiLayerPerceptron::loadStateFromFile(char* /*fileName*/)
{
  //....
}

//-------------------------------------------------------------------------------------------------
// inquiry:

rosic::Vector MultiLayerPerceptron::getWeightsAsVector()
{
  rosic::Vector wv(numWeights);
  vectorizeWeightMatrices(w, &wv);
  return wv;
}

//-------------------------------------------------------------------------------------------------
// computation:

void MultiLayerPerceptron::setInput(const rosic::Vector &xIn)
{
  // copy input vector into member x:
  for(int i=0; i<numInputs; i++)
  {
    x.v[i]      = xIn.v[i];
    z[0].v[i+1] = xIn.v[i];  // +1 to skip over bias node
  }
}

void MultiLayerPerceptron::forwardPropagate()
{
  int i;
  for(i=1; i<numWeightLayers; i++)
  {
    z[i] = w[i-1] * z[i-1];  // matrix times vector
    for(int j=1; j<z[i].dim; j++)
      z[i].v[j] = activationFunction( z[i].v[j] );
    z[i].v[0] = 1.0;  // restore bias node that has been corrputed
  }

  // final layer is linear - do the same procedure without the activation function:
  z[i]      = w[i-1] * z[i-1];  // matrix times vector
  z[i].v[0] = 1.0;              // restore (dummy) bias node that has been corrputed

  // copy output layer activations into y:
  for(int j=0; j<numOutputs; j++)
    y.v[j] = z[i].v[j+1];
}

rosic::Vector MultiLayerPerceptron::getOutput()
{
  return y;
}

void MultiLayerPerceptron::computeNetworkOutput(double *xIn, double *yOut)
{
  // todo: optimize this function - get rid of some overhead due to absorbing the biases into the
  // weight matrices
  setInput(rosic::Vector(numInputs, xIn));
  forwardPropagate();
  y.storeValuesInArray(yOut);
}

//-------------------------------------------------------------------------------------------------
// others:

void MultiLayerPerceptron::initializeWeightsToZeros()
{
  for(int i=0; i<numWeightLayers; i++)
    w[i].initWithZeros();
}

void MultiLayerPerceptron::initializeWeightsRandomly(double min, double max, int seed)
{
  randomUniform(min, max, seed); // init PNRG
  for(int i=0; i<numWeightLayers; i++)
  {
    w[i].randomizeElements(min, max);
    w[i] /= sqrt((double)z[i].dim);
  }
}

void MultiLayerPerceptron::printWeights()
{
  printf("%s", "Multilayer Perceptron - weight matrices:  \n");
  printf("%s", "\n");
  for(int i=0; i<numWeightLayers; i++)
  {
    w[i].print();
    printf("%s", "\n");
  }
}

void MultiLayerPerceptron::printActivations()
{
  for(int i=0; i<=numWeightLayers; i++)
    z[i].print();
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MultiLayerPerceptron::forEachWeight( double (*f) (double) )
{
  for(int i=0; i<numWeightLayers; i++)
    w[i].applyFunction(f);
}

void MultiLayerPerceptron::computeNumberOfWeights()
{
  numWeights = 0;
  for(int layer=0; layer < numWeightLayers; layer++)
  {
    Matrix *wm  = &(w[layer]);                        // weight matrix of layer 'layer'
    numWeights += (wm->numRows - 1) * wm->numColumns; // 1st row is irrelevant
  }
}

void MultiLayerPerceptron::vectorizeWeightMatrices(rosic::Matrix *weightMatrices,
  rosic::Vector *weightVector)
{
  weightVector->setDimensionality(numWeights);
  int i = 0;
  for(int layer=0; layer < numWeightLayers; layer++)
  {
    Matrix *wm  = &(weightMatrices[layer]);
    for(int r=1; r<wm->numRows; r++)
    {
      for(int c=0; c<wm->numColumns; c++)
      {
        weightVector->v[i] = wm->m[r][c];
        i++;
      }
    }
  }
}


/*

Ideas:

Fun with the "swish" function:
https://www.desmos.com/calculator/tfkbmkiq8p
https://www.desmos.com/calculator/soqbsrcmfh
the swish function itself could be interesting as an activation function


*/