using namespace RSLib;

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
  z    = new rsVectorDbl[numWeightLayers+1];
  z[0] = rsVectorDbl(numInputs+1);
  for(int i=1; i<numWeightLayers; i++)
  {
    z[i] = rsVectorDbl(numNeuronsInHiddenLayers[i-1]+1);
    rsVectorDbl* zDbg = &(z[i]);
  }
  z[numWeightLayers] = rsVectorDbl(numOutputs+1); // output layer has a dummy-bias, too
  for(int i=0; i<=numWeightLayers; i++)
    z[i].v[0] = 1.0;  // init bias nodes

  // create the weight layers:
  w = new rsMatrixDbl[numWeightLayers];
  for(int i=0; i<numWeightLayers; i++)
  {
    int numColumns = z[i].dim;        // number of input nodes for this weight matrix
    int numRows    = z[i+1].dim; // number of output nodes for this weight matrix
    w[i] = rsMatrixDbl(numRows, numColumns);
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

void MultiLayerPerceptron::setWeightVector(const rsVectorDbl &newWeightVector)
{
  int i = 0;
  for(int layer=0; layer < numWeightLayers; layer++)
  {
    rsMatrixDbl *wm  = &(w[layer]); 

    // assign the irrelevant weights to zeros:
    for(int c=0; c<wm->getNumColumns(); c++)
    {
      //wm->m[0][c] = 0.0;
      wm->set(0, c, 0.0);
    }

    // retrieve the relevant weights from the vector:
    for(int r=1; r<wm->getNumRows(); r++)
    {
      for(int c=0; c<wm->getNumColumns(); c++)
      {
        //wm->m[r][c] = newWeightVector.v[i];
        wm->set(r, c, newWeightVector.v[i]);
        i++;
      }
    }
  }
}

void MultiLayerPerceptron::saveStateToFile(char *fileName)
{
  //....
}

void MultiLayerPerceptron::loadStateFromFile(char *fileName)
{
  //....
}

//-------------------------------------------------------------------------------------------------
// inquiry:

rsVectorDbl MultiLayerPerceptron::getWeightsAsVector()
{
  rsVectorDbl wv(numWeights);
  vectorizeWeightMatrices(w, &wv);
  return wv;
}

//-------------------------------------------------------------------------------------------------
// computation:

void MultiLayerPerceptron::setInput(const rsVectorDbl &xIn)
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

rsVectorDbl MultiLayerPerceptron::getOutput()
{
  return y;
}

void MultiLayerPerceptron::computeNetworkOutput(double *xIn, double *yOut)
{
  // todo: optimize this function - get rid of some overhead due to absorbing the biases into the
  // weight matrices
  setInput(rsVectorDbl(numInputs, xIn));
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
  double dummy = rsRandomUniform(min, max, seed); // init PNRG
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
    rsMatrixDbl *wm  = &(w[layer]);                             // weight matrix of layer 'layer'
    numWeights += (wm->getNumRows() - 1) * wm->getNumColumns(); // 1st row is irrelevant 
  }
}

void MultiLayerPerceptron::vectorizeWeightMatrices(rsMatrixDbl *weightMatrices, rsVectorDbl *weightVector)
{
  weightVector->setDimensionality(numWeights);
  int i = 0;
  for(int layer=0; layer < numWeightLayers; layer++)
  { 
    rsMatrixDbl *wm  = &(weightMatrices[layer]); 
    for(int r=1; r<wm->getNumRows(); r++)
    {
      for(int c=0; c<wm->getNumColumns(); c++)
      {
        //weightVector->v[i] = wm->m[r][c];
        weightVector->v[i] = (*wm)(r, c);
        i++;
      }
    }
  }
}

//=================================================================================================

double zeroValue(double x)    { return 0.0; }

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiLayerPerceptronErrorFunction::MultiLayerPerceptronErrorFunction(MultiLayerPerceptron *mlpToTrain)
{
  mlp = mlpToTrain;

  deltas = new rsVectorDbl[mlp->numWeightLayers]; 
  for(int i=0; i < mlp->numWeightLayers; i++)
    deltas[i] = rsVectorDbl(mlp->z[i+1].dim);

  dwn = new rsMatrixDbl[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwn[i] = rsMatrixDbl(mlp->w[i].getNumRows(), mlp->w[i].getNumColumns());

  dwN = new rsMatrixDbl[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwN[i] = rsMatrixDbl(mlp->w[i].getNumRows(), mlp->w[i].getNumColumns());

  xTrain = NULL;
  yTrain = NULL;
  nTrain = 0;
}

MultiLayerPerceptronErrorFunction::~MultiLayerPerceptronErrorFunction()
{
  delete[] dwn;
  delete[] dwN;
  delete[] deltas;
}

//-------------------------------------------------------------------------------------------------
// setup:

void MultiLayerPerceptronErrorFunction::setTrainingData(rsVectorDbl *inputs, rsVectorDbl *targets, int numPatterns)
{
  xTrain = inputs;
  yTrain = targets;
  nTrain = numPatterns;
}

//-------------------------------------------------------------------------------------------------
// overrides:

double MultiLayerPerceptronErrorFunction::getValue(rsVectorDbl p)
{
  rsVectorDbl pTmp = mlp->getWeightsAsVector();
  mlp->setWeightVector(p);
  double result = getTrainingError();
  mlp->setWeightVector(pTmp);
  return result;
}

rsVectorDbl MultiLayerPerceptronErrorFunction::getGradient(rsVectorDbl p)
{
  // set passed weight vector (remember current state for restoring it afterwards):
  rsVectorDbl pTmp = mlp->getWeightsAsVector();
  mlp->setWeightVector(p);   

  // compute gradient of the (batch-) error function
  rsVectorDbl g(mlp->numWeights);
  g.initWithZeros();
  for(int i=0; i<nTrain; i++)
  {
    mlp->setInput(xTrain[i]); 
    mlp->forwardPropagate();
    computePatternGradient(yTrain[i]);
    g += getPatternGradient();
  }
  g /= nTrain;

  // restore old weight vector and return computed gradient
  mlp->setWeightVector(pTmp); 
  return g;
}

//-------------------------------------------------------------------------------------------------
// others:

void MultiLayerPerceptronErrorFunction::printPatternGradient()
{
  printf("%s", "Multilayer Perceptron Trainer - gradient matrices:  \n");
  printf("%s", "\n");
  for(int i=0; i<mlp->numWeightLayers; i++)
  {
    dwn[i].print();
    printf("%s", "\n");
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MultiLayerPerceptronErrorFunction::computePatternGradientByWeightPerturbation(const rsVectorDbl& yTarget)
{
  mlp->forwardPropagate();  // to obtain y, in case this has not already been done
  double eps   = 0.00001;   // epsilon
  double error = getPatternError(yTarget);
  double errorPerturbed, tmp;
  rsMatrixDbl *w, *dw1;
  for(int layer = 0; layer < mlp->numWeightLayers; layer++)
  {
    w   = &(mlp->w[layer]);
    dw1 = &(   dwn[layer]);
    for(int i=0; i<w->getNumRows(); i++)
    {
      for(int j=0; j<w->getNumColumns(); j++)
      {
        tmp             = (*w)(i,j);

        //(*w)(i,j)      += eps;
        w->set(i, j, tmp+eps);

        mlp->forwardPropagate();
        errorPerturbed  = getPatternError(yTarget);

        //(*dw1)(i,j)     = (errorPerturbed - error) / eps;  // (1),Eq.4.46
        dw1->set(i,j, (errorPerturbed-error)/eps);         // (1),Eq.4.46

        //(*w)(i,j)       = tmp;                             // restore original weight
        w->set(i, j, tmp);                                 // restore original weight
      }
    }
  }
  mlp->forwardPropagate();      // restore the y with unperturbed weights
}

void MultiLayerPerceptronErrorFunction::computePatternGradient(const rsVectorDbl& yTarget)
{
  computeDeltas(yTarget);
  for(int layer=0; layer<mlp->numWeightLayers; layer++)
  {
    rsVectorDbl *z   = &(mlp->z[layer]);          // current layer of input nodes
    rsVectorDbl *d   = &(deltas[layer]);          // current layer of deltas
    rsMatrixDbl *dw1 = &(   dwn[layer]);          // current matrix of weight derivatives
    for(int j=0; j<dw1->getNumRows(); j++)        // j is the output node index
    { 
      for(int i=0; i<dw1->getNumColumns(); i++)   // i is the input node index
      {
        //(*dw1)(j,i) = d->v[j] * z->v[i];
        dw1->set(j,i, d->v[j] * z->v[i]);
      }
    }
  }
}

void MultiLayerPerceptronErrorFunction::computeDeltas(const rsVectorDbl& yTarget)
{
  int numNeurons;     // number of neurons (excluding bias) of current layer
  int layer, j, k;    // indices
  rsVectorDbl *dc;         // vector of the deltas that are currently being computed
  rsVectorDbl *da;         // vector of the deltas that are already available and currently used
  rsMatrixDbl *w;          // layer of weights that is currently worked with
  rsVectorDbl *z;          // vector of nodes for which the deltas will be computed

  // compute deltas for the (linear) output-layer:
  layer      = mlp->numWeightLayers-1; // -1, because there are no deltas for the input nodes
  dc         = &deltas[layer];         // we now compute deltas for the output layer 
  numNeurons = dc->dim - 1;            // -1, because the 0th node is a (dummy) bias-node
  dc->v[0]   = 0.0;                    // delta for a bias node is zero
  for(k=0; k<numNeurons; k++)  
  {
    dc->v[k+1] = mlp->y.v[k] - yTarget.v[k];     // (1),Eq.4.41
      // this has to be modified for nonlinear output units and/or different error functions - the 
      // general expression would be (1),Eq.4.30: dE/dy[k] * g'(a[k])
  }

  // recursively compute the deltas for the hidden layers by means of back-propagating the deltas 
  // of the layer right to the layer in question - the delta for one unit with index j in some 
  // layer is the weighted sum over all deltas (indexed by k) in the subsequent layer multiplied 
  // by the derivative of the activation of unit j - the weights are given by the synaptic 
  // connections between the source neuron with index j and the target neuron with index k - this 
  // is the element w[k][j] in the weight-matrix between the successive layers
  int numSourceNodes, numTargetNodes;
  while( layer > 0 )
  {
    layer--;                              // layer for which the deltas will be computed now
    da             = dc;                  // available deltas most recently computed
    dc             = &deltas[layer];      // deltas to be computed
    z              = &(mlp->z[layer+1]);  // nodes for which we compute the delta
    w              = &(mlp->w[layer+1]);  // weights from z to the subsequent layer
    numSourceNodes = dc->dim;             // number of nodes for which to compute the delta
    numTargetNodes = da->dim;             // number of deltas to sum over
    dc->v[0]       = 0.0;                 // delta for a bias node is zero
    for(j=1; j<numSourceNodes; j++)  
    {
      // form the weighted sum, k runs over targets nodes, start at k=1 because 
      // da->v[k]==0 for k==0:
      double sum = 0.0;
      for(k=1; k<numTargetNodes; k++) 
        sum += (*w)(k,j) * da->v[k];

      // multiply the sum by the derivative of the activation to obtain the delta, 
      // this implements (1),Eq.4.36:
      dc->v[j] = activationDerivative(z->v[j]) * sum;   
    }
  }
}

double MultiLayerPerceptronErrorFunction::getPatternError(const rsVectorDbl& yTarget)
{
  // \todo: implement other error functions by using a switch statement
  return 0.5 * (mlp->y - yTarget).getSquaredNorm();
}

double MultiLayerPerceptronErrorFunction::getTrainingError()
{
  double sum = 0.0;
  for(int p=0; p<nTrain; p++)
  {
    mlp->setInput(xTrain[p]); 
    mlp->forwardPropagate();
    sum += getPatternError(yTrain[p]);
  }
  sum /= nTrain;
  return sum;
}

rsVectorDbl MultiLayerPerceptronErrorFunction::getPatternGradient()
{
  rsVectorDbl gv(mlp->numWeights);
  mlp->vectorizeWeightMatrices(dwn, &gv);
  return gv;
}

//=================================================================================================

// construction/destruction:

MultiLayerPerceptronTrainer::MultiLayerPerceptronTrainer(MultiLayerPerceptron *mlpToTrain)
{
  mlp = mlpToTrain;

  deltas = new rsVectorDbl[mlp->numWeightLayers]; 
  for(int i=0; i < mlp->numWeightLayers; i++)
    deltas[i] = rsVectorDbl(mlp->z[i+1].dim);

  dwn = new rsMatrixDbl[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwn[i] = rsMatrixDbl(mlp->w[i].getNumRows(), mlp->w[i].getNumColumns());

  dwN = new rsMatrixDbl[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwN[i] = rsMatrixDbl(mlp->w[i].getNumRows(), mlp->w[i].getNumColumns());

  xTrain = NULL;
  yTrain = NULL;
  nTrain = 0;
}

MultiLayerPerceptronTrainer::~MultiLayerPerceptronTrainer()
{
  delete[] dwn;
  delete[] dwN;
  delete[] deltas;
}

//-------------------------------------------------------------------------------------------------
// setup:

void MultiLayerPerceptronTrainer::setTrainingData(rsVectorDbl *inputs, rsVectorDbl *targets, int numPatterns)
{
  xTrain = inputs;
  yTrain = targets;
  nTrain = numPatterns;
}

//-------------------------------------------------------------------------------------------------
// inquiry:


//-------------------------------------------------------------------------------------------------
// others:

void MultiLayerPerceptronTrainer::initializeWeightsToZeros()
{
  
}

void MultiLayerPerceptronTrainer::initializeWeightsRandomly(double min, double max, int seed)
{

}

void MultiLayerPerceptronTrainer::printPatternGradient()
{
  printf("%s", "Multilayer Perceptron Trainer - gradient matrices:  \n");
  printf("%s", "\n");
  for(int i=0; i<mlp->numWeightLayers; i++)
  {
    dwn[i].print();
    printf("%s", "\n");
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

void MultiLayerPerceptronTrainer::computePatternGradientByWeightPerturbation(const rsVectorDbl& yTarget)
{
  mlp->forwardPropagate();  // to obtain y, in case this has not already been done
  double eps   = 0.000001;  // epsilon
  double error = computePatternError(yTarget);
  double errorPerturbed, tmp;
  rsMatrixDbl *w, *dw1;
  for(int layer = 0; layer < mlp->numWeightLayers; layer++)
  {
    w   = &(mlp->w[layer]);
    dw1 = &(   dwn[layer]);
    for(int i=0; i<w->getNumRows(); i++)
    {
      for(int j=0; j<w->getNumColumns(); j++)
      {
        tmp             = (*w)(i,j);

        //(*w)(i,j)      += eps;
        w->set(i,j, tmp+eps);

        mlp->forwardPropagate();
        errorPerturbed  = computePatternError(yTarget);

        //(*dw1)(i,j)     = (errorPerturbed - error) / eps;     // (1),Eq.4.46
        dw1->set(i,j, (errorPerturbed-error)/eps);            // (1),Eq.4.46

        //(*w)(i,j)       = tmp;   // restore original weight
        w->set(i,j, tmp);        // restore original weight
      }
    }
  }
  mlp->forwardPropagate();      // restore the y with unperturbed weights
}

void MultiLayerPerceptronTrainer::computePatternGradient(const rsVectorDbl& yTarget)
{
  computeDeltas(yTarget);
  for(int layer=0; layer<mlp->numWeightLayers; layer++)
  {
    rsVectorDbl *z   = &(mlp->z[layer]);          // current layer of input nodes
    rsVectorDbl *d   = &(deltas[layer]);          // current layer of deltas
    rsMatrixDbl *dw1 = &(   dwn[layer]);          // current matrix of weight derivatives
    for(int j=0; j<dw1->getNumRows(); j++)        // j is output node index
    { 
      for(int i=0; i<dw1->getNumColumns(); i++)   // i is input node index
      {
        //(*dw1)(j,i) = d->v[j] * z->v[i];
        dw1->set(j,i, d->v[j] * z->v[i]);
      }
    }
  }
}

void MultiLayerPerceptronTrainer::computeDeltas(const rsVectorDbl& yTarget)
{
  int numNeurons;     // number of neurons (excluding bias) of current layer
  int layer, /*i,*/ j, k; // indices
  rsVectorDbl *dc;    // vector of the deltas that are currently being computed
  rsVectorDbl *da;    // vector of the deltas that are already available and currently used
  rsMatrixDbl *w;     // layer of weights that is currently worked with
  rsVectorDbl *z;     // vector of nodes for which the deltas will be computed

  // compute deltas for the (linear) output-layer:
  layer      = mlp->numWeightLayers-1; // -1, because there are no deltas for the input nodes
  dc         = &deltas[layer];         // we now compute deltas for the output layer 
  numNeurons = dc->dim - 1;            // -1, because the 0th node is a (dummy) bias-node
  dc->v[0]   = 0.0;                    // delta for a bias node is zero
  for(k=0; k<numNeurons; k++)  
  {
    //double delta = mlp->y->v[k] - yTarget[k-1]; 
    //dc->v[k] = mlp->y->v[k] - yTarget[k-1];     // (1),Eq.4.41
    //double delta = mlp->y->v[k] - yTarget[k]; 
    dc->v[k+1] = mlp->y.v[k] - yTarget.v[k];     // (1),Eq.4.41
      // this has to be modified for nonlinear output units and/or different error functions - the 
      // general expression would be (1),Eq.4.30: dE/dy[k] * g'(a[k])
  }

  // recursively compute the deltas for the hidden layers by means of back-propagating the deltas 
  // of the layer right to the layer in question - the delta for one unit with index j in some 
  // layer is the weighted sum over all deltas (indexed by k) in the subsequent layer multiplied 
  // by the derivative of the activation of unit j - the weights are given by the synaptic 
  // connections between the source neuron with index j and the target neuron with index k - this 
  // is the element w[k][j] in the weight-matrix between the successive layers
  int numSourceNodes, numTargetNodes;
  while( layer > 0 )
  {
    da = dc;                         // available deltas most recently computed
    layer--;                         // layer for which the deltas will be computed now
    dc = &deltas[layer];             // deltas to be computed
    z  = &(mlp->z[layer+1]);         // nodes for which we compute the delta
    w  = &(mlp->w[layer+1]);         // weights from z to the subsequent layer
    numSourceNodes = dc->dim;        // number of nodes for which to compute the delta
    numTargetNodes = da->dim;        // number of deltas to sum over
    dc->v[0]   = 0.0;                // delta for a bias node is zero
    for(j=1; j<numSourceNodes; j++)  
    {
      // form the weighted sum:
      double sum = 0.0;
      for(k=1; k<numTargetNodes; k++)  // k runs over targets nodes, start at k=1 because da->v[k]==0 for k==0
        sum += (*w)(k,j) * da->v[k];

      // multiply the sum by the derivative of the activation to obtain the delta:
      double factor = activationDerivative(z->v[j]);
      dc->v[j]      = factor * sum;           // (1),Eq.4.36    
    }
  }
}

double MultiLayerPerceptronTrainer::computePatternError(const rsVectorDbl& yTarget)
{
  // todo: implement other error functions by using a switch statement
  double result = 0.0;
  double tmp;
  for(int i=0; i<mlp->numOutputs; i++)
  {
    tmp     = mlp->y.v[i] - yTarget.v[i];
    tmp    *= tmp;
    result += tmp;
  }
  return 0.5*result;
}

rsVectorDbl MultiLayerPerceptronTrainer::getPatternGradient()
{
  rsVectorDbl gv(mlp->numWeights);
  mlp->vectorizeWeightMatrices(dwn, &gv);
  return gv;
}

