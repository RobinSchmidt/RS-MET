//#include "rosic_MultiLayerPerceptronErrorFunction.h"
//using namespace rosic;

double zeroValue(double)    { return 0.0; }

//-------------------------------------------------------------------------------------------------
// construction/destruction:

MultiLayerPerceptronErrorFunction::MultiLayerPerceptronErrorFunction(MultiLayerPerceptron *mlpToTrain)
{
  mlp = mlpToTrain;

  deltas = new Vector[mlp->numWeightLayers];
  for(int i=0; i < mlp->numWeightLayers; i++)
    deltas[i] = Vector(mlp->z[i+1].dim);

  dwn = new Matrix[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwn[i] = Matrix(mlp->w[i].numRows, mlp->w[i].numColumns);

  dwN = new Matrix[mlp->numWeightLayers];
  for(int i=0; i<mlp->numWeightLayers; i++)
    dwN[i] = Matrix(mlp->w[i].numRows, mlp->w[i].numColumns);

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

void MultiLayerPerceptronErrorFunction::setTrainingData(rosic::Vector *inputs,
  rosic::Vector *targets, int numPatterns)
{
  xTrain = inputs;
  yTrain = targets;
  nTrain = numPatterns;
}

//-------------------------------------------------------------------------------------------------
// overrides:

double MultiLayerPerceptronErrorFunction::getValue(rosic::Vector p)
{
  rosic::Vector pTmp = mlp->getWeightsAsVector();
  mlp->setWeightVector(p);
  double result = getTrainingError();
  mlp->setWeightVector(pTmp);
  return result;
}

rosic::Vector MultiLayerPerceptronErrorFunction::getGradient(rosic::Vector p)
{
  // set passed weight vector (remember current state for restoring it afterwards):
  Vector pTmp = mlp->getWeightsAsVector();
  mlp->setWeightVector(p);

  // compute gradient of the (batch-) error function
  Vector g(mlp->numWeights);
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

void MultiLayerPerceptronErrorFunction::computePatternGradientByWeightPerturbation(
  const rosic::Vector& yTarget)
{
  mlp->forwardPropagate();  // to obtain y, in case this has not already been done
  double eps   = 0.00001;   // epsilon
  double error = getPatternError(yTarget);
  double errorPerturbed, tmp;
  Matrix *w, *dw1;
  for(int layer = 0; layer < mlp->numWeightLayers; layer++)
  {
    w   = &(mlp->w[layer]);
    dw1 = &(   dwn[layer]);
    for(int i=0; i<w->numRows; i++)
    {
      for(int j=0; j<w->numColumns; j++)
      {
        tmp             = w->m[i][j];
        w->m[i][j]     += eps;
        mlp->forwardPropagate();
        errorPerturbed  = getPatternError(yTarget);
        dw1->m[i][j]    = (errorPerturbed - error) / eps;  // (1),Eq.4.46
        w->m[i][j]      = tmp;                             // restore original weight
      }
    }
  }
  mlp->forwardPropagate();      // restore the y with unperturbed weights
}

void MultiLayerPerceptronErrorFunction::computePatternGradient(const rosic::Vector& yTarget)
{
  computeDeltas(yTarget);
  for(int layer=0; layer<mlp->numWeightLayers; layer++)
  {
    Vector *z   = &(mlp->z[layer]);          // current layer of input nodes
    Vector *d   = &(deltas[layer]);          // current layer of deltas
    Matrix *dw1 = &(   dwn[layer]);          // current matrix of weight derivatives
    for(int j=0; j<dw1->numRows; j++)        // j is the output node index
    {
      for(int i=0; i<dw1->numColumns; i++)   // i is the input node index
        dw1->m[j][i] = d->v[j] * z->v[i];
    }
  }
}

void MultiLayerPerceptronErrorFunction::computeDeltas(const rosic::Vector& yTarget)
{
  int numNeurons;     // number of neurons (excluding bias) of current layer
  int layer, j, k;    // indices
  Vector *dc;         // vector of the deltas that are currently being computed
  Vector *da;         // vector of the deltas that are already available and currently used
  Matrix *w;          // layer of weights that is currently worked with
  Vector *z;          // vector of nodes for which the deltas will be computed

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
        sum += w->m[k][j] * da->v[k];

      // multiply the sum by the derivative of the activation to obtain the delta,
      // this implements (1),Eq.4.36:
      dc->v[j] = activationDerivative(z->v[j]) * sum;
    }
  }
}

double MultiLayerPerceptronErrorFunction::getPatternError(const rosic::Vector& yTarget)
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

rosic::Vector MultiLayerPerceptronErrorFunction::getPatternGradient()
{
  rosic::Vector gv(mlp->numWeights);
  mlp->vectorizeWeightMatrices(dwn, &gv);
  return gv;
}
