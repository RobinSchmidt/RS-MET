//#include "rosic_RoutingMatrix.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

RoutingMatrix::RoutingMatrix(int numInputs, int numOutputs)
{
  this->numInputs  = numInputs;
  this->numOutputs = numOutputs;

  flatArray = new double[numInputs*numOutputs];
  mixMatrix = new double*[numInputs];
  for(int i=0; i<numInputs; i++)
    mixMatrix[i] = &flatArray[i*numOutputs];

  initializeMatrix();
}

RoutingMatrix::~RoutingMatrix()
{
  delete[] mixMatrix;
  delete[] flatArray;
}

//-------------------------------------------------------------------------------------------------
// parameter settings:

void RoutingMatrix::setMatrixEntry(int i, int o, double value)    
{ 
  if( i >=0 && i < numInputs && o >= 0 && o < numOutputs )
    mixMatrix[i][o] = value; 
}

//-------------------------------------------------------------------------------------------------
// inquiry:

double RoutingMatrix::getMatrixEntry(int i, int o)    
{ 
  if( i >=0 && i < numInputs && o >= 0 && o < numOutputs )
    return mixMatrix[i][o]; 
  else
    return 0.0;
}

//-------------------------------------------------------------------------------------------------
// others:

void RoutingMatrix::initializeMatrix()
{
  for(int i=0; i<numInputs; i++)
  {
    for(int o=0; o<numOutputs; o++)
      mixMatrix[i][o] = 0.0;
  }
}