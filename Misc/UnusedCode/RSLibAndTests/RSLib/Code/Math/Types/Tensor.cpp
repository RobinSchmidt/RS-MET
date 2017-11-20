using namespace RSLib;

// construction/destruction:

rsTensor::rsTensor()   
{
  N = 0;
  R = 0;
  c = NULL;
}

rsTensor::rsTensor(int numDimensions, int numIndices, bool initWithZeros)
{
  int C;

  // wrap these 4 lines into a function allocateMemory(N, R)
  N = numDimensions;
  R = numIndices;
  C = getNumComponents();
  c = new double[C];

  if( initWithZeros == true )
    rsFillWithZeros();
}
  
rsTensor::rsTensor(int numDimensions, int numIndices, double *components)
{
  int C;

  N = numDimensions;
  R = numIndices;
  C = getNumComponents();
  c = new double[C];

  memcpy(c, components, C*sizeof(double));
}

rsTensor::rsTensor(const rsTensor& other)
{
  int C;

  N = other.N;
  R = other.R;
  C = getNumComponents();
  c = new double[C];

  memcpy(c, other.c, C*sizeof(double));  
}
    
rsTensor::~rsTensor()
{
  delete[] c;
}

// inquiry:
   
int rsTensor::getNumComponents() const
{
  return (int) pow((double)N, R);  // \todo replace with a call to some powInt function
}