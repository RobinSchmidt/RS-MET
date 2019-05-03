//#include "rosic_Vector.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

Vector::Vector()
{
  dim = 0;
  v   = NULL;
}

Vector::Vector(int numElements)
{
  dim = numElements;
  v   = new double[dim];
}

Vector::Vector(int numElements, double *values)
{
  dim = numElements;
  v   = new double[dim];
  memcpy(v, values, dim*sizeof(double));
}

Vector::Vector(const Vector& v2)
{
  dim = v2.dim;
  v   = new double[dim];
 
  // copy the values:  
  memcpy(v, v2.v, dim*sizeof(double));
}

Vector::~Vector()
{
  delete[] v;
}

//-------------------------------------------------------------------------------------------------
// setup:


//-------------------------------------------------------------------------------------------------
// inquiry:

    
//-------------------------------------------------------------------------------------------------  
// others:

void Vector::print()
{
  printf("%s %d %s", "Vector - dimensionality: ", dim, "\n");
  for(int i=0; i<dim; i++)
    printf("%.4f %s", v[i], "  ");
  printf("%s", "\n");
}