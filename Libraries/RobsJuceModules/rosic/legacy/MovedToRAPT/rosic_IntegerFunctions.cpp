//#include "rosic_IntegerFunctions.h"
//using namespace rosic;

void rosic::getLineOfPascalTriangle(unsigned int *c, unsigned int n)
{
  RAPT::rsArrayTools::fillWithZeros(c, n+1);  
  c[0] = 1;
  for(unsigned int i = 0; i <= n; i++)
  {
    for(unsigned int j = i; j >= 1; j--)
      c[j] = c[j-1] + c[j];
  }
}

unsigned int rosic::binomialCoefficient(unsigned int n, unsigned int k)
{
  if( k == 0 || k == n )
    return 1;
  else
  {
    unsigned int *c = new unsigned int[n+1];
    getLineOfPascalTriangle(c, n);
    unsigned int result = c[k];
    delete[] c;
    return result;
  }
}

unsigned int rosic::binomialCoefficientUpTo20(unsigned int n, unsigned int k)
{
  rassert( n <= 20 );   // can only be used for n <= 20, otherwise internal overflow occurs
  unsigned long long nL  = n;
  unsigned long long kL  = k;
  return (unsigned int) (product(kL+1ULL, nL) / product(1ULL, nL-kL));
}


