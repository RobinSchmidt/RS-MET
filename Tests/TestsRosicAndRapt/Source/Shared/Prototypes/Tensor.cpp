#include "Tensor.h"


bool Tensor::isSymmetric(int index1, int index2)
{
  // we must run over i,j in 0...dim-1 of this(...,i,...,j,...) where i and j are at positions 
  // index1, index2 and compare all the elements - they must all be equal...how to express this?
  // oh - wait - we have to run over all elements where the dots (dummy indices) are..or not?

  return false;
}