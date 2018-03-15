#include "Tensor.h"

Tensor::Tensor(int numDimensions, int numSuperscripts, int numSubscripts)
{
  dim  = numDimensions;
  rank = numSuperscripts + numSubscripts;
  subscript.resize(rank);
  values.resize((int)pow(dim, rank)); // maybe use some powInt function
}

bool Tensor::isSymmetric(int index1, int index2) const
{
  // we must run over i,j in 0...dim-1 of this(...,i,...,j,...) where i and j are at positions 
  // index1, index2 and compare all the elements - they must all be equal...how to express this?
  // oh - wait - we have to run over all elements where the dots (dummy indices) are..or not?

  return false;
}