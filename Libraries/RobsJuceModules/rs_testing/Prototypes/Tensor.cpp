#include "Tensor.h"

Tensor::Tensor(int numDimensions, int numSuperscripts, int numSubscripts)
{
  dim  = numDimensions;
  rank = numSuperscripts + numSubscripts;
  subscript.resize(rank);
  values.resize((int)pow(dim, rank)); // maybe use some powInt function
}

bool Tensor::isSymmetric(int /*index1*/, int /*index2*/) const
{
  // we must run over i,j in 0...dim-1 of this(...,i,...,j,...) where i and j are at positions 
  // index1, index2 and compare all the elements - they must all be equal...how to express this?

  // maybe a regular nested loop (over i, j) but with somehow computed start and stride values for 
  // i and j

  // oh - wait - we have to run over all elements where the dots (dummy indices) are..or not?


  return false;
}

/*
Background:

Terms:
-scalar:
-vector:
-covector:
-tensor:
-heterogenous scalar product: sum-of-products of components of a vector and a covector
-homogenous scalar product: sum-of-products of components of a vector (or covector) and its dual 
 covector (or vector)
-dual-switch tensor:
-metric tensor


todo: check this:
https://github.com/QuantStack/xtensor/blob/master/notebooks/xtensor.ipynb


https://en.wikipedia.org/wiki/Tensor


Tensors for Beginneres by eigenchris:
https://www.youtube.com/watch?v=uDRzJIaN2qw&list=PLJHszsWbB6hrkmmq57lX8BV-o-YIOFsiG
-very good explanations, easy to follow


*/