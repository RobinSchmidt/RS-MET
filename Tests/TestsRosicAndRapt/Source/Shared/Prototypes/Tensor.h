#ifndef Tensor_h
#define Tensor_h

#include <vector>

/** Implements a simple proof-of-concept class for tensors */

class Tensor
{

public:


  int getRank() { return rank; }
  
  int getNumDimensions() { return dim; }
  
  bool isVector()   { return rank == 1 && subscript[0] == false; }
  
  bool isCovector() { return rank == 1 && subscript[0] == true;  }

  
  
  Tensor contract(const Tensor& A, int index);
  
  Tensor outerProduct(const Tensor& A, const Tensor& B);
  
  Tensor innerProduct(const Tensor& A, const Tensor& B);
  
  // operators: +,-,==,unary-, element access (i,j,k,...)
  

protected:

  int rank;                    // number of indices of the tensor
  int dim;                     // number of dimensions of the (co)vector space
  std::vector<double> v;       // holds the values as flat array (size = dim^rank)
  std::vector<bool> subscript; // indicates if a particular index is a subscript
  
};

#endif