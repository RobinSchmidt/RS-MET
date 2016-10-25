#ifndef RAPT_MATRIX_H
#define RAPT_MATRIX_H

namespace RAPT
{


/**

This is a class for treating C-arrays as matrices. It does not store/own the actual matrix data, 
it just acts as wrapper around an existing array for more conveniently accessing and manipulating 
matrix elements.

*/

template<class T>
class MatrixView
{

public:


  /** \name Operators */

  T& operator()(const int i, const int j) const; //< element access


  /** \name Data */

  size_t N, M;    // number of rows and columns
  T *d            // data pointer - factor out to baseclass

}


/**

This is a class for representing matrices and doing mathematical operations with them. 

*/

template<class T>
class Matrix : public MatrixView<T>, public vector<T>
{

public:

  /** \name Construction/Destruction */

  Matrix(const Matrix& other); // copy constructor
                               // move constructor
  ~Matrix();


    
  /** \name Manipulations */


  /** \name Inquiry */

    



  /** \name Decompositions */



  /** \name Data */

  //vector<T> v;    // data
  
}; 

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix

  // matrix/vector functions:

}

#endif
