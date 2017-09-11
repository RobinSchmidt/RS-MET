#ifndef RAPT_MATRIX_H
#define RAPT_MATRIX_H


/** This is a class for treating C-arrays as matrices. It does not store/own the actual matrix 
data, it just acts as wrapper around an existing array for more conveniently accessing and 
manipulating matrix elements.

*/

template<class T>
class rsMatrixView
{

public:


  /** \name Operators */

  /** Read and write access to matrix elements. */
  T& operator()(const int i, const int j);



protected:

  /** \name Data */

  size_t N, M;    // number of rows and columns
  T *d;           // data pointer - factor out to baseclass

};

//=================================================================================================

/** This is a class for representing matrices and doing mathematical operations with them. */

template<class T>
class rsMatrix : public rsMatrixView<T> /*, public std::vector<T>*/
{

public:

  /** \name Construction/Destruction */


  /** Standard constructor. You must pass the initial number of rows and columns */
  rsMatrix(size_t numRows, size_t numColumns);
 
  /** Copy constructor. */
  rsMatrix(const rsMatrix& other);
                              
  /** Move constructor. */
  rsMatrix(const rsMatrix&& other);


  /** Destructor. */
  ~rsMatrix();


    
  /** \name Manipulations */


  /** \name Inquiry */

    



  /** \name Decompositions */



  /** \name Data */

  std::vector<T> data;
  
}; 

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix

  // matrix/vector functions:

#endif
