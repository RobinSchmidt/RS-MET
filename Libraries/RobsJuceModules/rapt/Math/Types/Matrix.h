#ifndef RAPT_MATRIX_H
#define RAPT_MATRIX_H


/*
template<class T>
class rsMatrix2x2
{

};
*/


//=================================================================================================

/** This is a class for treating C-arrays as matrices. It does not store/own the actual matrix 
data, it just acts as wrapper around an existing array for more conveniently accessing and 
manipulating matrix elements.

*/

template<class T>
class rsMatrixView
{

public:

  /** \name Construction/Destruction */




  /**  */
  rsMatrixView(int numRows = 0, int numColumns = 0, T* data = nullptr)
  {
    this->numRows = numRows;
    this->numCols = numColumns;
    d = data;
  }

  /** \name Setup */

  inline void setAllValues(T value)
  {
    rsArray::fillWithValue(d, int(numRows * numCols), value);
  }

  inline void reshape(int newNumRows, int newNumColumns)
  {
    rsAssert(newNumRows*newNumColumns == numRows*numCols);
    numRows = newNumRows;
    numCols = newNumColumns;
  }

  // void setToIdentityMatrix(T scaler = 1);


  /** \name Operators */

  /** Read and write access to matrix elements. */
  inline T& operator()(const int i, const int j)
  {
    return d[numCols*i + j];
    // more general colStride*i + rowStride*j
    // even more general: colOffset + colStride*i + (rowOffset + rowStride)*j
  }



protected:

  /** \name Data */

  //size_t N, M;    // number of rows and columns

  int numRows, numCols;
  T *d;           // data pointer

};

//=================================================================================================

/** This is a class for representing matrices and doing mathematical operations with them. */

template<class T>
class rsMatrixNew : public rsMatrixView<T>
{

public:

  /** \name Construction/Destruction */


  /** Standard constructor. You must pass the initial number of rows and columns */
  rsMatrixNew(int numRows = 0, int numColumns = 0);
  //rsMatrix(size_t numRows = 1, size_t numColumns = 1); // leads to memory leaks

 
  /** Copy constructor. */
  //rsMatrix(const rsMatrix& other);
                              
  /** Move constructor. */
  //rsMatrix(const rsMatrix&& other);


  /** Destructor. */
  ~rsMatrixNew() {}


  /** \name Setup */

  /** Sets the number of rows and columns, this matrix should have. ToDo: provide a way to retain 
  the data (optionally) - what does std::vector's resize do? Does it retain data...but if it does,
  it would be useless anyway in case the number of columns changed. */
  void setSize(int numRows, int numColumns);

    
  /** \name Manipulations */


  /** \name Inquiry */

    



  /** \name Decompositions */



  /** \name Data */

  std::vector<T> data;
  
}; 

  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsMatrix

  // matrix/vector functions:

  // todo: make some special-case classes for 2x2, 3x3 matrices which can use simpler algorithms
  // for some of the computations

#endif
