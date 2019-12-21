#ifndef RAPT_MULTIARRAY_OLD_H
#define RAPT_MULTIARRAY_OLD_H

// this implementation is obsolete - there's a new implementation in the Data folder which uses 
// more modern programming techniques

/** This is a class for representing and manipulating multi-dimensional arrays. An 1D array may be
seen as a vector, a 2D array as a matrix, a 3D array as an array of matrices sliced one after
another in a 3rd dimension (constituting a 3D block), a 4D array as a set of such blocks, etc.

Elements are accessed by a C-style array of integer numbers representing the indexes in the
various dimensions.

The ranges for the different indices may be different. If they happen to be all equal, giving a
kind of hyper-cubic array (square-matrix in 2D case, cube in 3D case, etc.), the object may also
be used to represent tensors.  */

template<class ElementType>
class rsMultiArrayOld
{

public:

  /** \name Construction/Destruction */

  /** Default constructor. Creates a multiarray with 0 indices - that is: a scalar - and
  initializes this scalar to zero.
  \todo: make this initialization optional
  */
  rsMultiArrayOld();

  /** Constructor for constructing a multiarray of the desired dimensions with
  optional initialization of the elements with zero. */
  rsMultiArrayOld(const rsUint32 numIndices, const rsUint32 indexRanges[], bool initWithZeros = true);

  /** Copy constructor. */
  rsMultiArrayOld(const rsMultiArrayOld& other);


  /** Destructor. */
  ~rsMultiArrayOld();


  /** \name Operators */

  //rsMultiArray& operator=(const rsMultiArray& a);

  bool operator==(const rsMultiArrayOld& a) const;
  bool operator!=(const rsMultiArrayOld& a) const;

  /*
  rsMultiArray operator-();

  rsMultiArray& operator+=(const ElementType &x);
  rsMultiArray& operator-=(const ElementType &x);
  rsMultiArray& operator*=(const ElementType &x);
  rsMultiArray& operator/=(const ElementType &x);

  rsMultiArray& operator+=(const rsMultiArray &a);
  rsMultiArray& operator-=(const rsMultiArray &a);

  rsMultiArray operator+(const ElementType &x);
  rsMultiArray operator-(const ElementType &x);
  rsMultiArray operator*(const ElementType &x);
  rsMultiArray operator/(const ElementType &x);

  rsMultiArray operator+(const rsMultiArray &a);
  rsMultiArray operator-(const rsMultiArray &a);
  */


  /** \name Manipulations */

  /** Sets a new size [....]

  Should re-allocate only if necessary - if not necessary, it just re-indexes the already existing
  data.
  */
  void setSize(const rsUint32 newNumIndices, const rsUint32 newIndexRanges[]);

  /** Applies the passed function to each element of the multi-array. */
  void applyFunction(ElementType (*f) (ElementType));


  /** Sets all the values in the multi-array to the passed value (which is optional and defaults
  to zero). */
  void fillWithValue(const ElementType &value = ElementType(0));

  /** Sets the element at the given set of indices to the given value. */
  void setElement(const rsUint32 indices[], const ElementType &value);

  /** Copies the data from the passed flat array into our data-array member. The passed array
  should be at least as long as the number of elements in this object (given by the product of
  the index ranges). */
  void setDataFromFlatArray(const ElementType flatArray[]);


  /** \name Inquiry */

  /** Returns the offset in our flat data array given an array of N indices where N should be
  equal to the number of indices that this object has. No range-checking is done on the
  index-values themselves. */
  rsUint32 offsetFromIndices(const rsUint32 indices[]) const;

  /** Given an offset into our flat data array, this function fills the "indices" array
  (assumed to be of length given by our numIndices member) with the corresponding indices
  for the multidimensional array. */
  void indicesFromOffset(rsUint32 offset, rsUint32 indices[]) const;

  /** Returns the value at the given set of indices. */
  ElementType getElement(const rsUint32 indices[]) const;

  /** Returns the total number of elements in this multiarray */
  rsUint32 getNumElements() const;

  /** Fills the passed array with the data from this object. */
  void getDataAsFlatArray(ElementType flatArray[]) const;

  /** Returns true if this MultiArray has the same number of indices and the same index-ranges
  as the passed parameter. */
  bool isOfSameTypeAs(const rsMultiArrayOld<ElementType> &other) const;

  // getMin, getMax, getNormL1, getNormL2, etc.


  /** \name Algebraic Operations */

  /** Adds two multi-arrays (assumed to be of the same type) element-wise. */
  static void add(const rsMultiArrayOld<ElementType> &left,
    const rsMultiArrayOld<ElementType> &right, rsMultiArrayOld<ElementType> &result);

  /** Subtracts two multi-arrays (assumed to be of the same type) element-wise. */
  static void subtract(const rsMultiArrayOld<ElementType> &left,
    const rsMultiArrayOld<ElementType> &right, rsMultiArrayOld<ElementType> &result);

  /** Forms an outer product of the two given multiarrays. As an example, let A be an array with
  3 indices, B be an array with 2 indices. The outer product array C = A*B is an array with
  3 + 2 = 5 indices with elements given by:
  C_ijkpq = A_ijk * B_pq
  Likewise, the product C = B*A has also 2 + 3 = 5 elements given by:
  C_pqijk = B_pq  * A_ijk
  from which we see that the outer product is not commutative. A*B and B*A have the same elements
  but in a different order. */
  static void outerProduct(const rsMultiArrayOld<ElementType> &left,
    const rsMultiArrayOld<ElementType> &right, rsMultiArrayOld<ElementType> &result);

  /** Given a multi-array B ("rightFactor") and a multi-array C ("outerProduct") that is assumed
  to be is an outer product of A and B, such that C = A * B, this function retrieves the left
  factor A ("result"). */
  static void leftFactor(const rsMultiArrayOld<ElementType> &outerProduct,
    const rsMultiArrayOld<ElementType> &rightFactor, rsMultiArrayOld<ElementType> &result);

  /** Given a multi-array A ("leftFactor") and a multi-array C ("outerProduct") that is assumed
  to be is an outer product of A and B, such that C = A*B, this function retrieves the right
  factor B ("result"). */
  static void rightFactor(const rsMultiArrayOld<ElementType> &outerProduct,
    const rsMultiArrayOld<ElementType> &leftFactor, rsMultiArrayOld<ElementType> &result);

  /** Performs a contraction with respect to a pair of indices.

  B = contract(A_ipjqk, 1, 3)

  B_ijk = sum_p A_ipjpk
  */
  static rsMultiArrayOld<ElementType> contract(const rsMultiArrayOld<ElementType> &subject,
    rsUint32 index1, rsUint32 index2);

protected:

  /** \name Internal functions */

  void allocateAndInitIndexArray(const rsUint32 newNumIndices, const rsUint32 newIndexRanges[]);
  void allocateDataArray(bool initWithZeros);
  void freeIndexArray();
  void freeDataArray();


  /** \name Data */

  /** The number of indices required to access a single element. */
  rsUint32 numIndices;

  /** C-style array containing the range for each of the indices. For example, in a multiarray
  with 3 indices representing 5 rows, 7 columns and 3 slices, it would contain the 3 numbers
  5, 7, 3 and the object would be a 5x7x3 "hypermatrix". */
  rsUint32 *indexRanges;  // maybe use a std::vector

  /** The actual data as flat C-style array. */
  ElementType *data; // maybe use a std::vector

};


#endif
