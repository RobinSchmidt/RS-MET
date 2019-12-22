#ifndef RAPT_MULTIARRAY_H
#define RAPT_MULTIARRAY_H

// this is not yet finished

//=================================================================================================

/** This is a class for treating raw C-arrays as multidimensional arrays. It does not store/own the
actual data. It just acts as wrapper around an existing array for conveniently accessing and 
manipulating array elements via multiple indexes using the () operator, which has been overloaded
for as many indices as needed. Example code could look like:

float data[24];                             // flat C-array of data
rsMultiArrayView<float> A({2,4,3}, data);   // we want to interpret the data as 2x4x3 3D array
A(0,0,0) = 111.f;                           // set element at position 0,0,0 (the first one)
A(1,3,2) = 243.f;                           // set element at position 1,3,2 (the last one)

Note that the constructor has to perform two heap allocations (for two std::vectors), so wrapping
an array like this is not a free operation. 

todo: maybe reduce that to one allocation by using a single array for both shape and strides 
- maybe an array of struct DimInfo which has fields size, stride - ...as optimization later 

...or we should not own the shape and strides array - client code must allocate both and give 
pointers to this class - yes this will also allow the view to be used when the array has another
memory layout - we don't decide the memory layout in the view


*/

template<class T>
class rsMultiArrayView
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Default constructor. */
  rsMultiArrayView() {}

  /** Creates a multi-array view with the given shape for the given raw array of values in "data". 
  The view will *not* take ownership over the data. */
  rsMultiArrayView(const std::vector<int>& initialShape, T* data) : shape(initialShape)
  {
    updateStrides();
    updateSize();
    dataPointer = data;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up a new shape. The product of the extents (i.e. the total size) may be different from
  our former shape to facilitate using the same C-array with various sizes and shapes. */
  void setShape(const std::vector<int>& newShape) 
  { shape = newShape; updateStrides(); updateSize(); }
  // maybe rename to setExtents...but maybe not - shape seems to be common for that



  //-----------------------------------------------------------------------------------------------
  /** \name Manipulation */

  /** Sets all matrix elements to zero. */
  void setToZero() { rsArray::fillWithZeros(dataPointer, getSize()); }
    // same as in rsMatrixView...i think both classes really should have a common baseclass 
    // rsArrayView where we consolidate all these common functions - it needs to have the 
    // dataPointer and the size...which seems totally appropriate for an ArrayView class


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns true, iff the two arrays A and B have the same shape. */
  static bool areSameShape(const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B)
  { return A.shape == B.shape; }

  /** Returns the number of array elements. */
  int getSize() const { return size; }

  /** Returns the number of array dimensions, i.e. the number of indices. */
  int getNumIndices() const { return (int) shape.size(); }
  // maybe rename to getNumIndices - "dimension" is a bit ambiguous here - for example a vector in
  // 3D space is considered to be 3-dimensional, but has only one index
  // maybe cache that values - the conversion form size_t to int may be expensive, if this is done
  // often -> benchmark!

  /** Returns the one plus the maximum value that given index may assume */
  int getExtent(int index) const 
  { 
    rsAssert(index < getNumIndices());
    return shape[index]; 
  }
  // alternative names: getIndexRange, getIndexExtent, getIndexLimit ..."limit" is a bit bad because
  // it suggests that the limiting value itself is included, "range" is bad because it suggests to
  // consist of two values (lower and upper bound) where only the upper bound is meant...how about
  // getLengthAlong(int dimension)

  /** Returns a const reference to the shape array of the multidimensional array. The shape is 
  defined by the number of  values that each index may run over. For example, a 2x3 matrix has a 
  shape array of [2,3]. */
  //const std::vector<int>& getShape() const { return shape; }
  // we may not store the shape as vector<int> in an optimzed version, so i'm not sure, if i can 
  // provide that interface - maybe instead provide a function getIndexRange(int whichIndex) 
  // or getExtent(int index)

  // bool isSymmetricIn(const int i const int j);
  // returns true, iff array is symmetric with respect to the given pair of indices
  // maybe rename to isSymmetricRegarding(i, j)

  //

  // maybe have functions that return a pointer to a particular "row" - but the function should be 
  // a template, taking an arbitrary number of indices - for example A.getRowPointer(2, 3) would
  // return the 4th row in the 3rd slice ...or the other way around...i think so


  /** Returns a const pointer to the data for read access as a flat array. */
  const T* getDataPointerConst() const { return dataPointer; }

  /** Returns a pointer to the data for read and write access as a flat array. */
  T* getDataPointer() { return dataPointer; }


  //-----------------------------------------------------------------------------------------------
  /** \name Element Access */

  /** Read and write access to array elements. The syntax for accessing, for example, 3D array 
  elements is: A(i, j, k) = .... */
  template<typename... Rest>
  T& operator()(const int i, Rest... rest) { return dataPointer[flatIndex(0, i, rest...)]; }
  // maybe use const int as was doen in rsMatrix

  /** Read-only access to array elements. */
  template<typename... Rest>
  const T& operator()(const int i, Rest... rest) const 
  { return dataPointer[flatIndex(0, i, rest...)]; }


  /*
  const T& operator()(const int i, const int j) const
  {
    return dataPointer[flatIndex(i, j)];
  }
  */

  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic */

  /** Adds elements of A to corresponding elements in B and stores results in C. */
  static void add(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArray::add(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }
  // hmm...these functions add, subtrac, etc. are copy/pasted more or less exactly from rsMatrix - 
  // maybe we can factor out a common baseclass? maybe rsArrayView? it should contain the 
  // dataPointer and the size - then rsMatrix would also use a cached size - set up unit tests and 
  // then try it - maybe it could also handle application of functions in a uniform way
  // but wait - no - the implementation of areSameShape is different - on the other hand, that
  // function is only used in the assertion

  /** Subtracts elements of B from corresponding elements A in and stores results in C. */
  static void subtract(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArray::subtract(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Multiplies the two multiarrays element-wise. */
  static void multiply(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArray::multiply(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Divide the two multiarrays element-wise. */
  static void divide(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArray::divide(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }





protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Index Computation */
  // maybe move to public section

  /** Used for implementing the variadic template for the () operator. Takes a recursion depth (for
  recursive template instantiation) and a variable number of indices .... */
  template<typename... Rest>
  int flatIndex(const int depth, const int i, Rest... rest) const
  { return flatIndex(depth, i) + flatIndex(depth+1, rest...); }

  /** Base case for the variadic template. this version will be instatiated when, in addition to 
  the recursion depth, only one index is passed. */
  int flatIndex(const int depth, const int index) const
  { 
    rsAssert(index >= 0 && index < shape[depth], "invalid index"); // verify this!
    return index * strides[depth]; 
  }

  /** Converts a C-array (assumed to be of length getNumDimensions()) of indices to a flat 
  index. */
  int flatIndex(const int* indices) const
  {
    int fltIdx = 0;
    for(size_t i = 0; i < strides.size(); i++)
      fltIdx += indices[i] * strides[i];
    return fltIdx;
  }
  // needs test


  //-----------------------------------------------------------------------------------------------
  /** \name Member Updating */

  void updateStrides()
  {
    int rank = (int) shape.size();
    strides.resize(rank);
    int i = rank-1;         // last index has stride 1 -> row-major matrix storage
    int s = 1;
    while(i >= 0) {
      strides[i] = s;
      s *= shape[i];
      --i;
    }
  }
  // maybe move to cpp file

  /** Updates our size variable according to the values in the shape array. The total size is 
  (redundantly) cached in a member variable because it's used frequently. */
  void updateSize()
  {
    if(shape.size() == 0) { size = 0; return; }   // edge case
    size = 1;
    for(size_t i = 0; i < shape.size(); i++)
      size *= shape[i];
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Data */

  std::vector<int> shape;    // maybe rename to extents
  std::vector<int> strides;
  T* dataPointer = nullptr;
  int size = 0;

  // todo: get rid of strides, let shape be a non-owned pointer to int, store size of the shapes 
  // array - we want to avoid memory allocations when creating such a view object - creating a view
  // should be cheap! ...actually, we would only need two pointers: data and strides - to support 
  // the () syntax for accessing elements - but then we couldn't check for out-of-range indexes - 
  // for that, we need also the shape
  // or maybe make the number of indices a compile-time parameter, own strides and shape array (as
  // fixed arrays)...implement a constructor that takes an initializer list and copy its content 
  // into our members...but that implies that for each dimensionality, a template will be 
  // instantiated - so that makes 3 instantiations for 1D,2D,3D - whereas otherwise there would 
  // just be one instantiation...hmmm

};

//=================================================================================================

/** Implements a multidimensional array. Elements of an array A can be conveniently accessed via 
the natural syntax: 1D: A(i), 2D: A(i,j), 3D: A(i,j,k), etc. The data is stored in a std::vector. 
The implementation follows the same pattern as rsMatrix (which is in the Math folder). 

Note: This is still incomplete - all the required copy/move-constructors and -assignement operators
still need to be implemented. The class works already, but it's not yet return-value-optimized */

template<class T>
class rsMultiArray : public rsMultiArrayView<T>
{

public:


  rsMultiArray() {}


  rsMultiArray(const std::vector<int>& initialShape) : rsMultiArrayView<T>(initialShape, nullptr)
  {
    data.resize(this->size);
    updateDataPointer();
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up a new shape. */
  void setShape(const std::vector<int>& newShape) 
  { 
    rsMultiArrayView<T>::setShape(newShape); // updates extents, strides, size
    data.resize(this->size);
    updateDataPointer();
  }



  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two multiarrays element-wise. */
  rsMultiArray<T> operator+(const rsMultiArray<T>& B) const
  { rsMultiArray<T> C(this->shape); this->add(*this, B, &C); return C; }

  /** Subtracts two multiarrays element-wise */
  rsMultiArray<T> operator-(const rsMultiArray<T>& B) const
  { rsMultiArray<T> C(this->shape); this->subtract(*this, B, &C); return C; }

  /** Multiplies two multiarrays element-wise. */
  rsMultiArray<T> operator*(const rsMultiArray<T>& B) const
  { rsMultiArray<T> C(this->shape); this->multiply(*this, B, &C); return C; }

  /** Divides two multiarrays element-wise. */
  rsMultiArray<T> operator/(const rsMultiArray<T>& B) const
  { rsMultiArray<T> C(this->shape); this->divide(*this, B, &C); return C; }


  // todo: ==,!=


protected:

  /** Updates the data-pointer inherited from rsMultiArrayView to point to the begin of our 
  std::vector that holds the actual data. */
  void updateDataPointer()
  {
    if(data.size() > 0)
      this->dataPointer = &data[0];
    else
      this->dataPointer = nullptr;
  }

  std::vector<T> data;

};

#endif