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


  void fillRandomly(T min = T(0), T max = T(1), int seed = 0, bool roundToInt = false)
  {
    if(roundToInt)
      rsArrayTools::fillWithRandomIntegers(dataPointer, getSize(), int(min), int(max), seed);
    else
      rsArrayTools::fillWithRandomValues(dataPointer, getSize(), min, max, seed);
  }

  template<class T2>
  void setData(const std::vector<T2>& newData)
  {
    rsAssert((int)newData.size() == size, "Passed data vector must match size of our data array");
    rsArrayTools::convert(&newData[0], dataPointer, size);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Manipulation */

  /** Sets all matrix elements to zero. */
  void setToZero() { rsArrayTools::fillWithZeros(dataPointer, getSize()); }
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
  // maybe cache that values - the conversion from size_t to int may be expensive, if this is done
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
  const std::vector<int>& getShape() const { return shape; }
  // we may not store the shape as vector<int> in an optimized version, so i'm not sure, if i can 
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

  /** Returns true, iff all elements ar close to zero within a given tolerance, which is zero by 
  default. */
  bool isAllZeros(T tolerance = T(0)) const
  { return rsArrayTools::isAllZeros(dataPointer, getSize(), tolerance); }

  /** NOT YET READY FOR USE!!
  Returns true, if the corresponding elements of the arrays x and y are all close to each other
  up to some tolerance. Elements in one array that have no corresponding element in the other array
  must be (close to) zero. In this sense, the arrays are assumed to be zero padded as needed. */
  /*
  static bool areClosePadded(const rsMultiArrayView<T>& x, const rsMultiArrayView<T>& y, 
    T tol = T(0))
  {
    rsError("This seems to be still buggy - needs unit tests");
    rsWarning("rsMultiArrayView::areClosePadded not yet tested");

    int R = x.getNumIndices();
    if(y.getNumIndices() != R)  {
      rsError("x and y should have the same number of indices");
      return false;
    }
    int size = rsMax(x.getSize(), y.getSize());
    vector<int> ix(R), iy(R);                    // multiindices into x and y
    for(int i = 0; i < size; i++)
    {
      x.structuredIndices(i, &ix[0]);
      y.structuredIndices(i, &iy[0]);
      int jx = x.flatIndexSafe(&ix[0]);          // flat index into x
      int jy = y.flatIndexSafe(&iy[0]);          // flat index into y

      if(jx == -1 && jy == -1)                   // jx, jy both invalid: nothing to check
        continue;                                // ...can this actually happen?
      else if(jx == -1) {
        if(rsAbs(y.dataPointer[jy]) > tol)       // jx invalid: y[jy] must be zero
          return false; }
      else if(jy == -1) {
        if(rsAbs(x.dataPointer[jx]) > tol)       // jy invalid: x[jx] must be zero
          return false; }
      else {                                     // both jx and jy are valid... 
        T d = x.dataPointer[jx] - y.dataPointer[jy];
        if(rsAbs(d) > tol)  // ...so x[jx] must match y[jy]
          return false;  }
    }
    return true;
  }
  */
  // nope! that's not how it works! the whole idea with the index structuring and destructuring is
  // flawed!
  // needs more tests - compare a 2x3 with a 3x2 matrix

  /** NOT YET READY FOR USE!! */
  //bool isCloseToPadded(const rsMultiArrayView<T> y, T tolerance = T(0)) const
  //{ return areClosePadded(*this, y, tolerance); }




  //-----------------------------------------------------------------------------------------------
  /** \name Element Access */

  /** Read and write access to array elements. The syntax for accessing, for example, 3D array 
  elements is: A(i, j, k) = .... */
  template<typename... Rest>
  T& operator()(const int i, Rest... rest) { return dataPointer[flatIndex(0, i, rest...)]; }

  /** Read-only access to array elements. */
  template<typename... Rest>
  const T& operator()(const int i, Rest... rest) const 
  { return dataPointer[flatIndex(0, i, rest...)]; }


  /*
  template<typename... Rest>
  T getElementPadded(const int i, Rest... rest) const
  {
    int index = flatIndexSafe(0, i, rest...);
    if(index == -1)
      return T(0);
    return dataPointer[index];
  }
  */
  // needs test


  T getElementPadded3D(int i, int j, int k, T padding) const
  {
    rsAssert(shape.size() == 3, "Must be used only with 3D arrays");
    if(i < 0 || i >= shape[0] || j < 0 || j >= shape[1] || k < 0 || k >= shape[2]) 
      return T(0);
    return (*this)(i, j, k);
  }
  // todo: make 1D and 2D versions - these may be more efficient than using the recursive 
  // implementation for the general case (which is not yet ready)
  // -maybe make versions that leave out the check against < 0 - typically, this will be ensured 
  //  already by the caller (when a loop starts at 0), so the checks are superfluous
  //  maybe name it getElementPaddedRight3D or PaddedHigh or HighPadded or similar
  


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
    rsArrayTools::add(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
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
    rsArrayTools::subtract(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }


  static void weightedSum(const rsMultiArrayView<T>& A, T wA, const rsMultiArrayView<T>& B, T wB, 
    rsMultiArrayView<T>& C)
  {
    size_t numDims = C.shape.size();
    rsAssert(A.shape.size() == numDims && B.shape.size() == numDims);
    std::vector<int> indices(numDims);
    C.setToZero();
    for(int i = 0; i < C.size; i++)
    {
      C.structuredIndices(i, &indices[0]);
      int ia = A.flatIndexSafe(&indices[0]);
      int ib = B.flatIndexSafe(&indices[0]);
      if(ia != -1) C.dataPointer[i] += wA * A.dataPointer[ia];
      if(ib != -1) C.dataPointer[i] += wB * B.dataPointer[ib];
    }
  }
  // needs tests
  // ..can this be done more efficiently? the computation of the indices array is expensive



  /** Multiplies the two multiarrays element-wise. */
  static void multiply(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::multiply(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Divide the two multiarrays element-wise. */
  static void divide(
    const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArrayTools::divide(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }

  /** Scales all elements by a given factor. */
  void scale(T factor) { rsArrayTools::scale(dataPointer, getSize(), factor); }

  /** Negates all elements. */
  void negate() { rsArrayTools::negate(dataPointer, dataPointer, getSize()); }

  // maybe factor out common code with rsMatrixView into a class rsArrayView which serves as 
  // baseclass for both - a general "view" class for any sort of array, i.e. homogeneous data 
  // stored in a contiguous memory chunk. should have members dataPointer and size

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Computes the strides for a given shape. The arrays shape and strides muts both be numIndices 
  long. */
  static void computeStrides(int numIndices, int* shape, int* strides)
  {
    int i = numIndices-1;         // last index has stride 1 -> row-major matrix storage
    int s = 1;
    while(i >= 0) {
      strides[i] = s;
      s *= shape[i];
      --i;
    }
  }


  /** Resets all data fields to default vaules, indicating an empty object. */
  void reset()
  {
    shape.clear(); 
    strides.clear();
    dataPointer = nullptr;
    size = 0;
  }


protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Index Computation */
  // maybe move to public section

  /** Used for implementing the variadic template for the () operator. Takes a recursion depth (for
  recursive template instantiation) and a variable number of indices .... */
  template<typename... Rest>
  int flatIndex(const int depth, const int i, Rest... rest) const
  { 
    //int dbg = flatIndex(depth, i) + flatIndex(depth+1, rest...); // for debug
    return flatIndex(depth, i) + flatIndex(depth+1, rest...); 
  }

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
  // has this been tested?


  // Safe versions that check, if all the indices are within their proper range and returns -1, 
  // if any of them is not:
  /*
  template<typename... Rest>
  int flatIndexSafe(const int depth, const int i, Rest... rest) const
  {
    int i1 = flatIndex(depth,   i);
    int i2 = flatIndex(depth+1, rest...);
    if(i1 == -1 || i2 == -1)
      return -1;
    return i1 + i2;
  }
  // needs test

  int flatIndexSafe(const int depth, const int index) const
  {
    if(index < 0 || index >= shape[depth])
      return -1;
    return index * strides[depth]; 
  }
  // recursion base case
  */

  int flatIndexSafe(const int* indices) const
  {
    int fltIdx = 0;
    for(size_t i = 0; i < strides.size(); i++) {
      if(indices[i] < 0 || indices[i] >= shape[i])
        return -1;
      fltIdx += indices[i] * strides[i]; }
    return fltIdx;
  }




  /** Converts a flat index into an array of structured/hierarchical indices. */
  void structuredIndices(int flatIndex, int* indices) const
  {
    for(int i = 0; i < getNumIndices(); i++)
    {
      indices[i] = flatIndex / strides[i];
      flatIndex -= indices[i] * strides[i]; // remainder of previous division
      // maybe use divmod instead of div and mul
    }
    // there's a unit test in the research repo - maybe drag over into main repo
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Member Updating */

  void updateStrides()
  {
    int rank = (int) shape.size();
    strides.resize(rank);

    // call computeStrides instead of stuff below:

    int i = rank-1;         // last index has stride 1 -> row-major matrix storage
    int s = 1;
    while(i >= 0) {
      strides[i] = s;
      s *= shape[i];
      --i;
    }
  }
  // maybe move to cpp file, maybe have a static member function
  //   void computeStrides(int rank, int* shape, int* strides)
  // so we may defer the allocation of shape/strides arrays to client code - clent code may then 
  // call the stride-computation

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

  std::vector<int> shape;     // maybe rename to extents
  std::vector<int> strides;
  T* dataPointer = nullptr;   // rename to data or pData
  int size = 0;

  // todo: get rid of strides, let shape be a non-owned pointer to int, store size of the shapes 
  // array - we want to avoid memory allocations when creating such a view object - creating a view
  // should be cheap! ...actually, we would only need two pointers: data and strides - to support 
  // the () syntax for accessing elements - but then we couldn't check for out-of-range indexes - 
  // for that, we need also the shape - use: int numDims; int* shape; int* strides;


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

Note: This is still incomplete */

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

  // Copy/move / construction/assignment copy/paste/edited from rsMatrix - needs some unit tests

  /** Copy constructor. Copies data from B into this object.  */
  rsMultiArray(const rsMultiArray& B)
  {
    setShape(B.shape); // allocates the memory
    rsArrayTools::copy(B.dataPointer, this->dataPointer, this->getSize());
  }

  /** Move constructor. Takes over ownership of the data stored in B. */
  rsMultiArray(rsMultiArray&& B)
  {
    this->size    = B.size;
    this->data    = std::move(B.data);
    this->shape   = std::move(B.shape);
    this->strides = std::move(B.strides);
    rsAssert(B.data.size() == 0);
    rsAssert(B.shape.size() == 0);
    rsAssert(B.strides.size() == 0);
    updateDataPointer();
    B.reset();                         // invalidates pointer in B
  }

  /** Copy assignment operator. Copies data from rhs into this object. */
  rsMultiArray<T>& operator=(const rsMultiArray<T>& rhs)
  {
    if (this != &rhs) { // self-assignment check expected
      setShape(rhs.shape);
      rsArrayTools::copy(rhs.dataPointer, this->dataPointer, this->getSize());
    }
    return *this;
  }

  /** Move assignment operator. Takes over ownership of the data stored in rhs. */
  rsMultiArray<T>& operator=(rsMultiArray<T>&& rhs)
  {
    this->size    = rhs.size;
    this->data    = std::move(rhs.data);
    this->shape   = std::move(rhs.shape);
    this->strides = std::move(rhs.strides);
    rsAssert(rhs.data.size() == 0);
    rsAssert(rhs.shape.size() == 0);
    rsAssert(rhs.strides.size() == 0);
    updateDataPointer();
    rhs.reset();
    return *this;
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
  /** \name Inquiry */

  /** Returns true, iff this multiarry has the given shape. */
  bool hasShape(const std::vector<int>& shape)
  { return shape == this->shape; }


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

  bool operator==(const rsMultiArray<T>& B) const
  { return this->shape == B.shape && this->data == B.data; }

  bool operator!=(const rsMultiArray<T>& B) const
  { return !(*this == B); }


  /** Convolves the two arrays x and h which are assumed to be 3-dimensional (i.e. have 3 indices) 
  and stores the result in y. */
  static void convolve3D(const rsMultiArray<T>& x, const rsMultiArray<T>& h, rsMultiArray<T>& y);



protected:

  /** Updates the data-pointer inherited from rsMultiArrayView to point to the begin of our 
  std::vector that holds the actual data. */
  void updateDataPointer()
  {
    if(data.size() > 0)
      this->dataPointer = &data[0];
    else
      this->dataPointer = nullptr;
    // maybe assert that data.size == this->size
  }

  std::vector<T> data;

};

/** Multiplies a scalar and a multiarray. */
template<class T>
inline rsMultiArray<T> operator*(const T& s, const rsMultiArray<T>& A)
{ rsMultiArray<T> B(A); B.scale(s); return B; }

template<class T>
void rsMultiArray<T>::convolve3D(const rsMultiArray<T>& x, const rsMultiArray<T>& h, 
  rsMultiArray<T>& y)
{
  rsAssert(x.getNumIndices() == 3 && h.getNumIndices() == 3, "x and h must be 3-dimensional");
  int Lx = x.getExtent(0), Mx = x.getExtent(1), Nx = x.getExtent(2);
  int Lh = h.getExtent(0), Mh = h.getExtent(1), Nh = h.getExtent(2);
  int Ly = Lx + Lh - 1,    My = Mx + Mh - 1,    Ny = Nx + Nh - 1;
  y.setShape({Ly, My, Ny});
  for(int l = 0; l < Ly; l++) {
    for(int m = 0; m < My; m++) {
      for(int n = 0; n < Ny; n++) {
        T s = T(0);
        for(int i = rsMax(0, l-Lx+1); i <= rsMin(Lh-1, l); i++) {
          for(int j = rsMax(0, m-Mx+1); j <= rsMin(Mh-1, m); j++) {
            for(int k = rsMax(0, n-Nx+1); k <= rsMin(Nh-1, n); k++) {
              s += h(i, j, k) * x(l-i, m-j, n-k); }}}
        y(l, m, n) = s; }}}
}
// How can general nD convolution be implemented?


#endif