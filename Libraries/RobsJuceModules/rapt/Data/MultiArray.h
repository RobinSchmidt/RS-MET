#ifndef RAPT_MULTIARRAY_H
#define RAPT_MULTIARRAY_H

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
- maybe an array of struct DimInfo which has fields size, stride - ...as optimization later */

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

  // setShape, fill, 


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
  // alternative names: getIndexRange, getIndexExtent, getIndexLimit

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


  //-----------------------------------------------------------------------------------------------
  /** \name Element Access */

  /** Read and write access to array elements. The syntax for accessing, for example, 3D array 
  elements is: A(i, j, k) = .... */
  template<typename... Rest>
  T& operator()(int i, Rest... rest) { return dataPointer[flatIndex(0, i, rest...)]; }
  // maybe use const int as was doen in rsMatrix

  //-----------------------------------------------------------------------------------------------
  /** \name Arithmetic */

  /** Adds elements of A to corresponding elements in B and stores results in C. */
  static void add(const rsMultiArrayView<T>& A, const rsMultiArrayView<T>& B, 
    rsMultiArrayView<T>* C)
  {
    rsAssert(areSameShape(A, B) && areSameShape(A, *C), "arguments incompatible");
    rsArray::add(A.dataPointer, B.dataPointer, C->dataPointer, A.getSize());
  }



protected:

  //-----------------------------------------------------------------------------------------------
  /** \name Index Computation */
  // maybe move to public section

  /** Used for implementing the variadic template for the () operator. Takes a recursion depth (for
  recursive template instantiation) and a variable number of indices .... */
  template<typename... Rest>
  int flatIndex(int depth, int i, Rest... rest) 
  { return flatIndex(depth, i) + flatIndex(depth+1, rest...); }

  /** Base case for the variadic template. this version will be instatiated when, in addition to 
  the recursion depth, only one index is passed. */
  int flatIndex(int depth, int index) 
  { 
    rsAssert(index >= 0 && index < shape[depth], "invalid index"); // verify this!
    return index * strides[depth]; 
  }

  /** Converts a C-array (assumed to be of length getNumDimensions()) of indices to a flat 
  index. */
  int flatIndex(const int* indices)
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

  std::vector<int> shape;
  std::vector<int> strides;
  T* dataPointer = nullptr;
  int size = 0;

};

//=================================================================================================

/** Implements a multidimensional array. Elements of an array A can be conveniently accessed via 
the natural syntax: 1D: A(i), 2D: A(i,j), 3D: A(i,j,k), etc. The data is stored in a std::vector. 
The implementation follows the same pattern as rsMatrix (which is in the Math folder). */

template<class T>
class rsMultiArray : public rsMultiArrayView<T>
{

public:


  rsMultiArray(const std::vector<int>& initialShape) : rsMultiArrayView<T>(initialShape, nullptr)
  {
    data.resize(this->size);
    updateDataPointer();
  }




  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two arrays: C = A + B. */
  rsMultiArray<T> operator+(const rsMultiArray<T>& B) const
  { rsMultiArray<T> C(this->shape); this->add(*this, B, &C); return C; }

  // todo: -,*,/,==,!=




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



/*
// move to unit tests:
void testMultiArray()
{
  typedef std::vector<int> VecI;
  typedef std::vector<float> VecF;
  typedef rsMultiArray<float> MA;

  // 3D vector:
  MA a1 = MA(VecI{3});     
  a1(0) = 1;
  a1(1) = 2;
  a1(2) = 3;

  //a1 = a1 + a1;


  // 3x2 matrix:
  MA a2 = MA(VecI{3,2});  
  a2(0,0) = 11;
  a2(0,1) = 12;

  a2(1,0) = 21;
  a2(1,1) = 22;

  a2(2,0) = 31;
  a2(2,1) = 32;


  // 2x4x3 block/cuboid:
  MA a3 = MA(VecI{2,4,3}); 
  a3(0,0,0) = 111;
  a3(0,0,1) = 112;
  a3(0,0,2) = 113;

  a3(0,1,0) = 121;
  a3(0,1,1) = 122;
  a3(0,1,2) = 123;

  a3(0,2,0) = 131;
  a3(0,2,1) = 132;
  a3(0,2,2) = 133;

  a3(0,3,0) = 141;
  a3(0,3,1) = 142;
  a3(0,3,2) = 143;


  a3(1,0,0) = 211;
  a3(1,0,1) = 212;
  a3(1,0,2) = 213;

  a3(1,1,0) = 221;
  a3(1,1,1) = 222;
  a3(1,1,2) = 223;

  a3(1,2,0) = 231;
  a3(1,2,1) = 232;
  a3(1,2,2) = 233;

  a3(1,3,0) = 241;
  a3(1,3,1) = 242;
  a3(1,3,2) = 243;

  // move code over to RAPT and turn this into a unit test
  // ->retrieve the flat data pointer and check, if the data is inserted in the desired layout



  // allow the user to specify an allocator so we can unit-test the memory allocation avoidance


  int dummy = 0;
}
*/

#endif