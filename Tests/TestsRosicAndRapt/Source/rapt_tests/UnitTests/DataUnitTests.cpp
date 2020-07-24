//#include "DataUnitTests.h"
using namespace RAPT;

//#include "../../../../../Libraries/JUCE/modules/rapt/Data/Simd/Float64x2.h"
// needed when it's commented out in rapt -> reduce build time during tweaking the class

//double sum(double* a, size_t N)
//{
//  double accu = 0;
//  for(size_t i = 0; i < N; i++)
//    accu += a[i];
//  return accu;
//}

template <class T>
void movingAverage5pt(const T* x, int N, T* y)
{
  rsAssert(N >= 4);  // todo: handle other cases separately

  //rsAssert(N >= 0);
  //if(N == 0) return;
  //if(N == 1) { y[0] = x[0]; return; }

  T t1 = x[0];
  T t2 = x[1];
  T t3 = x[2];
  T t4 = x[3];
  y[0] = T(1/3.) * (t1 + t2 + t3);
  y[1] = T(1/4.) * (t1 + t2 + t3 + t4);
  for(int n = 2; n < N-2; n++) {
    y[n] = T(1/5.) * (t1 + t2 + t3 + t4 + y[n+2]);
    t1 = t2;
    t2 = t3;
    t3 = t4;
    t4 = y[n+2]; }
  y[N-2] = T(1/4.) * (t1 + t2 + t3 + t4);
  y[N-1] = T(1/3.) * (     t2 + t3 + t4);
}
// this is still under construction
// maybe make a fixed-ends version of that using: y[0] = x[0], y[1] = (x[0]+x[1]+x[2])/3 - so
// the first output equals the first input and the second output uses a symmetric 3-point average
// ...i think, that scheme generalizes more nicely to M-point smoothers without introducing 
// asymmetric averaging at the endpoints - also, fixed endpoints may themselves be a desirable
// feature - for example when smoothing parameter trajectories

// a function movingAverage that does not run in place and supports arbitrary odd order - useful
// for producing target output for unit tests
template <class T>
void movingAverage(const T* x, int N, T* y, int order)
{
  rsAssert(N >= 0);
  int M = order/2;
  int n, m;
  //T s = T(0);  // accumulator

  // initial section:

  // middle section:
  for(int n = M; n < N-M; n++)
  {
    y[n] = T(0);
    for(m = -M; m <= M; m++)
      y[n] += x[n+m];
    y[n] *= T(1.0/order);
  }

  // last section:

}
// under construction...


// this code here: 
// https://stackoverflow.com/questions/21128981/finding-gcd-of-array-code-c-language
// uses a signature like: int gcd_a(int n, int a[n])
// ...that's actually pretty nice - perhaps i should use this for array functions as well. does the n
// have to come before the a? ...hmm - well - it doesn't actually compile :-(
/*
template <class T>
T sum(int N, T a[N])
{
  T s = T(0);
  for(int n = 0; n < N; n++)
    s += a[n];
  return s;
}
*/


bool testArrayFiltering()
{
  bool r = true;      // test result
  typedef std::vector<double> Vec;
  typedef RAPT::rsArrayTools AR;
  Vec x,y;

  x = Vec({1,3,2,-2,3,5,1});
  y.resize(x.size());
  AR::movingAverage3pt(&x[0], (int)x.size(), &y[0], false);  // out-of-place
  r &= y == Vec({2,2,1,1,2,3,3});
  AR::movingAverage3pt(&x[0], (int)x.size(), &x[0], false);  // in-place
  r &= x == Vec({2,2,1,1,2,3,3});
  // todo: test edge cases (arrays of length 0,1,2), test with endsFixed condition true

  // test moving median:
  x = Vec({1,3,2,-2,3,5,1});
  AR::movingMedian3pt(&x[0], (int)x.size(), &y[0]);
  r &= y == Vec({ 2,2,2,2,3,3,3 });
  AR::movingMedian3pt(&x[0], (int)x.size(), &x[0]);
  r &= x == Vec({ 2,2,2,2,3,3,3 });


  x = Vec({60,120,-60,180,120,-120,60,240,120});
  y = x;
  movingAverage5pt(&y[0], (int)y.size(), &y[0]);
  r &= y == Vec({40,75,84,48,36,96,84,75,140});

  return r;
}

bool arrayUnitTest()
{
  bool r = true;      // test result

  typedef std::vector<int> Vec;
  Vec v({ 1,2,3 }), w({4,5});
  Vec u = v;
  rsAppend(u, w);
  r &= u == Vec({1,2,3,4,5});
  u = v;
  rsAppend(u, u);
  r &= u == Vec({1,2,3,1,2,3});


  Vec a = { 0,1,2,3,4,5,6,7,8,9 };
  rsRemoveRange(a, 4, 7);
  r &= a == Vec({ 0,1,2,3, 8,9 });


  r &= testArrayFiltering();

  // int s = sum(3, &u[0]); // sum function doesn't compile


  return r;
}


template<class T>
class rsBinaryHeapTest : public rsBinaryHeap<T>
{

public:

  using rsBinaryHeap::rsBinaryHeap;


  bool isMinHeap(int i = 0) const
  {
    if(i >= size)
      return true;
    bool result = true;
    int l = left(i);
    int r = right(i);
    if(l < size) result &= data[i] <= data[l] && isMinHeap(l);
    if(r < size) result &= data[i] <= data[r] && isMinHeap(r);
    return result;
  }

  bool isMaxHeap(int i = 0) const
  {
    if(i >= size)
      return true;
    bool result = true;
    int l = left(i);
    int r = right(i);
    if(l < size) result &= data[i] >= data[l] && isMaxHeap(l);
    if(r < size) result &= data[i] >= data[r] && isMaxHeap(r);
    return result;
  }


  // move these to baseclass when finished:


  /** Removes the element at given index i from the heap and re-orders the remaining elements to
  maintain the heap-property. */
  void remove(int i)
  {
    // remove it and replace the slot with either the left or right child, then replace the child 
    // that was promoted up with either it left or right child and so on, i.e. remove and promote
    // children
  }

  void sort()
  {

  }
  // todo: implement heap-sort in this class an test it with various random arrays

};

bool binaryHeapUnitTest()
{
  // maybe rename to binaryTreeUnitTest and integrate tests for rsBinarySearchTree

  bool r = true; 

  using Vec = std::vector<int>;

  Vec A = {2,8,14,16,4,1,7,9,10,3};
  int N = (int) A.size();


 
  rsBinaryHeapTest<int> H;
  r &= H.getSize() == 0;
  H.setData(&A[0], N, N);
  H.buildHeap();
  r &= H.getSize() == 10;
  r &= H.isMaxHeap();

  // test replacing:
  H.replace(7, 15);
  r &= H.isMaxHeap();
  int numTests = 100;
  rsNoiseGenerator<double> ng;
  for(int i = 1; i <= numTests; i++)
  {
    int newIndex = ng.getSampleRaw() % H.getSize();
    int newValue = ng.getSampleRaw() % 100;
    int k = H.replace(newIndex, newValue);
    r &= H.isMaxHeap();
  }
  
  // test inserting:
  A.resize(N + numTests);            // make space - re-allocates and fill up with zeros
  H.setData(&A[0], N, N + numTests); // size is N, capacity is N + numTests
  r &= H.getSize() == N;
  r &= H.isMaxHeap();
  for(int i = 1; i <= numTests; i++)
  {
    int newValue = ng.getSampleRaw() % 100;
    int k = H.insert(newValue);              // index where the value ended up
    r &= H.getSize() == N + i;
    r &= H.isMaxHeap();
  }

  // test removing:





  //N = (int) A.size();


  
  





  // test the double heap:

  rsDoubleHeap<int> D;
  Vec B;
  int i;
  A = Vec({5,2,3});
  B = Vec({6,7,8});
  D.setData(A, B);
  i = D.replace(1, 6);  // replace 2 by 6, float up to front of small, no exchange
  r &= i == 0 && A == Vec({6,5,3}) && B ==  Vec({6,7,8});
  i = D.replace(5, 4);  // replace 8 by 4, float to front of large, exchange
  r &= i == 1 && A == Vec({5,4,3}) && B ==  Vec({6,7,6});
  // do more tests, using larger heaps, maybe check property in a loop






  int dummy = 0;


  // test the binary search tree:

  // maybe test with the small array 1,2,3 in all permutations, test also 1,2,2 and 1,1,2 in all
  // possible permuations ...and 1,1,1
  rsBinarySearchTree<int> T;




  A = Vec({50,20,80,10,30,60,100,5,15,25,40,55,70});
  T.setData(A);
  r &= T.isSearchTree();
  T.replace(5, 35);   // the 60 becomes a 35
  r &= T.isSearchTree();  // 50,20,80,10,30,55,100,5,15,40,35,70
  T.replace(5, 65);
  r &= T.isSearchTree(); 
  A = Vec({50,20,80,10,30,60,100,5,15,25,40,55,70});
  T.setData(A);
  T.replace(1, 70);
  r &= T.isSearchTree(); 
  T.replace(4, 20);
  r &= T.isSearchTree(); 
  A = Vec({50,20,80,10,60,30,90});
  T.setData(A);
  r &= T.isSearchTree();

  A = Vec({50,20,80,10,30,60,100,5,15,25,40,55,70});
  T.setData(A);
  r &= T.isSearchTree();
  for(int i = 1; i <= numTests; i++)
  {
    int newIndex = ng.getSampleRaw() % H.getSize();
    int newValue = ng.getSampleRaw() % 100;
    int k = T.replace(newIndex, newValue);
    r &= T.isSearchTree();
    // does not work yet - i think the first thing that needs to be done after replacement is to
    // compare and possibly swap with the sibling? but no, if we replace a right node, and its 
    // now less than its sibling, it will be also less than its parent...ah - but nevertheless,
    // we may have to swap with the sibling..or maybe binary search trees are not so similar to 
    // heaps after all?
  }




  A = Vec({1,2,3}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });
  A = Vec({2,1,3}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });
  A = Vec({3,1,2}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });
  A = Vec({1,3,2}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });  
  A = Vec({2,3,1}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });
  A = Vec({3,2,1}); T.setData(A); T.buildSearchTree(); r &= A == Vec({ 2,1,3 });
  // The 3 cases where the left child of the root node is greater than the right child (1,3,2 (3>2), 
  // 2,3,1 (3>1), 3,2,1 (2>1))- previously failed when we were using a function similar to 
  // builHeap. When the left child of a node is greater than the right child, floating down the 
  // parent cannot fix the order properly - instead, we would have to swap the children (possibly 
  // in addition to swapping with the parent - maybe the general buildTree is misguided and we 
  // need to do something different to build the tree. However, to maintain the tree property when
  // we replace an element, we may assume that the left and right child are not in violation of the
  // property - the only node that may be in violation is the replaced node - so the data structure 
  // may still be useful for the moving quantile filter. However we cannot reuse the buildHeap 
  // function as general buildTree function - that simply doesn't work - we need a different 
  // implementation for the search trees (i think). Perhaps looping through all elements calling 
  // floatIntoPlace





  // ...oookay - that seems to work - at least in this very small case - now we need to test some
  // bigger cases...
  // hmm i think, etsblishing the property from an array of random values may be not so easy - but 
  // maintaining it, when a node is replaced should work exactly as in the heap


  /*
  int maxN = 100;
  A.resize(maxN);
  for(int i = 1; i <= numTests; i++)
  {
    N = ng.getSampleRaw() % maxN;
    A.resize(N);
    for(int n = 0; n < N; n++)
      A[n] = ng.getSampleRaw() % 100;
    T.setData(A);
    T.buildSearchTree();
    bool isTree = T.isSearchTree();
    r &= T.isSearchTree();
  }
  */




  /*
  T.setData(A);
  r &= T.getSize() == 3;
  r &= A == Vec({ 2,1,3 });      // desired array order is 2,1,3
  r &= T.isSearchTree();
  */








  /*
  rsBinarySearchTree<int> T;
  T.setData(&A[0], N, N);
  r &= T.getSize() == 10;
  r &= T.isSearchTree();
  */


  return r;
}



bool ringBufferUnitTest()
{
  bool r = true;

  rsRingBuffer<double> b(8);
  b.setLength(5);

  using Vec = std::vector<double>;
  Vec B(5);

  double y;

  // test getSample and copying the buffer content:
  y = b.getSample( 1); r &= y ==  0; b.copyTo(&B[0]); r &= B == Vec({  1,  0,  0,  0,  0 });
  y = b.getSample( 2); r &= y ==  0; b.copyTo(&B[0]); r &= B == Vec({  2,  1,  0,  0,  0 });
  y = b.getSample( 3); r &= y ==  0; b.copyTo(&B[0]); r &= B == Vec({  3,  2,  1,  0,  0 });
  y = b.getSample( 4); r &= y ==  0; b.copyTo(&B[0]); r &= B == Vec({  4,  3,  2,  1,  0 });
  y = b.getSample( 5); r &= y ==  0; b.copyTo(&B[0]); r &= B == Vec({  5,  4,  3,  2,  1 });
  y = b.getSample( 6); r &= y ==  1; b.copyTo(&B[0]); r &= B == Vec({  6,  5,  4,  3,  2 });
  y = b.getSample( 7); r &= y ==  2; b.copyTo(&B[0]); r &= B == Vec({  7,  6,  5,  4,  3 });
  y = b.getSample( 8); r &= y ==  3; b.copyTo(&B[0]); r &= B == Vec({  8,  7,  6,  5,  4 });
  y = b.getSample( 9); r &= y ==  4; b.copyTo(&B[0]); r &= B == Vec({  9,  8,  7,  6,  5 });
  y = b.getSample(10); r &= y ==  5; b.copyTo(&B[0]); r &= B == Vec({ 10,  9,  8,  7,  6 });
  y = b.getSample(11); r &= y ==  6; b.copyTo(&B[0]); r &= B == Vec({ 11, 10,  9,  8,  7 });
  y = b.getSample(12); r &= y ==  7; b.copyTo(&B[0]); r &= B == Vec({ 12, 11, 10,  9,  8 });
  y = b.getSample(13); r &= y ==  8; b.copyTo(&B[0]); r &= B == Vec({ 13, 12, 11, 10,  9 });
  y = b.getSample(14); r &= y ==  9; b.copyTo(&B[0]); r &= B == Vec({ 14, 13, 12, 11, 10 });
  y = b.getSample(15); r &= y == 10; b.copyTo(&B[0]); r &= B == Vec({ 15, 14, 13, 12, 11 });
  y = b.getSample(16); r &= y == 11; b.copyTo(&B[0]); r &= B == Vec({ 16, 15, 14, 13, 12 });


  b.reset();
  size_t i = 0, j;
  j = b.getIndexFromNewest(i);
  j = b.getIndexFromOldest(i);

  return r;
}

bool doubleEndedQueueUnitTest()
{
  bool r = true; 

  //rsDoubleEndedQueue<int> q(8);  // []

  rsDoubleEndedQueue<int> q(5);  // []

  r &= q.getMaxLength() == 6;     // a power of 2 minus 2, >= requested size
  //r &= q.getMaxLength() == 7;     // a power of 2 minus 1, >= requested size

  r &= q.isEmpty();
  r &= q.getLength() == 0;

  q.pushFront(5);                // [5]
  r &= q.getLength() == 1;
  r &= q.readHead()  == 5;
  r &= q.readTail()  == 5;

  q.pushFront(6);                // [5 6]
  r &= q.getLength() == 2;
  r &= q.readHead()  == 6;
  r &= q.readTail()  == 5;

  q.pushFront(7);                // [5 6 7]
  r &= q.getLength() == 3;
  r &= q.readHead()  == 7;
  r &= q.readTail()  == 5;

  q.pushBack(4);                 // [4 5 6 7]
  r &= q.getLength() == 4;
  r &= q.readHead()  == 7;
  r &= q.readTail()  == 4;

  q.pushBack(3);                 // [3 4 5 6 7]
  r &= q.getLength() == 5;
  r &= q.readHead()  == 7;
  r &= q.readTail()  == 3;

  q.popFront();                  // [3 4 5 6]
  r &= q.getLength() == 4;
  r &= q.readHead()  == 6;
  r &= q.readTail()  == 3;

  q.popBack();                   // [4 5 6]
  r &= q.getLength() == 3;
  r &= q.readHead()  == 6;
  r &= q.readTail()  == 4;

  q.clear();                     // []
  r &= q.isEmpty();
  r &= q.getLength() == 0;

  q.pushBack(1);                 // [1]
  r &= q.getLength() == 1;
  r &= q.readHead()  == 1;
  r &= q.readTail()  == 1;

  int i;
  for(i = 2; i <= 6; i++)        // [1 2 3 4 5 6]
    q.pushFront(i);

  r &= q.getLength() == 6;
  r &= q.readHead()  == 6;
  r &= q.readTail()  == 1;

  r &= q.isFull();


  //// this triggers an assertion but actually still works - however, the MovingMax filter does not
  //// work properly anymore when we fill the deque up to bufferSize-1
  //q.pushFront(7);                // [1 2 3 4 5 6 7]
  //r &= q.getLength() == 7;
  //r &= q.readHead()  == 7;
  //r &= q.readTail()  == 1;




  //// this does not work anymore due to overflow:
  //q.pushFront(8);                // [1 2 3 4 5 6 7 8]
  //r &= q.getLength() == 8;
  //r &= q.readHead()  == 8;
  //r &= q.readTail()  == 1;
  //// the best thing would be, if it would automatically resize

  return r;
}

bool float64x2UnitTest()
{
  bool r = true;      // test result

  // test constructors and getters:
  //rsFloat64x2 x00;              // default constructor, two zeros
  //r &= x00.get0() == 0.0;       // ...or actually no, we leave them uninitialized
  //r &= x00.get1() == 0.0;

  // construct from a double:
  rsFloat64x2 x11(1.0);  r &= x11.get0() == 1.0; r &= x11.get1() == 1.0;

  // construct from two doubles:
  rsFloat64x2 x12(1.0, 2.0); r &= x12.get0() == 1.0; r &= x12.get1() == 2.0;

  // construct from array of doubles:
  double arr[2] = { 3.0, 4.0 };
  rsFloat64x2 x34(arr); r &= x34.get0() == 3.0; r &= x34.get1() == 4.0;

  // construct from another instance:
  rsFloat64x2 y(x34); r &= y.get0() == 3.0; r &= y.get1() == 4.0;

  // setters:
  y.set0(5.0);     r &= y.get0() == 5.0; r &= y.get1() == 4.0;
  y.set1(6.0);     r &= y.get0() == 5.0; r &= y.get1() == 6.0;
  y.set(1.0, 2.0); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y.set(3.0);      r &= y.get0() == 3.0; r &= y.get1() == 3.0;

  // array access operator:
  y[0] = 1.0; y[1] = 2.0;  r &= y[0] == 1.0; r &= y[1] == 2.0;

  // getter:
  double v0, v1;
  y.get(v0, v1); r &= v0 == 1.0; r &= v1 == 2.0;

  // assignment and equality:
  y = x12; r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  //r &= y == x12;

  // binary arithmetic operators:
  y = x12 + x34; r &= y.get0() == 4.0; r &= y.get1() == 6.0;
  y = x34 - x12; r &= y.get0() == 2.0; r &= y.get1() == 2.0;
  y = x12 * x34; r &= y.get0() == 3.0; r &= y.get1() == 8.0;
  y = x34 / x12; r &= y.get0() == 3.0; r &= y.get1() == 2.0;

  // binary arithmetic operators with scalar lhs:
  y  = x12;
  y  = 1.0  + y; r &= y.get0() == 2.0; r &= y.get1() == 3.0;
  y  = 5.0  - y; r &= y.get0() == 3.0; r &= y.get1() == 2.0;
  y  = 2.0  * y; r &= y.get0() == 6.0; r &= y.get1() == 4.0;
  y  = 12.0 / y; r &= y.get0() == 2.0; r &= y.get1() == 3.0;

  // binary arithmetic operators with scalar rhs:
  y  = x34;
  y  = y + 2.0; r &= y.get0() == 5.0; r &= y.get1() == 6.0;
  y  = y - 2.0; r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y  = y * 2.0; r &= y.get0() == 6.0; r &= y.get1() == 8.0;
  y  = y / 2.0; r &= y.get0() == 3.0; r &= y.get1() == 4.0;

  // unary arithmetic operators -,+=,... :
  y  = x12;
  y += x34; r &= y.get0() ==  4.0; r &= y.get1() ==   6.0;
  y -= x12; r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;
  y *= x34; r &= y.get0() ==  9.0; r &= y.get1() ==  16.0;
  y /= x34; r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;
  y  = -y;  r &= y.get0() == -3.0; r &= y.get1() ==  -4.0;
  y  = -y;  r &= y.get0() ==  3.0; r &= y.get1() ==   4.0;

  // functions: sqrt, min, max, clip, abs, sign:
  y = x34 * x34;
  y = rsSqrt(y);       r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y = rsMin(x12, x34); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y = rsMin(x34, x12); r &= y.get0() == 1.0; r &= y.get1() == 2.0;
  y = rsMax(x12, x34); r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y = rsMax(x34, x12); r &= y.get0() == 3.0; r &= y.get1() == 4.0;
  y.set(1.0, 9.0);
  y = rsClip(y, 2.0, 7.0); r &= y.get0() == 2.0; r &= y.get1() == 7.0;
  y.set(-2.0, -3.0);
  y = rsAbs(y); r &= y.get0() == 2.0; r &= y.get1() == 3.0;
  y.set(-2.0, 3.0);
  y = rsSign(y); r &= y.get0() == -1.0; r &= y.get1() == +1.0;

  // reductions to scalar:
  double s;
  y.set(3.0, 5.0);
  s = y.getSum(); r &= s == 8.0;
  s = y.getMin(); r &= s == 3.0;
  s = y.getMax(); r &= s == 5.0;

  return r;
}

bool float32x4UnitTest()
{
  bool r = true;      // test result

  float a = 2.f, b = 3.f, c = 5.f, d = 7.f;

  // construct from a float and array access:
  rsFloat32x4 x1(a);  r &= x1[0] == a && x1[1] == a && x1[2] == a && x1[3] == a;

  // construct from 4 floats:
  rsFloat32x4 x2(a, b, c, d); 
  r &= x2[0] == a && x2[1] == b && x2[2] == c && x2[3] == d;

  // construct from array of floats:
  float arr[4] = { a, b, c, d };
  rsFloat32x4 x3(arr); r &= x3[0] == a && x3[1] == b && x3[2] == c && x3[3] == d;

  // construction from another instance and == operator:
  rsFloat32x4 y(x3); r &= y == x3;

  // setters:
  y[0] = d; r &= y[0] == d && y[1] == b && y[2] == c && y[3] == d;
  y[1] = c; r &= y[0] == d && y[1] == c && y[2] == c && y[3] == d;
  y.set(c); r &= y[0] == c && y[1] == c && y[2] == c && y[3] == c;
  y.set(a, b, c, d); r &= y[0] == a && y[1] == b && y[2] == c && y[3] == d;
  y.set(c, d, b, a); r &= y[0] == c && y[1] == d && y[2] == b && y[3] == a; // cdba

  // getters:
  //float y0, y1, y2, y3;
  //y.get(y0, y1, y2, y3); r &= y0 == a && y1 == b && y2 == c && y3 == d;
  y.get(arr); r &= arr[0] == c && arr[1] == d && arr[2] == b && arr[3] == a;

  // assignment and equality:
  y = x1; r &= y == x1;  // abcd again

  // binary arithmetic operators:
  rsFloat32x4 z, x;
  x.set(1.f, 2.f, 3.f, 4.f);
  y.set(5.f, 6.f, 7.f, 8.f);
  z = x + y; r &= z[0] == 6 && z[1] == 8 && z[2] == 10 && z[3] == 12;
  z = y - x; r &= z[0] == 4 && z[1] == 4 && z[2] == 4 && z[3] == 4;
  z = x * y; r &= z[0] == 5 && z[1] == 12 && z[2] == 21 && z[3] == 32;
  z = z / x; r &= z == y;

  // binary arithmetic operators with scalar lhs:
  z = 12.f + x; r &= z[0] == 13 && z[1] == 14 && z[2] == 15 && z[3] == 16;
  z = 12.f - x; r &= z[0] == 11 && z[1] == 10 && z[2] ==  9 && z[3] ==  8;
  z = 12.f * x; r &= z[0] == 12 && z[1] == 24 && z[2] == 36 && z[3] == 48;
  z = 12.f / x; r &= z[0] == 12 && z[1] ==  6 && z[2] ==  4 && z[3] ==  3;

  // binary arithmetic operators with scalar rhs:
  x.set(10, 20, 30, 40);
  z  = x + 2; r &= z[0] == 12 && z[1] == 22 && z[2] == 32 && z[3] == 42;
  z  = x - 2; r &= z[0] ==  8 && z[1] == 18 && z[2] == 28 && z[3] == 38;
  z  = x * 2; r &= z[0] == 20 && z[1] == 40 && z[2] == 60 && z[3] == 80;
  z  = x / 2; r &= z[0] ==  5 && z[1] == 10 && z[2] == 15 && z[3] == 20;

  // combined arithemtic/assignment
  z.set(10, 20, 30, 40);
  z += 2; r &= z[0] == 12 && z[1] == 22 && z[2] == 32 && z[3] == 42;
  z -= 2; r &= z[0] == 10 && z[1] == 20 && z[2] == 30 && z[3] == 40;
  z *= 2; r &= z[0] == 20 && z[1] == 40 && z[2] == 60 && z[3] == 80;
  z /= 2; r &= z[0] == 10 && z[1] == 20 && z[2] == 30 && z[3] == 40;



  // reductions to scalar:
  float s;
  x.set(1, 2, 3, 4);
  s = x.getSum(); r &= s == 10.f;

  return r;

  // when we later have a rsFloat64x4 class (SSE3), maybe templatize this unit test such that
  // it can be used with any kind of 4-vector
}

bool complexFloat64x2UnitTest()
{
  bool r = true;      // test result

  // we have 4 complex numbers z1[0] = 1 + 3i, z1[1] = 2 + 4i, z2[0] = 5 + 7i, z2[1] = 6 + 8i:
  std::complex<double> z10(1, 3), z11(2, 4), z20(5, 7), z21(6, 8), w0, w1;
  rsFloat64x2 re1(1, 2), im1(3, 4), re2(5, 6), im2(7, 8);
  std::complex<rsFloat64x2> z1(re1, im1), z2(re2, im2), w;
  //std::complex<rsFloat64x2> w; // for outputs
  //std::complex<double> w0, w1;

  // addition:
  w = z1 + z2;
  r &= w.real().get0() ==  6;
  r &= w.imag().get0() == 10;
  r &= w.real().get1() ==  8;
  r &= w.imag().get1() == 12;

  // division:
  w = z2 / z2;
  r &= w.real().get0() == 1;
  r &= w.imag().get0() == 0;
  r &= w.real().get1() == 1;
  r &= w.imag().get1() == 0;
  w   = z2;
  w  /= z2;
  r &= w.real().get0() == 1;
  r &= w.imag().get0() == 0;
  r &= w.real().get1() == 1;
  r &= w.imag().get1() == 0;

  //z = std::exp(z1); // this doesn't work - it doesn't try to invoke the exp for rsFloat64x2
  // i think, we need to implement explicit specializations for the math functions for
  // complex<rsFloat64x2>

  // exponential function:
  w0 = std::exp(z10);
  w1 = std::exp(z11);
  w  = rsExp(z1);
  r &= w0 == get0(w);
  r &= w1 == get1(w);
  //r &= w0 == rosic::get0(w);
  //r &= w1 == rosic::get1(w);

  return r;
}
