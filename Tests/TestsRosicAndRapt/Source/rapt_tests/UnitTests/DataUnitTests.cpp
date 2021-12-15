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

bool testContainerFuncs()
{
  bool ok = true;

  typedef std::vector<int> Vec;
  Vec u, v, w;

  // Test checking if an arbitrary number of vectors have all the same size:
  u = {1,2,3};
  v = {4,5,6};
  w = {7,8,9};
  ok &= rsAreSameSize(u, v);
  ok &= rsAreSameSize(u, v, w);
  ok &= rsAreSameSize(u, v, w, u);
  w = {7,8};
  ok &= !rsAreSameSize(u, v, w);
  ok &= !rsAreSameSize(u, v, u, w);

  // Check finding the smallest size of a bunch of containers
  u = {1,2};
  v = {3,4,5};
  w = {6,7,8,9};
  ok &= rsMinSize(u, v) == 2;
  ok &= rsMinSize(v, u) == 2;
  ok &= rsMinSize(u, v, w) == 2;
  ok &= rsMinSize(v, w, u) == 2;
  ok &= rsMinSize(w, u, v) == 2;

  return ok;
}

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

bool testArrayMisc()
{
  bool r = true;
  typedef std::vector<int> Vec;
  typedef RAPT::rsArrayTools AT;
  Vec x6, x7, x8, y6, y7, y8;

  // test reversal:
  x6 = Vec({1,2,3,4,5,6});
  x7 = Vec({1,2,3,4,5,6,7});
  x8 = Vec({1,2,3,4,5,6,7,8});
  y6 = x6; AT::reverse(&y6[0], 6); r &= y6 == Vec({6,5,4,3,2,1});
  y7 = x7; AT::reverse(&y7[0], 7); r &= y7 == Vec({7,6,5,4,3,2,1});

  // test circular shift:
  y7 = x7; rsCircularShift(&y7[0], 7, 3); r &= y7 == Vec({5,6,7,1,2,3,4});
  y8 = x8; rsCircularShift(&y8[0], 8, 3); r &= y8 == Vec({6,7,8,1,2,3,4,5});





  return r;
}

bool testTokenize()
{
  using namespace rosic; // ToDo: move into rosic tests!

  // Test splitting a string at a bunch of given split characters such as whitespaces and newlines.
  // The tokenizer should support tokens being seperated by any number of such split characters. We
  // have the unit test here because a string is basically also just an array (of characters). 
  // Maybe someday it can be moved elsewhere, when we have more string processing functions.

  bool ok = true;
  std::string str("ABC DE FG\nHI \nJKL\nMN   OPQRS \n TUVWXY\nZ \n");
  std::string sep(" \n");  // seperators

  int S =  0;  // start of the token
  int L = -1;  // length of token, initial value should be irrelevant

  rsFindToken(str, sep, &S, &L); ok &= S ==  0 && L == 3; S += L;  // ABC
  rsFindToken(str, sep, &S, &L); ok &= S ==  4 && L == 2; S += L;  // DE
  rsFindToken(str, sep, &S, &L); ok &= S ==  7 && L == 2; S += L;  // FG
  rsFindToken(str, sep, &S, &L); ok &= S == 10 && L == 2; S += L;  // HI
  rsFindToken(str, sep, &S, &L); ok &= S == 14 && L == 3; S += L;  // JKL
  rsFindToken(str, sep, &S, &L); ok &= S == 18 && L == 2; S += L;  // MN
  rsFindToken(str, sep, &S, &L); ok &= S == 23 && L == 5; S += L;  // OPQRS
  rsFindToken(str, sep, &S, &L); ok &= S == 31 && L == 6; S += L;  // TUVWXY
  rsFindToken(str, sep, &S, &L); ok &= S == 38 && L == 1; S += L;  // Z

  // Test some special cases:
  S = 0; L = -1; str = ("ABC");    // only 1 token, no separators
  rsFindToken(str, sep, &S, &L); ok &= S == 0 && L == 3; S += L;
  rsFindToken(str, sep, &S, &L); ok &= S == 3 && L == 0; S += L;
  // L == 0 indicates that the end of the string was reached. In practice, one could tokenize a 
  // string in a while(L != 0) loop..or maybe a while(s < str.length()) should also work

  S = 0; L = -1; str = ("   ");    // no token, only separators
  rsFindToken(str, sep, &S, &L); ok &= S == 3 && L == 0; S += L;
  rsFindToken(str, sep, &S, &L); ok &= S == 3 && L == 0; S += L;

  S = 0; L = -1; str = ("");       // empty string
  rsFindToken(str, sep, &S, &L); ok &= S == 0 && L == 0;
  rsFindToken(str, sep, &S, &L); ok &= S == 0 && L == 0;

  S = 0; L = -1; str = ("  ABC");  // starts with seperators
  rsFindToken(str, sep, &S, &L); ok &= S == 2 && L == 3; S += L;

  S = 0; L = -1; str = ("ABC  ");  // ends with seperators
  rsFindToken(str, sep, &S, &L); ok &= S == 0 && L == 3; S += L;
  rsFindToken(str, sep, &S, &L); ok &= S == 5 && L == 0; S += L;

  S = 0; L = -1; str = ("ABC DE"); // ends with token
  rsFindToken(str, sep, &S, &L); ok &= S == 0 && L == 3; S += L;
  rsFindToken(str, sep, &S, &L); ok &= S == 4 && L == 2; S += L;
  rsFindToken(str, sep, &S, &L); ok &= S == 6 && L == 0; S += L;

  return ok;
}

bool rsArrayViewTest() 
{
  bool ok = true;

  using Vec = std::vector<int>;
  using AV  = rsArrayView<int>;

  Vec vec     = {4,5,7,3,8,6,4,6,2,3};
  int raw[10] = {4,5,7,3,8,6,4,6,2,3};

  std::sort(vec.begin(), vec.end());

  // Now, let's try to apply std::sort ot our raw array by wrapping it into an rsArrayView which
  // gives it (partially) the interface of std::vector
  int N = 10;
  AV wrapped(raw, N);

  AV::iterator itBegin = wrapped.begin();
  AV::iterator itEnd   = wrapped.end();    // shall not be dereferenced
  ok &= wrapped[itBegin] == 4;

  AV::iterator it = itBegin;
  ok &= wrapped[it]   == 4; ++it;
  ok &= wrapped[it]   == 5; ++it;
  ok &= wrapped[it]   == 7; it++;
  ok &= wrapped[it]   == 3; it++;
  ok &= wrapped[it]   == 8;
  ok &= wrapped[it++] == 8;
  ok &= wrapped[it]   == 6;
  ok &= wrapped[++it] == 4;
  ok &= wrapped[--it] == 6;
  ok &= wrapped[it--] == 6;
  ok &= wrapped[it]   == 8;


  using AT = rsArrayTools;
  AT::stdSort(raw, N);
  ok &= AT::isSortedAscending(raw, N);


  // ToDo:
  // -Wrap and test some more algorithms from std::algorithm, see
  //  https://en.cppreference.com/w/cpp/algorithm
  //  https://docs.microsoft.com/en-us/cpp/standard-library/algorithm-functions?view=msvc-170
  //  https://www.youtube.com/watch?v=bFSnXNIsK4A
  // -What about algorithms that actually change the size of the container (like remove) or 
  //  produce a new container as output (like merge)? I think, remove returns an iterator to
  //  the new end and merge takes an iterator to the output maybe a bit like AT::convolve.
  //  -> figure out


  //std::sort(&raw[0], &raw[10]);  // OK - that works
  // can we wrap this into a call like sort(raw, 10) maybe in rsArrayTools, have a function
  // stdSort(T* x, int N) { std::sort(&x[0], &x[N]); }


  //std::sort(wrapped.begin(), wrapped.end());  // that doesn't yet
  // todo: implement binary -  for iterator

  return ok;
}



template<class T>
class rsNonReAllocArrayTest : public rsNonReAllocatingArray<T>
{

public:


  // Helper function to fill an array with given number of elements, linearly ascending, starting
  // at given startValue (like std::iota == "increment over the array")
  void fill(size_t numElems, T startValue = T(0))
  {
    for(size_t i = 0; i < numElems; i++)
      push_back(T(startValue + i));
  }

  static bool testIndexComputation()
  {
    bool ok = true;

    rsNonReAllocArrayTest<int> a;
    a.reserve(4);

    auto test = [&](size_t i, size_t jTarget, size_t kTarget)
    {
      size_t j, k;
      a.flatToChunkAndElemIndex(i, j, k); 
      return j == jTarget && k == kTarget;
    };
    ok &= test( 0, 0,  0);
    ok &= test( 1, 0,  1);
    ok &= test( 2, 0,  2);
    ok &= test( 3, 0,  3);

    ok &= test( 4, 1,  0);
    ok &= test( 5, 1,  1);
    ok &= test( 6, 1,  2);
    ok &= test( 7, 1,  3);

    ok &= test( 8, 2,  0);
    ok &= test( 9, 2,  1);
    ok &= test(10, 2,  2);
    ok &= test(11, 2,  3);
    // ...
    ok &= test(14, 2,  6);
    ok &= test(15, 2,  7);

    ok &= test(16, 3,  0);
    ok &= test(17, 3,  1);
    // ...
    ok &= test(30, 3, 14);
    ok &= test(31, 3, 15);

    ok &= test(32, 4,  0);
    ok &= test(33, 4,  1);
    // ...
    ok &= test(62, 4, 30);
    ok &= test(63, 4, 31);


    a.clear();
    ok &= a.chunks.size() == 0;
    //ok &= test( 0, 0,  0);      // should trigger assert - yep!
    a.reserve(1);
    ok &= test( 0, 0,  0);
    ok &= test( 1, 1,  0);
    ok &= test( 2, 2,  0);
    ok &= test( 3, 2,  1);
    ok &= test( 4, 3,  0);
    ok &= test( 5, 3,  1);
    ok &= test( 6, 3,  2);
    ok &= test( 7, 3,  3);
    ok &= test( 8, 4,  0);

    return ok;
  }
 
  static bool testPushBack()
  {
    bool ok = true;

    rsNonReAllocArrayTest<int> a;
    a.reserve(4);

    // The first 4 elements fit in the first chunk:
    a.push_back(0);
    a.push_back(1);
    a.push_back(2);
    a.push_back(3);
    ok &= a.chunks.size() == 1;

    // Now, growth should happen because the capacity of the initial chunk is exceeded:
    a.push_back(4);
    ok &= a.chunks.size() == 2;
    a.push_back(5);
    a.push_back(6);
    a.push_back(7);
    ok &= a.chunks.size() == 2;

    // Now, growth should happen for the 2nd time:
    a.push_back(8);
    ok &= a.chunks.size() == 3;

    // Test element random access:
    for(size_t i = 0; i <= 8; i++)
      ok &= a[i] == i;

    //a[9];  // should trigger assert - yep, does

    a.clear();                   // will put it back into initial state
    a.push_back(0);              // will reserve only 1 slot
    ok &= a.chunks.size() == 1;
    a.push_back(1);
    ok &= a.chunks.size() == 2;  // now we should have 2 slots of size 1
    a.push_back(2);
    ok &= a.chunks.size() == 3; 
    a.push_back(3);
    ok &= a.chunks.size() == 3; 
    a.push_back(4);
    ok &= a.chunks.size() == 4; 

    return ok;
  }

  static bool testIterator(int initialCapacity, int length)
  {
    bool ok = true;

    using Arr = rsNonReAllocArrayTest<int>;
    using It  = Arr::iterator;

    Arr a;
    a.reserve(initialCapacity);
    size_t N     = length;                // length of test array to generate
    size_t start = 1000;
    a.fill(N, 1000);
    ok &= a.size() == N;
    ok &= a.capacity() == RAPT::rsNextPowerOfTwo(length);  // fails!

    // Test manually incrementing and decrementing iterators:
    It it = a.begin();                    // get an iterator pointing to the begin
    for(size_t i = 0; i < a.size(); i++)  // iterate forward manually
    {
      ok &= *it == start + i; 
      it++; 
    }
    ok &= it == a.end();
    for(size_t i = N; i > 0; i--)
    {
      it--;
      ok &= *it == start + i-1;
    }
    ok &= it == a.begin(); 

    // Test a range-based loop over all values v in a. This also covers the pre-increment:
    size_t i = 0;
    for(const auto & v : a)
    {
      ok &= v == start + i;
      i++;
    }

    // ToDo: 
    // -Implement and test +,-,+=,-= operators -> figure out, how this iterator_difference
    //  thing is supposed to work

    // Notes:
    // -When we decrement a begin() iterator, we get an access violation. Maybe that should be
    //  expected because decrementing a begin() iterator is a bug anyway? Or should we somehow
    //  safeguard against that?

    // Ressources:
    // https://www.cplusplus.com/reference/iterator/

    return ok;
  }

  static bool testInsert()
  {
    bool ok = true;

    using Arr = rsNonReAllocArrayTest<int>;
    using It  = Arr::iterator;

    Arr a;
    a.reserve(4);
    size_t N = 6;
    size_t start = 1000;
    a.fill(N, 1000);
    ok &= a.size() == N;
    ok &= a.capacity() == 8;

    It it = a.begin();
    it++;
    it++;
    //a.insert(it, 2);


    // See:
    // https://en.cppreference.com/w/cpp/container/vector/insert

    return ok;
  }



  // template for copy-and-paste for adding a new test (nothing to do with c++ temaplates):
  static bool testTemplate()
  {
    bool ok = true;
    return ok;
  }

};




bool rsNonReAllocatingArrayTest()
{
  bool ok = true;

  using NAA = rsNonReAllocArrayTest<int>;

  NAA a;

  ok &= NAA::testIndexComputation();
  ok &= NAA::testPushBack();
  ok &= NAA::testIterator(4, 25);
  ok &= NAA::testIterator(2, 25);
  ok &= NAA::testIterator(1, 25);
  ok &= NAA::testIterator(8, 95);
  ok &= NAA::testInsert();

  // ToDo next:
  // -element insert (may move elements around, should use iterators and rsSwap)
  // -element removal (this may deallocate)


  // -try it also with initial capacities of 0,1 and 2

  return ok;
}


bool arrayUnitTest()  // maybe rename to stdVectorUnitTest
{
  bool ok = true;

  typedef std::vector<int> Vec;
  Vec v({ 1,2,3 }), w({4,5});
  Vec u = v;
  rsAppend(u, w);
  ok &= u == Vec({1,2,3,4,5});
  u = v;
  rsAppend(u, u);
  ok &= u == Vec({1,2,3,1,2,3});


  Vec a = { 0,1,2,3,4,5,6,7,8,9 };
  rsRemoveRange(a, 4, 7);
  ok &= a == Vec({ 0,1,2,3, 8,9 });


  ok &= testContainerFuncs();
  ok &= testArrayFiltering();
  ok &= testArrayMisc();
  ok &= testTokenize();
  ok &= rsArrayViewTest();
  ok &= rsNonReAllocatingArrayTest();

  // int s = sum(3, &u[0]); // sum function doesn't compile




  return ok;
}

template<class T>
class rsBinaryHeapTest : public rsBinaryHeap<T>
{

public:

  using rsBinaryHeap<T>::rsBinaryHeap;


  bool isMinHeap(int i = 0) const
  {
    if(i >= this->size)
      return true;
    bool result = true;
    int l = this->left(i);
    int r = this->right(i);
    if(l < this->size) result &= this->data[i] <= this->data[l] && isMinHeap(l);
    if(r < this->size) result &= this->data[i] <= this->data[r] && isMinHeap(r);
    return result;
  }

  bool isMaxHeap(int i = 0) const
  {
    if(i >= this->size)
      return true;
    bool result = true;
    int l = this->left(i);
    int r = this->right(i);
    if(l < this->size) result &= this->data[i] >= this->data[l] && isMaxHeap(l);
    if(r < this->size) result &= this->data[i] >= this->data[r] && isMaxHeap(r);
    return result;
  }

  int floatDownRec(int i)
  {
    int l = this->left(i);
    int r = this->right(i);
    int b = i;         // b for "big"
    if(l < this->size && this->less(this->data[i], this->data[l])) b = l;
    if(r < this->size && this->less(this->data[b], this->data[r])) b = r;
    if(b != i) { swap(this->data[i], this->data[b]); return floatDownRec(b); }
    return i;
  }
  // a.k.a. maxHeapify
  // That's the recursive implementation from (1) page 130. When the iterative version is ready,
  // move it to the rsBinaryHeapTest subclass - we don't need it anymore in production code, then.
  // But it may be interesting to figure out, if the recursion actually incurs an overhead since
  // it's tail recursion and smart compilers might be able to translate it to iteration themselves.
  // We also may want to keep it as reference for unit tests (to test, if the iterative version
  // really does the same thing).
  // maybe move to subclass rsBinaryHeapTest in the unit test section - we don't need it in
  // production code, when the iterative version works (which seems to be the case)

  // removing:
  // http://www.mathcs.emory.edu/~cheung/Courses/171/Syllabus/9-BinTree/heap-delete.html
  // https://www.geeksforgeeks.org/insertion-and-deletion-in-heaps/

  // -what if we have overflow in the l,r values - should we handle that? or maybe restrict the
  //  capacity to values which ensure that no overflow occurs? yes, that seems sensible

  void sort()
  {

  }
  // todo: implement heap-sort in this class an test it with various random arrays

};

bool binaryHeapUnitTest()
{

  //int a;
  //unsigned int b = 5;
  //int c = -10;
  //a = b + c;  // 5 - -10 = 15
  ////<signed int> = <unsigned int> <op> <signed int>
  //// https://github.com/fish-shell/fish-shell/issues/3493

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
  for(int i = 1; i <= numTests; i++)
  {
    int oldSize  = H.getSize();
    int remIndex = ng.getSampleRaw() % H.getSize();
    H.remove(remIndex);
    r &= H.getSize() == oldSize - 1;
    r &= H.isMaxHeap();
  }

  // test removing last element - it could be that we need to treat that case as special case - but
  // maybe not:
  N = H.getSize();
  H.remove(N-1);  // try to remove last element
  r &= H.getSize() == N-1;
  r &= H.isMaxHeap();
  // ...it seems to work - in this case - todo: test more cases

  // try to remove 0th element:
  H.remove(0);
  r &= H.getSize() == N-2;
  r &= H.isMaxHeap();





  // test the double heap:

  rsDoubleHeap<int> D;
  int min = 2147483648;
  Vec B;
  int i;
  A = Vec({5,2,3});
  B = Vec({6,7,8});
  D.setData(A, B);
  i = D.replace(1, 6);  // replace 2 by 6, float up to front of small, no exchange
  r &= i == 0 && A == Vec({6,5,3}) && B ==  Vec({6,7,8});
  i = D.replace(min+2, 4);  // replace 8 by 4, float to front of large, exchange
  r &= i == 1 && A == Vec({5,4,3}) && B ==  Vec({6,7,6});
  // do more tests, using larger heaps, maybe check property in a loop


  rsDoubleHeap<int> D2;
  D2.setData(A, B);

  int k;
  k = D2.indexToKey(0); r &= k == 0;
  i = D2.keyToIndex(0); r &= i == 0;
  k = D2.indexToKey(1); r &= k == 1;
  i = D2.keyToIndex(1); r &= i == 1;
  k = D2.indexToKey(2); r &= k == 2;
  i = D2.keyToIndex(2); r &= i == 2;
  k = D2.indexToKey(3); r &= k == min+0;
  i = D2.keyToIndex(k); r &= i == 3;
  k = D2.indexToKey(4); r &= k == min+1;
  i = D2.keyToIndex(k); r &= i == 4;
  k = D2.indexToKey(5); r &= k == min+2;
  i = D2.keyToIndex(k); r &= i == 5;

  int v;
  v = D2.atKey(0);
  v = D2.atKey(1);
  v = D2.atKey(2);
  //v = D2[3]; // these are access violations because rsDoubleHeap2 uses a different way to
  //v = D2[4]; // indicate a value from the large heap
  //v = D2[5];
  i = D2.replace(0, 9); r &= i < 1000;
  // should go into large heap - index should be a large neagtive number

  i = D2.replace(i, 2); // this brings it back to the small heap again








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

  rsDelayBuffer<double> b(8);
  b.setLength(5);

  using Vec = std::vector<double>;
  Vec B(5);

  double y, z;

  // test getSample and copying the buffer content:
  y = b.getSample( 1); r &= y ==  0; b.copyTo(&B[0], true ); r &= B == Vec({  1,  0,  0,  0,  0 });
  y = b.getSample( 2); r &= y ==  0; b.copyTo(&B[0], true ); r &= B == Vec({  2,  1,  0,  0,  0 });
  y = b.getSample( 3); r &= y ==  0; b.copyTo(&B[0], true ); r &= B == Vec({  3,  2,  1,  0,  0 });
  y = b.getSample( 4); r &= y ==  0; b.copyTo(&B[0], true ); r &= B == Vec({  4,  3,  2,  1,  0 });
  y = b.getSample( 5); r &= y ==  0; b.copyTo(&B[0], true ); r &= B == Vec({  5,  4,  3,  2,  1 });
  y = b.getSample( 6); r &= y ==  1; b.copyTo(&B[0], true ); r &= B == Vec({  6,  5,  4,  3,  2 });
  y = b.getSample( 7); r &= y ==  2; b.copyTo(&B[0], true ); r &= B == Vec({  7,  6,  5,  4,  3 });
  y = b.getSample( 8); r &= y ==  3; b.copyTo(&B[0], true ); r &= B == Vec({  8,  7,  6,  5,  4 });
  y = b.getSample( 9); r &= y ==  4; b.copyTo(&B[0], true ); r &= B == Vec({  9,  8,  7,  6,  5 });
  y = b.getSample(10); r &= y ==  5; b.copyTo(&B[0], true ); r &= B == Vec({ 10,  9,  8,  7,  6 });
  y = b.getSample(11); r &= y ==  6; b.copyTo(&B[0], true ); r &= B == Vec({ 11, 10,  9,  8,  7 });
  y = b.getSample(12); r &= y ==  7; b.copyTo(&B[0], true ); r &= B == Vec({ 12, 11, 10,  9,  8 });
                                     b.copyTo(&B[0], false); r &= B == Vec({  8,  9, 10, 11, 12 });
  y = b.getSample(13); r &= y ==  8; b.copyTo(&B[0], true ); r &= B == Vec({ 13, 12, 11, 10,  9 });
                                     b.copyTo(&B[0], false); r &= B == Vec({  9, 10, 11, 12, 13 });
  y = b.getSample(14); r &= y ==  9; b.copyTo(&B[0], true ); r &= B == Vec({ 14, 13, 12, 11, 10 });
                                     b.copyTo(&B[0], false); r &= B == Vec({ 10, 11, 12, 13, 14 });
  y = b.getSample(15); r &= y == 10; b.copyTo(&B[0], true ); r &= B == Vec({ 15, 14, 13, 12, 11 });
                                     b.copyTo(&B[0], false); r &= B == Vec({ 11, 12, 13, 14, 15 });
  y = b.getSample(16); r &= y == 11; b.copyTo(&B[0], true ); r &= B == Vec({ 16, 15, 14, 13, 12 });
                                     b.copyTo(&B[0], false); r &= B == Vec({ 12, 13, 14, 15, 16 });
  // todo: complete this - test forward and reverse copy (true/false) after each getSample

  z = b.getNewest(); y = b.fromNewest(0); r &= y == z;
  z = b.getOldest(); y = b.fromOldest(0); r &= y == z;

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

  //r &= q.getMaxLength() == 6;     // a power of 2 minus 2, >= requested size
  r &= q.getMaxLength() == 7;     // a power of 2 minus 1, >= requested size

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

  // this triggers an assertion but actually still works - however, the MovingMax filter does not
  // work properly anymore when we fill the deque up to bufferSize-1
  q.pushFront(7);                // [1 2 3 4 5 6 7]
  r &= q.getLength() == 7;
  r &= q.readHead()  == 7;
  r &= q.readTail()  == 1;

  r &= q.isFull();

  // maybe provide writeHead/writeTail function that overwrite the current data in head/tail


  // maybe getLength is flawed and returns the correct length only in certain situations?
  // -> do more tests!
  // i actually think, the maxLength is size-1 and not size-2!

  // i think, the getLength formula works only for h != t

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

template<class T, int N>
bool simdTemplateUnitTest()
{
  bool ok = true;

  // Some variables to work with:
  rsSimdVector<T, N> a, b, c;
  int i;

  // Test arithmetic operators for vector (op) vector:
  for(i=0; i<N; i++) { a[i] = T(2*i+1); b[i] = T(3)*a[i]; } // init operands
  c = a+b; for(i=0; i<N; i++) { ok &= c[i] == a[i]+b[i]; }
  c = a-b; for(i=0; i<N; i++) { ok &= c[i] == a[i]-b[i]; }
  c = a*b; for(i=0; i<N; i++) { ok &= c[i] == a[i]*b[i]; }
  c = b/a; for(i=0; i<N; i++) { ok &= c[i] == b[i]/a[i]; }

  // Test arithmetic operators for vector (op) scalar:
  T s = T(5);  // some scalar
  c = a+s; for(i=0; i<N; i++) { ok &= c[i] == a[i]+s; }
  c = a-s; for(i=0; i<N; i++) { ok &= c[i] == a[i]-s; }
  c = a*s; for(i=0; i<N; i++) { ok &= c[i] == a[i]*s; }
  c = a/s; for(i=0; i<N; i++) { ok &= c[i] == a[i]/s; }

  // Test arithmetic operators for scalar (op) vector:
  c = s+a; for(i=0; i<N; i++) { ok &= c[i] == s+a[i]; }
  c = s-a; for(i=0; i<N; i++) { ok &= c[i] == s-a[i]; }
  c = s*a; for(i=0; i<N; i++) { ok &= c[i] == s*a[i]; }
  c = s/a; for(i=0; i<N; i++) { ok &= c[i] == s/a[i]; }

  // Test unary plus/minus:
  c = +a; for(i=0; i<N; i++) { ok &= c[i] == +a[i]; }
  c = -a; for(i=0; i<N; i++) { ok &= c[i] == -a[i]; }

  // Test unary math functions:
  c = rsAbs( a); for(i=0; i<N; i++) { ok &= c[i] == rsAbs( a[i]); }
  c = rsSign(a); for(i=0; i<N; i++) { ok &= c[i] == rsSign(a[i]); }
  // ...

  // todo: test binary math functions, comparison operators, copy/conversion constructors, 
  // assignment operators

  return ok;
}

template<class T, int N>
bool simdFloatUnitTest()
{
  // Like simdTemplateUnitTest but with additional tests that apply only to floating point types, 
  // such as the transcendental functions (todo: test nan/inf stuff, too)

  bool ok = simdTemplateUnitTest<T, N>();

  rsSimdVector<T, N> a, b, c;
  int i;
  for(i=0; i<N; i++) { a[i] = T(2*i+1); b[i] = T(3)*a[i]; } // init operands

  // Test float-specific unary math functions:
  c = rsCos(  a); for(i=0; i<N; i++) { ok &= c[i] == rsCos(  a[i]); }
  c = rsExp(  a); for(i=0; i<N; i++) { ok &= c[i] == rsExp(  a[i]); }
  c = rsLog(  a); for(i=0; i<N; i++) { ok &= c[i] == rsLog(  a[i]); }
  c = rsSin(  a); for(i=0; i<N; i++) { ok &= c[i] == rsSin(  a[i]); }
  c = rsSqrt( a); for(i=0; i<N; i++) { ok &= c[i] == rsSqrt( a[i]); }
  c = rsTan(  a); for(i=0; i<N; i++) { ok &= c[i] == rsTan(  a[i]); }
  c = rsFloor(a); for(i=0; i<N; i++) { ok &= c[i] == rsFloor(a[i]); }
  c = rsCeil( a); for(i=0; i<N; i++) { ok &= c[i] == rsCeil( a[i]); }
  c = rsRound(a); for(i=0; i<N; i++) { ok &= c[i] == rsRound(a[i]); }
  c = rsCosh( a); for(i=0; i<N; i++) { ok &= c[i] == rsCosh( a[i]); }
  c = rsSinh( a); for(i=0; i<N; i++) { ok &= c[i] == rsSinh( a[i]); }
  c = rsTanh( a); for(i=0; i<N; i++) { ok &= c[i] == rsTanh( a[i]); }

  return ok;
}

bool simdInstantiationUnitTest()
{
  // Tests for specific instantitaions of the template

  bool ok = true;


  ok &= rsSimdVector<int,  1>::getEmulationLevel() == 0;
  ok &= rsSimdVector<int,  2>::getEmulationLevel() == 1;
  ok &= rsSimdVector<int,  4>::getEmulationLevel() == 2;
  ok &= rsSimdVector<int,  8>::getEmulationLevel() == 3;

  ok &= rsSimdVector<float,  1>::getEmulationLevel() == 0;
  ok &= rsSimdVector<float,  2>::getEmulationLevel() == 1;

#if !defined(RS_NO_SIMD_FLOAT32X4)
  ok &= rsSimdVector<float,  4>::getEmulationLevel() == 0;
  ok &= rsSimdVector<float,  8>::getEmulationLevel() == 1;
  ok &= rsSimdVector<float, 16>::getEmulationLevel() == 2;
#else
  ok &= rsSimdVector<float,  4>::getEmulationLevel() == 2;
  ok &= rsSimdVector<float,  8>::getEmulationLevel() == 3;
  ok &= rsSimdVector<float, 16>::getEmulationLevel() == 4;
#endif


  ok &= rsSimdVector<double, 1>::getEmulationLevel() == 0;
  ok &= rsSimdVector<double, 2>::getEmulationLevel() == 1;
  ok &= rsSimdVector<double, 4>::getEmulationLevel() == 2;
  // using real simd for float64x2 is not yet implemented

  // ToDo: 
  // -the expected outcomes of these tests should later depend on the compiler settings and
  //  macro definitions...not yet sure, how to implement this

  return ok;
}

bool simdUnitTest()
{
  bool ok = true;

  ok &= float64x2UnitTest();
  ok &= float32x4UnitTest();
  ok &= complexFloat64x2UnitTest();
  // fails on linux ("illegal instruction") ...seems that illegal instruction is our
  // rsAsserFalse debug-break

  //ok &= simdTemplateUnitTest<float, 8>();  // for debug

  // Test the new implementation:            // has explicit specialization
  ok &= simdTemplateUnitTest<int,     1>();  // no
  ok &= simdTemplateUnitTest<int,     2>();  // no
  ok &= simdTemplateUnitTest<int,     4>();  // no
  ok &= simdTemplateUnitTest<int,     8>();  // no
  ok &= simdTemplateUnitTest<int,    16>();  // no
  // try also with 16 chars

  ok &= simdFloatUnitTest<float,      1>();  // no
  ok &= simdFloatUnitTest<float,      2>();  // no
  ok &= simdFloatUnitTest<float,      4>();  // yes (but incomplete)
  ok &= simdFloatUnitTest<float,      8>();  // no
  ok &= simdFloatUnitTest<float,     16>();  // no
  ok &= simdFloatUnitTest<double,     1>();  // no
  ok &= simdFloatUnitTest<double,     2>();  // no
  ok &= simdFloatUnitTest<double,     4>();  // no
  ok &= simdFloatUnitTest<double,     8>();  // no
  ok &= simdFloatUnitTest<double,    16>();  // no

  ok &= simdInstantiationUnitTest();


  // Test:
  //using TestClass = rsOperatorTest<double>;
  //TestClass a(3.0), b(5.0);

  return ok;
}
