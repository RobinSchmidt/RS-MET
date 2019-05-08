#pragma once


/** A class for generating ratios of numbers that can be used - for example - as frequency ratios
for the various oscillators in an oscillator array (for "supersaw"/"superwave" stuff), 
delay-lengths in a feedback delay network, and so on. */

template<class T>
class rsRatioGenerator
{

public:

  /** Used to select the formula/algorithm by which a ratio is computed. Determines what kind of 
  ratios will be produced in our dispatching methods. */
  enum class RatioKind // maybe rename to RatioFormula
  {
    metallic,
    primeSqrt,
    primeSqrtDiff,
    // plastic,
    // intSqrt,
    rangeSplitSkewed,
    rangeSplitOdd,
    rangeSplitEven
  };

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** If you want to use the formulas based on prime numbers, you must pass a pointer to vector of
  prime numbers and ensure that the vector pointed is still alive when you call a function that 
  generates a frcation from prime numbers. The size of the table must be ...large enough...tbc.. */
  void setPrimeTable(std::vector<RAPT::rsUint32>* newTable)
  {
    primeTable = newTable;
  }

  /** Sets the kind of ratios that should be produced. */
  void setRatioKind(RatioKind newKind) { kind = newKind; }

  /** Sets the first parameter for the formula - the meaning of the parameter varies depending on
  which formula is selected. */
  void setParameter1(T newParam) { p1 = newParam; }
  //...


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Fills the given array with the ratios of the selected kind with selected parameters. */
  void fillRatioTable(T* ratios, int numRatios);

  /** Returns the so called metallic ratio of given index n. For n = 0, it's just 1, for n = 1 it
  is called the golden ratio, n = 2: silver, n = 3: bronze and beyond that, they don't have names.
  https://en.wikipedia.org/wiki/Metallic_mean
  The index is actually supposed to be a nonegative integer, but for extended flexibility, the 
  function allows you to pass a real number as well. The golden ratio has many interesting 
  mathematical properties, like being the "most irrational" number possible in the sense that it's
  hardest to approximate by a continued fraction expansion, see here:
  https://www.youtube.com/watch?v=CaasbfdJdJg  */
  static inline T metallic(T n) { return T(0.5) * (n + sqrt(n*n+T(4))); }

  /** Square root of n-th prime number. */
  inline T primeSqrt(int n) 
  { 
    rsAssert(primeTable != nullptr);
    rsAssert(n < primeTable->size());
    return sqrt(T(primeTable->at(n)));
  }

  /** Difference of the square roots of n+1-th and n-th prime number. */
  inline T primeSqrtDiff(int n)
  {
    return primeSqrt(n+1) - primeSqrt(n);
  }

  // what about plastic ratios? oh - there's only one such ratio - but maybe powers of that can 
  // be used? what about powers of some general base?
  // https://en.wikipedia.org/wiki/Plastic_number

protected:

  RatioKind kind = RatioKind::metallic;
  T p1; // p2, p3, ...

  std::vector<RAPT::rsUint32>* primeTable = nullptr;
  // table of prime numbers - we use pointer to share it among all existing instances

};