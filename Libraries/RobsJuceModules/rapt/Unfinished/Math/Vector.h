#ifndef RAPT_VECTOR_H_INCLUDED
#define RAPT_VECTOR_H_INCLUDED

// get rid of this - use the new implementation of rsMatrix

/** This is a class for representing n-dimensional vectors and doing mathematical operations with
them. To make the class fast, no consistency checking is done in the mathematical operators. You
must ensure yourself that the input arguments are compatible - for example don't try to add two
vectors of different dimensionality and such.

\todo: make an explicit specialization for Complex elements (norm, etc. has to be implemented a
little different) ...hmmm - but rsComplex is itself a template class. can we pass another
template class as template argument? ...check this out 

\todo: use std::vector for storing the data

*/

template<class ElementType>  // maybe replace ElementType with T
class rsVector
{

public:

  /** \name Public Data Members \todo: make them protected */

  /** Dimensionality. */
  int dim;

  /** The actual values of the vector. */
  ElementType *v;  // use std::vector





  /** \name Construction/Destruction */

  /** Default constructor - constructs a vector with zero dimensionality and a NULL pointer. */
  rsVector()
  {
    dim = 0;
    v   = NULL;
  }

  /** Constructor. You must pass the number of elements (i.e. the dimensionality) of the vector
  here. It will create a zero vector of given dimensionality. */
  rsVector(int numElements)
  {
    dim = numElements;
    v   = new ElementType[dim];
    initWithZeros();
  }

  /** Constructor. You must pass the number of elements (i.e. the dimensionality) of the vector
  and an array with the initial values here. */
  rsVector(int numElements, ElementType *values)
  {
    dim = numElements;
    v   = new ElementType[dim];
    memcpy(v, values, dim*sizeof(ElementType));
  }

  /** Copy constructor. */
  rsVector(const rsVector& other)
  {
    dim = other.dim;
    v   = new ElementType[dim];
    memcpy(v, other.v, dim*sizeof(ElementType));
  }

  /** Destructor. */
  ~rsVector()
  {
    delete[] v;
  }


  /** \name Operators */

  /** Accesses the element at given index for reading and writing.
  (shouldn't it be non-const for writing?) */
  ElementType& operator[](const int index) const
  {
    return v[index];
  }

  /** Assigns one vector with another one. */
  rsVector& operator=(const rsVector& v2)
  {
    // \todo - check hwo this responds to self-assignment
    if(dim != v2.dim)
    {
      delete[] v;
      dim = v2.dim;
      v   = new ElementType[dim];
    }
    memcpy(v, v2.v, dim*sizeof(ElementType));
    return *this;
  }

  /** Compares two vectors for equality. */
  bool operator==(const rsVector& v2) const
  {
    if(dim != v2.dim)
      return false;
    else
    {
      for(int i=0; i<dim; i++)
      {
        if(v[i] != v2.v[i])
          return false;
      }
    }
    return true;
  }

  /** Compares two vectors for inequality. */
  bool operator!=(const rsVector& v2) const
  {
    return !(*this == v2);
  }

  /** Defines the negative of a vector. */
  rsVector operator-()
  {
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = -v[i];
    return result;
  }

  /** Adds another vector to this vector and returns the result. */
  rsVector& operator+=(const rsVector &v2)
  {
    rsAssert(dim == v2.dim);   // vectors incompatible
    for(int i=0; i<dim; i++)
      v[i] += v2.v[i];
    return *this;
  }

  /** Adds a scalar to this vector and returns the result. */
  rsVector& operator+=(const ElementType &x)
  {
    for(int i=0; i<dim; i++)
      v[i] += x;
    return *this;
  }

  /** Subtracts another vector from this vector and returns the result. */
  rsVector& operator-=(const rsVector &v2)
  {
    rsAssert(dim == v2.dim);   // vectors incopatible
    for(int i=0; i<dim; i++)
      v[i] -= v2.v[i];
    return *this;
  }

  /** Subtracts a scalar from this vector and returns the result. */
  rsVector& operator-=(const ElementType &x)
  {
    for(int i=0; i<dim; i++)
      v[i] -= x;
    return *this;
  }

  /** Multiplies this vector by a scalar and returns the result. */
  rsVector& operator*=(const ElementType &x)
  {
    for(int i=0; i<dim; i++)
      v[i] *= x;
    return *this;
  }

  /** Divides this vector by a scalar and returns the result. */
  rsVector& operator/=(const ElementType &x)
  {
    ElementType scale = 1.0 / x;
    for(int i=0; i<dim; i++)
      v[i] *= scale;
    return *this;
  }

  /** Adds two vectors. */
  rsVector operator+(const rsVector &v2)
  {
    rsAssert(dim == v2.dim);   // vectors incompatible
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] + v2.v[i];
    return result;
  }

  /** Adds a vector and a scalar. */
  rsVector operator+(const ElementType &x)
  {
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] + x;
    return result;
  }

  /** Subtracts two vectors. */
  rsVector operator-(const rsVector &v2)
  {
    rsAssert(dim == v2.dim);   // vectors incopatible
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] - v2.v[i];
    return result;
  }

  /** Subtracts a scalar from a vector. */
  rsVector operator-(const ElementType &x)
  {
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] - x;
    return result;
  }

  /** Multiplies a vector and a scalar. */
  rsVector operator*(const ElementType &x)
  {
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] * x;
    return result;
  }

  /** Divides a vector by a scalar. */
  rsVector operator/(const ElementType &x)
  {
    ElementType scale = 1.0 / x;
    rsVector result(dim);
    for(int i=0; i<dim; i++)
      result.v[i] = v[i] * scale;
    return result;
  }


  /** \name Manipulations */

  // removeElement, addElement, setDimension, normalize, ...

  /** Sets a new dimensionality (number of elements) for this vector. */
  void setDimensionality(int newDimensionality)
  {
    // re-allocate memory, if necessary:
    if(dim != newDimensionality)
    {
      delete[] v;
      dim = newDimensionality;
      v   = new ElementType[dim];
    }
    // \todo retain old elements up to the new dimensionality and fill the rest with zeros
  }

  /** Applies the passed function to each element of the vector. */
  void applyFunction(ElementType (*f) (ElementType))
  {
    for(int i=0; i<dim; i++)
      v[i] = f(v[i]);
  }

  /** Assigns random values between 'min' and 'max' to each element. */
  void randomizeElements(ElementType min, ElementType max)
  {
    for(int i=0; i<dim; i++)
      v[i] = randomUniform(min, max);
  }

  /** Initializes the values of the vector with all zeros. */
  void initWithZeros()
  {
    for(int i=0; i<dim; i++)
      v[i] = 0.0;
  }

  /** Initializes the values of the vector with the indices. */
  void initWithIndex()
  {
    for(int i=0; i<dim; i++)
      v[i] = i;
  }


  /** \name Inquiry */

  /** Returns the squared Euclidean norm of the vector. This is the sum of the squares of the
  elements: x1^2 + x2^2 + x3^2 + .... */
  ElementType getSquaredNorm()
  {
    ElementType accu = 0.0;
    for(int i=0; i<dim; i++)
      accu += v[i] * v[i];
    return accu;
  }

  /** Returns the Euclidean norm of the vector defined as sqrt(x1^2 + x2^2 + x3^2 + ...). */
  ElementType getEuclideanNorm() { return rsSqrt(getSquaredNorm()); }


  /** \name Misc */

  /** Stores the values into the target array which has to be of length matching the
  dimensionality of this vector. */
  void storeValuesInArray(ElementType *targetArray)
  {
    memcpy(targetArray, v, dim*sizeof(ElementType));
  }

  /** Prints the values to the standard output - mainly for debugging purposes. */
  void print()
  {
    printf("%s %d %s", "rsVector - dimensionality: ", dim, "\n");
    for(int i=0; i<dim; i++)
      printf("%.4f %s", (double)v[i], "  ");
    printf("%s", "\n");
  }
  // maybe remove from class and move to the test suite


}; // end of class rsVector

// some binary operators are defined outside the class such that the left hand operand does
// not necesarrily need to be of class rsVector

/** Adds a scalar and a vector. */
template<class ElementType>
RS_INLINE rsVector<ElementType> operator+(const ElementType &x, const rsVector<ElementType> &v)
{
  rsVector<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = v.v[i] + x;
  return result;
}

/** Subtracts a vector from a scalar. */
template<class ElementType>
RS_INLINE rsVector<ElementType> operator-(const ElementType &x, const rsVector<ElementType> &v)
{
  rsVector<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = x - v.v[i];
  return result;
}

/** Multiplies a scalar and a vector. */
template<class ElementType>
RS_INLINE rsVector<ElementType> operator*(const ElementType &x, const rsVector<ElementType> &v)
{
  rsVector<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = x * v.v[i];
  return result;
}

/** Divides a scalar by a vector. */
template<class ElementType>
RS_INLINE rsVector<ElementType> operator/(const ElementType &x, const rsVector<ElementType> &v)
{
  rsVector<ElementType> result(v.dim);
  for(int i=0; i<v.dim; i++)
    result.v[i] = x / v.v[i];
  return result;
}


/** Computes the inner product between two vectors (which should have the same dimensionality) -
the result is a scalar. */
template<class ElementType>
RS_INLINE ElementType operator*(const rsVector<ElementType> &a, const rsVector<ElementType> &b)
{
  rsAssert(a.dim == b.dim);   // vectors incompatible
  ElementType result = ElementType(0);
  for(int i=0; i<a.dim; i++)
    result += a.v[i] * b.v[i];
  return result;
}

/** Applies the unary function to the values of v and returns the result. */
template<class ElementType>
rsVector<ElementType> rsApplyFunction(rsVector<ElementType> v, ElementType (*f) (ElementType))
{
  rsVector<ElementType> r(v.dim);
  for(int i = 0; i < r.dim; i++)
    r[i] = f(v[i]);
  return r;
}

/** Constructs a vector using the binary function f with first argument from vector v and second
argument x.
note: unfortunately, it is not possible to pass the function pointer as first argument to this
function - the (MS) compiler can't resolve the template parameter in this case
*/
template<class T>
rsVector<T> rsApplyFunction(rsVector<T> v, T x, T (*f) (T, T))
{
  rsVector<T> r(v.dim);
  for(int i = 0; i < r.dim; i++)
    r[i] = f(v[i], x);
  return r;
}

/** Constructs a vector using the binary function f with first argument x and second
argument from vector v. */
template<class T>
rsVector<T> rsApplyFunction(T x, rsVector<T> v, T (*f) (T, T))
{
  rsVector<T> r(v.dim);
  for(int i = 0; i < r.dim; i++)
    r[i] = f(x, v[i]);
  return r;
}

/** Constructs a vector using the binary function f with first argument from vector v1 and second
argument from vector v2.
*/
template<class T>
rsVector<T> rsApplyFunction(rsVector<T> v1, rsVector<T> v2, T (*f) (T, T))
{
  rsVector<T> r(v1.dim);
  for(int i = 0; i < r.dim; i++)
    r[i] = f(v1[i], v2[i]);
  return r;
}


// maybe define a vector product (generalized cross-product via the Levi-Civita symbol) - maybe
// this has an inverse operation (vector division?)


// typedefs for explicit instantiations:
typedef rsVector<double> rsVectorDbl;
//typedef rsVector<rsComplexDbl> RSLib_API rsVectorCmplxDbl;


#endif
