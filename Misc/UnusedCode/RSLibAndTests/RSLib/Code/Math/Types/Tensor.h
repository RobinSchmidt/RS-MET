#ifndef RS_TENSOR_H
#define RS_TENSOR_H

namespace RSLib
{

  /**

  This is a class for representing n-dimensional tensors and doing mathematical operations with
  them.

  !!! THIS CLASS IS STILL VERY PRELIMINARY. DO NOT YET USE IT. !!!

  \todo: check, if it makes sense to templatize it (is there such a thing as complex tensors?)
  \todo: maybe use MultiArray as baseclass

  */

  class rsTensor
  {

  public:

    /** \name Construction/Destruction */

    /** Default constructor - constructs a vector with zero rank and dimensionality and a NULL
    pointer to the components. */
    rsTensor();

    /** Constructor. You must pass the dimensionality and rank of the tensor here. It will create a
    tensor of the disired size. If the 3rd parameter is true, the components will be initialized
    with all zeros, otherwise they will have undefined values (this parameter is optional and
    defaults to "true"). */
    rsTensor(int numDimensions, int numIndices, bool initWithZeros = true);

    /** Constructor. ou must pass the dimensionality and rank of the tensor here and a flat array
    with the initial values for the components here. */
    rsTensor(int numDimensions, int numIndices, double *components);

    /** Copy constructor. */
    rsTensor(const rsTensor& other);

    /** Destructor. */
    ~rsTensor();


    /** \name Inquiry */

    /** Returns the number of dimensions of the space in which this tensor lives. */
    RS_INLINE int getNumDimensions() const
    {
      return N;
    }

    /** Returns the number indices which this tensor has (also known as rank or order of the
    tensor). */
    RS_INLINE int getNumIndices() const
    {
      return R;
    }

    /** Returns the number of (scalar) components which this tensor has. This is given by N^R,
    where N is the number of dimensions and R is the rank of the tensor. */
    int getNumComponents() const;

    /** Returns true if this tensor is of the same type as the passed tensor reference, false
    otherwise. Same type means that the number of indices and the number of dimensions must
    match. Tensors of the same type can be added and subtracted. */
    bool isOfSameTypeAs(const rsTensor &t2) const
    {
      return N == t2.N && R == t2.R;
    }


    /** \name Operators */
    // \todo move implementations into .cpp file

    /** Assigns one tensor with another one.
    // \todo - check how this responds to self-assignment
    */
    rsTensor& operator=(const rsTensor& other)
    {
      int C = other.getNumComponents();
      if( N != other.N || R != other.R )
      {
        delete[] c;
        N = other.N;
        R = other.R;
        c = new double[C];
      }
      memcpy(c, other.c, C*sizeof(double));
      return *this;
    }


    /** Compares two tensors for equality. */
    bool operator==(const rsTensor& t2) const
    {
      if( N != t2.N || R != t2.R )
        return false;
      else
      {
        for(int i = 0; i < getNumComponents(); i++)
        {
          if( c[i] != t2.c[i] )
            return false;
        }
        // \todo: maybe use memcmp
      }
      return true;
    }

    /** Compares two tensors for inequality. */
    bool operator!=(const rsTensor& t2) const
    {
      return !(*this == t2);
    }

    /** Defines the negative of a tensor. */
    rsTensor operator-()
    {
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = -c[i];
      return result;
    }

    /** Adds another tensor to this tensor and returns the result. */
    rsTensor& operator+=(const rsTensor &t2)
    {
      rsAssert(isOfSameTypeAs(t2));   // tensors incompatible
      for(int i = 0; i < getNumComponents(); i++)
        c[i] += t2.c[i];
      return *this;
    }

    /** Adds a scalar to this tensor and returns the result. */
    rsTensor& operator+=(const double &x)
    {
      for(int i = 0; i < getNumComponents(); i++)
        c[i] += x;
      return *this;
    }

    /** Subtracts another tensor from this tensor and returns the result. */
    rsTensor& operator-=(const rsTensor &t2)
    {
      rsAssert(isOfSameTypeAs(t2));   // tensors incompatible
      for(int i = 0; i < getNumComponents(); i++)
        c[i] -= t2.c[i];
      return *this;
    }

    /** Subtracts a scalar from this tensor and returns the result. */
    rsTensor& operator-=(const double &x)
    {
      for(int i = 0; i < getNumComponents(); i++)
        c[i] -= x;
      return *this;
    }

    /** Multiplies this tensor by a scalar and returns the result. */
    rsTensor& operator*=(const double &x)
    {
      for(int i = 0; i < getNumComponents(); i++)
        c[i] *= x;
      return *this;
    }

    /** Divides this tensor by a scalar and returns the result. */
    rsTensor& operator/=(const double &x)
    {
      double scale = 1.0 / x;
      for(int i = 0; i < getNumComponents(); i++)
        c[i] *= scale;
      return *this;
    }

    /** Adds two tensors. */
    rsTensor operator+(const rsTensor &t2)
    {
      rsAssert(isOfSameTypeAs(t2));   // tensors incompatible
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] + t2.c[i];
      return result;
    }

    /** Adds a tensor and a scalar. */
    rsTensor operator+(const double &x)
    {
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] + x;
      return result;
    }

    /** Subtracts two tensors. */
    rsTensor operator-(const rsTensor &t2)
    {
      rsAssert(isOfSameTypeAs(t2));   // tensors incompatible
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] - t2.c[i];
      return result;
    }

    /** Subtracts a scalar from a tensor. */
    rsTensor operator-(const double &x)
    {
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] - x;
      return result;
    }

    /** Multiplies a tensor and a scalar. */
    rsTensor operator*(const double &x)
    {
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] * x;
      return result;
    }

    /** Divides a tensor by a scalar. */
    rsTensor operator/(const double &x)
    {
      double scale = 1.0 / x;
      rsTensor result(N, R, false);
      for(int i = 0; i < getNumComponents(); i++)
        result.c[i] = c[i] * scale;
      return result;
    }


    /** \name Manipulations */

    /** Applies the passed function to each element of the tensor. */
    void applyFunction( double (*f) (double) )
    {
      for(int i = 0; i < getNumComponents(); i++)
        c[i] = f(c[i]);
    }

    /** Sets all components to zero. \todo: use memset instead of a loop */
    void rsFillWithZeros()
    {
      for(int i = 0; i < getNumComponents(); i++)
        c[i] = 0.0;
    }


    /** \name Misc */

    /** Prints the values to the standard output - mainly for debugging purposes. */
    /*
    void print()
    {
      printf("%s %d %s", "rsTensor - dimensionality: ", dim, "\n");
      for(int i=0; i<dim; i++)
        printf("%.4f %s", (double) v[i], "  ");
      printf("%s", "\n");
    }
    */

   protected:

    /** \name Data */

    /** Dimensionality of the space in which the tensor lives. */
    int N;

    /** Rank (aka order) of the tensor - that is: the number of indices. */
    int R;

    /** The components of the tensor. */
    double *c;

    // maybe use numDimensions, numIndices and components as self-explaining names, maybe have
    // another variable numComponents = numDimensions^numIndices (this would be redundant but maybe
    // worthwhile to increase efficiency nonetheless)
    // or maybe use nD, nI, nC
    // ...or numSuperscripts, numSubscripts, numDimensions, (numComponents)

    // maybe, we need some way of representing the basis tensors and a way to flag indices as
    // covariant or contravariant

    friend rsTensor operator+(const double &x, const rsTensor &t);
    friend rsTensor operator-(const double &x, const rsTensor &t);
    friend rsTensor operator*(const double &x, const rsTensor &t);

  }; // end of class rsTensor


  // some binary operators are defined outside the class such that the left hand operand does
  // not necesarrily need to be of class rsTensor

  /** Adds a scalar and a tensor. */
  RS_INLINE rsTensor operator+(const double &x, const rsTensor &t)
  {
    rsTensor result(t.N, t.R, false);
    for(int i = 0; i < t.getNumComponents(); i++)
      result.c[i] = t.c[i] + x;
    return result;
  }

  /** Subtracts a tensor from a scalar. */
  RS_INLINE rsTensor operator-(const double &x, const rsTensor &t)
  {
    rsTensor result(t.N, t.R, false);
    for(int i = 0; i < t.getNumComponents(); i++)
      result.c[i] = x - t.c[i];
    return result;
  }

  /** Multiplies a scalar and a tensor. */
  RS_INLINE rsTensor operator*(const double &x, const rsTensor &t)
  {
    rsTensor result(t.N, t.R, false);
    for(int i = 0; i < t.getNumComponents(); i++)
      result.c[i] = x * t.c[i];
    return result;
  }

  // \todo: write operator to compute the tensor product between two tensors
  // \todo: write a function to contract a tensor

}

#endif
