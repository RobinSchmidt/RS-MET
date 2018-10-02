#ifndef rosic_Vector_h
#define rosic_Vector_h

//// rosic-indcludes:
//#include "rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a class for representing n-dimensional vectors and doing mathematical operations with 
  them. To make the class fast, no consistency checking is done in the mathematical operators. You
  must ensure yourself that the input arguments are compatible - for example don't try to add two
  vectors of different dimensionality and such.

  \todo: templatize the class, make an explicit specialization for Complex elements (norm, etc.
  has to be implemented a little different)

  */

  class Vector  
  {

  public:

    //---------------------------------------------------------------------------------------------
    // public member variables:

    /** Dimensionality. */
    int dim;   

    /** The actual values of the vector. */
    double *v;   

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Default constructor - constructs a vector with zero dimensionality and a NULL pointer. */
    Vector();
    // check

    /** Constructor. You must pass the number of elements (i.e. the dimensionality) of the vector 
    here. */
    Vector(int numElements); 
    // check

    /** Constructor. You must pass the number of elements (i.e. the dimensionality) of the vector 
    and an array with the initial values here. */
    Vector(int numElements, double *values); 
    // check

    /** Copy constructor. */
    Vector(const Vector& other); 
    // check

    /** Destructor. */
    ~Vector(); 

    //---------------------------------------------------------------------------------------------
    // overloaded operators:

    /** Assigns one vector with another one. */
    Vector& operator=(const Vector& v2)
    {
      // re-allocate memory, if necessary:
      if( dim != v2.dim )
      {
        delete[] v;
        dim = v2.dim;
        v   = new double[dim];
      }

      // copy the values:
      memcpy(v, v2.v, dim*sizeof(double));

      return *this;
    }
    // check

    /** Compares two vectors of equality. */
    bool operator==(const Vector& v2) const  
    {
      if( dim != v2.dim )
        return false;
      else
      {
        for(int i=0; i<dim; i++)
        {
          if( v[i] != v2.v[i] )
            return false;
        }
      }
      return true;
    }
    // check

    /** Compares two vectors of inequality. */
    bool operator!=(const Vector& v2) const  
    {
      return !(*this == v2);
    }
    // check

    /** Defines the negative of a vector. */
    Vector operator-()
    { 
      Vector result(dim);
      for(int i=0; i<dim; i++)
        result.v[i] = -v[i];
      return result; 
    }
    // check

    /** Adds another vector to this vector and returns the result. */
    Vector& operator+=(const Vector &v2)
    {
      rassert( dim == v2.dim );   // vectors incompatible
      for(int i=0; i<dim; i++)
        v[i] += v2.v[i];
      return *this;
    }
    // check

    /** Adds a scalar to this vector and returns the result. */
    Vector& operator+=(const double &x)
    {
      for(int i=0; i<dim; i++)
        v[i] += x;
      return *this;
    }
    // check

    /** Subtracts another vector from this vector and returns the result. */
    Vector& operator-=(const Vector &v2)
    {
      rassert( dim == v2.dim );   // vectors incopatible
      for(int i=0; i<dim; i++)
        v[i] -= v2.v[i];
      return *this;
    }
    // check

    /** Subtracts a scalar from this vector and returns the result. */
    Vector& operator-=(const double &x)
    {
      for(int i=0; i<dim; i++)
        v[i] -= x;
      return *this;
    }
    // check

    /** Multiplies this vector by a scalar and returns the result. */
    Vector& operator*=(const double &x)
    {
      for(int i=0; i<dim; i++)
        v[i] *= x;
      return *this;
    }
    // check

    /** Divides this vector by a scalar and returns the result. */
    Vector& operator/=(const double &x)
    {
      double scale = 1.0 / x;
      for(int i=0; i<dim; i++)
        v[i] *= scale;
      return *this;
    }
    // check

    /** Adds two vectors. */
    Vector operator+(const Vector &v2)
    { 
      rassert( dim == v2.dim );   // vectors incompatible
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] + v2.v[i];
      return result; 
    }
    // check

    /** Adds a vector and a scalar. */
    Vector operator+(const double &x)
    { 
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] + x;
      return result; 
    }
    // check

    /** Subtracts two vectors. */
    Vector operator-(const Vector &v2)
    { 
      rassert( dim == v2.dim );   // vectors incopatible
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] - v2.v[i];
      return result; 
    }
    // check

    /** Subtracts a scalar from a vector. */
    Vector operator-(const double &x)
    { 
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] - x;
      return result; 
    }
    // check

    /** Multiplies a vector and a scalar. */
    Vector operator*(const double &x)
    { 
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] * x;
      return result; 
    }
    // check

    /** Divides a vector by a scalar. */
    Vector operator/(const double &x)  
    {
      double scale = 1.0 / x;
      Vector result(dim);
      for(int i=0; i<dim; i++)    
        result.v[i] = v[i] * scale;
      return result; 
    }
    // check

    //---------------------------------------------------------------------------------------------
    // manipulations:

    // removeElement, addElement, setDimension, normalize, ...

    /** Sets a new dimensionality (number of elements) for this vector. */
    void setDimensionality(int newDimensionality)
    {
      // re-allocate memory, if necessary:
      if( dim != newDimensionality )
      {
        delete[] v;
        dim = newDimensionality;
        v   = new double[dim];
      }
    }

    /** Applies the passed function to each element of the vector. */
    void applyFunction( double (*f) (double) )
    {
      for(int i=0; i<dim; i++)    
        v[i] = f(v[i]);
    }
    // check

    /** Assigns random values between 'min' and 'max' to each element. */
    void randomizeElements(double min, double max)
    {
      for(int i=0; i<dim; i++)    
        v[i] = RAPT::rsRandomUniform(min, max);
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

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the squared Euclidean norm of the vector. This is the sum of the squares of the 
    elements: x1^2 + x2^2 + x3^2 + .... */
    double getSquaredNorm()
    {
      double accu = 0.0;
      for(int i=0; i<dim; i++)    
        accu += v[i] * v[i];
      return accu;
    }

    /** Returns the Euclidean norm of the vector defined as sqrt(x1^2 + x2^2 + x3^2 + ...). */
    double getEuclideanNorm() { return sqrt( getSquaredNorm() ); }

    //---------------------------------------------------------------------------------------------
    // others:

    /** Stores the values into the target array which has to be of length matching the 
    dimensionality of this vector. */
    void storeValuesInArray(double *targetArray) { memcpy(targetArray, v, dim*sizeof(double)); }

    /** Prints the values to the standard output - mainly for debugging purposes. */
    void print();

  }; // end of class Vector

  // some binary operators are defined outside the class such that the left hand operand does 
  // not necesarrily need to be of class Vector

  /** Adds a scalar and a vector. */
  inline Vector operator+(const double &x, const Vector &v)
  { 
    Vector result(v.dim);
    for(int i=0; i<v.dim; i++)    
      result.v[i] = v.v[i] + x;
    return result; 
  }
  // check

  /** Subtracts a vector from a scalar. */
  inline Vector operator-(const double &x, const Vector &v)
  { 
    Vector result(v.dim);
    for(int i=0; i<v.dim; i++)    
      result.v[i] = x - v.v[i];
    return result; 
  }
  // check

  /** Multiplies a scalar and a vector. */
  inline Vector operator*(const double &x, const Vector &v)
  { 
    Vector result(v.dim);
    for(int i=0; i<v.dim; i++)    
      result.v[i] = x * v.v[i];
    return result; 
  }
  // check

  /** Computes the inner product between two vectors (which should have the same dimensionality) - 
  the result is a scalar. */
  inline double operator*(const Vector &a, const Vector &b)
  {
    rassert( a.dim == b.dim );   // vectors incopatible
    double result = 0.0;
    for(int i=0; i<a.dim; i++)
      result += a.v[i] * b.v[i];
    return result;
  }
  // check

}  // end namespace rosic

#endif // rosic_Vector_h
