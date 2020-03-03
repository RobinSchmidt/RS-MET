#ifndef RAPT_VECTOR_H
#define RAPT_VECTOR_H

/** Class for representing 2-dimensional vectors. */

template<class T>
class rsVector2D  // maybe it should be a struct?
{

public:

  /** The 2 cartesian coordinate values. */
  T x, y;

  /** Constructor. Initializes coordinates with the passed values. */
  rsVector2D(T _x = 0, T _y = 0) : x(_x), y(_y) {}
  // for optimization, make a constructor without initialization

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y; }
  // rename to squaredNorm or getSquaredLength or getSquaredNorm

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }
  // rename to getLength or getNorm

  /** Tests, if another vector v is close to this vector within a given tolerance (both components
  of the difference must be <= tolerance). */
  bool isCloseTo(const rsVector2D<T>& v, T tol)
  {
    if(rsAbs(v.y - y) <= tol && rsAbs(v.x - x) <= tol)
      return true;
    return false;
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Manipulations */

  /** Normalizes this vector to unit length.  */
  void normalize()
  {
    T s = T(1) / getEuclideanNorm();
    x *= s;
    y *= s;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds two vectors. */
  rsVector2D<T> operator+(const rsVector2D<T> &v) const { return rsVector2D<T>(x+v.x, y+v.y); }

  /** Subtracts two vectors. */
  rsVector2D<T> operator-(const rsVector2D<T> &v) const { return rsVector2D<T>(x-v.x, y-v.y); }

  /** Multiplies a vector by a number. */
  rsVector2D<T> operator*(const T &r) const { return rsVector2D<T>(x*r, y*r); }

  /** Divides a vector by a number. */
  rsVector2D<T> operator/(const T &r) const { T s = T(1) / r; return rsVector2D<T>(x*s, y*s); }

  /** Adds vector v to this vector. */
  rsVector2D<T>& operator+=(const rsVector2D<T>& v) { x += v.x; y += v.y; return *this; }

  /** Subtracts vector v from this vector. */
  rsVector2D<T>& operator-=(const rsVector2D<T>& v) { x -= v.x; y -= v.y; return *this; }

  /** Multiplies this vector by scalar s. */
  rsVector2D<T>& operator*=(const T& s) { x *= s; y *= s; return *this; }

  /** Divides this vector by scalar s. */
  rsVector2D<T>& operator/=(T s) { s = 1/s; x *= s; y *= s; return *this; }

  /** Defines the negative of a vector. */
  rsVector2D<T> operator-() const { return rsVector2D(-x, -y); }

  /** Compares two vectors for equality. */
  bool operator==(const rsVector2D<T>& p) const 
  {
    if(x == p.x && y == p.y)
      return true;
    else
      return false;
  }

  /** Compares two vectors for inequality. */
  bool operator!=(const rsVector2D<T>& p) const
  {
    if(x != p.x || y != p.y)
      return true;
    else
      return false;
  }

  //-----------------------------------------------------------------------------------------------

  /** \name Static Member Functions */

  /** Returns the cross-product of two vectors defined as x1*y2-x2*y1. The cross-product can be
  interpreted as the signed area of the parallelogram formed by the vectors (0,0),p1,p2,p1+p2
  where the sign is positive when p1 is clockwise from p2 with respect to the origin, negative
  when it's counterclockwise and zero when they are collinear (pointing in the same or opposite
  directions).  */
  static T crossProduct(const rsVector2D<T> &p1, const rsVector2D<T> &p2)
  {
    return p1.x*p2.y - p2.x*p1.y;
  }
  // rename to cross

  /** Returns the dot-product (aka inner product or scalar product) of two vectors defined as
  x1*x2+y1*y2. When the two vectors are of unit length, it can be interpreted as the cosine of
  the angle between the two vectors. Because the cosine of a 90 degree angle is zero, this
  implies that a zero scalar-product indicates orthogonal vectors. If only p1=(x1,y1) is a unit
  vector, then the scalar-product can be interpreted as the length of the orthogonal projection
  of p2 onto p1. */
  static T scalarProduct(const rsVector2D<T> &p1, const rsVector2D<T> &p2)
  {
    return p1.x*p2.x + p1.y*p2.y;
  }
  // rename to dot

};

/** Multiplies a number and a vector. We need to define this operator outside the class because the 
left operand is not of class rsVector2D (but of some real number type). We use the fact that 
number*vector == vector*number to compute the result.  */
template<class T>
inline rsVector2D<T> operator*(const T &r, const rsVector2D<T> &p)
{
  rsVector2D<T> tmp(p);
  tmp *= r;
  return tmp;
}

/** Returns the determinant of the matrix that results from writing the 2 given vectors as columns
into a 2x2 matrix. If this determinant is 0, the 2 vectors are collinear, i.e. one is a scalar 
multiple of the other, i.e. they are both on the same line through the origin. */
// not yet tested
template<class T>
T rsDet(const rsVector2D<T>& a, const rsVector2D<T>& b)
{
  return a.x * b.y - a.y * b.x;  // verify formula
}
// maybe move as static function into class

/** Returns the dot product of vectors a and b. */
template<class T>
T rsDot(const rsVector2D<T>& a, const rsVector2D<T>& b)
{
  return a.x * b.x + a.y * b.y;
}

/** Returns the Euclidean norm of vector a. */
template<class T>
T rsNorm(const rsVector2D<T>& a)
{
  return sqrt(rsDot(a, a));
}



// todo: dot, cross

//=================================================================================================

/** Class for representing 3-dimensional vectors. 

\todo:
 -complete the set of operators
*/

template<class T>
class rsVector3D
{

public:

  /** The 3 cartesian coordinate values. */
  T x, y, z;

  /** Constructor. Initializes coordinates with the passed values. */
  rsVector3D(T _x, T _y, T _z) : x(_x), y(_y), z(_z) {}

  /** Standard constructor. Leaves components uninitialized. */
  rsVector3D() {}



  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y + z*z; }
    // rename to squaredNorm or getSquaredNorm

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }

  // maybe just define "norm", "normSquared" functions outside the class (like dot, det, etc)

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds vector v to this vector. */
  rsVector3D<T>& operator+=(const rsVector3D<T>& v) { x += v.x; y += v.y; z += v.z; return *this; }

  /** Subtracts vector v from this vector. */
  rsVector3D<T>& operator-=(const rsVector3D<T>& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }

  /** Multiplies this vector by scalar s. */
  rsVector3D<T>& operator*=(const T& s) { x *= s; y *= s; z *= s; return *this; }

  /** Divides this vector by scalar s. */
  rsVector3D<T>& operator/=(T s) { s = 1/s; x *= s; y *= s; z *= s; return *this; }


  // implement dot- and

};

/** Adds two vectors and returns result as new vector. */
template<class T>
rsVector3D<T> operator+(const rsVector3D<T>& v, const rsVector3D<T>& w)
{
  return rsVector3D<T>(v.x+w.x, v.y+w.y, v.z+w.z);
}
// todo: v+s, s+v (vector+scalar, scalar+vector)

/** Subtracts two vectors and returns result as new vector. */
template<class T>
rsVector3D<T> operator-(const rsVector3D<T>& v, const rsVector3D<T>& w)
{
  return rsVector3D<T>(v.x-w.x, v.y-w.y, v.z-w.z);
}
// todo: v-s, s-v

/** Multiplies a vector and a scalar and returns result as new vector. */
template<class T>
rsVector3D<T> operator*(const rsVector3D<T>& v, const T& s)
{
  return rsVector3D<T>(v.x*s, v.y*s, v.z*s);
}

/** Multiplies a scalar and a vector and returns result as new vector. */
template<class T>
rsVector3D<T> operator*(const T& s, const rsVector3D<T>& v)
{
  return rsVector3D<T>(v.x*s, v.y*s, v.z*s);
}

/** Divides a vector by a scalar and returns result as new vector. */
template<class T>
rsVector3D<T> operator/(const rsVector3D<T>& v, T s)
{
  s = 1/s;
  return rsVector3D<T>(v.x*s, v.y*s, v.z*s);
}

/** Computes the dot-product between two 3D vectors a and b. */
template<class T>
T dot(const rsVector3D<T>& a, const rsVector3D<T>& b)
{
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

/** Computes the cross-product between two 3D vectors a and b. Here, a is the left operand. That's 
important, because the cross-product is not commutative. Instead, we have (a x b) = -(b x a) where 
the x symbol is used here to denote the cross-product. */
template<class T>
rsVector3D<T> cross(const rsVector3D<T>& a, const rsVector3D<T>& b)
{
  return rsVector3D<T>(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x); 
}

/** Returns the triple product of the 3 vectors a,b,c, defined as: 
V = (a x b) * c = a * (b x c) = (b x c) * a = (c x a) * b
Its a scalar that gives the volume V of the parallelepiped spanned by the 3 vectors. */
// not yet tested
template<class T>
rsVector3D<T> triple(const rsVector3D<T>& a, const rsVector3D<T>& b, const rsVector3D<T>& c)
{
  dot(cross(a, b), c); 
}

/** Returns the determinant of the matrix that results from writing the 3 given vectors as columns
into a 3x3 matrix. If this determinant is 0, the 3 vectors are linearly dependent, i.e. one can be 
expressed in terms of the other two, i.e. they are all in the same plane. */
// not yet tested
template<class T>
T det(const rsVector3D<T>& a, const rsVector3D<T>& b, const rsVector3D<T>& c)
{
  return a.x*b.y*c.z + b.x*c.y*a.z + c.x*a.y*b.z - c.x*b.y*a.z - b.x*a.y*c.z - a.x*c.y*b.z;
}

/** Returns the angle between vectors a and b. */
// not yet tested
template<class T>
T angle(const rsVector3D<T>& a, const rsVector3D<T>& b)
{
  return acos(dot(a, b) / sqrt(a.getSquaredEuclideanNorm() * b.getSquaredEuclideanNorm()));
}

// maybe return references from operators - avoid copying

/*
Vector identities:
(a x b) x c = (ac)b - (bc)a  "double cross"
b = (1/a^2) * (ab)a + (1/a^2) (a x b) x a = parallel + normal component of b with respect to a


*/

//=================================================================================================

/*
// maybe, we don't really need a general n-dimensional rsVector class - we may represent vectors as 
column matrices ...or: maybe have rsRowVector, rsColumnVector classes as subclasses (special cases) 
of rsMatrix - on the other hand, it may be convenient to have a vector class for doing 
matrix-vector computations...but then we need a lot of multiplication operators: 
matrix * columnVector = columnVector, rowVector * matrix = rowVector, 
rowVector * columnVector = scalar, columnVector * rowVector = matrix - and these are actually all 
already subsumed under a general matrix * matrix = matrix multiplication as special cases..hmmm...
// ...or maybe call the rsRowMatrix, rsColumnMatrix
template<class T>
class rsVectorView  // rename to rsVector and the current rsVector class to rsVectorOld
{

public:

  int N;   // dimensionality
  T*  v;   // the actual data

protected:



};

template<class T>
class rsVectorNew  // rename to rsVector and the current rsVector class to rsVectorOld
{

public:

protected:

};
*/


#endif
