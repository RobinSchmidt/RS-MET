#ifndef RAPT_VECTOR_H
#define RAPT_VECTOR_H

/** Class for representing 2-dimensional vectors. */

template<class T>
class rsVector2D
{

public:

  /** The 2 cartesian coordinate values. */
  T x, y;

  /** Constructor. Initializes coordinates with the passed values. */
  rsVector2D(T _x = 0, T _y = 0) : x(_x), y(_y) {}
  // for optimization, make a constructor without initialization

  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y; }
  // rename to squaredNorm

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }

  //-----------------------------------------------------------------------------------------------
  /** \name Operators */

  /** Adds vector v to this vector. */
  rsVector2D<T>& operator+=(const rsVector2D<T>& v) { x += v.x; y += v.y; return *this; }

  /** Subtracts vector v from this vector. */
  rsVector2D<T>& operator-=(const rsVector2D<T>& v) { x -= v.x; y -= v.y; return *this; }

  /** Multiplies this vector by scalar s. */
  rsVector2D<T>& operator*=(const T& s) { x *= s; y *= s; return *this; }

  /** Divides this vector by scalar s. */
  rsVector2D<T>& operator/=(T s) { s = 1/s; x *= s; y *= s; return *this; }

};


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
  rsVector3D(T _x = 0, T _y = 0, T _z = 0) : x(_x), y(_y), z(_z) {}
    // for optimization, make a constructor without initialization

  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y + z*z; }
    // rename to squaredNorm

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
Its a scalar that gives the volume V of parallelepiped spanned by the 3 vectors. */
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




#endif
