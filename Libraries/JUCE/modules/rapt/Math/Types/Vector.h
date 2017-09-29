#ifndef RAPT_VECTOR_H
#define RAPT_VECTOR_H


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

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }

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

/** Computes the dot-product between two 3D vectors v and w. */
template<class T>
T dot(const rsVector3D<T>& v, const rsVector3D<T>& w)
{
  return v.x*w.x + v.y*w.y + v.z*w.z;
}

/** Computes the cross-product between two 3D vectors v and w. Here, v is the left operand. That's 
important, because the cross-product is not commutative. Instead, we have (v*w) = -(w*v) where the
* symbol is used here to denote the cross-product. */
template<class T>
rsVector3D<T> cross(const rsVector3D<T>& v, const rsVector3D<T>& w)
{
  return rsVector3D<T>(v.y*w.z-v.z*w.y, v.z*w.x-v.x*w.z, v.x*w.y-v.y*w.x); 
}

// todo: triple(u, v, w); - compute triple-product
// maybe return references from operators - avoid copying




#endif
