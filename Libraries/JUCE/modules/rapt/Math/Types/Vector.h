#ifndef RAPT_VECTOR_H
#define RAPT_VECTOR_H


/** Class for representing 3-dimensional vectors. */

template<class T>
class rsVector3D
{

public:

  /** Constructor. Initializes coordinates with the passed values. */
  rsVector3D(T _x = 0, T _y = 0, T _z = 0) : x(_x), y(_y), z(_z) {}

  /** Returns the squared Euclidean norm of this vector. */
  T getSquaredEuclideanNorm() { return x*x + y*y + z*z; }

  /** Returns the Euclidean norm of this vector. */
  T getEuclideanNorm() { return sqrt(getSquaredEuclideanNorm()); }

  // we need operators +,-,+=,-= that accept vectors and scalars as right operands and
  // *,/,*=,/= that accept scalars. we also need a * operator that accepts scalars on the
  // left


  /** The 3 cartesian coordinate values. */
  T x, y, z;

};

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





#endif
