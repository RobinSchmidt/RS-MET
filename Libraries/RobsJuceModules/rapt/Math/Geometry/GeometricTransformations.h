#ifndef RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED
#define RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED

/** Collection of functions for geometric transforms. */

template<class T>
class rsGeometricTransforms
{

public:

  /** Computes the perspective projection matrix in 4D homogeneous coordinates. Corresponds to 
  vmath::frustum in OpenGL, but we use row-major indexing. see
  https://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/opengl-perspective-projection-matrix 
  */
  static void perspectiveProjection(T A[4][4], T left, T right, T bottom, T top, T near, T far);
  // not yet tested
  // rename to perspectiveFrustum - corresponds to OpenGL vmath::frustum

  /** Computes the orthographic projection matrix in 4D homogeneous coordinates.  */
  static void orthographicProjection(T A[4][4], T left, T right, T bottom, T top, T near, T far);
  // not yet tested
  // rename to orthographicFrustum - corresponds to OpenGL vmath::ortho

  // todo: see also vmath::ortho in OpenGL - lets user specify the projection in terms of 3 vectors
  // eye, center, up

  static void rotationAroundAxis(T A[3][3], T angle, T x, T y, T z);
  // not yet tested


  /** Creates a rotation matrix that rotates vector u into vector v. */
  static void rotationMatrixFromTo(rsVector3D<T> u, rsVector3D<T> v, T A[3][3]);

};
// class needs tests
// make the output matrix consistently the first or last argument - first might be better because 
// it allows to have optional input parameters


//=================================================================================================

/** A class for 2D rotations. */

template<class T>
class rsRotationXY
{

public:

  rsRotationXY(T angle = 0)
  {
    r = angle;
    updateCoeffs();
  }

  void setAngle(T newAngle)
  {
    r = newAngle;
    updateCoeffs();
  }

  /** Applies the rotation matrix to the coordinate values. */
  void apply(T* x, T* y);


protected:

  /** Updates the rotation matrix coefficients. */
  void updateCoeffs();

  T r;              // rotation angle (around z-axis)
  T xx, xy, yx, yy; // matrix coeffs

};

//=================================================================================================

/** A class for 3D rotations. You can set up the rotation angles around the 3 coordinate axes x, y 
and z and apply the rotation to a vector via the apply function. The rotations are applied in the 
order: x y z. The order is important because rotations are not commutative. */

template<class T>
class rsRotationXYZ
{

public:

  rsRotationXYZ(T angleX = 0, T angleY = 0, T angleZ = 0)
  {
    rx = angleX;
    ry = angleY;
    rz = angleZ;
    updateCoeffs();
  }

  void setAngleX(T newAngle)
  {
    rx = newAngle;
    updateCoeffs();
  }

  void setAngleY(T newAngle)
  {
    ry = newAngle;
    updateCoeffs();
  }

  void setAngleZ(T newAngle)
  {
    rz = newAngle;
    updateCoeffs();
  }

  void setAngles(T newAngleX, T newAngleY, T newAngleZ)
  {
    rx = newAngleX;
    ry = newAngleY;
    rz = newAngleZ;
    updateCoeffs();
  }

  /** Applies the rotation matrix to the coordinate values. */
  void apply(T* x, T* y, T* z);

protected:

  /** Updates the rotation matrix coefficients. */
  void updateCoeffs();

  T rx, ry, rz; // rotation angles around x, y and z-axis
  T xx, xy, xz, yx, yy, yz, zx, zy, zz; // matrix coeffs

  // 

};

#endif