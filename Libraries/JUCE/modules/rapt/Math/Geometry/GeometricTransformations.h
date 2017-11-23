#ifndef RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED
#define RAPT_GEOMETRICTRANSFORMATIONS_H_INCLUDED

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

  /** Applies the rotation matrix to the coordinate values. */
  void apply(T* x, T* y, T* z);

protected:

  /** Updates the rotation matrix coefficients. */
  void updateCoeffs();

  T rx, ry, rz; // rotation angles around x, y and z-axis
  T xx, xy, xz, yx, yy, yz, zx, zy, zz; // matrix coeffs

};

#endif