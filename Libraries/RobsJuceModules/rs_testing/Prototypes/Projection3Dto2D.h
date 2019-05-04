#pragma once


/**

References:
1: Computer Graphics and Geometric Modeling (David Salomon)

todo: rename to rsPerspectiveProjection and make another class rsParallelProjection (maybe later 
also others like rsFishEyeProjection, etc.)
*/

template<class T>
class rsProjection3Dto2D
{

public:

  rsProjection3Dto2D();

  /** (xc,yc,zc): camera coordinates, rot: camera rotation, zoom: camera zoom setting, 
  (xt,yt,zt): target point viewed/aimed at. */
  void setup(T xc, T yc, T zc, T rot, T zoom, T xt, T yt, T zt);

  // maybe have two separate setPerspectiveProjection, setParallelProjection functions
  // and maybe some functions for special cases (standardized parallel projections - 
  // isometric, oblique, etc.)
  // maybe rename to rsPerspectiveProjection - i think, we will need the homogeneous coordinate
  // approach - the resulting equations involve a division by w ("perspective division"), so i 
  // think, they are nonlinear and can't be boiled down to a 3D affine transform - but for an 
  // orthographic/parallel projection, it may still be possible, so maybe have another class for 
  // that

  /** Applies the transformation to an input vector (x,y,z) giving an output vector (X,Y). */
  void apply(T x, T y, T z, T* X, T* Y)
  {
    *X = xx*x + xy*y + xz*z + dx;
    *Y = yx*x + yy*y + yz*z + dy;
  }

  /** Applies the transformation to an input vector (x,y,z) giving an output vector (X,Y,Z). The
  actual projection onto 2D will be obtained by just discarding Z, but this function transforms and
  retains Z, so you can use it for distance effects. */
  void apply(T x, T y, T z, T* X, T* Y, T* Z)
  {
    *X = xx*x + xy*y + xz*z + dx;
    *Y = yx*x + yy*y + yz*z + dy;
    *Z = zx*x + zy*y + zz*z + dz;
  }


  // functions for creating homogenous coordinate matrices for basic transformations:
  void translation(T M[4][4], T tx, T ty, T tz);
  void rotationX(  T M[4][4], T rx);
  void rotationY(  T M[4][4], T ry);
  void rotationZ(  T M[4][4], T rz);

  /** Composes transformation matrices A and B into C = A*B. */
  void compose(T A[4][4], T B[4][4], T C[4][4]);

protected:

  // matrix and shift coeffs (maybe factor out into rsAffineTransform3D):
  T xx, xy, xz,
    yx, yy, yz,
    zx, zy, zz,
    dx, dy, dz;
};