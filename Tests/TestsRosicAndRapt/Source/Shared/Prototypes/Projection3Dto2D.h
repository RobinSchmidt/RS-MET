#pragma once

template<class T>
class rsProjection2Dto3D
{

public:

  /** (xc,yc,zc): camera coordinates, rot: camera rotation, zoom: camera zoom setting, 
  (xv,yv,zv): point, viewed/aimed at. */
  void setup(T xc, T yc, T zc, T rot, T zoom, T xv, T yv, T zv);

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


protected:

  // matrix and shift coeffs (maybe factor out into rsAffineTransform3D):
  T xx, xy, xz,
    yx, yy, yz,
    zx, zy, zz,
    dx, dy, dz;
};