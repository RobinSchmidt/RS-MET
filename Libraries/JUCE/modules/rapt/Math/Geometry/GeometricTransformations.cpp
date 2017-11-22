template<class T>
void rsRotateXYZ(T* x, T* y, T* z, T rx, T ry, T rz)
{
  // temporaries:
  T X = *x; 
  T Y = *y; 
  T Z = *z;

  // sines/cosines:
  T sx = sin(rx); T cx = cos(rx);
  T sy = sin(ry); T cy = cos(ry);
  T sz = sin(rz); T cz = cos(rz);

  /*
  // old rotation matrix coeffs (i think, they are actually for a ZYX-rotation):
  T xx =  cy*cz; 
  T xy = -cy*sz; 
  T xz =  sy;
  T yx =  sx*sy*cz + cx*sz;
  T yy = -sx*sy*sz + cx*cz;
  T yz = -sx*cy;
  T zx = -cx*sy*cz + sx*sz;
  T zy =  cx*sy*sz + sx*cz;
  T zz =  cx*cy;
  // maybe we can use them for a reverse rotation (we would have to set the angles to teit negative
  // values too).
  */

  // rotation matrix coeffs:
  T xx =  cz*cy;
  T xy = -sz*cx + cz*sy*sx;
  T xz =  sz*sx + cz*sy*cx;
  T yx =  sz*cy;
  T yy =  cz*cx + sz*sy*sx;
  T yz = -cz*sx + sz*sy*cx;
  T zx = -sy;
  T zy =  cy*sx;
  T zz =  cy*cx;

  // new vector is given by matrix-vector product:
  *x = xx*X + xy*Y + xz*Z;    // |x|   |xx xy xz|   |X|
  *y = yx*X + yy*Y + yz*Z;    // |y| = |yx yy yz| * |Y|
  *z = zx*X + zy*Y + zz*Z;    // |z|   |zx zy zz|   |Z|
}