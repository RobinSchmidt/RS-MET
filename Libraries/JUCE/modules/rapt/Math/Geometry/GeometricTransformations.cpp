template<class T>
void rsRotationXY<T>::apply(T* x, T* y)
{
  // temporaries:
  T X = *x; 
  T Y = *y; 

  // new vector is given by matrix-vector product:
  *x = xx*X + xy*Y;    // |x| = |xx xy| * |X|
  *y = yx*X + yy*Y;    // |y|   |yx yy|   |Y|
}

template<class T>
void rsRotationXY<T>::updateCoeffs()
{
  T s = sin(r); 
  T c = cos(r);
  xx =  c;
  xy = -s;
  yx =  s;
  yy =  c;
}

//=================================================================================================

template<class T>
void rsRotationXYZ<T>::apply(T* x, T* y, T* z)
{
  // temporaries:
  T X = *x; 
  T Y = *y; 
  T Z = *z;

  // new vector is given by matrix-vector product:
  *x = xx*X + xy*Y + xz*Z;    // |x|   |xx xy xz|   |X|
  *y = yx*X + yy*Y + yz*Z;    // |y| = |yx yy yz| * |Y|
  *z = zx*X + zy*Y + zz*Z;    // |z|   |zx zy zz|   |Z|
}

template<class T>
void rsRotationXYZ<T>::updateCoeffs()
{
  // sines/cosines:
  T sx = sin(rx); T cx = cos(rx);
  T sy = sin(ry); T cy = cos(ry);
  T sz = sin(rz); T cz = cos(rz);

  // rotation matrix coeffs:
  xx =  cz*cy;
  xy = -sz*cx + cz*sy*sx;
  xz =  sz*sx + cz*sy*cx;
  yx =  sz*cy;
  yy =  cz*cx + sz*sy*sx;
  yz = -cz*sx + sz*sy*cx;
  zx = -sy;
  zy =  cy*sx;
  zz =  cy*cx;
}

// for a rotation around the z-axis first, then the y-axis then the x axis, the matrix would be:
// xx =  cy*cz; 
// xy = -cy*sz; 
// xz =  sy;
// yx =  sx*sy*cz + cx*sz;
// yy = -sx*sy*sz + cx*cz;
// yz = -sx*cy;
// zx = -cx*sy*cz + sx*sz;
// zy =  cx*sy*sz + sx*cz;
// zz =  cx*cy;
// maybe we can use them for a reverse rotation (we would have to set the angles to their negative
// values too).