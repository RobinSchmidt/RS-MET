
template<class T>
void rsGeometricTransforms<T>::perspectiveProjection(T* A[4][4], T l, T r, T b, T t, T n, T f)
{
  // todo: precompute 1/(r-l), 1/(t-b), 1/(f-n) -> replace 3 divisions by multiplications

  A[0][0] = 2*n/(r-l);
  A[0][1] = 0;
  A[0][2] = (r+l)/(r-l);
  A[0][3] = 0;

  A[1][0] = 0;
  A[1][1] = (2*n)/(t-b);
  A[1][2] = (t+b)/(t-b);
  A[1][4] = 0;

  A[2][0] = 0;
  A[2][1] = 0;
  A[2][2] = -(f+n)/(f-n);
  A[2][4] = -(2*f*n)/(f-n);

  A[3][0] =  0;
  A[3][1] =  0;
  A[3][2] = -1;
  A[3][4] =  0;
}

//=================================================================================================

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
// maybe inline these two functions

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

/*
ToDo:
-make a rotation class that lets the user set up the rotation axis and angle
 (see OpenGL Programming Guide, p. 852 for formulas - maybe make it possible to decompose the 
 rotation into an x-, y-, z-roation (by successively dividing out the respective inverse
 rotation-matrices, recover also the rotation angles (acos)

*/