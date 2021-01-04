
template<class T>
void rsGeometricTransforms<T>::perspectiveProjection(T A[4][4], T l, T r, T b, T t, T n, T f)
{
  // todo: precompute 1/(r-l), 1/(t-b), 1/(f-n) -> replace 3 divisions by multiplications

  // compute the 4 columns of the matrix (we order the assignments by last index because that 
  // better reveals the structure of the matrix)
  A[0][0] = 2*n/(r-l);
  A[1][0] = 0;
  A[2][0] = 0;
  A[3][0] = 0;

  A[0][1] = 0;
  A[1][1] = (2*n)/(t-b);
  A[2][1] = 0;
  A[3][1] = 0;

  A[0][2] =  (r+l)/(r-l);
  A[1][2] =  (t+b)/(t-b);
  A[2][2] = -(f+n)/(f-n);
  A[3][2] = -1;

  A[0][3] = 0;
  A[1][3] = 0;
  A[2][3] = -(2*f*n)/(f-n);
  A[3][3] = 0;
}
// not yet tested
// formula from OpenGL Programming Guide, 9th Edition, page 853 (Appendix E), i think, it 
// corresponds to vmath::frustum in OpenGL (see page 219)

// also for the function below:

template<class T>
void rsGeometricTransforms<T>::orthographicProjection(T A[4][4], T l, T r, T b, T t, T n, T f)
{
  A[0][0] = 2/(r-l);
  A[1][0] = 0;
  A[2][0] = 0;
  A[3][0] = 0;

  A[0][1] = 0;
  A[1][1] = 2/(t-b);
  A[2][1] = 0;
  A[3][1] = 0;

  A[0][2] = 0;
  A[1][2] = 0;
  A[2][2] = -2/(f-n);
  A[3][2] = 0;

  A[0][3] = -(r+l)/(r-l);
  A[1][3] = -(t+b)/(t-b);
  A[2][3] = -(f+n)/(f-n);
  A[3][3] = 1;
}
// not yet tested
// formula from OpenGL Programming Guide, 9th Edition, page 853 (Appendix E)

template<class T>
void rsGeometricTransforms<T>::rotationAroundAxis(T A[3][3], T a, T x, T y, T z)
{
  T xx, xy, xz, yy, yz, zz, k;

  // naive (1 sqrt, 1 div, 12 mul, 2 add):
  k = 1 / sqrt(x*x + y*y + z*z);
  x *= k;
  y *= k;
  z *= k;
  xx = x*x;
  xy = x*y;
  xz = x*z;
  yy = y*y;
  yz = y*z;
  zz = z*z;

  /*
  // optimized (no sqrt, 1 div, 12 mul, 2 add);
  xx  = x*x;
  yy  = y*y;
  zz  = z*z;
  k   = T(1) / (xx + yy + zz);
  xx *= k;
  yy *= k;
  zz *= k;
  xy  = k*x*y;
  xz  = k*x*z;
  yz  = k*y*z;
  */

  T s = sin(a);
  T c = cos(a);

  A[0][0] = -(xx - 1)*c + xx;
  A[0][1] = -c*x*y + x*y - s*z;
  A[0][2] = -c*x*z + s*y + x*z;

  A[1][0] = -c*x*y + x*y + s*z;
  A[1][1] = -(yy - 1)*c + yy;
  A[1][2] = -c*y*z - s*x + y*z;

  A[2][0] = -c*x*z - s*y + x*z;
  A[2][1] = -c*y*z + s*x + y*z;
  A[2][2] = -(zz - 1)*c + zz;
}
// not yet tested

/*
formula derived from OpenGL Programming Guide, page 852 with the sage code:
var("x y z c s")
u  = matrix([[x], [y], [z]]) # u is column matrix (i.e. vector)
u2 = u * u.transpose()       # u2 is 3x3 matrix u * u^T
S  = matrix([[0,-z,y], [z,0,-x], [-y,x,0]])
I  = matrix.identity(3)
M  = u2 + c*(I-u2) + s*S
M
*/

template<class T>
void rsGeometricTransforms<T>::rotationMatrixFromTo(rsVector3D<T> u, rsVector3D<T> v, T A[3][3])
{
  u.normalize();
  v.normalize();
  rsVector3D<T> a = cross(v, u); // exchanged arguments with respect to Weitz's code
  T c  = dot(u, v);              // cos(u,v)
  T s  = rsNorm(a);              // sin(u,v)
  T c1 = T(1)-c;                 // 1-cos(u,v)
  a.normalize();

  A[0][0] = a.x*a.x*c1 + c;
  A[0][1] = a.x*a.y*c1 + s*a.z;
  A[0][2] = a.x*a.z*c1 - s*a.y;

  A[1][0] = a.x*a.y*c1 - s*a.z;
  A[1][1] = a.y*a.y*c1 + c;
  A[1][2] = a.y*a.z*c1 + s*a.x;

  A[2][0] = a.x*a.z*c1 + s*a.y;
  A[2][1] = a.y*a.z*c1 - s*a.x;
  A[2][2] = a.z*a.z*c1 + c; 
}
// has been tested to rotate every canonical basis vector (1,0,0),(0,1,0),(0,0,1) into every other
// -> works as it should
// Adapted from Weitz - Differentialgeometrie, page 103, but i had to change the argument order in 
// the call to the cross-product with respect to Weitz's code - taking the code as is would have 
// resulted in an opposite rotation - apparently some different conventions are in use (todo: 
// figure out what exactly is going on - maybe it's a convention about how transformations work in 
// applyMatrix in p5.js?)



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
 see also: https://www.youtube.com/watch?v=PDgG2Z6T1ho for more on 3D rotations

*/