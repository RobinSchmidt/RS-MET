template<class T>
rsProjection3Dto2D<T>::rsProjection3Dto2D()
{
  setup(0, 0, -1, 0, 1, 0, 0, 0); 
  // camera at (0,0,-1) looking at (0,0,0) with no rotation and unit zoom (standard position)

  // maybe late setup a parallel projection by default
}

template<class T>
void rsProjection3Dto2D<T>::setup(T xc, T yc, T zc, T rot, T zoom, T xt, T yt, T zt)
{
  // compute viewing direction (aka view plane normal, VPN):
  T vx = xc-xt;
  T vy = yc-yt;
  T vz = zc-zt;

  // normalize view direction vector:
  T s = 1 / sqrt(vx*vx + vy*vy + vz*vz); // scaler
  vx *= s;
  vy *= s;
  vz *= s;


  //T T1[4][4], 

  // create homogenous matrix via formula (1), Eq 3.44 (factor out):
  T a = xc, b = yc, c = zc;
  T d = vx, e = vy, f = vz;
  T r = 1; // preliminary - later 1/k where k is the distance of the viewer from the screen
  s = 1 / (1+f);  // later: catch special case where f = -1 (div-by-zero)
  T Tg[4][4];
  Tg[0][0] = s*(e*e + f + f*f);  // 1st row
  Tg[0][1] = s*(-d*e);
  Tg[0][2] = 0;
  Tg[0][3] = d*r;
  Tg[1][0] = s*(-d*e);           // 2nd row
  Tg[1][1] = s*(d*d + f + f*f);
  Tg[1][2] = 0;
  Tg[1][3] = e*r;
  Tg[2][0] = -d;                 // 3rd row
  Tg[2][1] = -e;
  Tg[2][2] = 0;
  Tg[2][3] = f*r;
  Tg[3][0] = s*(c*d + b*d*e - a*e*e - a*f + c*d*f - a*f*f);  // 4th row
  Tg[3][1] = s*(-b*d*d + c*e + a*d*e - b*f + c*e*f - b*f*f);
  Tg[3][2] = 0;
  Tg[3][3] = -(a*d + b*e + c*f)*r;
  // actually, it would be nice to absorb the final rotation as well
  // damn - this matrix sets the z-coordinate to 0 because it includes the final projection matrix 
  // it would be better to have one that pre

  // maybe make a class rsAffineTransform3D
  // it should be able to create the unique rotation that rotates one vector into (the direction 
  // of) another vector which can then be used to rotate viewing direction into (0,0,1) (catch 
  // singular case of original vector being (0,0,-1))
  // -compose a translation of viewer to the origin with this rotation (or maybe use simpler 
  //  formulas taking into account the simplicity of the target vector)
  // -a perspective projection can use the same transform as a parallel projection and just apply
  //  a diffent final step (instead of X = x, Y = y, do: X = x/(1+z/k), Y = y/(1+z/k), Eq. 33.8)
  // -to preserve depth in a perspective projection, use Z / z/(1+z/k) - this is proprtionla to 
  //  what is given on page 137)
  // -maybe it's best to re-derive 3.44 but without the final projection step, i.e. only the matrix
  //  product of the translation*rotation*translation (maybe even the last translation (by k)) is 
  //  not necessarry - k is irrelevant in parallel projections and for perspective, it can be 
  //  included in the final step)
  //  -so, (re)-derive a formula for T1*T2 (maybe using column-vectors and pre-multiplication) and
  //   for parallel projection, just strip off z and for perspective use X = x/(1+z/k), etc.
  //  -the parallel projection can (i think) indeed be expressed a 3D affine transform - so the
  //   perspective can be represented by that same affine transform followed by multiplying the
  //   end results by 1/(1+z/k) -> all of that can be done in one class with functions
  //   projectParallel(x,y,z, X,Y,Z); projectPerspective(x,y,z, X,Y,Z) where the latter calls the 
  //   former and multiplies results by 1/(1+z/k) (note how that factor goes to one when k goes to 
  //   inf as it should). of course, we can have versions opf the "project" functions that don't 
  //   compute z when it's not needed  - for optimzation

  // now we have the 4D homogeneous matrix Tg - next step: translate to 3D affine transform








  //...


  // idea: to derive the transformation coeffs, consider a plane which has as its normal the vector
  // that connects the camera to the viewing point. Each point is then projected onto this plane
  // ...dunno if this will work out...
  // transformation should be in such a way that the view plane normal (vpn) becomes the new 
  // z-axis

  // maybe we need to cast the problem into a sequence of affine transformations:
  // 1: translate camera to the origin
  // 2: apply a rotation that turns the viewing direction into the z-axis
  // 3: 
  // and the view up
  //    vector (VUP) into the y-axis

  // ...hmm...maybe figure it out in 2D-to-1D first and get the 3D-to-2D version by analogy

  // in the transformations and projections, page 109, there is a matrix given in 4D homogeneous
  // coordinates, see also some comments on 111-112 for a special case that needs to be catched
}


/*
terminology from Computer Graphics (Foley, ...)
VRP: view reference point (on the view plane)
VPN: view plane normal
VRC: viewing-reference coordinates (origin at VRP, n-axis/forward-axis is given by the VPN)
VUP: view up vector (v-axis, upward axis in VRC), u-axis is defined such that u,v,n for a 
     right-handed system
CW:  center of window, sometimes identical to VRP but not necessarily
WC:  window coordinates
COP: center of projection
DOP: direction of projection
PRP: projection reference point (specified in the VRC system)
NPC: normalized projection coordinates


VPD: view plane distance (from VRP? but this is supposed to be in the plane?)

References:
http://www.scratchapixel.com/lessons/mathematics-physics-for-computer-graphics/lookat-function

https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-generating-camera-rays/generating-camera-rays

http://www.scratchapixel.com/lessons/3d-basic-rendering/perspective-and-orthographic-projection-matrix/opengl-perspective-projection-matrix



*/