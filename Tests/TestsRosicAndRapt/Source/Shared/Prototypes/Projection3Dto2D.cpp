template<class T>
void rsProjection2Dto3D<T>::setup(T xc, T yc, T zc, T rot, T zoom, T xv, T yv, T zv)
{
  // compute view plane normal:
  T nx = xc-xv;
  T ny = yc-yv;
  T nz = zc-zv;

  //...


  // idea: to derive the transformation coeffs, consider a plane which has as its normal the vector
  // that connects the camera to the viewing point. Each point is then projected onto this plane
  // ...dunno if this will work out...
  // transformation should be in such a way that the view plane normal (vpn) becomes the new 
  // z-axis

  // maybe we need to cast the problem into a sequence of affine transformations:
  // 1: translate model to the origin, apply same translation to camera (i.e. subtract aiming 
  //    point from both)
  // 2: apply a rotation that turns the view plane (VPN) normal into the z-axis and the view up
  //    vector (VUP) into the y-axis

  // ...hmm...maybe figure it out in 2D-to-1D first and get the 3D-to-2D version by analogy
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

*/