template<class T>
void rsRotationOscillator<T>::processSampleFrame(T* x, T* y, T* z)
{
  // update matrices, if necessarry:
  if(oscMatrixNeedsUpdate)
    updateOscillationMatrix();
  if(trafoMatrixNeedsUpdate)
    updateTransformMatrix();

  // obtain new vector on unit sphere:
  oscillateRotation.apply(&X, &Y, &Z);

  // transform to ellipsoid:
  T tx = (X + shiftX) * scaleX;
  T ty = (Y + shiftY) * scaleY;
  T tz = (Z + shiftZ) * scaleZ;

  // rotate the ellipsoid:
  ellipsoidRotation.apply(&tx, &ty, &tz);

  // renormalize length:
  T s = 1 / sqrt(tx*tx + ty*ty + tz*tz);
  tx *= s;
  ty *= s;
  tz *= s;

  // apply the output rotation and assign outputs:
  //outputRotation.apply(&tx, &ty, &tz);
  *x = tx;
  *y = ty;
  *z = tz;


  // todo: maybe apply decay, maybe inject and input signal in state update -> turns it into a 
  // filter
}

template<class T>
void rsRotationOscillator<T>::updateOscillationMatrix()
{
  T w  = 2*T(PI)*frequency/sampleRate;
  oscillateRotation.setAngles(w*freqScaleX, w*freqScaleY, w*freqScaleZ);
}

template<class T>
void rsRotationOscillator<T>::updateTransformMatrix()
{
  transformRotation.setAngles(rotX, rotY, rotZ);
}