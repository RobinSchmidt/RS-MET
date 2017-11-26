template<class T>
void rsRotationOscillator<T>::processSampleFrame(T* x, T* y, T* z)
{
  // hmm...no - this is boring - it just rotates with one frequency
  // instead we should use 3 sine-oscillators with freq and phase
  // ...make a 3D Lissajous oscillator - it reduces to xoxos osc when
  // wz = 0, wx = wy, px = 0°, py = 90°

  // update matrices, if necessarry:
  if(oscMatrixNeedsUpdate)
    updateOscillationMatrix();
  if(trafoMatrixNeedsUpdate)
    updateTransformMatrix();
  if(outputMatrixNeedsUpdate)
    updateOutputMatrix();

  // obtain new vector on unit sphere:
  oscillateRotation.apply(&X, &Y, &Z);

  // transform to ellipsoid:
  T tx = (X + shiftX) * scaleX;
  T ty = (Y + shiftY) * scaleY;
  T tz = (Z + shiftZ) * scaleZ;

  // rotate the ellipsoid:
  transformRotation.apply(&tx, &ty, &tz);

  // renormalize length:
  T s = 1 / sqrt(tx*tx + ty*ty + tz*tz); // what about div-by-zero?
  tx *= s;
  ty *= s;
  tz *= s;

  //...maybe apply spherical (soft) clipping?

  // apply the output rotation and assign outputs:
  outputRotation.apply(&tx, &ty, &tz);
  *x = tx;
  *y = ty;
  *z = tz;


  // todo: maybe apply decay, maybe inject and input signal in state update -> turns it into a 
  // filter
}

template<class T>
void rsRotationOscillator<T>::reset()
{
  X = 1;
  Y = 0;
  Z = 0;
}

template<class T>
void rsRotationOscillator<T>::updateOscillationMatrix()
{
  T w  = 2*T(PI)*freq/sampleRate;
  oscillateRotation.setAngles(w*freqScaleX, w*freqScaleY, w*freqScaleZ);
  oscMatrixNeedsUpdate = false;
}

template<class T>
void rsRotationOscillator<T>::updateTransformMatrix()
{
  transformRotation.setAngles(trafoRotX, trafoRotY, trafoRotZ);
  trafoMatrixNeedsUpdate = false;
}

template<class T>
void rsRotationOscillator<T>::updateOutputMatrix()
{
  outputRotation.setAngles(outRotX, outRotY, outRotZ);
  outputMatrixNeedsUpdate = false;
}