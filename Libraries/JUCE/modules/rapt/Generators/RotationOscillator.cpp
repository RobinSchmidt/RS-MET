

template<class T>
void rsLissajousOscillator3D<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;
  updatePhaseIncrements();
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequency(T newFrequency)
{
  freq = newFrequency;
  updatePhaseIncrements();
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyScalerX(T newScaler)
{
  freqScaleX = newScaler;
  incX = freqScaleX * 2*T(PI)*freq/sampleRate;
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyScalerY(T newScaler)
{
  freqScaleY = newScaler;
  incY = freqScaleY * 2*T(PI)*freq/sampleRate;
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyScalerZ(T newScaler)
{
  freqScaleZ = newScaler;
  incZ = freqScaleZ * 2*T(PI)*freq/sampleRate;
}



template<class T>
void rsLissajousOscillator3D<T>::processSampleFrame(T* x, T* y, T* z)
{
  // hmm...no - this is boring - it just rotates with one frequency
  // instead we should use 3 sine-oscillators with freq and phase
  // ...make a 3D Lissajous oscillator - it reduces to xoxos osc when
  // wz = 0, wx = wy, px = 0°, py = 90°

  // update matrices, if necessarry:
  if(trafoMatrixNeedsUpdate)
    updateTransformMatrix();
  if(outputMatrixNeedsUpdate)
    updateOutputMatrix();

  // obtain new vector on unit sphere:
  //oscillateRotation.apply(&X, &Y, &Z);

  T X = sin(posX);
  T Y = sin(posY);
  T Z = sin(posZ);

  // transform to ellipsoid:
  T tx = (X + shiftX) * scaleX;
  T ty = (Y + shiftY) * scaleY;
  T tz = (Z + shiftZ) * scaleZ;

  // rotate the ellipsoid:
  transformRotation.apply(&tx, &ty, &tz);

  // renormalize length:
  T s = 1 / sqrt(tx*tx + ty*ty + tz*tz); // what about div-by-zero?
  s   = pow(s, renormExponent);
  tx *= s;
  ty *= s;
  tz *= s;

  //...maybe apply spherical (soft) clipping?

  // apply the output rotation and assign outputs:
  outputRotation.apply(&tx, &ty, &tz);
  *x = tx;
  *y = ty;
  *z = tz;

  posX += incX;
  posY += incY;
  posZ += incZ;

  posX = rsWrapAround(posX, 2*PI);
}

template<class T>
void rsLissajousOscillator3D<T>::reset()
{
  posX = phaseX;
  posY = phaseY;
  posZ = phaseZ;
}

template<class T>
void rsLissajousOscillator3D<T>::updatePhaseIncrements()
{
  T w  = 2*T(PI)*freq/sampleRate;
  incX = w * freqScaleX;
  incY = w * freqScaleY;
  incZ = w * freqScaleZ;
}

template<class T>
void rsLissajousOscillator3D<T>::updateTransformMatrix()
{
  transformRotation.setAngles(trafoRotX, trafoRotY, trafoRotZ);
  trafoMatrixNeedsUpdate = false;
}

template<class T>
void rsLissajousOscillator3D<T>::updateOutputMatrix()
{
  outputRotation.setAngles(outRotX, outRotY, outRotZ);
  outputMatrixNeedsUpdate = false;
}