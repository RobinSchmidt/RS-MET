template<class T>
void rsLissajousOscillator3D<T>::setSampleRate(T newSampleRate)
{
  sampleRate = newSampleRate;
  omegaFactor = 2*PI / sampleRate;
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
  freqScalerX = newScaler;
  incX = omegaFactor * (freqScalerX * freq + freqOffsetX);
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyScalerY(T newScaler)
{
  freqScalerY = newScaler;
  incY = omegaFactor * (freqScalerY * freq + freqOffsetY);
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyScalerZ(T newScaler)
{
  freqScalerZ = newScaler;
  incZ = omegaFactor * (freqScalerZ * freq + freqOffsetZ);
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyOffsetX(T newOffset)
{
  freqOffsetX = newOffset;
  incX = omegaFactor * (freqOffsetX * freq + freqOffsetX);
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyOffsetY(T newOffset)
{
  freqOffsetY = newOffset;
  incY = omegaFactor * (freqOffsetY * freq + freqOffsetY);
}

template<class T>
void rsLissajousOscillator3D<T>::setFrequencyOffsetZ(T newOffset)
{
  freqOffsetZ = newOffset;
  incZ = omegaFactor * (freqOffsetZ * freq + freqOffsetZ);
}

template<class T>
void rsLissajousOscillator3D<T>::processSampleFrame(T* x, T* y, T* z)
{
  // i think, it reduces to xoxos osc when
  // wz = 0, wx = wy, px = 0°, py = 90°

  // update matrices, if necessarry:
  if(trafoMatrixNeedsUpdate)
    updateTransformMatrix();
  if(outputMatrixNeedsUpdate)
    updateOutputMatrix();

  // obtain new vector on unit sphere:
  //oscillateRotation.apply(&X, &Y, &Z);

  T X = sin(posX + phaseX);
  T Y = sin(posY + phaseY);
  T Z = sin(posZ + phaseZ);

  // transform to ellipsoid:
  T tx = (X + shiftX) * scaleX;
  T ty = (Y + shiftY) * scaleY;
  T tz = (Z + shiftZ) * scaleZ;

  // rotate the ellipsoid:
  transformRotation.apply(&tx, &ty, &tz);

  // renormalize length:
  T c = 0.001; // to avoid div-by-zero
  T s = 1 / sqrt(c + tx*tx + ty*ty + tz*tz); // what about div-by-zero?
  s   = pow(s, renormExponent);
  tx *= s;
  ty *= s;
  tz *= s;

  //...maybe apply spherical (soft) clipping?
  T clipInv = 1 / clip; // make member
  tx = clip * rsNormalizedSigmoids<T>::softClipHexic(clipInv*tx);
  ty = clip * rsNormalizedSigmoids<T>::softClipHexic(clipInv*ty);
  tz = clip * rsNormalizedSigmoids<T>::softClipHexic(clipInv*tz);
  //ty /= 1 / (1+ty*ty);
  //tz /= 1 / (1+tz*tz);


  // apply the output rotation and assign outputs:
  outputRotation.apply(&tx, &ty, &tz);
  *x = tx;
  *y = ty;
  *z = tz;

  posX += incX;
  posY += incY;
  posZ += incZ;

  posX = rsWrapAround(posX, 2*PI);
  posY = rsWrapAround(posY, 2*PI);
  posZ = rsWrapAround(posZ, 2*PI);
}

template<class T>
void rsLissajousOscillator3D<T>::reset()
{
  posX = 0;
  posY = 0;
  posZ = 0;
}

template<class T>
void rsLissajousOscillator3D<T>::updatePhaseIncrements()
{
  incX = omegaFactor * (freqScalerX * freq + freqOffsetX);
  incY = omegaFactor * (freqScalerY * freq + freqOffsetY);
  incZ = omegaFactor * (freqScalerZ * freq + freqOffsetZ);
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