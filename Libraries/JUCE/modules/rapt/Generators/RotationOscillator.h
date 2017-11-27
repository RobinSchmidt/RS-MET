#ifndef RAPT_ROTATIONOSCILLATOR_H_INCLUDED
#define RAPT_ROTATIONOSCILLATOR_H_INCLUDED

/** This is a sound generator based on 3-dimensional Lissajous figures. The x, y and z coordinates
all oscillate with their own frequency....

todo: 
-Level (dB), LevelByFreq (dB/oct), shift, scale, rot, waveshape parameters to warp sines into saws and/or
 squares (phase-shaping)
-rename to rsOscillator3D


*/

template<class T>
class rsLissajousOscillator3D
{

public:

  void setSampleRate(T newSampleRate);
  void setFrequency(T newFrequency);
  void setFrequencyScalerX(T newScaler);
  void setFrequencyScalerY(T newScaler);
  void setFrequencyScalerZ(T newScaler);
  void setFrequencyOffsetX(T newOffset);
  void setFrequencyOffsetY(T newOffset);
  void setFrequencyOffsetZ(T newOffset);


  void setPhaseX(T newPhase) { phaseX = rsDegreeToRadiant(newPhase); }
  void setPhaseY(T newPhase) { phaseY = rsDegreeToRadiant(newPhase); }
  void setPhaseZ(T newPhase) { phaseZ = rsDegreeToRadiant(newPhase); }

  void setOutputRotationX(T newAngle)
  {
    outRotX = rsDegreeToRadiant(newAngle);
    outputMatrixNeedsUpdate = true;
  }

  void setOutputRotationY(T newAngle)
  {
    outRotY = rsDegreeToRadiant(newAngle);
    outputMatrixNeedsUpdate = true;
  }

  void setOutputRotationZ(T newAngle)
  {
    outRotZ = rsDegreeToRadiant(newAngle);
    outputMatrixNeedsUpdate = true;
  }

  /** Sets the maount by which the length of the x,y,z will be renormalized to unit length. It 
  works as follows: the length of the vector is computes, then the reciprocal r = 1/length of this 
  length is obtained and x,y,z value are multiplied by r^p where "p" is a power that is adjusted
  by this function.  */
  void setRenormalizationAmount(T newAmount) { renormExponent = newAmount; }

  /** Sets the amplitude at which x,y,z values are clipped. Th nonlinar range begins at half of 
  this value, so when you pass 2, all amplitudes between -1 and +1 will be passed without 
  distortion. */
  void setClipAmplitude(T newAmplitude) { clip = newAmplitude; clipInv = 1/clip; }

  void processSampleFrame(T* x, T* y, T* z);




  void getSampleFrameStereo(T* left, T* right)
  {
    T x, y, z;
    processSampleFrame(&x, &y, &z);
    *left  = x;
    *right = y;
  }

  T getSample()
  {
    T x, y, z;
    processSampleFrame(&x, &y, &z);
    return x + y + z;
  }


  /** Resets the phases. */
  void reset();


protected:


  void updatePhaseIncrements();
  void updateTransformMatrix();
  void updateOutputMatrix();

  rsRotationXYZ<T> transformRotation;
  rsRotationXYZ<T> outputRotation; // ..or maybe use 2x3 projection matrix

  T sampleRate  = 44100; 
  T omegaFactor = 2*PI/44100;
  T freq = 100;

  T freqScalerX = 1, freqScalerY = 1, freqScalerZ = 1;
  T freqOffsetX = 1, freqOffsetY = 1, freqOffsetZ = 1;

  T phaseX = 0;
  //T phaseY = T(0.5*PI);
  T phaseY = 0;
  T phaseZ = 0;

  T incX = 0, incY = 0, incZ = 0;  // phase increments
  T posX = 0, posY = 0, posZ = 0;  // current position/phase

  T renormExponent = 0;
  T clip = 2, clipInv = T(0.5);


  T shiftX = 0, shiftY = 0, shiftZ = 0;
  T scaleX = 1, scaleY = 1, scaleZ = 1;

  T trafoRotX = 0, trafoRotY = 0, trafoRotZ = 0;
  T outRotX   = 0, outRotY   = 0, outRotZ   = 0;

  bool trafoMatrixNeedsUpdate  = true;
  bool outputMatrixNeedsUpdate = true;
};

#endif