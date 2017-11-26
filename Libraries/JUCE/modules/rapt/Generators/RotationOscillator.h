#ifndef RAPT_ROTATIONOSCILLATOR_H_INCLUDED
#define RAPT_ROTATIONOSCILLATOR_H_INCLUDED

/** An oscillator based on a vector that rotates in 3D space. The vector produced by this basic 
rotation produces values that lie on the unit sphere. this unit-sphere vector is then shifted and 
scaled along the 3 coordinate axes - this produces a vector on an ellipsoid whose axes are aligned
with the coordinate axes. Then, this ellipsoid is rotated around the 3 axes to produce a general 
ellipsoid. After that, the length of the vector is renormalized (so it will end up on the unit 
sphere again). 

...Then, an output signal is produced by projecting.. 

maybe rename to SphereProjectionOsc or EllispoidOsc

This oscillator is a 3D generalization of an algorithm proposed by xoxos on KVR:

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

  void setRenormalizationAmount(T newAmount)
  {
    renormExponent = newAmount;
  }

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

  T sampleRate = 44100; 
  T freq = 100;
  T freqScaleX = 1;
  T freqScaleY = 1;
  T freqScaleZ = 0;

  T phaseX = 0;
  //T phaseY = T(0.5*PI);
  T phaseY = 0;
  T phaseZ = 0;

  T incX = 0, incY = 0, incZ = 0;  // phase increments
  T posX = 0, posY = 0, posZ = 0;  // current position/phase

  T renormExponent = 0;


  T shiftX = 0;
  T shiftY = 0;
  T shiftZ = 0;

  T scaleX = 1;
  T scaleY = 1;
  T scaleZ = 1;

  T trafoRotX = 0;
  T trafoRotY = 0;
  T trafoRotZ = 0;

  T outRotX = 0;
  T outRotY = 0;
  T outRotZ = 0;

  bool trafoMatrixNeedsUpdate  = true;
  bool outputMatrixNeedsUpdate = true;
};

#endif