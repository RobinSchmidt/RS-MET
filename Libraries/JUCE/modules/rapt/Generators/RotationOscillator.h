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
class rsRotationOscillator
{

public:



  void setSampleRate(T newSampleRate)
  {
    sampleRate = newSampleRate;
    oscMatrixNeedsUpdate = true;
  }

  void setFrequency(T newFrequency)
  {
    freq = newFrequency;
    oscMatrixNeedsUpdate = true;
  }

  void setFrequencyScalerX(T newScaler)
  {
    freqScaleX = newScaler;
    oscMatrixNeedsUpdate = true;
  }

  void setFrequencyScalerY(T newScaler)
  {
    freqScaleY = newScaler;
    oscMatrixNeedsUpdate = true;
  }

  void setFrequencyScalerZ(T newScaler)
  {
    freqScaleZ = newScaler;
    oscMatrixNeedsUpdate = true;
  }



  void processSampleFrame(T* x, T* y, T* z);


  T getSample()
  {
    T x, y, z;
    processSampleFrame(&x, &y, &z);
    return x + y + z;
  }


  /** Resets the state vactor to (x,y,z) = (1,0,0). */
  void reset();


protected:

  void updateOscillationMatrix();

  void updateTransformMatrix();



  rsRotationXYZ<T> oscillateRotation;
  rsRotationXYZ<T> transformRotation;
  //rsRotationXYZ<T> outputRotation; // ..or maybe use 2x3 projection matrix

  T sampleRate = 44100; 
  T freq = 1000;
  T freqScaleX = 1;
  T freqScaleY = 0;
  T freqScaleZ = 0;
  T shiftX = 0;
  T shiftY = 0;
  T shiftZ = 0;
  T scaleX = 0;
  T scaleY = 0;
  T scaleZ = 0;
  T rotX = 0;
  T rotY = 0;
  T rotZ = 0;

  T X = 1, Y = 0, Z = 0;

  bool oscMatrixNeedsUpdate = true;
  bool trafoMatrixNeedsUpdate = true;

};

#endif