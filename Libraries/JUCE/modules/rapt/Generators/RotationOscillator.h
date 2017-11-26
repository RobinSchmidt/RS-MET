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

  void setOutputRotationX(T newAngle)
  {
    outRotX = PI * newAngle / 180;
    outputMatrixNeedsUpdate = true;
  }

  void setOutputRotationY(T newAngle)
  {
    outRotY = PI * newAngle / 180;
    outputMatrixNeedsUpdate = true;
  }

  void setOutputRotationZ(T newAngle)
  {
    outRotZ = PI * newAngle / 180;
    outputMatrixNeedsUpdate = true;
  }



  void processSampleFrame(T* x, T* y, T* z);


  void getSampleFrameStereo(T* left, T* right)
  {
    T x, y, z;
    processSampleFrame(&x, &y, &z);
    //*left  = x+z;  // preliminary - use projection matrix later
    //*right = y-z;

    *left  = x;  // preliminary - use projection matrix later
    *right = y;
  }

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

  void updateOutputMatrix();



  rsRotationXYZ<T> oscillateRotation;
  rsRotationXYZ<T> transformRotation;
  rsRotationXYZ<T> outputRotation; // ..or maybe use 2x3 projection matrix

  T sampleRate = 44100; 
  T freq = 1000;
  T freqScaleX = 0;
  T freqScaleY = 0;
  T freqScaleZ = 1;
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


  T X = 1, Y = 0, Z = 0;

  bool oscMatrixNeedsUpdate    = true;
  bool trafoMatrixNeedsUpdate  = true;
  bool outputMatrixNeedsUpdate = true;

};

#endif