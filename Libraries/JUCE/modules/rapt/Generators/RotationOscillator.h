#ifndef RAPT_ROTATIONOSCILLATOR_H_INCLUDED
#define RAPT_ROTATIONOSCILLATOR_H_INCLUDED

template<class T>
class rsRotationOscillator
{

public:

  void setFrequency(T newFrequency);

  void setSampleRate(T newSampleRate);

  void setFrequencyScalerX(T newScaler);

  void setFrequencyScalerY(T newScaler);

  void setFrequencyScalerZ(T newScaler);


  void processSampleFrame(T* x, T* y, T* z);


  /** Resets the state vactor to (x,y,z) = (1,0,0). */
  void reset()
  {
    X = 1;
    Y = 0;
    Z = 0;
  }


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