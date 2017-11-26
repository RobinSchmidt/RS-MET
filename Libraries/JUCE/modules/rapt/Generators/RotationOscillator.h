#ifndef RAPT_ROTATIONOSCILLATOR_H_INCLUDED
#define RAPT_ROTATIONOSCILLATOR_H_INCLUDED

template<class T>
class rsRotationOscillator
{

public:

  void setFrequency(T newFrequency);

  void setSampleRate(T newSampleRate);

  void setFrequencyScaler1(T newScaler);

  void setFrequencyScaler2(T newScaler);

  void setFrequencyScaler3(T newScaler);


  void processSampleFrame(T* x, T* y, T* z);


  /** Resets the state vactor to (x,y,z) = (1,0,0). */
  void reset()
  {
    X = 1;
    Y = 0;
    Z = 0;
  }


protected:

  void updateRotationMatrix();

  rsRotationXYZ<T> oscillateRotation;
  rsRotationXYZ<T> ellipsoidRotation;
  //rsRotationXYZ<T> outputRotation;

  T sampleRate = 44100; 
  T freq = 1000;
  T freqScale1 = 1;
  T freqScale2 = 0;
  T freqScale3 = 0;
  T shiftX = 0;
  T shiftY = 0;
  T shiftZ = 0;
  T scaleX = 0;
  T scaleY = 0;
  T scaleZ = 0;

  T X = 1, Y = 0, Z = 0;

  bool matrixNeedsUpdate = true;

};

#endif