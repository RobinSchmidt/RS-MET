#ifndef RAPT_ROTATIONOSCILLATOR_H_INCLUDED
#define RAPT_ROTATIONOSCILLATOR_H_INCLUDED

template<class T>
class rsRotationOscillator
{

public:

protected:

  rsRotationXYZ<T> oscRotation;

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

};

#endif