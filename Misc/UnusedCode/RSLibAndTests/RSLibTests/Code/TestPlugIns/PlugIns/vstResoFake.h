#ifndef VST_RESOFAKE_H
#define VST_RESOFAKE_H

#include "../Common/Utilities/vstPlugIn.h"

class vstResoFake : public vstPlugIn
{

public:
	
  vstResoFake(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void  updateCoreParameter(VstInt32 index, float value);
  virtual void  getParameterLabel  (VstInt32 index, char* label);
  virtual void  getParameterDisplay(VstInt32 index, char* text);
  virtual void  getParameterName   (VstInt32 index, char* text);

protected:

  // parameter indices:
  enum parameters
  {
    FREQUENCY = 0,   // lowpass cutoff frequency
    RESOSHIFT,       // resonance pitch-shift with respect to cutoff frequency
    RESOGAIN,        // resonance gain
    RESOPHASE,       // resonance phase
    RESODECAY,       // resonance decay time
    DECAYBYFREQ,     // ...scaling of decaytime by frequency
    RESOATTACK,      // resonance attack time
    RESODELAY,       // delay of resonance path

    NUM_PARAMETERS
  };

  // mapped parameters (maybe use an array of doubles):
  double freq, shift, gain, phase, decay, decayByFreq, attack, delay;

  // audio processing core:
  rsFakeResonanceFilter filter;

};

#endif
