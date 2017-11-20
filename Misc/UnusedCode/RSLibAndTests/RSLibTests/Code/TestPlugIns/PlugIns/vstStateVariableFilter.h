#ifndef VST_STATEVARIABLEFILTER_H
#define VST_STATEVARIABLEFILTER_H

#include "../Common/Utilities/vstPlugIn.h"

class vstStateVariableFilter : public vstPlugIn
{

public:
	
  vstStateVariableFilter(audioMasterCallback audioMaster);	
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
    MODE = 0,
    FREQUENCY, 
    GAIN,
    BANDWIDTH,

    // todo: MORPH

    NUM_PARAMETERS
  };

  // mapped parameters:
  int    mode;
  double frequency, gain, bandwidth;

  // audio processing core:
  rsStateVariableFilter filterL, filterR; 

};

#endif
