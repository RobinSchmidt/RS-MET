#ifndef VST_PHONOFILTER_H
#define VST_PHONOFILTER_H

#include "../Common/Utilities/vstPlugIn.h"

class vstPhonoFilter : public vstPlugIn
{

public:
	
  vstPhonoFilter(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void  updateCoreParameter(VstInt32 index, float value);
  virtual void  getParameterLabel  (VstInt32 index,  char* label);
  virtual void  getParameterDisplay(VstInt32 index,  char* text);
  virtual void  getParameterName   (VstInt32 index,  char* text);

protected:

  // parameter indices:
  enum parameters
  {
    MODE = 0,

    NUM_PARAMETERS
  };

  // mapped parameters:
  int mode;

  // audio processing core:
  rsPhonoFilter filterL, filterR;  // later - use one object of type rsPhonoFilterStereo

};

#endif
