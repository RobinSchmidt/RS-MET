#ifndef VST_PLUGIN_H
#define VST_PLUGIN_H

#include "../../../../../RSLib/Code/RSLib.h"
#include "public.sdk/source/vst2.x/audioeffectx.h"

// common baseclass for the various test-plugins:
class vstPlugIn : public AudioEffectX 
{

public:
	
  vstPlugIn(audioMasterCallback audioMaster, VstInt32 numPrograms, VstInt32 numParams);	

  ~vstPlugIn();	


  virtual void processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames);
  virtual void processDoubleReplacing(double** inputs, double** outputs, VstInt32 sampleFrames);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR) = 0;

  virtual void  setParameter(VstInt32 index, float value);
  virtual void  updateCoreParameter(VstInt32 index, float value) = 0;
  virtual float getParameter(VstInt32 index);

  virtual void setProgramName(char* name);
  virtual void getProgramName(char* name);

  virtual bool getVendorString(char* text);
  virtual bool getProductString(char* text);
  virtual VstInt32 getVendorVersion ();

protected:

  char programName[kVstMaxProgNameLen+1];
  static const VstInt32 defaultNumPrograms = 1;

  float *params;   // normalized vst parameter values

  //rsMutex mutex; // use later - for thread-safe parameter updates

};

// make a class rsParameter

#endif
