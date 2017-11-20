#include "vstPhonoFilter.h"

// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstPhonoFilter(audioMaster);
}

vstPhonoFilter::vstPhonoFilter(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('PhFl');
  setParameter(0, rsPhonoFilter::DE_EMPHASIS);
}

bool vstPhonoFilter::getEffectName(char* name)
{
  vst_strncpy(name, "PhonoFilter", kVstMaxEffectNameLen);
  return true;
}

void vstPhonoFilter::processStereoFrame(double *inL, double *inR, double *outL, double *outR)
{
  *outL = filterL.getSample(*inL);
  *outR = filterR.getSample(*inR);
}

void vstPhonoFilter::resume()
{
  filterL.reset();
  filterR.reset();
}

void vstPhonoFilter::setSampleRate(float sampleRate)
{
  filterL.setSampleRate(sampleRate);
  filterR.setSampleRate(sampleRate);
}

// Parameters:

void vstPhonoFilter::updateCoreParameter(VstInt32 index, float value)
{
  mode = rsRoundToInt(value);
  filterL.setMode(mode);
  filterR.setMode(mode);
}

void vstPhonoFilter::getParameterName(VstInt32 index, char* label)
{
  vst_strncpy(label, "Mode", kVstMaxParamStrLen);
}

void vstPhonoFilter::getParameterDisplay(VstInt32 index, char* text)
{
  switch( mode )
  {
  case rsPhonoFilter::PRE_EMPHASIS: vst_strncpy(text, "PreEmphasis", kVstMaxParamStrLen); break;
  case rsPhonoFilter::DE_EMPHASIS:  vst_strncpy(text, "DeEmphasis",  kVstMaxParamStrLen); break;
  }
}

void vstPhonoFilter::getParameterLabel(VstInt32 index, char* label)
{
  vst_strncpy(label, "", kVstMaxParamStrLen);
}
