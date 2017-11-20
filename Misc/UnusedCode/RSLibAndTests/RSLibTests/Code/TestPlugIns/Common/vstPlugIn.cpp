#include "vstPlugIn.h"

vstPlugIn::vstPlugIn(audioMasterCallback audioMaster, VstInt32 numPrograms, VstInt32 numParams) 
  : AudioEffectX(audioMaster, numPrograms, numParams)
{
  setNumInputs(2);
  setNumOutputs(2);
  canProcessReplacing();
  canDoubleReplacing();
  programsAreChunks(false); // check, if this is needed
  vst_strncpy(programName, "Default", kVstMaxProgNameLen);
  params = new float[numParams];
}

vstPlugIn::~vstPlugIn()
{
  delete[] params;
}

void vstPlugIn::setProgramName (char* name)
{
  vst_strncpy(programName, name, kVstMaxProgNameLen);
}
void vstPlugIn::getProgramName (char* name)
{
  vst_strncpy(name, programName, kVstMaxProgNameLen);
}

bool vstPlugIn::getProductString(char* text)
{
  return getEffectName(text);
}
bool vstPlugIn::getVendorString (char* text)
{
  vst_strncpy(text, "RS-MET", kVstMaxVendorStrLen);
  return true;
}
VstInt32 vstPlugIn::getVendorVersion ()
{ 
  return 1000; 
}

void vstPlugIn::processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames)
{
  //mutex.lock();

  // maybe implement parameter-smoothing here later

  for(int n = 0; n < sampleFrames; n++)
  {
    double inL = (double) inputs[0][n];
    double inR = (double) inputs[1][n];
    processStereoFrame(&inL, &inR, &inL, &inR);
    outputs[0][n] = (float) inL;
    outputs[1][n] = (float) inR;
  }
  //mutex.unlock();
}
void vstPlugIn::processDoubleReplacing(double** inputs, double** outputs, VstInt32 sampleFrames)
{
  //mutex.lock();
  for(int n = 0; n < sampleFrames; n++)
    processStereoFrame(&inputs[0][n], &inputs[1][n], &outputs[0][n], &outputs[1][n]);
  //mutex.unlock();
}

void vstPlugIn::setParameter(VstInt32 index, float value)
{
  params[index] = value;
  //mutex.lock();
  updateCoreParameter(index, value);
  //mutex.unlock();
}
 
float vstPlugIn::getParameter(VstInt32 index)
{
  return params[index];
}
