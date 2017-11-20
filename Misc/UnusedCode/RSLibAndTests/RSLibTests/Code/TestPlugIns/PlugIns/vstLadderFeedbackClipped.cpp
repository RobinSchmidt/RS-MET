#include "vstLadderFeedbackClipped.h"


// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstLadderFeedbackClipped(audioMaster);
}

vstLadderFeedbackClipped::vstLadderFeedbackClipped(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('LdFC');

  // set up intial settings:
  setParameter(FREQUENCY,   0.5);
  //setParameter(RESODECAY,   0.5);
  // setParameter(FB_LPF,      0.5);
  setParameter(RESO,        0.0);     // 0
  setParameter(LOLIMIT,     0.5);     // 1 
  setParameter(HILIMIT,     0.5);
  setParameter(SHAPE,       1.0);     // hard clip
  setParameter(MODE,        0.0);
}

bool vstLadderFeedbackClipped::getEffectName(char* name)
{
  vst_strncpy(name, "LadderClippedFB", kVstMaxEffectNameLen);
  return true;
}

void vstLadderFeedbackClipped::processStereoFrame(double *inL, double *inR, double *outL, 
  double *outR)
{
  // preliminary - mix input to mono, process, broadcast output to both channels (later introduce
  // stereo processing):
  double in  = 0.5 * (*inL + *inR);
  double out = filter.getSample(in);
  *outL = out;
  *outR = out;
}

void vstLadderFeedbackClipped::resume()
{
  filter.reset();
}

void vstLadderFeedbackClipped::setSampleRate(float sampleRate)
{
  filter.setSampleRate(sampleRate);
}

// Parameters:

void vstLadderFeedbackClipped::updateCoreParameter(VstInt32 index, float value)
{
  switch(index)
  {
  case FREQUENCY: {
    freq = rsLinToExp(params[FREQUENCY], 0, 1, 20, 20000);
    filter.setCutoff(freq);
  } break;

  //case FB_LPF: {
  //  filter.fbLpCutoff = rsLinToExp(value, 0, 1, 0.1, 10.0);  // temporary
  //} break;

  case RESO: {
    reso = rsLinToLin(params[RESO], 0, 1, 0, 10);
    //reso = rsLinToExpWithOffset(params[RESO], 0, 1, 0.01, 100, -0.01);
    filter.setResonance(reso);
  } break;
  case LOLIMIT: {
    loLimit = -rsLinToExp(params[LOLIMIT], 0, 1, 0.01, 10.0);
    filter.setLowerFeedbackLimit(loLimit);
  } break;
  case HILIMIT: {
    hiLimit = rsLinToExp(params[HILIMIT], 0, 1, 0.01, 10.0);
    filter.setUpperFeedbackLimit(hiLimit);
  } break;
  case SHAPE: {
    shape = rsLinToLin(params[SHAPE], 0, 1, 0.5, 1.0);
    filter.setSaturationGainAt1(shape);
  } break;
  case MODE: {
    mode = rsRoundToInt(rsLinToLin(value, 0, 1, 1, rsLadderFilterFeedbackSaturated::NUM_SATPLACES-1));
    filter.setSaturationMode(mode);
  } break;

  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only space for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstLadderFeedbackClipped::getParameterName(VstInt32 index, char* label)
{
  switch( index )
  {
  case FREQUENCY:   vst_strncpy(label, "Cutoff",    kVstMaxParamStrLen); break;
  //case FB_LPF:      vst_strncpy(label, "FBLowpass", kVstMaxParamStrLen); break;
  case RESO:        vst_strncpy(label, "Reso",      kVstMaxParamStrLen); break;
  case LOLIMIT:     vst_strncpy(label, "LoLimit",   kVstMaxParamStrLen); break;
  case HILIMIT:     vst_strncpy(label, "HiLimit",   kVstMaxParamStrLen); break;
  case SHAPE:       vst_strncpy(label, "SatShape",  kVstMaxParamStrLen); break;
  case MODE:        vst_strncpy(label, "SatMode",   kVstMaxParamStrLen); break;
  } 
}

void vstLadderFeedbackClipped::getParameterDisplay(VstInt32 index, char* text)
{
  switch( index )
  {
  case FREQUENCY: sprintf(text, "%.2f", freq);    break;
  //case FB_LPF:    sprintf(text, "%.2f", filter.fbLpCutoff);  break;
  case RESO:      sprintf(text, "%.2f", reso);    break;
  case LOLIMIT:   sprintf(text, "%.2f", loLimit); break;
  case HILIMIT:   sprintf(text, "%.2f", hiLimit); break;
  case SHAPE:     sprintf(text, "%.2f", shape);   break;
  case MODE:      sprintf(text, "%d",   mode);    break;
  }
}

void vstLadderFeedbackClipped::getParameterLabel(VstInt32 index, char* label)
{
  switch( index )
  {
  case FREQUENCY: vst_strncpy(label, "Hz",  kVstMaxParamStrLen);  break;

  //case RESODECAY: vst_strncpy(label, "ms",  kVstMaxParamStrLen);  break;
  //case RESODRIVE: vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;
  //case LOLIMIT:   vst_strncpy(label, "",    kVstMaxParamStrLen);  break;
  //case HILIMIT:   vst_strncpy(label, "",    kVstMaxParamStrLen);  break;

  default:        vst_strncpy(label, "",    kVstMaxParamStrLen);  break;
  } 
}

//void vstLadderFeedbackClipped::setWaveShape(VstInt32 index)
//{
//  switch( index )
//  {
//  case CLIP:    filter.setSaturationFunction(&rsNormalizedSigmoids::clip);    break;
//  case CUBIC:   filter.setSaturationFunction(&rsNormalizedSigmoids::cubic);   break;
//  case QUARTIC: filter.setSaturationFunction(&rsNormalizedSigmoids::quartic); break;
//  case SIXTIC:  filter.setSaturationFunction(&rsNormalizedSigmoids::sixtic);  break;
//  case TANH:    filter.setSaturationFunction(&rsNormalizedSigmoids::tanh);    break;
//  case ATAN:    filter.setSaturationFunction(&rsNormalizedSigmoids::atan);    break;
//  } 
//}
//
//void vstLadderFeedbackClipped::getWaveShapeName(VstInt32 index, char* text)
//{
//  switch( index )
//  {
//  case CLIP:    vst_strncpy(text, "Clip",    kVstMaxParamStrLen);  break;
//  case CUBIC:   vst_strncpy(text, "Cubic",   kVstMaxParamStrLen);  break;
//  case QUARTIC: vst_strncpy(text, "Quartic", kVstMaxParamStrLen);  break;
//  case SIXTIC:  vst_strncpy(text, "Sixtic",  kVstMaxParamStrLen);  break;
//  case TANH:    vst_strncpy(text, "Tanh",    kVstMaxParamStrLen);  break;
//  case ATAN:    vst_strncpy(text, "Atan",    kVstMaxParamStrLen);  break;
//  } 
//}
