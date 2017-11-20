#include "vstStateVariableFilter.h"

// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstStateVariableFilter(audioMaster);
}

vstStateVariableFilter::vstStateVariableFilter(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('SvFl');
  setParameter(MODE,      0.0);
  setParameter(FREQUENCY, 0.5);
  setParameter(GAIN,      0.1875);
  setParameter(BANDWIDTH, 0.5);
}

bool vstStateVariableFilter::getEffectName(char* name)
{
  vst_strncpy(name, "StateVariableFilter", kVstMaxEffectNameLen);
  return true;
}

void vstStateVariableFilter::processStereoFrame(double *inL, double *inR, 
                                                double *outL, double *outR)
{
  *outL = filterL.getSample(*inL);
  *outR = filterR.getSample(*inR);
}

void vstStateVariableFilter::resume()
{
  filterL.reset();
  filterR.reset();
}

void vstStateVariableFilter::setSampleRate(float sampleRate)
{
  filterL.setSampleRate(sampleRate);
  filterR.setSampleRate(sampleRate);
}

// Parameters:

void vstStateVariableFilter::updateCoreParameter(VstInt32 index, float value)
{
  switch( index )
  {
  case MODE:
    {
      mode  = rsRoundToInt(rsLinToLin(params[MODE], 0, 1, 0, rsStateVariableFilter::NUM_MODES-1));
      filterL.setMode(mode);
      filterR.setMode(mode);
    }
    break;
  case FREQUENCY:
    {
      frequency  = rsLinToExp(params[FREQUENCY], 0, 1, 20, 20000);
      filterL.setFrequency(frequency);
      filterR.setFrequency(frequency);
    }
    break;
  case GAIN:
    {
      gain  = rsLinToLin(params[GAIN], 0, 1, -12, 36);
      filterL.setGain(rsDB2amp(gain));
      filterR.setGain(rsDB2amp(gain));
    }
    break;
  case BANDWIDTH:
    {
      bandwidth  = rsLinToExp(params[BANDWIDTH], 0, 1, 0.25, 4);
      filterL.setBandwidth(bandwidth);
      filterR.setBandwidth(bandwidth);
    }
    break;
  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only spec for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstStateVariableFilter::getParameterName(VstInt32 index, char* label)
{
  switch( index )
  {
  case MODE:      vst_strncpy(label, "Mode",  kVstMaxParamStrLen); break;
  case FREQUENCY: vst_strncpy(label, "Freq",  kVstMaxParamStrLen); break;
  case GAIN:      vst_strncpy(label, "Gain",  kVstMaxParamStrLen); break;
  case BANDWIDTH: vst_strncpy(label, "Width", kVstMaxParamStrLen); break;
  } 
}

void vstStateVariableFilter::getParameterDisplay(VstInt32 index, char* text)
{
  switch( index )
  {
  case MODE:
    {
      switch( mode )
      {
      case rsStateVariableFilter::BYPASS:         sprintf(text, "%s", "Bypass");  break;
      case rsStateVariableFilter::LOWPASS:        sprintf(text, "%s", "LP");      break;
      case rsStateVariableFilter::HIGHPASS:       sprintf(text, "%s", "HP");      break;
      case rsStateVariableFilter::BANDPASS_SKIRT: sprintf(text, "%s", "BP-S");    break;
      case rsStateVariableFilter::BANDPASS_PEAK:  sprintf(text, "%s", "BP-P");    break;
      case rsStateVariableFilter::BANDREJECT:     sprintf(text, "%s", "BR");      break;
      case rsStateVariableFilter::BELL:           sprintf(text, "%s", "PK");      break;
      case rsStateVariableFilter::LOWSHELF:       sprintf(text, "%s", "LS");      break;
      case rsStateVariableFilter::HIGHSHELF:      sprintf(text, "%s", "HS");      break;
      case rsStateVariableFilter::ALLPASS:        sprintf(text, "%s", "AP");      break;
      case rsStateVariableFilter::MORPH_LP_BP_HP: sprintf(text, "%s", "M-LBH");   break;
      default:                                    sprintf(text, "%s", "Unknown");
      }
    }
    break;
  case FREQUENCY: sprintf(text, "%.2f", frequency); break;
  case GAIN:      sprintf(text, "%.2f", gain);      break;
  case BANDWIDTH: sprintf(text, "%.2f", bandwidth); break;
  }
}

void vstStateVariableFilter::getParameterLabel(VstInt32 index, char* label)
{
  switch( index )
  {
  case MODE:      vst_strncpy(label, "",    kVstMaxParamStrLen); break;
  case FREQUENCY: vst_strncpy(label, "Hz",  kVstMaxParamStrLen); break;
  case GAIN:      vst_strncpy(label, "dB",  kVstMaxParamStrLen); break;
  case BANDWIDTH: vst_strncpy(label, "oct", kVstMaxParamStrLen); break;
  } 
}
