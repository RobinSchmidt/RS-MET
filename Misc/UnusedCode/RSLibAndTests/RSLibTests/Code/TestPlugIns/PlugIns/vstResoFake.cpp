#include "vstResoFake.h"

// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstResoFake(audioMaster);
}

vstResoFake::vstResoFake(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('RsFk');

  // set up intial settings:
  setParameter(FREQUENCY,   0.5);
  setParameter(RESOSHIFT,   0.5);     // 0 semitones
  setParameter(RESOGAIN,    0.5);     // 0 amplitude
  setParameter(RESOPHASE,   0.0);     // PI/2
  setParameter(RESODECAY,   0.25);
  setParameter(DECAYBYFREQ, 2/3.f);   // 100% -> constant Q
  setParameter(RESOATTACK,  0.0);     // (almost) instant attack
  setParameter(RESODELAY,   0.5);     // 100% -> optimally aligned
}

bool vstResoFake::getEffectName(char* name)
{
  vst_strncpy(name, "ResoFake", kVstMaxEffectNameLen);
  return true;
}

void vstResoFake::processStereoFrame(double *inL, double *inR, double *outL, double *outR)
{
  // preliminary - mix input to mono, process, broadcast output to both channels (later introduce
  // stereo processing):
  double in  = 0.5 * (*inL + *inR);
  double out = filter.getSample(in);
  *outL = out;
  *outR = out;
}

void vstResoFake::resume()
{
  filter.reset();
}

void vstResoFake::setSampleRate(float sampleRate)
{
  filter.setSampleRate(sampleRate);
}

// Parameters:

void vstResoFake::updateCoreParameter(VstInt32 index, float value)
{
  switch(index)
  {
  case FREQUENCY: {
    freq = rsLinToExp(params[FREQUENCY], 0, 1, 20, 20000);
    filter.setLowpassCutoff(freq);
  } break;
  case RESOSHIFT: {
    shift = rsLinToLin(params[RESOSHIFT], 0, 1, -24, 24);
    filter.setResonanceShift(shift);
  } break;
  case RESOGAIN: {
    gain = rsLinToLin(params[RESOGAIN], 0, 1, -2, 2);
    filter.setResonanceAmplitude(gain);
  } break;
  case RESOPHASE: {
    phase = rsLinToLin(params[RESOPHASE], 0, 1, 0, 360);            // degrees here
    filter.setResonancePhase(phase*PI/180);                         // radians there
  } break;
  case RESODECAY: {
    decay = rsLinToExp(params[RESODECAY], 0, 1, 0.1, 1000.0);       // milliseconds here
    filter.setResonanceDecay(0.001*decay);                          // seconds there
  } break;
  case DECAYBYFREQ: {
    decayByFreq = rsLinToLin(params[DECAYBYFREQ], 0, 1, -100, 200);  // % here
    filter.setDecayByFrequency(0.01*decayByFreq);                    // 0..1 there
  } break;
  case RESOATTACK: {
    attack = rsLinToLin(params[RESOATTACK], 0, 1, 1, 99);            // % here
    filter.setResonanceAttack(0.01*attack);                          // 0..1 there
  } break;
  case RESODELAY: {
    delay = rsLinToLin(params[RESODELAY], 0, 1, 0, 200);             // % here
    filter.setResonanceDelay(0.01*delay);                            // 0..1 there
  } break;

  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only space for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstResoFake::getParameterName(VstInt32 index, char* label)
{
  switch( index )
  {
  case FREQUENCY:   vst_strncpy(label, "Cutoff",   kVstMaxParamStrLen); break;
  case RESOSHIFT:   vst_strncpy(label, "ResShift", kVstMaxParamStrLen); break;
  case RESOGAIN:    vst_strncpy(label, "ResGain",  kVstMaxParamStrLen); break;
  case RESOPHASE:   vst_strncpy(label, "ResPhs",   kVstMaxParamStrLen); break;
  case RESODECAY:   vst_strncpy(label, "ResDec",   kVstMaxParamStrLen); break;
  case DECAYBYFREQ: vst_strncpy(label, "-byFreq",  kVstMaxParamStrLen); break;
  case RESOATTACK:  vst_strncpy(label, "ResAtt",   kVstMaxParamStrLen); break;
  case RESODELAY:   vst_strncpy(label, "ResDly",   kVstMaxParamStrLen); break;
  } 
}

void vstResoFake::getParameterDisplay(VstInt32 index, char* text)
{
  switch( index )
  {
  case FREQUENCY:   sprintf(text, "%.2f", freq);         break;
  case RESOSHIFT:   sprintf(text, "%.2f", shift);        break;
  case RESOGAIN:    sprintf(text, "%.2f", gain);         break;
  case RESOPHASE:   sprintf(text, "%.2f", phase);        break;
  case RESODECAY:   sprintf(text, "%.2f", decay);        break;
  case DECAYBYFREQ: sprintf(text, "%.2f", decayByFreq);  break;
  case RESOATTACK:  sprintf(text, "%.2f", attack);       break;
  case RESODELAY:   sprintf(text, "%.2f", delay);        break;
  }
}

void vstResoFake::getParameterLabel(VstInt32 index, char* label)
{
  switch( index )
  {
  case FREQUENCY:   vst_strncpy(label, "Hz",  kVstMaxParamStrLen);  break;
  case RESOSHIFT:   vst_strncpy(label, "st",  kVstMaxParamStrLen);  break;
  case RESOGAIN:    vst_strncpy(label, "",    kVstMaxParamStrLen);  break;
  case RESOPHASE:   vst_strncpy(label, "deg", kVstMaxParamStrLen);  break;
  case RESODECAY:   vst_strncpy(label, "ms",  kVstMaxParamStrLen);  break;
  case DECAYBYFREQ: vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  case RESOATTACK:  vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  case RESODELAY:   vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  } 
}
