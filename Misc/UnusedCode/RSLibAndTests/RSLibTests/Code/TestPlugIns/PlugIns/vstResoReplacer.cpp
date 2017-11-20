#include "vstResoReplacer.h"

// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstResoReplacer(audioMaster);
}

vstResoReplacer::vstResoReplacer(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('RsRp');

  // set up intial settings:

  setParameter(INGAIN,       1./3.);      // +- 0 dB
  setParameter(LEAK,         0.5);        // +- 0
  setParameter(FREQUENCY,    0.5);
  setParameter(RESODECAY,    0.25);
  setParameter(DECAYBYFREQ,  2/3.f);      // 100% -> constant Q
  setParameter(RESOATTACK,   0.0);        // instant attack
  setParameter(RESOGAIN,     12./72.);    // +- 0 dB
  setParameter(AMPOFFSET,    0.0);        // 0
  setParameter(RESOLIMIT,    1.0);        // 10
  setParameter(RESORANGE,    1.0);        // 10
  setParameter(RESOPHASE,    0.0);        // no phase shift
  setParameter(RESOWAVE,     0.0);        // sine
  setParameter(WAVECUTOFF,   1.0);        // maximally open
  setParameter(WAVEMOD,      0.0);        // no modulation
  setParameter(CHAOSAMOUNT,  0.5);        // no chaos
  setParameter(CHAOSCUTOFF,  0.5);        // 100 Hz
  setParameter(CHAOSPARAM,   0.5);        // 0

  setParameter(SELFEXCITE,   0.0);        // 0

  setParameter(FBDRIVE,      0.0);        // 0 dB
  setParameter(FBLOLIMIT,    0.5);        // 1 
  setParameter(FBHILIMIT,    0.5);
  setParameter(FBSATGAINAT1, 1.0);        // hardclip
  setParameter(FBSATPLACE,   0.0);
}

bool vstResoReplacer::getEffectName(char* name)
{
  vst_strncpy(name, "ResoReplacer", kVstMaxEffectNameLen);
  return true;
}

void vstResoReplacer::processStereoFrame(double *inL, double *inR, double *outL, double *outR)
{
  // preliminary - mix input to mono, process, broadcast output to both channels (later introduce
  // stereo processing):
  double in  = 0.5 * (*inL + *inR);
  double out = filter.getSample(in);
  *outL = out;
  *outR = out;
}

void vstResoReplacer::resume()
{
  filter.reset();
}

void vstResoReplacer::setSampleRate(float sampleRate)
{
  filter.setSampleRate(sampleRate);
}

// Parameters:

void vstResoReplacer::updateCoreParameter(VstInt32 index, float value)
{
  switch(index)
  {
  case INGAIN: {
    inGain = rsLinToLin(params[INGAIN], 0, 1, -30, +60);
    filter.setInputGain(rsDB2amp(inGain));
  } break;
  case LEAK: {
    leak = rsLinToLin(params[LEAK], 0, 1, -4, +4);
    filter.setInputLeak(leak);
  } break;
  case FREQUENCY: {
    //freq = rsLinToExp(params[FREQUENCY], 0, 1, 20, 20000);
    freq = rsLinToExp(params[FREQUENCY], 0, 1, 1.0, 20000);
      //maybe use rsLinToExpWithOffset here, such that the low frequency range is squashed
    filter.setCutoff(freq);
  } break;
  case RESODECAY: {
    decay = rsLinToExp(params[RESODECAY], 0, 1, 0.1, 10000.0);       // milliseconds here
    filter.setResonanceDecay(0.001*decay);                           // seconds there
  } break;
  case DECAYBYFREQ: {
    decayByFreq = rsLinToLin(params[DECAYBYFREQ], 0, 1, -100, 200);  // % here
    filter.setDecayByFrequency(0.01*decayByFreq);                    // 0..1 there
  } break;
  case RESOATTACK: {
    attack = rsLinToLin(params[RESOATTACK], 0, 1, 0, 100);           // % here
    filter.setResonanceAttack(0.01*attack);                          // 0..1 there
      // maybe a mapping function should be applied - it feels nonuniform - in the lower range,
      // we need more precision and less in the upper...
  } break;
  case RESOGAIN: {
    gain = rsLinToLin(params[RESOGAIN], 0, 1, -12, 60);
    filter.setResonanceGain(rsDB2amp(gain));
  } break;
  case AMPOFFSET: {
    ampOffset = rsLinToLin(params[AMPOFFSET], 0, 1, 0, 1);
    filter.setAmplitudeOffset(ampOffset);
  } break;
  case RESOLIMIT: {
    limit = rsLinToExp(params[RESOLIMIT], 0, 1, 0.01, 10.0);
    filter.setAmplitudeLimit(limit);
  } break;
  case RESORANGE: {
    range = rsLinToExp(params[RESORANGE], 0, 1, 0.01, 10.0);
    filter.setInputRange(range);
  } break;
  case RESOPHASE: {
    phase = rsLinToLin(params[RESOPHASE], 0, 1, 0, 360);            // degrees here
    filter.setResonancePhase(phase*PI/180);                         // radians there
  } break;
  case RESOWAVE: {
    waveForm = (int)round(4*params[RESOWAVE]);
    filter.setResonanceWaveform(waveForm);
  } break;  
  case WAVECUTOFF: {
    waveCut = rsLinToExp(params[WAVECUTOFF], 0, 1, 0.1, 100.0);
    filter.setResoCutoffMultiplier(waveCut);
  } break;
  case WAVEMOD: {
    waveMod = rsLinToExp(params[WAVEMOD], 0, 1, 0.1, 100.0);
    filter.setResoCutoffModulation(waveMod);
  } break;
  case CHAOSAMOUNT: {
    chaosAmt = rsLinToLin(params[CHAOSAMOUNT], 0, 1, -100.0, 100.0);
    filter.setBumpFactor(chaosAmt);
  } break;
  case CHAOSCUTOFF: {
    chaosCut = rsLinToExp(params[CHAOSCUTOFF], 0, 1, 10.0, 1000.0);
    filter.setBumpCutoffs(chaosCut, 0.2*chaosCut); // the 2nd is ad-hoc
  } break;
  case CHAOSPARAM: {
    chaosPar = rsLinToLin(params[CHAOSPARAM], 0, 1, -2.0, 2.0);
    filter.setChaosParameter(chaosPar);
  } break;


  case SELFEXCITE: {
    //selfEx = rsLinToLin(params[SELFEXCITE], 0, 1, 0.0, 1.0);
    selfEx = rsLinToLin(params[SELFEXCITE], 0, 1, -100.0, 0.0);
    filter.setSelfExitation(rsDB2amp(selfEx));
  } break;

  case FBDRIVE: {
    fbDrive = rsLinToLin(params[FBDRIVE], 0, 1, 0, 24);
    filter.setFeedbackDrive(rsDB2amp(fbDrive));
  } break;
  case FBLOLIMIT: {
    fbLoLimit = -rsLinToExp(params[FBLOLIMIT], 0, 1, 0.01, 10.0);
    filter.setFeedbackLowerLimit(fbLoLimit);
  } break;
  case FBHILIMIT: {
    fbHiLimit = rsLinToExp(params[FBHILIMIT], 0, 1, 0.01, 10.0);
    filter.setFeedbackUpperLimit(fbHiLimit);
  } break;
  case FBSATGAINAT1: {
    fbGainAt1 = rsLinToLin(params[FBSATGAINAT1], 0, 1, 0.5, 1.0);
    filter.setFeedbackSaturationGainAt1(fbGainAt1);
  } break;
  case FBSATPLACE: {
    fbSatPlace = (int)round((rsLadderFilterFeedbackSaturated::NUM_SATPLACES-1)*params[FBSATPLACE]);
    filter.setFeedbackSaturationPlace(fbSatPlace);
  } break;

  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only space for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstResoReplacer::getParameterName(VstInt32 index, char* label)
{
  switch( index )
  {
  case INGAIN:       vst_strncpy(label, "Gain",     kVstMaxParamStrLen); break;
  case LEAK:         vst_strncpy(label, "Leak",     kVstMaxParamStrLen); break;
  case FREQUENCY:    vst_strncpy(label, "Cutoff",   kVstMaxParamStrLen); break;
  case RESODECAY:    vst_strncpy(label, "ResDec",   kVstMaxParamStrLen); break;
  case DECAYBYFREQ:  vst_strncpy(label, "-byFreq",  kVstMaxParamStrLen); break;
  case RESOATTACK:   vst_strncpy(label, "ResAtt",   kVstMaxParamStrLen); break;
  case RESOGAIN:     vst_strncpy(label, "ResGain",  kVstMaxParamStrLen); break;
  case AMPOFFSET:    vst_strncpy(label, " -Add",    kVstMaxParamStrLen); break;
  case RESOLIMIT:    vst_strncpy(label, "ResLimit", kVstMaxParamStrLen); break;
  case RESORANGE:    vst_strncpy(label, "ResRange", kVstMaxParamStrLen); break;
  case RESOPHASE:    vst_strncpy(label, "ResPhs",   kVstMaxParamStrLen); break;
  case RESOWAVE:     vst_strncpy(label, "ResWave",  kVstMaxParamStrLen); break;
  case WAVECUTOFF:   vst_strncpy(label, "-Cutoff",  kVstMaxParamStrLen); break;
  case WAVEMOD:      vst_strncpy(label, "--Mod",    kVstMaxParamStrLen); break;
  case CHAOSAMOUNT:  vst_strncpy(label, "Chaos",    kVstMaxParamStrLen); break;
  case CHAOSCUTOFF:  vst_strncpy(label, "-Cutoff",  kVstMaxParamStrLen); break;
  case CHAOSPARAM:   vst_strncpy(label, "-Param",   kVstMaxParamStrLen); break;

  case SELFEXCITE:   vst_strncpy(label, "Excite",   kVstMaxParamStrLen); break;

  case FBDRIVE:      vst_strncpy(label, "FBDrive",  kVstMaxParamStrLen); break;
  case FBLOLIMIT:    vst_strncpy(label, "FBLoLim",  kVstMaxParamStrLen); break;
  case FBHILIMIT:    vst_strncpy(label, "FBHiLim",  kVstMaxParamStrLen); break;
  case FBSATGAINAT1: vst_strncpy(label, "FBValAt1", kVstMaxParamStrLen); break;
  case FBSATPLACE:   vst_strncpy(label, "FBSatPos", kVstMaxParamStrLen); break;
  } 
}

void vstResoReplacer::getParameterDisplay(VstInt32 index, char* text)
{
  switch( index )
  {
  case INGAIN:       sprintf(text, "%.2f", inGain);             break;
  case LEAK:         sprintf(text, "%.2f", leak);               break;
  case FREQUENCY:    sprintf(text, "%.2f", freq);               break;
  case RESODECAY:    sprintf(text, "%.2f", decay);              break;
  case DECAYBYFREQ:  sprintf(text, "%.2f", decayByFreq);        break;
  case RESOATTACK:   sprintf(text, "%.2f", attack);             break;
  case RESOGAIN:     sprintf(text, "%.2f", gain);               break;
  case AMPOFFSET:    sprintf(text, "%.2f", ampOffset);          break;
  case RESOLIMIT:    sprintf(text, "%.2f", limit);              break;
  case RESORANGE:    sprintf(text, "%.2f", range);              break;
  case RESOPHASE:    sprintf(text, "%.2f", phase);              break;
  case RESOWAVE:     getWaveformName(waveForm, text);           break;
  case WAVECUTOFF:   sprintf(text, "%.2f", waveCut);            break;
  case WAVEMOD:      sprintf(text, "%.2f", waveMod);            break;
  case CHAOSAMOUNT:  sprintf(text, "%.2f", chaosAmt);           break;
  case CHAOSCUTOFF:  sprintf(text, "%.2f", chaosCut);           break;
  case CHAOSPARAM:   sprintf(text, "%.2f", chaosPar);           break;

  case SELFEXCITE:   sprintf(text, "%.2f", selfEx);             break;

  case FBDRIVE:      sprintf(text, "%.2f", fbDrive);            break;
  case FBLOLIMIT:    sprintf(text, "%.2f", fbLoLimit);          break;
  case FBHILIMIT:    sprintf(text, "%.2f", fbHiLimit);          break;
  case FBSATGAINAT1: sprintf(text, "%.2f", fbGainAt1);          break;
  case FBSATPLACE:   getFeedbackSatPlaceName(fbSatPlace, text); break;
  }
}

void vstResoReplacer::getParameterLabel(VstInt32 index, char* label)
{
  switch( index )
  {
  case INGAIN:      vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;
  case FREQUENCY:   vst_strncpy(label, "Hz",  kVstMaxParamStrLen);  break;
  case RESODECAY:   vst_strncpy(label, "ms",  kVstMaxParamStrLen);  break;
  case DECAYBYFREQ: vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  case RESOGAIN:    vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;
  case RESOATTACK:  vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  case RESOPHASE:   vst_strncpy(label, "deg", kVstMaxParamStrLen);  break;
  case CHAOSCUTOFF: vst_strncpy(label, "Hz",  kVstMaxParamStrLen);  break;

  case SELFEXCITE:  vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;

  case FBDRIVE:     vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;

  default:          vst_strncpy(label, "",    kVstMaxParamStrLen);  break;
  } 
}

void vstResoReplacer::getWaveformName(VstInt32 index, char* text)
{
  typedef rsLadderFilterFeedbackSaturated Flt;
  switch( index )
  {
  case 0:  vst_strncpy(text, "Sine",     kVstMaxParamStrLen);  break;
  case 1:  vst_strncpy(text, "Triangle", kVstMaxParamStrLen);  break;
  case 2:  vst_strncpy(text, "Square",   kVstMaxParamStrLen);  break;
  case 3:  vst_strncpy(text, "SawUp",    kVstMaxParamStrLen);  break;
  case 4:  vst_strncpy(text, "SawDown",  kVstMaxParamStrLen);  break;
  default: vst_strncpy(text, "Unknown",  kVstMaxParamStrLen);  break;
  }
}

void vstResoReplacer::getFeedbackSatPlaceName(VstInt32 index, char* text)
{
  typedef rsLadderFilterFeedbackSaturated Flt;
  switch( index )
  {
  case Flt::NOWHERE:           vst_strncpy(text, "Bypass",   kVstMaxParamStrLen);  break;
  case Flt::PRE_FB_GAIN:       vst_strncpy(text, "PreK",     kVstMaxParamStrLen);  break;
  case Flt::POST_FB_GAIN:      vst_strncpy(text, "PostK",    kVstMaxParamStrLen);  break;
  case Flt::POST_INPUT_ADD:    vst_strncpy(text, "PostAdd",  kVstMaxParamStrLen);  break;
  case Flt::POST_GAIN_AND_ADD: vst_strncpy(text, "kAndAdd",  kVstMaxParamStrLen);  break;

  case Flt::POST_1ST_STAGE:    vst_strncpy(text, "Post1",    kVstMaxParamStrLen);  break;
  case Flt::POST_2ND_STAGE:    vst_strncpy(text, "Post2",    kVstMaxParamStrLen);  break;
  case Flt::POST_3RD_STAGE:    vst_strncpy(text, "Post3",    kVstMaxParamStrLen);  break;
  case Flt::POST_4TH_STAGE:    vst_strncpy(text, "Post4",    kVstMaxParamStrLen);  break;
  case Flt::POST_EACH_STAGE:   vst_strncpy(text, "PostEach", kVstMaxParamStrLen);  break;

  case Flt::IN_1ST_STAGE:      vst_strncpy(text, "In1",      kVstMaxParamStrLen);  break;
  case Flt::IN_2ND_STAGE:      vst_strncpy(text, "In2",      kVstMaxParamStrLen);  break;
  case Flt::IN_3RD_STAGE:      vst_strncpy(text, "In3",      kVstMaxParamStrLen);  break;
  case Flt::IN_4TH_STAGE:      vst_strncpy(text, "In4",      kVstMaxParamStrLen);  break;
  case Flt::IN_EACH_STAGE:     vst_strncpy(text, "InEach",   kVstMaxParamStrLen);  break;
  }
}