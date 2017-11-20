#include "vstResoShaper.h"

//// Waveshaping functions (maybe move to RSLib - or maybe remove):
//
//double hardClip(double x)
//{
//  if(x > 1.0)
//    return 1.0;
//  if(x < -1.0)
//    return -1.0;
//  return x;
//}
//double softClip(double x)
//{
//  // this function is an identity function for 0.5 < x < 0.5, equal to +-1 for abs(x) > 1
//  // and a third order order polynomial in between, such that at the juction points, the
//  // function values and 1st derivatives are matched
//  double s = rsSign(x); // maybe use a more efficient sign-function
//  double a = rsAbs(x);
//  if(a < 0.5)
//    return x;
//  if(a > 1.0)
//    return s;
//  double a2 = a*a;
//  return s * (1 + 8*a2 - 4*(1+a2)*a);
//}
// move toRSLib::NormalizedSigmoids



// Initialization, Processing, Setup:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstResoShaper(audioMaster);
}

vstResoShaper::vstResoShaper(audioMasterCallback audioMaster) 
  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('RsSh');

  // set up intial settings:

  setParameter(INGAIN,       0.5);        // +- 0 dB
  setParameter(LEAK,         0.5);        // +- 0
  setParameter(FREQUENCY,    0.5);
  setParameter(RESODECAY,    0.25);
  setParameter(DECAYBYFREQ,  2/3.f);      // 100% -> constant Q
  setParameter(RESOATTACK,   0.0);        // instant attack
  setParameter(RESOGAIN,     0.5);        // +- 0 dB
  setParameter(RESOPHASE,    0.0);        // no phase shift
  setParameter(FBDRIVE,      0.0);        // 0 dB
  setParameter(FBLOLIMIT,    0.5);        // 1 
  setParameter(FBHILIMIT,    0.5);
  setParameter(GATE,         0.0);        // no gating
  setParameter(GATEMIX,      0.0);        // use input as gate-signal
  setParameter(FBSATGAINAT1, 1.0);        // hardclip
  setParameter(FBSATPLACE,   0.0);
  setParameter(SATGAINAT1,   1.0);        // hardclip
  setParameter(DRIVE,        0.25);       // +- 0 dB
  setParameter(DRIVE_COMP,   1/3.f);      // 0%
  setParameter(SATMODE,      0.0);        // bypass
  setParameter(SATADDCONST,  0.5);        // +- 0
  setParameter(SATADDIN,     0.5);        // +- 0
  setParameter(SATADDFLT,    0.5);        // +- 0
}

bool vstResoShaper::getEffectName(char* name)
{
  vst_strncpy(name, "ResoShaper", kVstMaxEffectNameLen);
  return true;
}

void vstResoShaper::processStereoFrame(double *inL, double *inR, double *outL, double *outR)
{
  // preliminary - mix input to mono, process, broadcast output to both channels (later introduce
  // stereo processing):
  double in  = 0.5 * (*inL + *inR);
  double out = filter.getSample(in);
  *outL = out;
  *outR = out;
}

void vstResoShaper::resume()
{
  filter.reset();
}

void vstResoShaper::setSampleRate(float sampleRate)
{
  filter.setSampleRate(sampleRate);
}

// Parameters:

void vstResoShaper::updateCoreParameter(VstInt32 index, float value)
{
  switch(index)
  {
  case INGAIN: {
    inGain = rsLinToLin(params[INGAIN], 0, 1, -30, +30);
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
    decay = rsLinToExp(params[RESODECAY], 0, 1, 0.1, 1000.0);       // milliseconds here
    filter.setResonanceDecay(0.001*decay);                          // seconds there
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
    gain = rsLinToLin(params[RESOGAIN], 0, 1, -12, 12);
    filter.setResonanceGain(rsDB2amp(gain));
  } break;
  case RESOPHASE: {
    phase = rsLinToLin(params[RESOPHASE], 0, 1, 0, 360);            // degrees here
    filter.setResonancePhase(phase*PI/180);                         // radians there
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


  case GATE: {
    gate = rsLinToLin(params[GATE], 0, 1, 0, 1);
    filter.setGateSensitivity(gate);
  } break;
  case GATEMIX: {
    gateMix = rsLinToLin(params[GATEMIX], 0, 1, 0, 1);
    filter.setGateMix(gateMix);
  } break;

  case FBSATGAINAT1: {
    fbGainAt1 = rsLinToLin(params[FBSATGAINAT1], 0, 1, 0.5, 1.0);
    filter.setFeedbackSaturationGainAt1(fbGainAt1);
  } break;
  case FBSATPLACE: {
    fbSatPlace = (int)round((rsLadderFilterFeedbackSaturated::NUM_SATPLACES-1)*params[FBSATPLACE]);
    filter.setFeedbackSaturationPlace(fbSatPlace);
  } break;
  case SATGAINAT1: {
    gainAt1 = rsLinToLin(params[SATGAINAT1], 0, 1, 0.5, 1.0);
    filter.setSaturationGainAt1(gainAt1);
  } break;
  case DRIVE: {
    drive = rsLinToLin(params[DRIVE], 0, 1, -12, 36);
    filter.setDrive(rsDB2amp(drive));
  } break;
  case DRIVE_COMP: {
    driveComp = rsLinToLin(params[DRIVE_COMP], 0, 1, -100, 200);     // % here
    filter.setDriveCompensation(0.01*driveComp);                     // 0..1 there
  } break;
  case SATMODE: {
    satMode = (int)round((rsSidechainSaturator::NUM_MODES-1) * params[SATMODE]);
    filter.setSaturationMode(satMode);
  } break;  
  case SATADDCONST: {
    addConst = rsLinToLin(params[SATADDCONST], 0, 1, -5, 5);
    filter.setSaturationAddConstant(addConst);
  } break;
  case SATADDIN: {
    addIn = rsLinToLin(params[SATADDIN], 0, 1, -5, 5);
    filter.setSaturationAddInput(addIn);
  } break;
  case SATADDFLT: {
    addFlt = rsLinToLin(params[SATADDFLT], 0, 1, -5, 5);
    filter.setSaturationAddFiltered(addFlt);
  } break;
  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only space for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstResoShaper::getParameterName(VstInt32 index, char* label)
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
  case RESOPHASE:    vst_strncpy(label, "ResPhs",   kVstMaxParamStrLen); break;
  case FBDRIVE:      vst_strncpy(label, "FBDrive",  kVstMaxParamStrLen); break;
  case FBLOLIMIT:    vst_strncpy(label, "FBLoLim",  kVstMaxParamStrLen); break;
  case FBHILIMIT:    vst_strncpy(label, "FBHiLim",  kVstMaxParamStrLen); break;
  case GATE:         vst_strncpy(label, "Gate",     kVstMaxParamStrLen); break;
  case GATEMIX:      vst_strncpy(label, "GateMix",  kVstMaxParamStrLen); break;
  case FBSATGAINAT1: vst_strncpy(label, "FBValAt1", kVstMaxParamStrLen); break;
  case FBSATPLACE:   vst_strncpy(label, "FBSatPos", kVstMaxParamStrLen); break;
  case SATGAINAT1:   vst_strncpy(label, "ValAt1",   kVstMaxParamStrLen); break;
  case DRIVE:        vst_strncpy(label, "Drive",    kVstMaxParamStrLen); break;
  case DRIVE_COMP:   vst_strncpy(label, "-Comp",    kVstMaxParamStrLen); break;
  case SATMODE:      vst_strncpy(label, "SatMode",  kVstMaxParamStrLen); break;
  case SATADDCONST:  vst_strncpy(label, "-Offset",  kVstMaxParamStrLen); break;
  case SATADDIN:     vst_strncpy(label, "-AddIn",   kVstMaxParamStrLen); break;
  case SATADDFLT:    vst_strncpy(label, "-AddFlt",  kVstMaxParamStrLen); break;
  } 
}

void vstResoShaper::getParameterDisplay(VstInt32 index, char* text)
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
  case RESOPHASE:    sprintf(text, "%.2f", phase);              break;
  case FBDRIVE:      sprintf(text, "%.2f", fbDrive);            break;
  case FBLOLIMIT:    sprintf(text, "%.2f", fbLoLimit);          break;
  case FBHILIMIT:    sprintf(text, "%.2f", fbHiLimit);          break;
  case GATE:         sprintf(text, "%.2f", gate);               break;
  case GATEMIX:      sprintf(text, "%.2f", gateMix);            break;
  case FBSATGAINAT1: sprintf(text, "%.2f", fbGainAt1);          break;
  case FBSATPLACE:   getFeedbackSatPlaceName(fbSatPlace, text); break;
  case SATGAINAT1:   sprintf(text, "%.2f", gainAt1);            break;
  case DRIVE:        sprintf(text, "%.2f", drive);              break;
  case DRIVE_COMP:   sprintf(text, "%.2f", driveComp);          break;
  case SATMODE:      getResoSatModeName(satMode, text);         break;
  case SATADDCONST:  sprintf(text, "%.2f", addConst);           break;
  case SATADDIN:     sprintf(text, "%.2f", addIn);              break;
  case SATADDFLT:    sprintf(text, "%.2f", addFlt);             break;
  }
}

void vstResoShaper::getParameterLabel(VstInt32 index, char* label)
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
  case FBDRIVE:     vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;
  case DRIVE:       vst_strncpy(label, "dB",  kVstMaxParamStrLen);  break;
  case DRIVE_COMP:  vst_strncpy(label, "%",   kVstMaxParamStrLen);  break;
  default:          vst_strncpy(label, "",    kVstMaxParamStrLen);  break;
  } 
}

void vstResoShaper::getResoSatModeName(VstInt32 index, char* text)
{
  typedef rsSidechainSaturator Sat;
  switch( index )
  {
  case Sat::BYPASS:         vst_strncpy(text, "Bypass",   kVstMaxParamStrLen);  break;
  case Sat::DIRECT1:        vst_strncpy(text, "Direct1",  kVstMaxParamStrLen);  break;
  case Sat::DIRECT2:        vst_strncpy(text, "Direct2",  kVstMaxParamStrLen);  break;
  case Sat::ADDITIVE:       vst_strncpy(text, "Add",      kVstMaxParamStrLen);  break;
  case Sat::MULTIPLICATIVE: vst_strncpy(text, "Mul",      kVstMaxParamStrLen);  break;
  }
}

void vstResoShaper::getFeedbackSatPlaceName(VstInt32 index, char* text)
{
  typedef rsLadderFilterFeedbackSaturated Flt;
  switch( index )
  {
  case Flt::NOWHERE:             vst_strncpy(text, "Bypass",    kVstMaxParamStrLen);  break;
  case Flt::PRE_FB_GAIN:         vst_strncpy(text, "PreK",      kVstMaxParamStrLen);  break;
  case Flt::POST_FB_GAIN:        vst_strncpy(text, "PostK",     kVstMaxParamStrLen);  break;
  case Flt::INPUT_AND_POST_GAIN: vst_strncpy(text, "InAndPost", kVstMaxParamStrLen);  break;

  case Flt::POST_INPUT_ADD:      vst_strncpy(text, "PostAdd",  kVstMaxParamStrLen);  break;
  case Flt::POST_GAIN_AND_ADD:   vst_strncpy(text, "kAndAdd",  kVstMaxParamStrLen);  break;

  case Flt::POST_1ST_STAGE:      vst_strncpy(text, "Post1",    kVstMaxParamStrLen);  break;
  case Flt::POST_2ND_STAGE:      vst_strncpy(text, "Post2",    kVstMaxParamStrLen);  break;
  case Flt::POST_3RD_STAGE:      vst_strncpy(text, "Post3",    kVstMaxParamStrLen);  break;
  case Flt::POST_4TH_STAGE:      vst_strncpy(text, "Post4",    kVstMaxParamStrLen);  break;
  case Flt::POST_EACH_STAGE:     vst_strncpy(text, "PostEach", kVstMaxParamStrLen);  break;

  case Flt::IN_1ST_STAGE:        vst_strncpy(text, "In1",      kVstMaxParamStrLen);  break;
  case Flt::IN_2ND_STAGE:        vst_strncpy(text, "In2",      kVstMaxParamStrLen);  break;
  case Flt::IN_3RD_STAGE:        vst_strncpy(text, "In3",      kVstMaxParamStrLen);  break;
  case Flt::IN_4TH_STAGE:        vst_strncpy(text, "In4",      kVstMaxParamStrLen);  break;
  case Flt::IN_EACH_STAGE:       vst_strncpy(text, "InEach",   kVstMaxParamStrLen);  break;
  }
}
