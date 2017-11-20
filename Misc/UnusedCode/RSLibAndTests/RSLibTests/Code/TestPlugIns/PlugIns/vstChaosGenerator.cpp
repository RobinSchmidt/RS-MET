#include "vstChaosGenerator.h"


//=================================================================================================
// MonoSynth:

MonoSynth::MonoSynth()
{
  sampleRate            = 44100.0;
  currentNote           = 0;
  glideSamples          = 0;
  remainingGlideSamples = 0; 
  glideTime             = 0.0;
  targetFrequency       = 1000.0;
  currentFrequency      = 1000.0;
  freqFactorPerSample   = 1.0; 
  detuneFactor          = 1.0;
}

void MonoSynth::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  setGlideTime(glideTime);
}

void MonoSynth::setGlideTime(double newTime)
{
  glideTime    = newTime;
  glideSamples = rsRoundToInt(glideTime * sampleRate);
}

void MonoSynth::setDetune(double newDetune)
{
  detuneFactor = rsPitchOffsetToFreqFactor(newDetune);
}

void MonoSynth::noteOn(int key, int velocity)
{
  if(velocity > 0)
  {
    currentNote = key;
    if( isSilent() )
      jumpTo((double)key);
    else
      glideTo((double)key); 
    triggerAttack(key, velocity);
  }
  else if( key == currentNote )
    triggerRelease();
}

void MonoSynth::glideTo(double pitch)
{
  targetFrequency       = rsPitchToFreq(pitch);
  remainingGlideSamples = glideSamples;
  if( remainingGlideSamples > 0 )
    freqFactorPerSample = pow(targetFrequency/currentFrequency, 1.0/remainingGlideSamples);
  else
  {
    freqFactorPerSample = 1.0;
    currentFrequency    = targetFrequency;
  }
}

void MonoSynth::jumpTo(double pitch)
{
  targetFrequency       = rsPitchToFreq(pitch);
  currentFrequency      = targetFrequency;
  freqFactorPerSample   = 1.0;
  remainingGlideSamples = 0;
}

double MonoSynth::getCurrentFrequency()
{
  if(remainingGlideSamples > 0)
  {
    currentFrequency *= freqFactorPerSample;
    remainingGlideSamples--;
  }
  else
    currentFrequency = targetFrequency; // get rid of accumulated rounding errors
  return detuneFactor * currentFrequency;
}

//=================================================================================================
// ChaosGenerator:

ChaosGenerator::ChaosGenerator()
{
  fbFilter.setMode(rsStateVariableFilter::BELL);

  //inFilter.setMode(rsLadderFilter::LP_6);
  inFilter.setMode(rsLadderFilter::LP_24);
  //inFilter.setMode(rsLadderFilter::HP_24);

  setSampleRate(44100.0);
  setFrequency(1000.0);
  setInputMix(0.0);
  setInputFreqScaler(1.0);
  setFrequencyRatio(1.0);
  setModulationDepthScaler(1.0);
  setModDepth12(0.0);
  setModDepth21(0.0);
  setFreqVsPhaseMod12(0.0);
  setFreqVsPhaseMod21(0.0);
  setModEnvDepth(0.0);

  setFilter1FreqScaler(0.1);
  setFilter2FreqScaler(10.0);
  setFilter1Mode(rsLadderFilter::HP_6);
  setFilter2Mode(rsLadderFilter::LP_6);

  setFeedbackFreqScaler(1.0);       // 1.0 * master-frequency
  setFeedbackBoostAmount(0.0);      // 0 dB
  setFeedbackBoostWidth(1.0);       // 1 octave wide
  setFreqOffset1(0.0);
  setFreqOffset2(0.0);
  setWaveform1(0);
  setWaveform2(0);
  setClipLevel1(2.0);
  setClipLevel2(2.0);
  setPhase1(0.0);
  setPhase2(0.0);

  // make sure that later breakpoints are shifted in time when attack/decay gets modified:
  ampEnv.setEditMode(rsBreakpointModulator::EDIT_WITH_SHIFT);
  modEnv.setEditMode(rsBreakpointModulator::EDIT_WITH_SHIFT);

  // adjust the shape for both envelopes:
  ampEnv.setAllBreakpointsShape(rsModBreakpoint::ANALOG);
  ampEnv.setAllBreakpointsShapeAmount(1.0);

  modEnv.setAllBreakpointsShape(rsModBreakpoint::ANALOG);
  //modEnv.setAllBreakpointsShape(rsModBreakpoint::LINEAR);
  modEnv.setAllBreakpointsShapeAmount(1.0);
    // maybe for the modulation, a linear shape is better suited? we may experiment with this

  setSyncMode(0);
  setSyncThreshold(1.0);
  setSyncAmount(1.0);
  setFilterResetOnSync(true);
  setOutputClipLevel(2.0);

  reset();
}

void ChaosGenerator::setSampleRate(double newSampleRate)
{
  MonoSynth::setSampleRate(newSampleRate);
  hpf.setSampleRate(sampleRate);
  lpf.setSampleRate(sampleRate);
  fbFilter.setSampleRate(sampleRate);
  inFilter.setSampleRate(sampleRate);
  ampEnv.setSampleRate(sampleRate);
  modEnv.setSampleRate(sampleRate);
}

void ChaosGenerator::setFrequency(double newFrequency)
{
  frequency = newFrequency;

  // update filter cutoffs (they track the master frequency):
  setInputFreqScaler(   inScl);
  setFilter1FreqScaler( hpScl);
  setFilter2FreqScaler( lpScl);
  setFeedbackFreqScaler(fbScl);
}

void ChaosGenerator::setInputMix(double newMixFactor)
{
  inMix = newMixFactor;
}

void ChaosGenerator::setInputFreqScaler(double newScaler)
{
  inScl = newScaler;
  inFilter.setCutoff(rsLimitToRange(inScl * frequency, 0.0, 20000.0));
}

void ChaosGenerator::setInputFilterDetune(double newDetune)
{
  setInputFreqScaler(rsPitchOffsetToFreqFactor(newDetune));
}

void ChaosGenerator::setFrequencyRatio(double newRatio)
{
  freqRatio = newRatio;
}

void ChaosGenerator::setModulationDepthScaler(double newDepth)
{
  modScale = newDepth;
}
void ChaosGenerator::setModDepth12(double newDepth)
{
  mod12 = newDepth;
}
void ChaosGenerator::setModDepth21(double newDepth)
{
  mod21 = newDepth;
}
void ChaosGenerator::setFreqVsPhaseMod12(double newValue)
{
  fmpm12 = newValue;
}
void ChaosGenerator::setFreqVsPhaseMod21(double newValue)
{
  fmpm21 = newValue;
}

void ChaosGenerator::setModEnvDepth(double newDepth)
{
  modEnvDepth = newDepth;
}
void ChaosGenerator::setModEnvAttack(double newAttack)
{
  modEnv.setBreakpointTime(1, 0.001*newAttack);
}
void ChaosGenerator::setModEnvDecay(double newDecay)
{
  double time = modEnv.getBreakpointTime(1) + 0.001*newDecay; 
  modEnv.setBreakpointTime(2, time);
}
void ChaosGenerator::setModEnvSustain(double newSustain)
{
  modEnv.setBreakpointLevel(2, newSustain);
  modEnv.setBreakpointLevel(3, newSustain);
}
void ChaosGenerator::setModEnvRelease(double newRelease)
{
  double time = modEnv.getBreakpointTime(3) + 0.001*newRelease;
  modEnv.setBreakpointTime(4, time);
}

void ChaosGenerator::setFilter1FreqScaler(double newScaler)
{
  hpScl = newScaler;
  hpf.setCutoff(rsLimitToRange(hpScl * frequency, 0.0, 20000.0));
}
void ChaosGenerator::setFilter2FreqScaler(double newScaler)
{
  lpScl = newScaler;
  lpf.setCutoff(rsLimitToRange(lpScl * frequency, 0.0, 20000.0));
}

void ChaosGenerator::setFilter1Detune(double newDetune)
{
  setFilter1FreqScaler(rsPitchOffsetToFreqFactor(newDetune));
}
void ChaosGenerator::setFilter2Detune(double newDetune)
{
  setFilter2FreqScaler(rsPitchOffsetToFreqFactor(newDetune));
}

void ChaosGenerator::setFilter1Resonance(double newReso)
{
  hpf.setResonance(newReso);
}
void ChaosGenerator::setFilter2Resonance(double newReso)
{
  lpf.setResonance(newReso);
}

void ChaosGenerator::setFilter1Mode(int newMode)
{
  hpf.setMode(newMode);
}
void ChaosGenerator::setFilter2Mode(int newMode)
{
  lpf.setMode(newMode);
}

void ChaosGenerator::setFeedbackFreqScaler(double newScaler)
{
  fbScl = newScaler;
  fbFilter.setFrequency(rsLimitToRange(fbScl * frequency, 0.0, 20000.0));
}

void ChaosGenerator::setFeedbackDetune(double newDetune)
{
  setFeedbackFreqScaler(rsPitchOffsetToFreqFactor(newDetune));
}

void ChaosGenerator::setFeedbackBoostAmount(double newAmount)
{
  fbFilter.setGain(rsDB2amp(newAmount));
}

void ChaosGenerator::setFeedbackBoostWidth(double newWidth)
{
  fbFilter.setBandwidth(newWidth);
}

void ChaosGenerator::setGain(double newGain)
{
  gain = newGain;
}


void ChaosGenerator::setFreqOffset1(double newOffset)
{
  freqOffset1 = newOffset;
}
void ChaosGenerator::setFreqOffset2(double newOffset)
{
  freqOffset2 = newOffset;
}
void ChaosGenerator::setWaveform1(int newWaveform)
{
  shape1 = newWaveform;
}
void ChaosGenerator::setWaveform2(int newWaveform)
{
  shape2 = newWaveform;
}
void ChaosGenerator::setClipLevel1(double newLevel)
{
  clip1 = newLevel;
}
void ChaosGenerator::setClipLevel2(double newLevel)
{
  clip2 = newLevel;
}
void ChaosGenerator::setPhase1(double newPhase)
{
  phase1 = newPhase;
}
void ChaosGenerator::setPhase2(double newPhase)
{
  phase2 = newPhase;
}

void ChaosGenerator::setAmpEnvAttack(double newAttack)
{
  ampEnv.setBreakpointTime(1, 0.001*newAttack);
}
void ChaosGenerator::setAmpEnvDecay(double newDecay)
{
  double time = ampEnv.getBreakpointTime(1) + 0.001*newDecay; 
  ampEnv.setBreakpointTime(2, time);
}
void ChaosGenerator::setAmpEnvSustain(double newSustain)
{
  ampEnv.setBreakpointLevel(2, newSustain);
  ampEnv.setBreakpointLevel(3, newSustain);
}
void ChaosGenerator::setAmpEnvRelease(double newRelease)
{
  double time = ampEnv.getBreakpointTime(3) + 0.001*newRelease;
  ampEnv.setBreakpointTime(4, time);
}

void ChaosGenerator::setSyncMode(int newMode)
{
  syncMode = newMode;
}
void ChaosGenerator::setSyncThreshold(double newThreshold)
{
  syncThresh = newThreshold;
}
void ChaosGenerator::setSyncAmount(double newAmount)
{
  syncAmount = newAmount;
}
void ChaosGenerator::setFilterResetOnSync(bool shouldReset)
{
  syncFlt = shouldReset;
}

void ChaosGenerator::setDuckingRange(double newRange)
{
  duckWidth = 2*newRange;
  duckLo = duckCenter - 0.5 * duckWidth;
  duckHi = duckCenter + 0.5 * duckWidth;
  ducker.setWidth(duckWidth);
}

void ChaosGenerator::setDuckingFlatness(double newFlatness)
{
  ducker.setFlatTopWidth(newFlatness);
}

void ChaosGenerator::setDuckingCenter(double newCenter)
{
  duckCenter = newCenter;
  duckLo = duckCenter - 0.5 * duckWidth;
  duckHi = duckCenter + 0.5 * duckWidth;
  ducker.setCenter(duckCenter);
}

void ChaosGenerator::setDuckingShape(int newShape)
{
  switch(newShape)
  {
  case 0: ducker.setPrototypeBell(rsPositiveBellFunctions::linear);  break;
  case 1: ducker.setPrototypeBell(rsPositiveBellFunctions::cubic);   break;
  case 2: ducker.setPrototypeBell(rsPositiveBellFunctions::quintic); break;
  case 3: ducker.setPrototypeBell(rsPositiveBellFunctions::heptic);  break;
  }
}

void ChaosGenerator::setOutputClipLevel(double newLevel)
{
  clipOut = newLevel;
}

void ChaosGenerator::processSampleFrame(double *inL, double *inR, double *outL, double *outR)
{
  double freq, inc, tmp, env, fm, pm, mod, modScaler, flt, duck;

  // pitch glide:
  tmp = getCurrentFrequency();
  if(tmp != frequency)
    setFrequency(tmp); // costly call, therefore the conditional

  // mix/process input:
  tmp  = 0.5 * (*inL + *inR);              // mix input to mono
  flt  = inFilter.getSample(tmp);          // lowpass filter input
  duck = ducker.getValue(tmp);             // gain multiplier from unfiltered input   

  // sync:
  if(needsSyncTrigger(tmp, inOld)) 
    triggerSync();
  inOld = tmp;                             // remember current input for next call
    // if we want to change back to use the filtered input for ducking, we should pass the current
    // and previous filtered input to the needsSyncTrigger() function, too


  // old:
  //tmp = 0.5 * (*inL + *inR);        // mix input to mono
  //if(fabs(tmp-inOld) > syncThresh)  // sync to input edges
  //  triggerSync();
  //inOld = tmp;                      // remember current input for next call
  //flt   = inFilter.getSample(tmp);  // lowpass filter input

  ////duck  = ducker.getValue(flt);     // gain multiplier from lowpassed input
  //duck  = ducker.getValue(tmp);     // gain multiplier from unfiltered input
  //if(duck == 0.0)                   // sync to microenvelope start
  //  triggerSync();

  // read modulation envelope:
  env = modEnv.getSample();
  modScaler = modScale * (1 + modEnvDepth*env);

  // establish feedback signal:
  tmp = fbFilter.getSample(out2);

  // oscillator 1:
  mod   = modScaler  * mod21;                                // scaled, enveloped modulation
  fm    = cos(0.5*PI*fmpm21) * mod;                          // frequency modulation
  pm    = sin(0.5*PI*fmpm21) * mod;                          // phase modulation
  freq  = freqOffset1 + frequency*freqRatio * (1 + fm*tmp);  // frequency in Hz
  inc   = direction * 2*PI*freq / sampleRate;                // phase increment
  pos1 += inc;
  while(pos1 > 2*PI)
    pos1 -= 2*PI;
  while(pos1 < 0.0)
    pos1 += 2*PI;
  out1 = getWaveform(pos1+phase1+pm*tmp, shape1); 
  out1 = clip(out1, clip1); // can also be optimized

  // filters:
  tmp = out1;
  tmp = hpf.getSample(out1);               // highpass 
  tmp = lpf.getSample(tmp);                // lowpass

  // oscillator 2:
  mod   = modScaler  * mod12;
  fm    = cos(0.5*PI*fmpm12) * mod;
  pm    = sin(0.5*PI*fmpm12) * mod;
  freq  = freqOffset2 + frequency * (1 + fm*tmp);
  inc   = direction * 2*PI*freq / sampleRate;
  pos2 += inc;
  while(pos2 > 2*PI)
    pos2 -= 2*PI;
  while(pos2 < 0.0)
    pos2 += 2*PI;
  out2 = getWaveform(pos2+phase2+pm*tmp, shape2);
  out2 = clip(out2, clip2); 

  // \todo:
  // -precompute cos(0.5*PI*fmpm21), sin(0.5*PI*fmpm21), cos(0.5*PI*fmpm12), sin(0.5*PI*fmpm12)
  //  whenever a parameter that affects their value changes - we don't need to recompute them
  //  at each sample instant
  // -maybe factor out the increment/wraparound code (avoid code duplication)

  // create filtered input, read amplitude envelope and assign output slots:
  env   = duck * ampEnv.getSample();
  *outL = gain * clip(env * out1, clipOut) + inMix * flt;
  *outR = gain * clip(env * out2, clipOut) + inMix * flt;
    // maybe have a mono output that is an adjustable mix between out1 and out2
    // compute duck*env*gain, inMix*tmp just once


  // a strange observation: when just dialing in some 1->2 modulation, there's a kind of
  // attack phase (when the filters are reset on note-on) i think, this comes from the highpass
  // filter. the attack is longer, the lower the cutoff frequency of the highpass - i think
  // this is because the impulse response becomes longer with lower cutoffs
}

void ChaosGenerator::reset()
{
  resetOscillators();
  resetFilters();
  direction = +1.0;   // forward
  //currentNote = 0;  // interferes with sync?
}

void ChaosGenerator::resetOscillators()
{
  pos1 = pos2 = out1 = out2 = 0.0;
}

void ChaosGenerator::resetFilters()
{
  hpf.reset();
  lpf.reset();
  fbFilter.reset();
}

bool ChaosGenerator::needsSyncTrigger(double in, double inOld)
{       
  // sync to input edges:
  if(fabs(in-inOld) > syncThresh)
    return true;

  // sync, when input enters allowed amplitude range for ducker:
  if(inOld < duckLo && in >= duckLo)   // ...from below
    return true;
  if(inOld > duckHi && in <= duckHi)   // ...from above
    return true;

  return false;
}

void ChaosGenerator::triggerSync()
{
  syncOscillators();
  if(syncFlt)
    syncFilters();
}

void ChaosGenerator::syncOscillators()
{
  double s = syncAmount;
  double t = 1-s;

  //// old:
  //pos1 = t * pos1;            // s * 0.0 + t * pos1
  //pos2 = t * pos2;
  //out1 = t * out1;
  //out2 = t * out2;

  switch(syncMode)
  {
  case SYNC_RESET:
  {
    pos1 = t * pos1;            // s * 0.0 + t * pos1
    pos2 = t * pos2;
    //out1 = t * out1;            // does this make sense?
    //out2 = t * out2;
  }; 
  break;
  case SYNC_JUMP:
  {
    pos1 += s * 2*PI;
    pos2 += s * 2*PI;
  }; 
  break;
  case SYNC_REVERSE:
  {
    direction *= -1;
  }; 
  break;

  };
}

void ChaosGenerator::syncFilters()
{
   resetFilters(); // preliminary - maybe later we can have a partial reset, too
}

double ChaosGenerator::clip(double x, double level)
{
  return level * rsNormalizedSigmoids::softClipHexic(x / level);
}

double ChaosGenerator::getWaveform(double phase, int shape)
{
  // we need this here too because the phase-modulation +pm*tmp may bump it out of range in the 
  // call getWaveform(pos1+phase1+pm*tmp, shape1);
  while(phase > 2*PI)
    phase -= 2*PI;
  while(phase < 0.0)
    phase += 2*PI;

  switch(shape)
  {
  case 1:  return  rsTriWave(phase);
  case 2:  return  rsSqrWave(phase);
  case 3:  return  rsSawWave(phase);
  case 4:  return -rsSawWave(phase);
  default: return  sin(phase);
  }
}

bool ChaosGenerator::isSilent()
{
  return ampEnv.endIsReached;
}

void ChaosGenerator::triggerAttack(int key, int velocity)
{
  ampEnv.noteOn(true,  key, velocity);
  modEnv.noteOn(false, key, velocity);
}

void ChaosGenerator::triggerRelease()
{
  ampEnv.noteOff();
  modEnv.noteOff(); 
}

//=================================================================================================
// VST PlugIn:

AudioEffect* createEffectInstance(audioMasterCallback audioMaster)  
{
  return new vstChaosGenerator(audioMaster);
}

//vstChaosGenerator::vstChaosGenerator(audioMasterCallback audioMaster) 
//  : vstPlugIn(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
vstChaosGenerator::vstChaosGenerator(audioMasterCallback audioMaster) 
  : vstInstrument(audioMaster, defaultNumPrograms, NUM_PARAMETERS)
{
  setUniqueID('ChGn');

  setParameter(TUNE,            0.5);
  setParameter(GLIDE_TIME,      0.0);
  setParameter(IN_MIX,          0.5);
  setParameter(IN_DETUNE,       1.0);      // filter totally open
  setParameter(FREQ_RATIO,      0.5);

  setParameter(MOD_SCALE,       1.0);
  setParameter(MOD_12,          0.5);
  setParameter(FMPM_12,         0.0);
  setParameter(MOD_21,          0.5);
  setParameter(FMPM_21,         0.0);
  setParameter(MOD_ENV,         0.5);
  setParameter(MOD_ATT,         0.1f);
  setParameter(MOD_DEC,         0.3f);
  setParameter(MOD_SUS,         0.5f);
  setParameter(MOD_REL,         1.0);

  setParameter(FLT1_SCALE, 0.0);
  //setParameter(FLT1_SCALE, 0.1);
  setParameter(FLT1_RESO,  0.0);
  setParameter(FLT1_MODE, (rsLadderFilter::HP_6+0.5f) / rsLadderFilter::NUM_MODES);
  setParameter(FLT2_SCALE, 1.0);
  //setParameter(FLT2_SCALE, 0.9);
  setParameter(FLT2_RESO,  0.0);
  setParameter(FLT2_MODE, (rsLadderFilter::LP_6+0.5f) / rsLadderFilter::NUM_MODES);

  setParameter(FB_DETUNE,       0.5);
  setParameter(FB_BOOST,        0.5);
  setParameter(FB_WIDTH,        0.5);

  setParameter(FREQ_OFFSET1,    0.5);
  setParameter(FREQ_OFFSET2,    0.5);
  setParameter(WAVE1,           0.0);
  setParameter(WAVE2,           0.0);
  setParameter(CLIP1,           1.0);
  setParameter(CLIP2,           1.0);
  setParameter(PHASE1,          0.0);
  setParameter(PHASE2,          0.0);

  setParameter(VOLUME,          0.5);
  setParameter(AMP_ATT,         0.1f);
  setParameter(AMP_DEC,         0.1f);
  setParameter(AMP_SUS,         1.0);
  setParameter(AMP_REL,         0.1f);
  setParameter(CLIP_OUT,        1.0);

  setParameter(SYNC_MODE,       0.0);
  setParameter(SYNC_THRESH,     1.0);
  setParameter(SYNC_AMOUNT,     1.0);
  setParameter(SYNC_FILTERS,    1.0);

  setParameter(DUCK_RANGE,      1.0);
  setParameter(DUCK_FLATNESS,   0.5);    // 0.5
  setParameter(DUCK_CENTER,     0.5);    // 0
  setParameter(DUCK_SHAPE,      0.0);    // linear

  //setParameter(CLIP_LEVEL,      1.0);
  //setParameter(CLIP_POSITION,   0.0);

  resume();
}

bool vstChaosGenerator::getEffectName(char* name)
{
  vst_strncpy(name, "ChaosGenerator(2016-09-21)", kVstMaxEffectNameLen);
  return true;
}

void vstChaosGenerator::processStereoFrame(double *inL, double *inR, double *outL, double *outR)
{
  chaosGen.processSampleFrame(inL, inR, outL, outR);
}

void vstChaosGenerator::resume()
{
  chaosGen.reset();
}

void vstChaosGenerator::setSampleRate(float sampleRate)
{
  chaosGen.setSampleRate(sampleRate);
}

// Parameters:
void vstChaosGenerator::updateCoreParameter(VstInt32 index, float value)
{
  switch(index)
  {
  case TUNE:
  {
    tune = rsLinToLin(value, 0, 1, -36, 36);
    chaosGen.setDetune(tune);
  } break;
  case GLIDE_TIME:
  {
    glide_time = rsLinToLin(value, 0, 1, 0.0, 500.0);
    chaosGen.setGlideTime(0.001*glide_time);
  } break;
  case IN_MIX:
  {
    in_mix = rsLinToLin(value, 0, 1, -8.0, 8.0);
    chaosGen.setInputMix(in_mix);
  } break;
  case IN_DETUNE:
  {
    in_det = rsLinToLin(value, 0, 1, -24, 120.0);
    chaosGen.setInputFilterDetune(in_det);
  } break;
  case FREQ_RATIO:
  {
    freq_ratio = rsLinToExp(value, 0, 1, 1./64, 64);
    chaosGen.setFrequencyRatio(freq_ratio);
  } break;

  case MOD_SCALE:
  {
    mod_scale = rsLinToLin(value, 0, 1, -1.0, 1.0);
    chaosGen.setModulationDepthScaler(mod_scale);
  } break;
  case MOD_12:
  {
    mod_12 = rsLinToLin(value, 0, 1, -20.0, 20.0);
    chaosGen.setModDepth12(mod_12);
  } break;
  case FMPM_12:
  {
    fmpm12 = value;
    chaosGen.setFreqVsPhaseMod12(fmpm12);
  } break;
  case MOD_21:
  {
    mod_21 = rsLinToLin(value, 0, 1, -20.0, 20.0);
    chaosGen.setModDepth21(mod_21);
  } break;
  case FMPM_21:
  {
    fmpm21 = value;
    chaosGen.setFreqVsPhaseMod21(fmpm21);
  } break;
  case MOD_ENV:
  {
    mod_env = rsLinToLin(value, 0, 1, -400.0, 400.0);
    chaosGen.setModEnvDepth(0.01*mod_env);
  } break;
  case MOD_ATT:
  {
    mod_att = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setModEnvAttack(mod_att);
  } break;
  case MOD_DEC:
  {
    mod_dec = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setModEnvDecay(mod_dec);
  } break;
  case MOD_SUS:
  {
    mod_sus = rsLinToLin(value, 0, 1, 0.0, 1.0);
    chaosGen.setModEnvSustain(mod_sus);
  } break;
  case MOD_REL:
  {
    mod_rel = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setModEnvRelease(mod_rel);
  } break;

  case FLT1_SCALE:
  {
    fltScl1 = rsLinToExp(value, 0, 1, 0.01, 100.0);
    //fltScl1 = rsLinToExpWithOffset(value, 0, 1, 0.01, 100.0, -0.01);
    chaosGen.setFilter1FreqScaler(fltScl1); 
    // allow for true 0 cutoff frequency scaler which effectively turns the highpass off
    // ..but maybe that's not good because it will produce DC at 0 cutoff (the lowpass stages will 
    // just output their stored state value again and again (the feedback coefficient becomes 
    // unity)
  } break;
  case FLT1_RESO:
  {
    fltRes1 = value;
    chaosGen.setFilter1Resonance(fltRes1);
  } break;
  case FLT1_MODE:
  {
    fltMode1 =  rsRoundToInt(value * (rsLadderFilter::NUM_MODES-1));
    chaosGen.setFilter1Mode(fltMode1);
  } break;
  case FLT2_SCALE:
  {
    fltScl2 = rsLinToExp(          value, 0, 1, 0.01, 100.0);
    //fltScl2 = rsLinToExpWithOffset(value, 0, 1, 0.01, 100.0, -0.01);
    chaosGen.setFilter2FreqScaler(fltScl2);
  } break;
  case FLT2_RESO:
  {
    fltRes2 = value;
    chaosGen.setFilter2Resonance(fltRes2);
  } break;
  case FLT2_MODE:
  {
    fltMode2 =  rsRoundToInt(value * (rsLadderFilter::NUM_MODES-1));
    chaosGen.setFilter2Mode(fltMode2);
  } break;

  case FB_DETUNE:
  {
    fb_det = rsLinToLin(value, 0, 1, -48, 48);
    chaosGen.setFeedbackDetune(fb_det);
  } break;
  case FB_BOOST:
  {
    fb_boost = rsLinToLin(value, 0, 1, -100.0, 100.0);
    chaosGen.setFeedbackBoostAmount(fb_boost);
  } break;
  case FB_WIDTH:
  {
    fb_width = rsLinToExp(value, 0, 1, 0.25, 4.0);
    chaosGen.setFeedbackBoostWidth(fb_width);
  } break;


  case FREQ_OFFSET1:
  {
    freq_offset1 = rsLinToLin(value, 0, 1, -2000.0, 2000.0);
    chaosGen.setFreqOffset1(freq_offset1);
  } break;
  case FREQ_OFFSET2:
  {
    freq_offset2 = rsLinToLin(value, 0, 1, -2000.0, 2000.0);
    chaosGen.setFreqOffset2(freq_offset2);
  } break;
  case WAVE1:
  {
    wave1 = rsRoundToInt(value * (ChaosGenerator::NUM_WAVEFORMS-1));
    chaosGen.setWaveform1(wave1);
  } break;
  case WAVE2:
  {
    wave2 = rsRoundToInt(value * (ChaosGenerator::NUM_WAVEFORMS-1));
    chaosGen.setWaveform2(wave2);
  } break;
  case CLIP1:
  {
    clip1 = rsLinToLin(value, 0, 1, 0.0, 2.0);
    chaosGen.setClipLevel1(clip1);
  } break;
  case CLIP2:
  {
    clip2 = rsLinToLin(value, 0, 1, 0.0, 2.0);
    chaosGen.setClipLevel2(clip2);
  } break;
  case PHASE1:
  {
    phase1 = rsLinToLin(value, 0, 1, 0.0, 360.0);
    chaosGen.setPhase1(PI*phase1/180.0);
  } break;
  case PHASE2:
  {
    phase2 = rsLinToLin(value, 0, 1, 0.0, 360.0);
    chaosGen.setPhase2(PI*phase2/180.0);
  } break;

  case VOLUME:
  {
    volume = value;
    chaosGen.setGain(volume);
  } break;
  case AMP_ATT:
  {
    amp_att = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setAmpEnvAttack(amp_att);
  } break;
  case AMP_DEC:
  {
    amp_dec = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setAmpEnvDecay(amp_dec);
  } break;
  case AMP_SUS:
  {
    amp_sus = rsLinToLin(value, 0, 1, 0.0, 1.0);
    chaosGen.setAmpEnvSustain(amp_sus);
  } break;
  case AMP_REL:
  {
    amp_rel = rsLinToLin(value, 0, 1, 0.0, 3000.0);
    chaosGen.setAmpEnvRelease(amp_rel);
  } break;
  case CLIP_OUT:
  {
    clipOut = rsLinToLin(value, 0, 1, 0.0, 2.0);
    chaosGen.setOutputClipLevel(clipOut);
  } break;

  case SYNC_MODE:
  {
    syncMode = rsRoundToInt(value * (ChaosGenerator::NUM_SYNCMODES-1));
    chaosGen.setSyncMode(syncMode);
  } break;
  case SYNC_THRESH:
  {
    sync_thresh = rsLinToExp(value, 0, 1, 0.001, 2.0);
    chaosGen.setSyncThreshold(sync_thresh);
  } break;
  case SYNC_AMOUNT:
  {
    syncAmount = value;
    chaosGen.setSyncAmount(syncAmount);
  } break;
  case SYNC_FILTERS:
  {
    sync_filters = value >= 0.5;
    chaosGen.setFilterResetOnSync(sync_filters);
  } break;

  case DUCK_RANGE:
  {
    duckRange = rsLinToExp(value, 0, 1, 0.001, 2.0);
    chaosGen.setDuckingRange(duckRange);
  } break;
  case DUCK_FLATNESS:
  {
    duckFlat = value;
    chaosGen.setDuckingFlatness(duckFlat);
  } break;
  case DUCK_CENTER:
  {
    duckCenter = rsLinToLin(value, 0, 1, -1.0, 1.0);
    chaosGen.setDuckingCenter(duckCenter);
  } break;
  case DUCK_SHAPE:
  {
    duckShape = rsRoundToInt(value * 3);
    chaosGen.setDuckingShape(duckShape);
  } break;

  }
}

// beware - according to the VST-spec, the char-arrays passed to the following functions may have
// only spec for 8 characters (...including or excluding the terminating zero?)
// - write some utility functions vstPrint, vstFloat2Str, etc. that take care of this

void vstChaosGenerator::getParameterName(VstInt32 index, char* label)
{
  switch( index )
  {
  case TUNE:           vst_strncpy(label, "Tune",      kVstMaxParamStrLen); break;
  case GLIDE_TIME:     vst_strncpy(label, "Glide",     kVstMaxParamStrLen); break;
  case IN_MIX:         vst_strncpy(label, "InMix",     kVstMaxParamStrLen); break;
  case IN_DETUNE:      vst_strncpy(label, "-Detune",   kVstMaxParamStrLen); break;
  case FREQ_RATIO:     vst_strncpy(label, "Ratio",     kVstMaxParamStrLen); break;

  case MOD_SCALE:      vst_strncpy(label, "Mod Scale", kVstMaxParamStrLen); break;
  case MOD_12:         vst_strncpy(label, "Mod 1->2",  kVstMaxParamStrLen); break;
  case FMPM_12:        vst_strncpy(label, "-FM/PM",    kVstMaxParamStrLen); break;
  case MOD_21:         vst_strncpy(label, "Mod 2->1",  kVstMaxParamStrLen); break;
  case FMPM_21:        vst_strncpy(label, "-FM/PM",    kVstMaxParamStrLen); break;
  case MOD_ENV:        vst_strncpy(label, "ModEnv",    kVstMaxParamStrLen); break;
  case MOD_ATT:        vst_strncpy(label, "-Att",      kVstMaxParamStrLen); break;
  case MOD_DEC:        vst_strncpy(label, "-Dec",      kVstMaxParamStrLen); break;
  case MOD_SUS:        vst_strncpy(label, "-Sus",      kVstMaxParamStrLen); break;
  case MOD_REL:        vst_strncpy(label, "-Rel",      kVstMaxParamStrLen); break;

  case FLT1_SCALE:     vst_strncpy(label, "Cutoff1",   kVstMaxParamStrLen); break;
  case FLT1_RESO:      vst_strncpy(label, "-Reso",     kVstMaxParamStrLen); break;
  case FLT1_MODE:      vst_strncpy(label, "-Mode",     kVstMaxParamStrLen); break;
  case FLT2_SCALE:     vst_strncpy(label, "Cutoff2",   kVstMaxParamStrLen); break;
  case FLT2_RESO:      vst_strncpy(label, "-Reso",     kVstMaxParamStrLen); break;
  case FLT2_MODE:      vst_strncpy(label, "-Mode",     kVstMaxParamStrLen); break;
  case FB_DETUNE:      vst_strncpy(label, "FB Det",    kVstMaxParamStrLen); break;
  case FB_BOOST:       vst_strncpy(label, "FB Boost",  kVstMaxParamStrLen); break;
  case FB_WIDTH:       vst_strncpy(label, "FB Width",  kVstMaxParamStrLen); break;

  case FREQ_OFFSET1:  vst_strncpy(label, "dFreq1",    kVstMaxParamStrLen); break;
  case FREQ_OFFSET2:  vst_strncpy(label, "dFreq2",    kVstMaxParamStrLen); break;
  case WAVE1:         vst_strncpy(label, "Wave1",     kVstMaxParamStrLen); break;
  case WAVE2:         vst_strncpy(label, "Wave2",     kVstMaxParamStrLen); break;
  case CLIP1:         vst_strncpy(label, "Clip1",     kVstMaxParamStrLen); break;
  case CLIP2:         vst_strncpy(label, "Clip2",     kVstMaxParamStrLen); break;
  case PHASE1:        vst_strncpy(label, "Phase1",    kVstMaxParamStrLen); break;
  case PHASE2:        vst_strncpy(label, "Phase2",    kVstMaxParamStrLen); break;

  case VOLUME:        vst_strncpy(label, "Volume",    kVstMaxParamStrLen); break;
  case AMP_ATT:       vst_strncpy(label, "-Att",      kVstMaxParamStrLen); break;
  case AMP_DEC:       vst_strncpy(label, "-Dec",      kVstMaxParamStrLen); break;
  case AMP_SUS:       vst_strncpy(label, "-Sus",      kVstMaxParamStrLen); break;
  case AMP_REL:       vst_strncpy(label, "-Rel",      kVstMaxParamStrLen); break;
  case CLIP_OUT:      vst_strncpy(label, "-Clip",     kVstMaxParamStrLen); break;

  case SYNC_MODE:     vst_strncpy(label, "SyncMode",  kVstMaxParamStrLen); break;
  case SYNC_THRESH:   vst_strncpy(label, "-Thr",      kVstMaxParamStrLen); break;
  case SYNC_AMOUNT:   vst_strncpy(label, "-Amt",      kVstMaxParamStrLen); break;
  case SYNC_FILTERS:  vst_strncpy(label, "-Flt",      kVstMaxParamStrLen); break;

  case DUCK_RANGE:    vst_strncpy(label, "DuckLim",   kVstMaxParamStrLen); break;
  case DUCK_FLATNESS: vst_strncpy(label, "-flat",     kVstMaxParamStrLen); break;
  case DUCK_CENTER:   vst_strncpy(label, "-center",   kVstMaxParamStrLen); break;
  case DUCK_SHAPE:    vst_strncpy(label, "-smooth",   kVstMaxParamStrLen); break;

  } 
}

void vstChaosGenerator::getParameterDisplay(VstInt32 index, char* text)
{
  switch( index )
  {
  case TUNE:           sprintf(text, "%.2f", tune);             break;
  case GLIDE_TIME:     sprintf(text, "%.2f", glide_time);       break;
  case IN_MIX:         sprintf(text, "%.2f", in_mix);           break;
  case IN_DETUNE:      sprintf(text, "%.2f", in_det);           break;
  case FREQ_RATIO:     sprintf(text, "%.3f", freq_ratio);       break;

  case MOD_SCALE:      sprintf(text, "%.2f", mod_scale);        break;
  case MOD_12:         sprintf(text, "%.2f", mod_12);           break;
  case FMPM_12:        sprintf(text, "%.2f", fmpm12);           break;
  case MOD_21:         sprintf(text, "%.2f", mod_21);           break;
  case FMPM_21:        sprintf(text, "%.2f", fmpm21);           break;
  case MOD_ENV:        sprintf(text, "%.2f", mod_env);          break;
  case MOD_ATT:        sprintf(text, "%.2f", mod_att);          break;
  case MOD_DEC:        sprintf(text, "%.2f", mod_dec);          break;
  case MOD_SUS:        sprintf(text, "%.2f", mod_sus);          break;
  case MOD_REL:        sprintf(text, "%.2f", mod_rel);          break;

  //case LOWPASS_DETUNE: sprintf(text, "%.2f", lp_det);          break;
  //case LP_ORDER:       sprintf(text, "%i",   lp_order);         break;

  case FLT1_SCALE:     sprintf(text, "%.3f", fltScl1);          break;
  case FLT1_RESO:      sprintf(text, "%.3f", fltRes1);          break;
  case FLT1_MODE:      filterModeName(fltMode1, text);          break;
  case FLT2_SCALE:     sprintf(text, "%.3f", fltScl2);          break;
  case FLT2_RESO:      sprintf(text, "%.3f", fltRes2);          break;
  case FLT2_MODE:      filterModeName(fltMode2, text);          break;
  case FB_DETUNE:      sprintf(text, "%.2f", fb_det);           break;
  case FB_BOOST:       sprintf(text, "%.2f", fb_boost);         break;
  case FB_WIDTH:       sprintf(text, "%.2f", fb_width);         break;

  case FREQ_OFFSET1:   sprintf(text, "%.2f", freq_offset1);     break;
  case FREQ_OFFSET2:   sprintf(text, "%.2f", freq_offset2);     break;
  case WAVE1:          sprintf(text, "%i",   wave1);            break;
  case WAVE2:          sprintf(text, "%i",   wave2);            break;
  case CLIP1:          sprintf(text, "%.3f", clip1);            break;
  case CLIP2:          sprintf(text, "%.3f", clip2);            break; 
  case PHASE1:         sprintf(text, "%.2f", phase1);           break;
  case PHASE2:         sprintf(text, "%.2f", phase2);           break;

  case VOLUME:         sprintf(text, "%.2f", rsAmp2dB(volume)); break;
  case AMP_ATT:        sprintf(text, "%.2f", amp_att);          break;
  case AMP_DEC:        sprintf(text, "%.2f", amp_dec);          break;
  case AMP_SUS:        sprintf(text, "%.2f", amp_sus);          break;
  case AMP_REL:        sprintf(text, "%.2f", amp_rel);          break;
  case CLIP_OUT:       sprintf(text, "%.3f", clipOut);          break;

  case SYNC_MODE:      sprintf(text, "%i",   syncMode);         break;
  case SYNC_THRESH:    sprintf(text, "%.3f", sync_thresh);      break;
  case SYNC_AMOUNT:    sprintf(text, "%.3f", syncAmount);       break;
  case SYNC_FILTERS:   sprintf(text, "%i", (int) sync_filters); break;

  case DUCK_RANGE:     sprintf(text, "%.3f", duckRange);        break;
  case DUCK_FLATNESS:  sprintf(text, "%.3f", duckFlat);         break;
  case DUCK_CENTER:    sprintf(text, "%.3f", duckCenter);       break;
  case DUCK_SHAPE:     sprintf(text, "%i",   duckShape);        break;
  }
}

void vstChaosGenerator::getParameterLabel(VstInt32 index, char* label)
{
  switch( index )
  {
  case TUNE:           vst_strncpy(label, "st",  kVstMaxParamStrLen); break;
  case GLIDE_TIME:     vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case IN_DETUNE:      vst_strncpy(label, "st",  kVstMaxParamStrLen); break;
  case MOD_ENV:        vst_strncpy(label, "%",   kVstMaxParamStrLen); break;
  case MOD_ATT:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case MOD_DEC:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case MOD_REL:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case FB_DETUNE:      vst_strncpy(label, "st",  kVstMaxParamStrLen); break;
  case FB_BOOST:       vst_strncpy(label, "dB",  kVstMaxParamStrLen); break;
  case FB_WIDTH:       vst_strncpy(label, "oct", kVstMaxParamStrLen); break;
  case FREQ_OFFSET1:   vst_strncpy(label, "Hz",  kVstMaxParamStrLen); break;
  case FREQ_OFFSET2:   vst_strncpy(label, "Hz",  kVstMaxParamStrLen); break;
  case VOLUME:         vst_strncpy(label, "dB",  kVstMaxParamStrLen); break;
  case AMP_ATT:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case AMP_DEC:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  case AMP_REL:        vst_strncpy(label, "ms",  kVstMaxParamStrLen); break;
  default:             vst_strncpy(label, "",    kVstMaxParamStrLen); break;
  } 
}

void vstChaosGenerator::filterModeName(int index, char* text)
{
  switch(index)
  {
  case rsLadderFilter::FLAT:     vst_strncpy(text, "Flat",    kVstMaxParamStrLen); break;
  case rsLadderFilter::LP_6:     vst_strncpy(text, "LP6",     kVstMaxParamStrLen); break;
  case rsLadderFilter::LP_12:    vst_strncpy(text, "LP12",    kVstMaxParamStrLen); break;
  case rsLadderFilter::LP_18:    vst_strncpy(text, "LP18",    kVstMaxParamStrLen); break;
  case rsLadderFilter::LP_24:    vst_strncpy(text, "LP24",    kVstMaxParamStrLen); break;
  case rsLadderFilter::HP_6:     vst_strncpy(text, "HP6",     kVstMaxParamStrLen); break;
  case rsLadderFilter::HP_12:    vst_strncpy(text, "HP12",    kVstMaxParamStrLen); break;
  case rsLadderFilter::HP_18:    vst_strncpy(text, "HP18",    kVstMaxParamStrLen); break;
  case rsLadderFilter::HP_24:    vst_strncpy(text, "HP24",    kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_6_6:   vst_strncpy(text, "BP6/6",   kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_6_12:  vst_strncpy(text, "BP6/12",  kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_6_18:  vst_strncpy(text, "BP6/18",  kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_12_6:  vst_strncpy(text, "BP12/6",  kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_12_12: vst_strncpy(text, "BP12/12", kVstMaxParamStrLen); break;
  case rsLadderFilter::BP_18_6:  vst_strncpy(text, "BP18/6",  kVstMaxParamStrLen); break;
  default:                       vst_strncpy(text, "Unknown", kVstMaxParamStrLen); break;
  }
}

void vstChaosGenerator::onNoteOn(int note, int velocity, int detune)
{
  chaosGen.noteOn(note, velocity);
}

void vstChaosGenerator::onControlChange(int index, int value)
{
  float val = value / 127.f;
  switch( index )
  {
  case   7: setParameterAutomated(VOLUME,   val); break;
  case  74: setParameterAutomated(TUNE,     val); break;
    // more to come....
  }
}

void vstChaosGenerator::onPitchWheel(int value)
{
  // something to do...
  int dummy = 0;
}