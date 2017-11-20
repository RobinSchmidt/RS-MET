#ifndef VST_RESOREPLACER_H
#define VST_RESOREPLACER_H

#include "../Common/Utilities/vstPlugIn.h"

class vstResoReplacer : public vstPlugIn
{

public:
	
  vstResoReplacer(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void updateCoreParameter(VstInt32 index, float value);
  virtual void getParameterLabel(  VstInt32 index, char* label);
  virtual void getParameterDisplay(VstInt32 index, char* text);
  virtual void getParameterName(   VstInt32 index, char* text);

  virtual void getWaveformName(     VstInt32 index, char* text);
  virtual void getFeedbackSatPlaceName(VstInt32 index, char* text);

protected:

  // parameter indices:
  enum parameters
  {
    INGAIN = 0,      // input gain
    LEAK,            // input leakage
    FREQUENCY,       // lowpass cutoff frequency
    RESODECAY,       // resonance decay time
    DECAYBYFREQ,     // ...scaling of decaytime by frequency
    RESOATTACK,      // resonance attack time
    RESOGAIN,        // resonance gain
    AMPOFFSET,       // value added to instantaneous amplitude
    RESOLIMIT,       // amplitude limit for resonance
    RESORANGE,       // resonance range based on input
    RESOPHASE,       // resonance phase
    RESOWAVE,        // resonance waveform
    WAVECUTOFF,      // frequency multiplier for resonance waveform filter
    WAVEMOD,         // ...modulation by instantaneous amplitude
    CHAOSAMOUNT,     // amount of phase-bumping chaos - maybe rename to growl
    CHAOSCUTOFF,     // cutoff frequency for the smoother
    CHAOSPARAM,      // parameter that determines type of chaos

    SELFEXCITE,      // self-excitation (experimental - maybe remove)

    // feedback saturation parameters (new):
    FBDRIVE,         // feedback drive
    FBLOLIMIT,       // feedback low limit
    FBHILIMIT,       // feedback high limit
    FBSATGAINAT1,    // feedback saturation gain at x=1
    FBSATPLACE,      // feedback saturation place


    NUM_PARAMETERS
  };

  // mapped parameters:
  double inGain, leak, freq, decay, decayByFreq, attack, gain, ampOffset, limit, range, phase, 
    waveCut, waveMod, chaosAmt, chaosCut, chaosPar, selfEx, fbDrive, fbLoLimit, fbHiLimit, 
    fbGainAt1;
  int waveForm, fbSatPlace;

  // audio processing core:
  //rsResoReplacer filter;

  rsResoReplacerPhaseBumped filter;

};

// todo: 
// remove gate stuff from resoShape
// add saturation parameters here
// maybe use a more interesting filter instead of a simple lowpass to post-filter the resonance 
// waveform
//

#endif