#ifndef VST_RESOSHAPER_H
#define VST_RESOSHAPER_H

#include "../Common/Utilities/vstPlugIn.h"

class vstResoShaper : public vstPlugIn
{

public:
	
  vstResoShaper(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void  updateCoreParameter(VstInt32 index, float value);
  virtual void  getParameterLabel(  VstInt32 index, char* label);
  virtual void  getParameterDisplay(VstInt32 index, char* text);
  virtual void  getParameterName(   VstInt32 index, char* text);

  virtual void getResoSatModeName(     VstInt32 index, char* text);
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
    RESOPHASE,       // resonance phase
    FBDRIVE,         // feedback drive
    FBLOLIMIT,       // feedback low limit
    FBHILIMIT,       // feedback high limit
    GATE,            // gate sensitivity
    GATEMIX,         // mix between lowpass-signal and input-signal for gate-signal
    FBSATGAINAT1,    // feedback saturation gain at x=1
    FBSATPLACE,      // feedback saturation place
    SATGAINAT1,      // saturation gain at x=1 - 0.5: x/(|x|+1), 1: hardclip
    DRIVE,           // saturation drive in dB
    DRIVE_COMP,      // post-saturation drive-compensation amount
    SATMODE,         // saturation mode
    SATADDCONST,     // constant added to saturator input
    SATADDIN,        // amount of input signal added to saturator input
    SATADDFLT,       // amount of nonresonant filtered signal added to saturator input

    NUM_PARAMETERS
  };

  // mapped parameters (maybe use an array of doubles):
  double inGain, leak, freq, decay, decayByFreq, attack, gain, phase;
  double fbDrive, fbLoLimit, fbHiLimit, gate, gateMix, fbGainAt1;
  double gainAt1, drive, driveComp, addConst, addIn, addFlt;
  int fbSatPlace;
  int satMode;

  // audio processing core:
  rsLadderResoShaped2 filter;

};

#endif

/*
ToDo:
-There are two places of saturation: in the feedback loop of the resonant filter and in the
 resonance post-processor. Maybe we could place another saturator on the filtered output
-filter-fm by input 
-resonance-phase-modulation by input, resonance-signal (the latter is a kind of phase-distortion,
 i think)
-implement thresholding/ducking of the resonance-signal by the input signal
 see mails around 2016/02/08
-the "smash" leads to artifacts - some kind of buzzing leaks through
-make parameter ranges larger than practical
*/
