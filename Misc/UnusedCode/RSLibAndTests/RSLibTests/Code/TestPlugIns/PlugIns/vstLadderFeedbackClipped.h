#ifndef VST_LADDERFEEDBACKCLIPPED_H
#define VST_LADDERFEEDBACKCLIPPED_H

#include "../Common/Utilities/vstPlugIn.h"

class vstLadderFeedbackClipped : public vstPlugIn
{

public:
	
  vstLadderFeedbackClipped(audioMasterCallback audioMaster);	
  virtual bool getEffectName(char* name);
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR);
  virtual void resume();
  virtual void setSampleRate(float sampleRate);

  virtual void  updateCoreParameter(VstInt32 index, float value);
  virtual void  getParameterLabel  (VstInt32 index, char* label);
  virtual void  getParameterDisplay(VstInt32 index, char* text);
  virtual void  getParameterName   (VstInt32 index, char* text);

  //virtual void setWaveShape(VstInt32 index);
  //virtual void getWaveShapeName(VstInt32 index, char* text);

protected:

  // parameter indices:
  enum parameters
  {
    FREQUENCY,       // lowpass cutoff frequency

    //RESODECAY,       // resonance decay time
    //RESODRIVE,       // resonance drive in dB
    //FB_LPF,         // cutoff for feedback lowpass (temporary);

    RESO,
    LOLIMIT,         // resonance low limit
    HILIMIT,         // resonance high limit
    SHAPE,           // saturation waveshape - determines the gain at 1 of prototype sigmoid
    MODE,            // saturation mode (temporary)

    NUM_PARAMETERS
  };

  //// waveshape indices:
  //enum shapes
  //{
  //  CLIP = 0,
  //  CUBIC,
  //  QUARTIC,
  //  SIXTIC,
  //  TANH,
  //  ATAN,

  //  NUM_SHAPES
  //};

  // mapped parameters (maybe use an array of doubles):
  double freq, reso, loLimit, hiLimit, shape;

  //int shape;  // selects the waveshape of the saturator
  int mode;   // selects, where the saturator is applied

  // audio processing core:
  rsLadderFilterFeedbackSaturated filter;

};

#endif


/*
ToDo:

-the actual amplitude of the self-oscillation signal seems to depend on the cutoff frequency
 -has this to with the k-factor? it should depend on the cutoff...
 -try a fixed feedback gain

-figure out what causes this gain increase with some settings - and try get rid of it
 -example settings: cutoff: 2kHz, decay: 1ms, hi: 10, lo: -0.2, input: saw at -6dB
 -it looks like a linear decay
 -try other waveforms as input - see if the shape of the decay depends on the input shape
 -check, if the a/symmetry of limits affects it
 -try to add an offset pre- and/or post clipper
 -visualize the feedback signal
-maybe introduce a "Mode" parameter to select, where the nonlineraity is applied
 -pre/post multiplication by k, maybe also to the difference in - k*y[4], etc.
-use center/width instead of loLimit/hiLimit internally: l=c-w/2, h=c+w/2; c=(l+h)/2, w=h-l
-try different saturation functions: clip, tanh, atan, asinh, x/(1+abs(x)), x/(1+x^2)
 -see, how the nonmonotonicity of some of them affects output
-try a zero-delay feedback version



*/