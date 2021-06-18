#ifndef VST_PLUGIN_H
#define VST_PLUGIN_H

#include "../../../../../RSLib/Code/RSLib.h"
#include "public.sdk/source/vst2.x/audioeffectx.h" // replace with vestige.h
#include <vector>  // for the event buffer

/** This is the baseclass used for the test plugins. It simplifies the interface of the original   
VST baseclass. When you derive your plugin from this class instead of the original AudioEffectX
from the VST-SDK, you no longer have to override processReplacing, processDoubleReplacing - 
instead, you have to override processStereoFrame which is supposed to produce one stereo output
sampl-frae at a time (and is called here in our processReplacing override in a loop). This 
simplifies your implementation (but makes it also potenetially less efficient, but we don't care
about optimization here). Also, you don't need to ovveride setParameter anymore, instead 
override updateCoreParameter, which is called here in our overriden setParameter. This baseclass
keeps track of the normalized parameter in the array params - but you will want to keep track of
the corresponding unnormalized, mapped parameter in your subclass. */

class vstPlugIn : public AudioEffectX 
{

public:
	
  vstPlugIn(audioMasterCallback audioMaster, VstInt32 numPrograms, VstInt32 numParams);	
  ~vstPlugIn();	

  /** Your subclass has to override this in order to produce one stereo output frame at a time. */
  virtual void processStereoFrame(double *inL, double *inR, double *outL, double *outR) = 0;

  /** Your subclass has to override this in order to update its internal mapped parameter 
  representation for the parameter with giev index and possibly set up some variables in the core 
  DSP code. */
  virtual void  updateCoreParameter(VstInt32 index, float value) = 0;

  // functions inherited from AudioEffectX (in the VST-SDK) that are already overriden here for you
  // to simplify your implementation:
  virtual void processReplacing(float** inputs, float** outputs, VstInt32 sampleFrames);
  virtual void processDoubleReplacing(double** inputs, double** outputs, VstInt32 sampleFrames);
  virtual void  setParameter(VstInt32 index, float value);
  virtual float getParameter(VstInt32 index);
  virtual void setProgramName(char* name);
  virtual void getProgramName(char* name);
  virtual bool getVendorString(char* text);
  virtual bool getProductString(char* text);
  virtual VstInt32 getVendorVersion ();


protected:

  char programName[kVstMaxProgNameLen+1];
  static const VstInt32 defaultNumPrograms = 1;

  float *params; // normalized vst parameter values

  // buffering (used in the single-precision processReplacing function):
  void allocateBuffers(int size);
  void freeBuffers();
  int bufSize;
  double *buf[2];  // left and right channel buffer

};

//=================================================================================================

/** This class is supposed to be used as baseclass for VST instruments, i.e. VST plugins that
handle MIDI inputs. It simplifies the original VST-SDK baseclass interface by not requiring you
to override processEvents anymore - a function that is supposed to collect an array of time-stamped
events which are supposed to be used subsequently in the process call. Instead, this baseclass here
collects the events and let's you handle one event at a time by overriding handleEvent or 
alternatively overriding onNoteOn/onControlChange/onPitchWheel. Overriding the latter 3 is less 
flexible but more convenient and often just what you need. Just override the event handlers and
processStereoFrame (inherited from vstPlugIn), and this baseclass takes over responsibility for 
sample accurate event handling. */

class vstInstrument : public vstPlugIn 
{

public:

  vstInstrument(audioMasterCallback audioMaster, VstInt32 numPrograms, VstInt32 numParams);	
  ~vstInstrument();	

  /** Callback method that you may override, if you want to handle note events. The note, velocity
  and detune parameters are between 0 and 127. The detune parameter is an addition by the VST
  specification that is not present in the original MIDI specification and is supposed to be used
  for per note microtuning. */
  virtual void onNoteOn(int note, int velocity, int detune) {}

  /** Callback method that you may override, if you want to handle controller events. The index
  and value is between 0 and 127. */
  virtual void onControlChange(int index, int value) {}

  /** Callback method that you may override, if you want to handle pitch wheel events. The value
  is between -8192 and +8191. */
  virtual void onPitchWheel(int value) {}

  /** Your subclass may override this function in order to handle one MIDI event at a time. If you 
  don't override this, the baseclass implementation will do a dispatch to onNoteOn, 
  onControlChange, onPitchWheel which are more specific event handling callbacks that you may
  override alternatively for more convenience. This baseclass makes sure, that the MIDI events 
  arrive with sample accurate timing right before the corresponding subsequent call to your 
  overriden processStereoFrame function. */
  virtual void handleEvent(VstMidiEvent midiEvent);


  // functions from AudioEffectXthat are already overriden here for you:
  virtual VstInt32 processEvents(VstEvents* events);
  virtual void processDoubleReplacing(double** inputs, double** outputs, VstInt32 sampleFrames);
  virtual VstInt32 canDo(char* text);
  //virtual VstInt32 getNumMidiInputChannels () { return 1; }


protected:

  vector<VstEvent> events; // collected in processEvents, used in processDoubleReplacing

};

#endif
