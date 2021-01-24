#ifndef rosic_PolyModule_h
#define rosic_PolyModule_h

namespace rosic
{


class rsMusicalEvent
{

public:

  enum class EventType
  {
    unknown,
    noteOn,
    noteOff,
    controller,
    pitchBend
    // allNotesOff, 
    // aftertouch,
  };

protected:

  EventType type = EventType::unknown;
  void* data = nullptr;
  //int channel;
  //int voice;

};

//=================================================================================================


/*
// moved to jura - it fits better there
class rsVoiceManager
{

public:

  rsVoiceManager()
  {
    setMaxNumVoices(16);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setMaxNumVoices(int newNumber)
  {
    maxNumVoices    = newNumber;
    numVoices       = RAPT::rsMin(numVoices,       maxNumVoices);
    numActiveVoices = RAPT::rsMin(numActiveVoices, maxNumVoices);
    playingVoices.resize(maxNumVoices);
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getMaxNumVoices() const { return maxNumVoices; }

  int getNumVoices() const { return numVoices; }

  int getNumActiveVoices() const { return numActiveVoices; } 


protected:

  int maxNumVoices    = 16;  // maximum number of voices
  int numVoices       =  8;  // number of available voices
  int numActiveVoices =  0;  // number of currently playing voices

  std::vector<int> playingVoices;

};
*/







//=================================================================================================

class rsPolyModule;

/** 

*/

class rsPolyModuleManager
{

public:

  void addManagedModule(rsPolyModule* newModule)
  { RAPT::rsAppendIfNotAlreadyThere(managedModules, newModule); }

  void removeManagedModule(rsPolyModule* moduleToRemove)
  { RAPT::rsRemoveFirstOccurrence(managedModules, moduleToRemove); }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of voices that are currently playing. */
  int getNumActiveVoices() const { return 0; } // preliminary

  /** Returns the index of the i-th voice that is currently playing. */
  int getVoiceIndex(int i) const { return 0; } // preliminary


protected:

  std::vector<rsPolyModule*> managedModules;   // managed modules

};

//=================================================================================================

class rsPolyModuleObserver
{

public:

  enum class MessageType
  {
    unknown,
    numInputsChanged,
    numOutputsChanged,
    willBeDeleted
  };


  virtual void handleModuleMessage(rsPolyModule* messengerModule, MessageType* msg) = 0;

};

//=================================================================================================

/** Baseclass for polyphonic audio modules. */

class rsPolyModule
{

public:



  virtual ~rsPolyModule()
  {
    if(polyManager != nullptr)
      polyManager->removeManagedModule(this);

    // todo: send observers willBeDeleted message so they can invalidate their pointers
  }

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setPolyModuleManager(rsPolyModuleManager* newManager) 
  {
    if(polyManager != nullptr)
      polyManager->removeManagedModule(this);  // de-register from old manager
    polyManager = newManager;                  // set new manager
    if(polyManager != nullptr)
      polyManager->addManagedModule(this);     // de-register at old manager
  }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the number of signal (audio or modulation) inputs of this module. */
  //int getNumInputs() const { return numInputs; }

  /** Returns the number of signal (audio or modulation) outputs of this module. */
  //int getNumOutputs() const { return numOutputs; }

  //-----------------------------------------------------------------------------------------------
  // \name Processing

  //virtual void processFrame(const double* in, int numIns, double* out, int numOuts, int voice)
  //{ RAPT::rsArrayTools::clear(out, numOuts); }
  // maybe make purely virtual
  // maybe the user should pass numIn/Outputs - if the number disaggrees with what the module 
  // naturally produces, apply some rules, like:
  // -if more outputs are requested than produced, duplicate the existing ones
  // -maybe use getSample instead - it's less flexible but simpler and more efficient - maybe 
  //  handle multichannel stuff later

  /** Must be overriden by subclasses to produce one sample at a time for the given voice. */
  virtual rsFloat64x2 getSample(const rsFloat64x2& in, int voice) { return rsFloat64x2(0, 0); }

  /** Block processing function. Can be overriden optionally for optimization. */
  virtual void processBlock(rsFloat64x2* inOut, int numSamples, int voice)
  {
    for(int n = 0; n < numSamples; n++)
      inOut[n] = getSample(inOut[n], voice);
  }

  /** Event handler. Must be overriden, if the module needs to respond to events. */
  virtual void handleEvent(rsMusicalEvent* event) {};

  /** under construction */
  void updateModulatedParameters(int voice);


  /** Must be overriden by subclasses to produce one sample at a time. */
  //virtual rsFloat64x2 getSample(const rsFloat64x2& in) = 0;



  // get rid - use processFrame instead

  // maybe have a function processFrame(double* in, double *out); the number of in/out channels may
  // vary from module to module






protected:

  //int numInputs  = 0;   // number of input signals
  //int numOutputs = 0;   // number of output signals

  rsPolyModuleManager *polyManager = nullptr;  
  // needs to be passed to the constructor, should contain info like (max) number of voices, etc.

  //std::vector<ModulatableParameter*> params;
  // maybe move into baseclass rsModulatable..

  std::string name;      // must be unique to identify modulation targets
  std::string typeName;  // Ladder, WaveOsc, etc.

};

// todo: 
// -maybe don't distinguish between audio-sources and modulators at this code level - modulators 
//  are also just rsPolyModules - that means they must also produce stereo output - maybe in most
//  cases, just copy the output into both channels


}

#endif