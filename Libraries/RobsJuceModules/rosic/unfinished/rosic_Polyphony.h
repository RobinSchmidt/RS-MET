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

  virtual void processFrame(const double* in, int numIns, double* out, int numOuts, int voice)
  { RAPT::rsArrayTools::clear(out, numOuts); }
  // maybe make purely virtual
  // maybe the user should pass numIn/Outputs - if the number disaggrees with what the module 
  // naturally produces, apply some rules, like:
  // -if more outputs are requested than produced, duplicate the existing ones

  virtual void handleEvent(rsMusicalEvent* event) {};


  /** Must be overriden by subclasses to produce one sample at a time. */
  //virtual rsFloat64x2 getSample(const rsFloat64x2& in) = 0;


  /** Must be overriden by subclasses to produce one sample at a time for the given voice. */
  //virtual rsFloat64x2 getSample(const rsFloat64x2& in, int voice) = 0;
  virtual rsFloat64x2 getSample(const rsFloat64x2& in, int voice) { return rsFloat64x2(0, 0); }
  // get rid - use processFrame instead

  // maybe have a function processFrame(double* in, double *out); the number of in/out channels may
  // vary from module to module


  /*
  virtual void processBlockStereo(rsFloat64x2* inOut, int numSamples)
  {
    for(int n = 0; n < numSamples; n++)
      inOut[n] = getSample(inOut[n]);
  }
  */

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
// -don't distinguish between audio-sources and modulators at this level


}

#endif