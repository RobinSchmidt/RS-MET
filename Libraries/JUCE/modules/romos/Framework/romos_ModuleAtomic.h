#ifndef romos_ModuleAtomic_h
#define romos_ModuleAtomic_h

#include "romos_Module.h"

namespace romos
{


  // maybe rename file to ModuleBaseClasses - we have more specific baseclasses here too


  //=======================================================================================================================================
  // class ModuleAtomic:

  /**

  This is the baseclass for atomic modules (except for some very special ones). Atomic modules differ from composite modules in the way 
  they deal with connections. Atomic modules deal with their connections directly while composite modules delegate their connections to 
  their I/O child-modules.

  */

  class ModuleAtomic : public Module
  {

    friend class ModuleContainer;

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup of pins:



    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry about pins:

    /** Returns the name of one of our pins. */
    virtual rosic::String getPinName(int kind, int direction, int pinIndex) const; 


    // these two are already defined in Module:
    /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured from our first input memory location 
    which is returned by getAudioInputAddress(). */
    /*
    INLINE int getInputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
    {
      return (frameIndex * ProcessingStatus::getSoleInstance()->getNumAllocatedVoices() + voiceIndex) * getNumAudioInputs() + pinIndex;
    }
    */
    /** Returns the offset of the memory-location for the given sample-frame/voice/input-pin measured from our first output memory location 
    which is returned by getAudioOutputAddress(). */
    /*
    INLINE int getOutputPinMemoryOffset(int frameIndex, int voiceIndex, int pinIndex)
    {
      return (frameIndex * ProcessingStatus::getSoleInstance()->getNumAllocatedVoices() + voiceIndex) * getNumAudioOutputs() + pinIndex;
    }
    */

    //=====================================================================================================================================

  protected:

    /** Creates an initial number of input pins. To be called from initialize only. When more pins should be added later, use 
    addAudioInput. */
    virtual void initInputPins(int numberOfPins, ...);

    /** Creates an initial number of output pins. To be called from initialize only. When more pins should be added later, use 
    addAudioOutput. */
    virtual void initOutputPins(int numberOfPins, ...);

    /** Adds an audio input. A name for the pin can optionally be passed. */
    virtual void addAudioInput(const rosic::String &pinName = rosic::String());

    /** Adds an audio output. A name for the pin can optionally be passed. */
    virtual void addAudioOutput(const rosic::String &pinName = rosic::String());

    /** Deletes the audio input with the given index. */
    virtual void deleteAudioInput(int index);

    /** Deletes the audio output with the given index. */
    virtual void deleteAudioOutput(int index);

    /** Constructor. Protected because instances should be created only through the ModuleFactory. */
    ModuleAtomic(const rosic::String &name = rosic::String(), int x = 0, int y = 0, bool polyphonic = false);

    /** Destructor. Protected because instances should be deleted only through the ModuleFactory. */
    virtual ~ModuleAtomic();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // data members:
 
    std::vector<rosic::String> audioInputNames;
    std::vector<rosic::String> audioOutputNames;

  private:

    // copy and assignment - not possible:
    ModuleAtomic(const ModuleAtomic &other) {}                          
    ModuleAtomic& operator=(const ModuleAtomic &other) { return *this; }    

  };




  // factor out the parameter-stuff into a mix-in class so we can derive - for example - ParameterModule from ConstantModule and 
  // ParameterMixIn






  //=======================================================================================================================================
  // class ParameterMixIn:

  /**

  A mix-in class that can be used for modules that need parameters that should be presented on the GUI editor. We use a mix-in approach to
  make it possible to derive one specific modue type with GUI-parameters from another specific module type that doesn't have GUI-parameters
  without resorting to multiple virtual inheritance.

  */

  class ParameterMixIn
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets the parameter with given name to the given new value. If no parameter with the given name is found, it returns false, 
    otherwise true. */
    bool setParameter(const rosic::String &parameterName, const rosic::String &newValue, bool callInternalCallback = true);

    /** Sets the parameter with given index to the given new value. When the "callInternalCallback" parameter is true, it will call 
    the virtual parameterChanged() callback function after setting the value. You may override this function in order to do any internal 
    variable updates that may be necessary when a parameter changes. Normally, you want this but in certain circumstances it may be necessary
    to supress the callback by passing false to avoid endless (indirect) recursions. */
    void setParameter(int index, const rosic::String &newValue, bool callInternalCallback = true);

    /** Appends a parameter with given name and default value to this module - the current value will also be initialized with
    the default value. */
    void addParameter(const rosic::String &parameterName, const rosic::String &defaultValue);

    /** Internal callback that is triggered from setParameter - you may override it when you need to re-compute some internal variables 
    when a parameter was changed. */
    virtual void parameterChanged(int index) { }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of parameters that this module has. */
    int getNumParameters() const { return parameters.size(); }

    /** Returns the name of the parameter with given index. */
    rosic::String getParameterName(int index) const;

    /** Returns the value of the parameter with given index. */
    rosic::String getParameterValue(int index) const;

    /** Returns the value of the parameter with the given name. If no parameter with the given name is found, it return 0.0. */
    rosic::String getParameterValue(const rosic::String &parameterName) const;

    /** Returns the default value of the parameter with given index. */
    rosic::String getParameterDefaultValue(int index) const;

    /** Returns the index of the parameter with given name. Returns -1, if no matching name is found. */
    int findIndexOfParameterWithName(const rosic::String &nameToFind) const;

    //=====================================================================================================================================

  protected:

    // data:
    struct Parameter
    {
      rosic::String name;
      rosic::String value;         
      rosic::String defaultValue; 
    };

    std::vector<Parameter> parameters;

  };









  //=======================================================================================================================================
  // class ModuleWithParameters:

  /**

  Baseclass for modules that have some parameters that are adjusted from the GUI. This class is suitable as baseclass for modules with 
  parameters that do not need to derive from other modules (without parameters).

  */

  class ModuleWithParameters : public ModuleAtomic, public ParameterMixIn
  {

  protected:

    ModuleWithParameters(const rosic::String &name = rosic::String(), int x = 0, int y = 0, bool polyphonic = false)
      : ModuleAtomic(name, x, y, polyphonic)
    {

    }

    virtual ~ModuleWithParameters() 
    { 
    
    }

  private:

    // copy and assignment - not possible:
    ModuleWithParameters(const ModuleWithParameters &other) {}                          
    ModuleWithParameters& operator=(const ModuleWithParameters &other) { return *this; }   

  };











  //=======================================================================================================================================
  // class ModuleProxy:

  /**

  Baseclass for modules that must redirect/delegate the pointers (and related variables) for the input signals to their connected 
  source-modules and notify target modules when their source-module has changed. Input- and and output modules fall into this class.

  */

  class ModuleProxy : public ModuleAtomic
  {

  public:



    virtual void mapApparentSourceToProcessingSource(Module * &sourceModule, int &sourceOutputPinIndex);
    virtual void mapProcessingSourceToSourceApparent(Module * &sourceModule, int &sourceOutputPinIndex);
    virtual void initialize();


    /** Assigns dummy functions that do nothing. */
    virtual void assignProcessingFunctions();

    static void processFrameDummy(Module *module, int voiceIndex)
    {
      //DEBUG_BREAK; // dummy function - should not be called
    }
    static void processBlockDummy(Module *module, int voiceIndex, int blockSize)
    {
      //DEBUG_BREAK; // dummy function - should not be called
    }


  protected:

        
    ModuleProxy(const rosic::String &name = rosic::String(), int x = 0, int y = 0, bool polyphonic = false);
    virtual ~ModuleProxy();


  private:

    // copy and assignment - not possible:
    ModuleProxy(const ModuleProxy &other) {}                          
    ModuleProxy& operator=(const ModuleProxy &other) { return *this; }   

  };

} // end namespace romos

#endif 
