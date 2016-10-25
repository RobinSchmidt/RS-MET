#ifndef rosof_AudioModule_h
#define rosof_AudioModule_h

#include "../../rosic/infrastructure/rosic_PolyphonicInstrument.h"
using namespace rosic;

#include "../../rojue/components/widgets/rojue_AutomatableModule.h"
#include "../../rojue/misc/rojue_StringTools.h"
#include "../../rojue/file_management/rojue_StateFileManager.h"
using namespace rojue;

#include "../xml_tools/rosof_XmlToolsForRosic.h"

namespace rosof
{



  class AudioModule;

  /**

  A class that can be informed (via a callback method) when an AudioModule object is going to be deleted. Mainly intended as baseclass for 
  GUI elements that keep a pointer to an AudioModule that is being edited

  */

  class AudioModuleDeletionWatcher
  {

  public:

    /** Destructor. */
    virtual ~AudioModuleDeletionWatcher();

    /** Callback function that subclasses should override in order to invalidate theri references/pointers to the Audiomodule in 
    question. */
    virtual void audioModuleWillBeDeleted(AudioModule *moduleToBeDeleted) = 0;

    /** Registers ourselves as deletion-watcher with the passed AudioModule such that we will get callbacks to audioModuleWillBeDeleted 
    when the module in question is going to be deletd. */
    virtual void addWatchedAudioModule(AudioModule *moduleToBeWatched);

    /** De-registers ourselves from and AudioModule to which we presumably had previously registered via addWatchedAudioModule - if we 
    didn't, this function will do nothing. */
    virtual void removeWatchedAudioModule(AudioModule *moduleNotToBeWatchedAnymore);

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    juce::Array<AudioModule*> watchedAudioModules;

  };


    



  /**

  This class is the base class for all audio modules.

  */

  class AudioModule : public AutomatableModule, public StateFileManager
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. You must pass a CriticalSection object to be used for wrapping accesses to the underlying core dsp object. */
    AudioModule(CriticalSection *newPlugInLock);   

    /** Destructor. */
    virtual ~AudioModule();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Override this to set the sample-rate for this module. */
    virtual void setSampleRate(double newSampleRate); 

    /** Override this to set the tempo in bpm for synced modules. */
    virtual void setBeatsPerMinute(double newBpm); 

    /** Sets up the name for this AudioModule. */
    virtual void setModuleName(const juce::String& newName);

    /** Sets up an appendix (like "Demo Version") for the name of this AudioModule. */
    virtual void setModuleNameAppendix(const juce::String& newAppendix);

    /** Adds a child AudioModule to this one. */
    virtual void addChildAudioModule(AudioModule* moduleToAdd);

    /** Removes a child AudioModule from this one and optionally deletes the object (it will only delete it if it is actually in the array 
    of childAudioModules). */
    virtual void removeChildAudioModule(AudioModule* moduleToRemove, bool deleteOject);

    /** Use this function to turn on/off save/recall of the state - the idea is that some module may have child modules which are currently 
    inactive and therefore don't need to save and recall their state. This makes preset files more economical. */
    virtual void setStateSaveAndRecall(bool shouldSaveAndRecall)
    { saveAndRecallState = shouldSaveAndRecall; }

    /** Checks, if this is a cracked version and if so, it sets up the appendix for the headline accordingly. Return value informs also 
    whether or not a cracked version was detected. */
    virtual bool checkForCrack();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the name of this module. */
    virtual juce::String getModuleName() const { return moduleName; }

    /** When we have several child-modules with the same name (member "moduleName"), this function can be used to find the index of the
    passed child-module among them. It will return 0 when the passed AudioModule is the first (or only one) with that name, 1 for the 
    second and so on. When the passed module is not one of our child modules, it will return -1. */
    virtual int getIndexAmongNameSakes(AudioModule *child);

    /** Returns a string that is to be used as headline for a GUI-editor for this module - this will be typically the name of the module, 
    perhaps appended by some qualifier like "Demo Version", "Cracked By bLaH" or something. */
    virtual juce::String getModuleHeadlineString();

    /** Returns the interval at which the module wants to receive callbacks to trigger(). */
    virtual double getTriggerInterval() const { return triggerInterval; }

    /** Returns, whether this mdoule wants to save/recall its state - the idea is that some module may have child modules which are 
    currently inactive and therefore don't need to save and recall their state. This makes preset files more economical. */
    bool wantsSaveAndRecallState() const { return saveAndRecallState; }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // automation and state management:

    /** Callback to indicate that a parameter has changed - subclasses should override this and 
    update their signal processing accordingly. */
    virtual void parameterChanged(rojue::Parameter* parameterThatHasChanged);

    /** Calls a parameterChanged for each of the observed parameters - this should trigger the 
    appropriate updating of the signal processing core in the subclasses. */
    virtual void updateCoreObjectAccordingToParameters();

    /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
    virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, bool markAsClean);

    /** Returns the state (i.e. the settings of all relevant parameters) in form of an XmlElement. */
    virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean);
        
    /** Converts a state which might possibly be from an older version to the current patch-format. The baseclass implementation just 
    returns the state as is, but will trigger a debug-break if the patchFormatIndex of the state and the module don't match. Override this
    function in your subclass to do the actual conversion. */
    virtual XmlElement convertXmlStateIfNecessary(const XmlElement& xmlState);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Processes a block of sample-frames that is a sub-block of the passed 'buffer' starting at sample-index 'startSample' and having a 
    length of 'length' (thus, ending at sample-index startSample+length-1). */
    virtual void processBlock(AudioSampleBuffer& buffer, int startSample, int length)
    {
      if( buffer.getNumChannels() == 2 )
      {
        float *left  = buffer.getSampleData(0, startSample);
        float *right = buffer.getSampleData(1, startSample);
        processBlockStereo(left, right, length);
      }
      else
        jassertfalse; // if your subclass is not 2 in / 2 out, you should override this function
    }

    /** Processes a block of stereo sample-frames. The baseclass implementation just clears the buffers. For 2 in / 2 out modules, it is 
    easiest to override this function, for other I/O configurations override processBlock. */
    virtual void processBlockStereo(float *left, float *right, int numSamples)
    {
      for(int n=0; n<numSamples; n++)
      {
        left[n]  = 0.f;
        right[n] = 0.f;
      }
    }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // event processing:

    /** Handles a generic MidiMessage. */
    virtual void handleMidiMessage(MidiMessage message);

    /** Triggers a note-on event. */
    virtual void noteOn(int noteNumber, int velocity);

    /** Triggers a note-off event. */
    virtual void noteOff(int noteNumber);

    /** Triggers an all-notes-off event. */
    virtual void allNotesOff();

    /** Overrides setMidiController which is inherited from both base-classes - and we simply we pass through the function call to both of 
    them here. */
    virtual void setMidiController(int controllerNumber, int controllerValue); 

    /** Triggers a pitch-bend event. */
    virtual void setPitchBend(int pitchBendValue);

    /** Override this and set triggerInterval to some nonzero value if you need to re-trigger something at regular intervals (like LFOs, 
    for example). This function will be called from the process-function at the given triggerInterval (if nonzero) - this value has to be 
    specified in beats. */
    virtual void trigger() {}

    /** Override this to reset this module (audio buffers and such). */
    virtual void reset() {}

    /** Override this to reset the state of this module to defaults (user parameters). */
    virtual void setStateToDefaults() {}

    /** Flag to indicate that this module needs tempo sync information (current BPM). */
    bool wantsTempoSyncInfo;

    rosic::PolyphonicInstrument *underlyingRosicInstrument;

    //=====================================================================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** Must be overriden by subclasses to fill the inherited array of observed parameters. */
    virtual void initializeAutomatableParameters();

    /** Our child modules to which we will distribute MIDI-events and of which we manage the 
    states. */
    juce::Array<AudioModule*, CriticalSection> childModules;

    CriticalSection *plugInLock;            // mutex to access the wrapped core dsp object 

    juce::String moduleName;                // name of this AudioModule
    juce::String moduleNameAppendix;        // string to be appended to the name on the GUI (such as Demo-Version, etc.)

    int    patchFormatIndex;                // the version number of this module
    double triggerInterval;                 // interval (in beats) for calls to trigger()
    bool   saveAndRecallState;              // flag to indicate that this module wants to save/recall its state
    
  private:

    /** Registers an AudioModuleDeletionWatcher that will be called back when this object is deleted. */
    virtual void registerDeletionWatcher(AudioModuleDeletionWatcher *watcher);
        
    /** De-registers a previously registered AudioModuleDeletionWatcher. */
    virtual void deRegisterDeletionWatcher(AudioModuleDeletionWatcher *watcher);

    juce::Array<AudioModuleDeletionWatcher*> deletionWatchers;
    friend class AudioModuleDeletionWatcher;

  };

}

#endif 