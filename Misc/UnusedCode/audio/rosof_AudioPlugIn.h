#ifndef rosof_AudioPlugIn_h
#define rosof_AudioPlugIn_h

#include "rosof_AudioModule.h"

#include "../../rojue/components/misc/rojue_RectangleComponent.h"

namespace rosof
{

  class AudioPlugIn : public AudioProcessor, public ChangeBroadcaster
  {

  public:


    AudioPlugIn();

    virtual ~AudioPlugIn();

    virtual void prepareToPlay(double sampleRate, int samplesPerBlock);

    virtual void releaseResources();

    virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

    // by default, plugins have an editor:
    virtual bool hasEditor() const { return true; }

    virtual AudioProcessorEditor* createEditor();

    virtual const juce::String  getName() const;
    virtual int                 getNumParameters();
    virtual float               getParameter(int index);
    virtual void                setParameter(int index, float newValue);
    virtual const juce::String  getParameterName (int index);
    virtual const juce::String  getParameterText (int index);
    virtual const juce::String  getInputChannelName(int channelIndex) const;
    virtual const juce::String  getOutputChannelName(int channelIndex) const;
    virtual bool                isInputChannelStereoPair (int index) const;
    virtual bool                isOutputChannelStereoPair(int index) const;  
    virtual bool                acceptsMidi() const;  
    virtual bool                producesMidi() const;
    virtual int                 getNumPrograms()                                { return 0; }
    virtual int                 getCurrentProgram()                             { return 0; }
    virtual void                setCurrentProgram (int index)                   { }
    virtual const juce::String  getProgramName (int index)               { return juce::String::empty; }
    virtual void                changeProgramName (int index, const juce::String& newName) { }
    virtual void                getStateInformation(MemoryBlock& destData);
    virtual void                setStateInformation(const void* data, int sizeInBytes);

    /** Implements the response to a midi-message. */
    virtual void handleMidiMessage(MidiMessage message);

    /** Recalls a state (i.e. the settings of all relevant parameters) from an XmlElement. */
    virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& name);

    /** Returns the state (i.e. the settings of all relevant parameters) in form of an XmlElement. */
    virtual XmlElement* getStateAsXml(XmlElement* xmlElementToStartFrom = NULL) const;

    
    AudioModule *underlyingAudioModule;  // the actual dsp work engine

    CriticalSection plugInLock; // mutex-lock for all accesses to the underlyingAudioModule's member functions - a pointer to the lock is 
                                // passed to the embedded AudioModule and should be used there also and the AudioModule should also pass
                                // this lock on to the GUI Editors 


    //===============================================================================================
    juce_UseDebuggingNewOperator;

  protected:

    /** This must be overriden in subclasses to calculate one stereo sample frame - the slots serve 
    as inputs and outputs at the same time. */
    //virtual void getSampleFrameStereo(double *left, double *right);

    /** Calculates a sample frame for arbitrary number of channels. */
    //virtual void getSampleFrame(double* inOutSamples, int numChannels); 



    juce::String plugInName;    // name of the plugIn - assign this in the constructor of your subclass



  };

}

#endif  
