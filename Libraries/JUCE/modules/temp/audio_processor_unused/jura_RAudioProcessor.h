#ifndef jura_RAudioProcessor_h
#define jura_RAudioProcessor_h

class JUCE_API RAudioProcessor : public AudioProcessor, public ChangeBroadcaster, public Timer
{

public:

  RAudioProcessor();

  virtual ~RAudioProcessor();

  virtual void prepareToPlay(double sampleRate, int samplesPerBlock);

  virtual void releaseResources();

  virtual void processBlock(AudioSampleBuffer& buffer, MidiBuffer& midiMessages);

  virtual AudioProcessorEditor* createEditor() = 0;

  virtual const String  getName() const;
  virtual int           getNumParameters();
  virtual float         getParameter(int index);
  virtual void          setParameter(int index, float newValue);
  virtual const         String getParameterName (int index);
  virtual const         String getParameterText (int index);
  //virtual const         String getInputChannelName (const int channelIndex) const;
  //virtual const         String getOutputChannelName (const int channelIndex) const;
  //virtual bool          isInputChannelStereoPair (int index) const;
  //virtual bool          isOutputChannelStereoPair(int index) const;  
  virtual bool          acceptsMidi() const;  
  virtual bool          producesMidi() const;

  virtual int  getNumPrograms()                                { return 0; }
  virtual int  getCurrentProgram()                             { return 0; }
  virtual void setCurrentProgram (int index)                   { }
  virtual const String getProgramName (int index)              { return String::empty; }
  virtual void changeProgramName (int index, const String& newName)  { }

  virtual void getStateInformation (MemoryBlock& destData);
  virtual void setStateInformation (const void* data, int sizeInBytes);

  /** Implements the response to a midi-message. */
  virtual void respondToMidiMessage(MidiMessage theMessage);

  /** Returns the current state as pointer to an XmlElement - the object must be deleted by the
  caller when it's not needed anymore. */
  virtual XmlElement* getStateAsXml() = 0;

  /** Recalls a state according to an to an XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState) = 0;

  /** This is the callback which triggers the gain to be set to zero after 30 minutes for demo 
  versions. */
  virtual void timerCallback();

  /** Informs whether this is a demo version of not. */
  virtual bool isDemoVersion();

  //// the audio-engine:
  //rosic::PlugInEngine *plugInEngine;

protected:

  /** This method is supposed to be used for initialization which have to be done on loading the 
  plugIn such as checking whether this is a demo version and imposing the required 
  demo-restrictions (initialize the timer, etc.). Make sure to override it in your subclass and 
  call it in the constructor of the subclass. Can be ignored for freeware plugins. */
  virtual void initializePlugIn();

  /** This function checks for the existence and validity of a keyfile and returns true if a valid 
  keyfile is there and false otherwise. */
  virtual void checkForKeyFile();

  /** Check if a keyfile is valid or not. */
  virtual bool checkIfKeyFileIsValid(const File& keyFileToCheck);

  /** This must be overriden in subclasses to calculate one stereo sample frame - the slots serve 
  as inputs and outputs at the same time. */
  virtual void getSampleFrameStereo(double *left, double *right);

  /** Can be overriden by subclasses which need tempo syncronisation. The default implementation 
  is empty. */
  virtual void setBeatsPerMinute(double newBpm);

  /** An index for the product - this is needed to check for the keyfile. */
  int productIndex;

  /** Flag to indicate whether this is a commercial product proteccted by a keayfile or not - we 
  use a double to obfuscate the meaning for crackers. */
  double isKeyFileProtected;

  /** Flag to indicate whether this is a demo-version or not - we use a double to obfuscate the
  meaning for crackers. */
  double isInDemoMode;

  /** 'Demo-Gain' - a gain factor to multiply the signal with. The factor will be 1.0 for plugIns 
  which are either freeware, licensed payware or the 30 minute demo period is still running, 0.0 
  otherwise (demo-versions running longer than 30 minutes). */
  double demoGain;

  /** The name of the plugIn - assign this in the constructor of your subclass. */  
  String plugInName;

  // temporarily deactivated - this comes from the old codebase:
  ///** A Smoother for the 'Demo-Gain' to avoid clicks after the 30 minute demo period. */
  //rosic::OnePoleFilter dgSmoother;

  ///** The object which will be used to verify the validity of the key. */
  //rosic::KeyGenerator keyValidator;

  juce_UseDebuggingNewOperator;
};


#endif  
