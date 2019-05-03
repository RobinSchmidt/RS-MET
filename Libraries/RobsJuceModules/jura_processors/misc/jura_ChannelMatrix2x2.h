#ifndef jura_ChannelMatrix2x2_h
#define jura_ChannelMatrix2x2_h

class ChannelMatrix2x2AudioModule : public AudioModule
{

  friend class ChannelMatrix2x2ModuleEditor;

public:

  ChannelMatrix2x2AudioModule(CriticalSection *newPlugInLock, 
    rosic::ChannelMatrix2x2 *channelMatrix2x2ToWrap = nullptr);

  virtual ~ChannelMatrix2x2AudioModule();

  AudioModuleEditor* createEditor(int type) override;



  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
    bool markAsClean) override;

  virtual XmlElement* getStateAsXml(const juce::String& stateName, bool markAsClean) override;

  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedChannelMatrix2x2->getSampleFrameStereo(inOutL, inOutR);
  }

protected:

  rosic::ChannelMatrix2x2 *wrappedChannelMatrix2x2;
  bool wrappedChannelMatrix2x2IsOwned = false;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

class ChannelMatrix2x2ModuleEditor : public AudioModuleEditor, public RTextEntryFieldObserver
{

public:

  ChannelMatrix2x2ModuleEditor(CriticalSection *newPlugInLock, 
    ChannelMatrix2x2AudioModule* newChannelMatrix2x2AudioModule);

  virtual void textChanged(RTextEntryField *rTextEntryFieldThatHasChanged);

  virtual void resized();

protected:

  virtual void updateWidgetsAccordingToState();

  ChannelMatrix2x2AudioModule *channelMatrix2x2AudioModule;

  // the widgets:
  RTextField      *labelLeftToLeft, *labelRightToLeft, *labelLeftToRight, *labelRightToRight, 
    *leftEquationLabel, *rightEquationLabel;  
  RTextEntryField *editLabelLeftToLeft, *editLabelRightToLeft, *editLabelLeftToRight, 
    *editLabelRightToRight;

  juce_UseDebuggingNewOperator;
};

#endif 
