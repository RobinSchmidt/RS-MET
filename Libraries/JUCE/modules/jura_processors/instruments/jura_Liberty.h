#ifndef jura_Liberty_h
#define jura_Liberty_h

//#include "../../../romos/Source/romos.h"
//using namespace romos;

//=================================================================================================

/** This is a class for keeping track of the state of the user interface of the modular synth 
(i.e., which panel is open, what is the scroll/zoom-setting, etc.).  */

class LibertyInterfaceState
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LibertyInterfaceState();

  //-----------------------------------------------------------------------------------------------
  // persistence:

  //  ...getAsXml, setFromXml

  enum panels
  {
    STRUCTURE_PANEL = 0,
    INTERFACE_PANEL,
    HYBRID_PANEL
  };

  int activePanel;

};

//=================================================================================================

/** This class wraps romos::Liberty into a jura::AudioModule to facilitate its use as plugIn. */

//class LibertyAudioModule : public PolyphonicInstrumentAudioModule  
class LibertyAudioModule : public AudioModule
{

  friend class LibertyEditor;
  friend class LibertyInterfaceMediator;
  friend class ModularBlockDiagramPanel;

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  LibertyAudioModule(CriticalSection *newPlugInLock, romos::Liberty *modularSynthToWrap);

  //-----------------------------------------------------------------------------------------------
  // parameter settings:

  virtual void setSampleRate(double newSampleRate)
  {
    wrappedLiberty->setSampleRate(newSampleRate);
  }

  //-----------------------------------------------------------------------------------------------
  // persistence:

  /** Adds module-type specific state data to an existing xml-element, if any. For example, a 
  constant-module may have its value stored here, or a multi-stage equalizer may store the number 
  of biquad stages and settings per stage, etc. */
  static void writeModuleTypeSpecificStateDataToXml(romos::Module *module, XmlElement* xmlState);

  /** Restores module-type specific state data from a xml-element, if any. 
  @see writeModuleTypeSpecificStateDataToXml */
  static void restoreModuleTypeSpecificStateDataFromXml(romos::Module *module, 
    const XmlElement& xmlState);

  /** Returns the state of the passed module as a (pointer to) an XmlElement. The caller is 
  responsible for eventually deleting the XmlElement or adding it as child-element to some other 
  XmlElement in which case it will be deleted together the parent-element. The second argument 
  specifies whether or not the external connections of the module should be stored - when the user 
  stores a container from the diagram-panel, we don't want this - otherwise, normally yes. */
  static XmlElement* getModuleStateAsXml(romos::Module *module, bool withExternalConnections);

  /** Sets up the state of the passed module from the passed XmlElement. First, it calls 
  createAndSetupEmbeddedModulesFromXml, then it calls createConnectionsFromXml making the 
  restoration of the state a two-pass process through the XmlElement. */
  static void setModuleStateFromXml(const XmlElement& xmlState, romos::Module *module);

  /** Creates all the child modules (recursively) for the passed module and sets theri state from 
  the passed XmlElement. This function  does not wire up the connections. */
  static void createAndSetupEmbeddedModulesFromXml(const XmlElement& xmlState, 
    romos::Module *module);

  /** Assuming that the child modules have been already properly created, it creates all the 
  connections in the passed module (and recursively inside its children, if any. */
  static void createConnectionsFromXml(const XmlElement& xmlState, romos::Module *module);

  /** Returns the state of the whole instrument as XmlElement. */
  virtual XmlElement* getStateAsXml(const juce::String &stateName, bool markAsClean);

  /** Restores the state of the whole instrument from the passed XmlElement. */
  virtual void setStateFromXml(const XmlElement& xmlState, const juce::String& stateName, 
    bool markAsClean);

  /** Calling this functions inside setStateFromXml is a quick and dirty ad-hoc solution to restore 
  the names and positions of the top-level I/O modules. */
  virtual void restoreTopLevelInOutStates(const XmlElement& xmlState);

  //  \todo override the xml get/set methods to include the GUI state

  /** Returns a pointer to the interface state object that is maintained here. This pointer can be 
  used by the GUI to set and retrieve interface settings from the LibertyAudioModule object. We 
  maintain the interface state here for storing it in the presets. ...nah- we don't need that 
  function - the GUI can access the member directly */
  //LibertyInterfaceState* getInterfaceState() { return &interfaceState; }

  //-----------------------------------------------------------------------------------------------
  // event-handling

  /** Triggers a note-on event. */
  virtual void noteOn(int noteNumber, int velocity);

  /** Triggers a note-off event. */
  virtual void noteOff(int noteNumber);


  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates a stereo-ouput frame. */
  virtual void getSampleFrameStereo(double* inOutL, double* inOutR)
  {
    wrappedLiberty->getSampleFrameStereo(inOutL, inOutR);
  }

  virtual void processBlockStereo(float *left, float *right, int numSamples)
  {
    if(wrappedLiberty->isSilent())
    {
      rosic::fillWithZeros(left, numSamples);
      rosic::fillWithZeros(right, numSamples);
    }
    else
    {
      wrappedLiberty->getBlockOfSampleFramesStereo(left, right, numSamples);


      /*
      for(int n = 0; n < numSamples; n++)
      {
        double dL = (double) left[n];
        double dR = (double) right[n];
        wrappedLiberty->getSampleFrameStereo(&dL, &dR);
        left[n]  = (float) dL;
        right[n] = (float) dR;
      }
      */

    }
  }

  //-----------------------------------------------------------------------------------------------
  // others:

  virtual void reset()
  {
    wrappedLiberty->reset();
  }

protected:

  /** Pointer to the underlying object which is wrapped. */
  romos::Liberty *wrappedLiberty;

  LibertyInterfaceState interfaceState; // maintains info about open panels, scroll-positions, etc.

  juce::File macroDirectory;

  juce_UseDebuggingNewOperator;
};

//=================================================================================================


/** A mix-in class to mix into the basic rojue widget classes. */
class LibertyModuleWidget
{

public:

  LibertyModuleWidget(const juce::String &name) 
  { 
    widgetParameterName = name; 
  }

  juce::String getWidgetParameterName() const   
  { 
    return widgetParameterName; 
  }

  juce_UseDebuggingNewOperator;

protected:

  juce::String widgetParameterName;

};

/** RSlider, augmented by the LibertyModuleWidget mix-in class. */
class LibertySlider : public RSlider, public LibertyModuleWidget
{
public:
  LibertySlider(const juce::String &name) : LibertyModuleWidget(name), RSlider(name) { }
  juce_UseDebuggingNewOperator;
};

/** RComboBox, augmented by the LibertyModuleWidget mix-in class. */
class LibertyComboBox : public RComboBox, public LibertyModuleWidget
{
public:
  LibertyComboBox(const juce::String &name) : LibertyModuleWidget(name), RComboBox(name) { }
  juce_UseDebuggingNewOperator;
};

/** RNamedComboBox, augmented by the LibertyModuleWidget mix-in class. */
class LibertyNamedComboBox : public RNamedComboBox, public LibertyModuleWidget
{
public:
  LibertyNamedComboBox(const juce::String &name) : LibertyModuleWidget(name), 
    RNamedComboBox(name, name) { }
  juce_UseDebuggingNewOperator;
};

/** RTextEntryField, augmented by the LibertyModuleWidget mix-in class. */
class LibertyTextEntryField : public RTextEntryField, public LibertyModuleWidget
{
public:
  LibertyTextEntryField(const juce::String &name) 
    : LibertyModuleWidget(name), RTextEntryField(name) { }
  juce_UseDebuggingNewOperator;
};


/** RLabeledTextEntryField, augmented by the LibertyModuleWidget mix-in class. */
class LibertyLabeledTextEntryField : public RLabeledTextEntryField, public LibertyModuleWidget
{
public:
  LibertyLabeledTextEntryField(const juce::String &name) 
    : LibertyModuleWidget(name), RLabeledTextEntryField(name) 
  { 
    removeChildWidget(entryField, true, true);
    entryField = new LibertyTextEntryField(name);
    addChildWidget(entryField);
  }

  juce_UseDebuggingNewOperator;
};

// todo: make classes LibertyButton, etc.
// we also need to override setValueFromString in RButton, RComboBox, RTextEntryField, etc. then

//=================================================================================================
// class ModulePropertiesEditor:

/** This class is used to edit some properties of modules in the modular synthesizer such as the 
name, the polyphony setting, etc. */


class ModulePropertiesEditor : public Editor, public RSliderListener, public RComboBoxObserver, 
  public RTextEntryFieldObserver, public RButtonListener
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor.  */  
  ModulePropertiesEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);

  /** Destructor. */
  //virtual ~ModulePropertiesEditor();


  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rSliderValueChanged(RSlider         *rSlider);
  virtual void rComboBoxChanged(   RComboBox       *comboBoxThatHasChanged);
  virtual void textChanged(        RTextEntryField *rTextEntryFieldThatHasChanged);
  virtual void rButtonClicked(     RButton         *buttonThatWasClicked);
  virtual void resized();

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Updates the widgets from the moduleToEdit member. */
  virtual void updateWidgetsFromModuleState();


protected:


  /** Called internally by the callbacks for the specific widgets rSliderValueChanged, rComboBoxChanged, etc.. The function tries to cast
  the widget to a LibertyWidget, reads out the parameter-name and sets the parameter with given in the moduleToEdit (provided, that this 
  moduleToEdit has parameters, which it should when this function is called). */
  virtual void widgetChanged(RWidget *widgetThatHasChanged);


  CriticalSection *plugInLock;   
  romos::Module   *moduleToEdit;

  RTextField *moduleTypeLabel, *moduleTypeField;
  RButton    *polyButton;


  // preset load/save section

  // non-editable fields (just for info):
  // String moduleKind;
  // int x, y;
  // numAudioInputs, numChildModules, numAudioOutputs, numEventInputs, numEventOutputs

  juce_UseDebuggingNewOperator;
};

//=================================================================================================
// subclasses for concrete module types:

class ParameterModuleEditor : public ModulePropertiesEditor 
{
public:
  ParameterModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit); 
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void resized();
  virtual void updateWidgetsFromModuleState();
  juce_UseDebuggingNewOperator;
protected:

  juce::Rectangle<int> topSectionRect, setupRect, controlSetupRect;

  RTextField *parameterSetupLabel;
  RTextField *helpTextLabel;
  LibertyTextEntryField        *helpTextField;
  RTextField                   *minValueLabel, *maxValueLabel;
  LibertyTextEntryField        *minValueField, *maxValueField;
  LibertySlider                *valueSlider;
  LibertyLabeledTextEntryField *valueField, *defaultField, *unitField;
  LibertyNamedComboBox         *scalingComboBox; 
  RClickButton                 *enterValueButton, *setToDefaultButton;

  // current, default, scaling, unit, quantization
  // ControllerNumber, ControlRangeMin, ControlRangeMax, Smoothing

  /*
  RTextField            *minValueLabel,  *defaultValueLabel,  *maxValueLabel;
  LibertyTextEntryField *minValueField,  *defaultValueField,  *maxValueField;
  RClickButton          *setToMinButton, *setToDefaultButton, *setToMaxButton;

  */

  // mappingBox
  // metaSilder, metaMinField, metaMaxField
  // smoothingSlider
};

class ContainerModuleEditor : public ModulePropertiesEditor 
{
public:
  ContainerModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit); 
  juce_UseDebuggingNewOperator;

protected:

};

class TopLevelModuleEditor : public ModulePropertiesEditor 
{
public:
  TopLevelModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit); 
  juce_UseDebuggingNewOperator;

protected:

};

class VoiceKillerModuleEditor : public ModulePropertiesEditor 
{
public:
  VoiceKillerModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);   
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  //RLabeledTextEntryField *thresholdField, *timeOutField;
  LibertySlider *thresholdSlider, *timeOutSlider;
};


class WhiteNoiseModuleEditor : public ModulePropertiesEditor 
{
public:
  WhiteNoiseModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);   
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  LibertySlider *seedSlider;
};


class BiquadDesignerModuleEditor : public ModulePropertiesEditor 
{
public:
  BiquadDesignerModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);   
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  LibertyNamedComboBox *modeComboBox;
};


class LibertyLadderFilterModuleEditor : public ModulePropertiesEditor 
{
public:
  LibertyLadderFilterModuleEditor(CriticalSection *newPlugInLock, romos::Module* newModuleToEdit);   
  virtual void resized();
  juce_UseDebuggingNewOperator;
protected:
  LibertyNamedComboBox *filterModeComboBox, *saturationModeComboBox;
};






#endif
