#ifndef jura_AudioModuleEditor_h
#define jura_AudioModuleEditor_h

//=================================================================================================

/** Baseclass for GUI editors for AudioModule objects. */

class AudioModuleEditor : public jura::Editor, public ChangeListener, public RDialogBoxListener, 
  public RButtonListener
{

public:

  enum positions
  {
    INVISIBLE,
    RIGHT_TO_HEADLINE,
    BELOW_HEADLINE,
    RIGHT_TO_INFOLINE
  };

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AudioModuleEditor(AudioModule* newModuleToEdit);

  /** Destructor. */
  virtual ~AudioModuleEditor();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Passes a new AudioModule objcet to be edited. This should be used when the same editor object 
  should be re-used for editing another AudioModule. */
  virtual void setModuleToEdit(AudioModule* newModuleToEdit);

  /** Sets the pointer to the moduleToEdit member to NULL without doing anything else. This should 
  be called whenever the underlying AudioModule was deleted. */
  virtual void invalidateModulePointer();

  /** Makes this a top-level editor meaning that some additional widgets (global preferences 
  button, infoline etc.) should be drawn. */
  virtual void setAsTopLevelEditor(bool isTopLevel) { isTopLevelEditor = isTopLevel; }

  /** Sets the position of the link to the website. @see positions */
  virtual void setLinkPosition(int newPosition) { linkPosition = newPosition; }

  /** Sets the position of preset load/saev section. @see positions */
  virtual void setPresetSectionPosition(int newPosition) { presetSectionPosition = newPosition; }

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the bottom (in pixels) of the preset section. */
  virtual int getPresetSectionBottom();

  //-----------------------------------------------------------------------------------------------
  // callbacks:

  virtual void rDialogBoxChanged(RDialogBox* dialogBoxThatHasChanged);
  virtual void rDialogBoxOKClicked(RDialogBox* dialogBoxThatWantsToAcceptAndLeave);
  virtual void rDialogBoxCancelClicked(RDialogBox* dialogBoxThatWantsToBeCanceled);
  virtual void rButtonClicked(RButton *buttonThatWasClicked);
  virtual void changeListenerCallback(juce::ChangeBroadcaster *objectThatHasChanged);
  virtual void resized();

  /** Updates the widgets according to the state of the assignedParameter (if any) and updates the 
  state-widget set. calls updateWidgetEnablement(). */
  virtual void updateWidgetsAccordingToState();

  /** Override this if you want to update the enablement of some widgets according to the state
  of the module. Will be called from updateWidgetsAccordingToState(). */
  virtual void updateWidgetEnablement() {}

  //-----------------------------------------------------------------------------------------------
  // public data members:

  StateLoadSaveWidgetSet* stateWidgetSet;  // \todo check, why we have this in the public area?

protected:

  /** Automatically generates a slider for each parameter in the module which is being edited. */
  //virtual void autoGenerateSliders();

  /** Returns a poiner to an RSlider object with the given name or NULL if no such slider exists
  (in our array automatableSliders) */
  //virtual RSlider* getSliderByName(const juce::String& sliderName);

  /** Opens a dialog to adjust the global preferences like the colour-scheme, preset paths etc.
  If your subclass needs some special settings (like, for example, a sample-path), you may override
  this an open a custom dialog in your class. */
  virtual void openPreferencesDialog();

  /** Loads the current colorscheme into a file. @see aveColorSchemeToFile(). */
  virtual void loadPreferencesFromFile();

  /** Saves the current colorscheme into a file. The filename will be given by the name of the 
  underlying AudioModule concatenated with 'Preferences'. Later we may want to store other settings
  there as well (such as preset- and sample-paths etc.) - we may then have to move the function 
  into AudioModule. */
  virtual void savePreferencesToFile();

  // todo: replace loadColorSchemeFromFile()/saveColorSchemeToFile() with loadPreferencesFromFile()/savePreferencesToFile(),
  // introduce methods getPreferencesAsXml/setPreferencesFromXml - these can then be overrided by subclasses

  /** Returns the xml tag-name that should be used for storing the preferences. */
  virtual juce::String getPreferencesTagName();

  /** Returns the xml filename that should be used for storing the preferences. */
  virtual juce::String getPreferencesFileName();

  RClickButton            *setupButton;
  RHyperlinkButton        *webLink;
  RTextField              *infoField;
  ColourSchemeSetupDialog *setupDialog;

  //CriticalSection moduleToEditLock;  // replace with a pointer to the global plugInLock
  CriticalSection *plugInLock;         // mutex to access the moduleToEdit member
  AudioModule     *moduleToEdit;

  int presetSectionPosition, linkPosition;

  /** This is an array of the automatable sliders - if add RSlider objects here, they will be
  updated in AudioModuleEditor::updateWidgetsAccordingToState via calls to their
  updateWidgetFromAssignedParameter() methods. In the destructor, this array is cleared first
  without deleting the objects, such that it does not interfere with the deleteAllChildren-function
  (which is supposed to be called in the destructor). */
  //OwnedArray<RSlider,   CriticalSection> sliders;
  //OwnedArray<RButton,   CriticalSection> buttons;
  //OwnedArray<RComboBox, CriticalSection> comboBoxes;

  // factor out into a class TopLevelEditor (or something like that):
  bool   drawGradientsBasedOnOutlines, isTopLevelEditor;
  Colour gradientMidColour;
  int    numHueOffsets; 

  juce_UseDebuggingNewOperator;
};

//=================================================================================================

/** A wrapper class that wraps an object of a subclass of AudioModuleEditor into a
juce::AudioProcessorEditor. We cannot derive AudioModuleEditor directly from AudioProcessorEditor
because we already derive it from jura::Editor to get all the basic look-and-feel stuff 
(colorscheme, etc.) from that superclass. ...and Editor already derives from juce::Component - as 
does AudioProcessorEditor - so we get in all kinds of trouble if we try to do such a messy 
branch-and-reunite multiple inheritance. I tried it - it makes problems. We get all the inherited 
Component member variables twice and only one set is correctly set up with our setters (but some 
getters may return the wrong ones...things like that...it's a mess. Don't do it. */

class AudioModuleEditorWrapper : public juce::AudioProcessorEditor
{

public:

  AudioModuleEditorWrapper(AudioModuleEditor *newContentComponent, AudioModule* newModuleToEdit) 
    : AudioProcessorEditor(newModuleToEdit)
  {
    contentComponent = newContentComponent;
    setSize(contentComponent->getWidth(), contentComponent->getHeight());
    addAndMakeVisible(contentComponent);
  }

  AudioModuleEditorWrapper::~AudioModuleEditorWrapper()
  {
    delete contentComponent;  
      // we need to delete it here because the baseclass destructor does not delete it's child 
      // components - is this a change with respect to the old juce?
  }


  virtual void paint(Graphics &g) override {} // we hit a breakpoint if we don't override this

protected:

  AudioModuleEditor *contentComponent;

};

#endif