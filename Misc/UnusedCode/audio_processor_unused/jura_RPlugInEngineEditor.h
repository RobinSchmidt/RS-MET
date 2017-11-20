#ifndef jura_RPlugInEngineEditor_h
#define jura_RPlugInEngineEditor_h

class JUCE_API RPlugInEngineEditor /*: public PresetRemembererEditor*/
{

public:

  /** Constructor. */
  RPlugInEngineEditor();
  //RPlugInEngineEditor(rosic::PlugInEngine* plugInEngineToEdit = NULL);

  /** Destructor. */
  ~RPlugInEngineEditor();

  //-----------------------------------------------------------------------------------------------
  // appearance:

  /** Triggers loading of a color-scheme from a file ColorScheme.xml - if the file doesn't exist,
  it will fall back to a default color-scheme. */
  virtual void loadColorScheme();

  //-----------------------------------------------------------------------------------------------
  // state-management:

  /** Overrides the method inherited from RobsEditorBase. */
  virtual XmlElement* getStateAsXml(
    const String& stateName = String("RAudioProcessorState")) const;

  /** Overrides the method inherited from RobsEditorBase. */
  virtual bool setStateFromXml(const XmlElement& xmlState);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Attaches this editor to the actual plugin which is to be edited. */
  //virtual void setPlugInEngineToEdit(rosic::PlugInEngine* newPlugInToEdit);

  /** Overrides the resized-method of the RobsEditorBase base-class. */
  virtual void resized();

protected:

  /** Overrides the method inherited from RobsEditorBase. */
  virtual void updateWidgetsAccordingToState();

  HyperlinkButton     *webLink;
  RLabel              *infoLabel;
  RLabel              *infoField;
  //rosic::PlugInEngine *plugInEngine;
  //RAudioProcessor     *plugIn;

  //Colour webLinkColour;

  /** This is an array of the automatable sliders - if add RSlider objects here, they will be 
  updated in RPlugInEngineEditor::updateWidgetsAccordingToState via calls to their
  updateWidgetFromAssignedParameter() methods. In the destructor, this array is cleared first 
  without deleting the objects, such that it does not interfere with the deleteAllChildren-function  
  (which is supposed to be called in the destructor). */
  OwnedArray<RSlider> automatableSliders;

};

#endif
