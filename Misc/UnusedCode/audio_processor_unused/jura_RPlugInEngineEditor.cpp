#include "RPlugInEngineEditor.h"

RPlugInEngineEditor::RPlugInEngineEditor(rosic::PlugInEngine* plugInEngineToEdit) 
: PresetRemembererEditor(String(T("RPlugInEngineEditor")))
{
  // assign the pointer to the plugins audio engine:
  //plugIn       = plugInToEdit;
  plugInEngine = plugInEngineToEdit;

  // create the info-label and info-field:
  addAndMakeVisible( infoLabel = new RLabel(String(T("InfoField")), String(T("Info:"))) );
  infoLabel->setColour(Label::backgroundColourId, Colours::transparentWhite);
  infoLabel->setColour(Label::outlineColourId, Colours::transparentWhite);
  addAndMakeVisible( infoField = new RLabel(String(T("InfoField")), String::empty) );
  infoField->setColour(Label::backgroundColourId, Colours::transparentWhite);
  infoField->setColour(Label::outlineColourId, Colours::transparentWhite);

  // set up the inherited PresetRemembererEditor-object:
  drawWithEnclosingRectangle = false;
  enablePresetFileManagement(true);
  presetLabel->setFont(Font(16));
  presetLabel->setJustificationType(Justification::centredLeft);
  setPresetRemembererToEdit(plugInEngine);
  presetLabel->setDescriptionField(infoField);
  presetFileNameLabel->setDescriptionField(infoField);
  presetSaveButton->setDescriptionField(infoField);
  presetLoadButton->setDescriptionField(infoField);
  presetPlusButton->setDescriptionField(infoField);
  presetMinusButton->setDescriptionField(infoField);
  updatePresetField();

  // set the plugIn-headline:
  /*
  String headline = plugInEngine->getName();
  if( plugInEngine->isDemoVersion() )
    headline += String(T(" (Demo Version)"));
  setHeadline( headline );
  */
  headline->setFont(Font(18));
  headline->setJustificationType(Justification::centred);
  presetLabel->setFont(Font(16));

  // create the link:
  webLink = new HyperlinkButton(T("www.rs-met.com"), 
    URL(T("http://www.rs-met.com")));
  webLink->setFont(Font(14), false, Justification::horizontallyCentred);
  addAndMakeVisible(webLink);

  // set up the inherited RobsEditorBase-object:
  enablePresetFileManagement(true);

  // register as listener to change-messages from the plugIn (to reflect automation on the GUI):
  //plugInToEdit->addChangeListener(this);

  //setSize (944, 756);

  loadColorScheme();

  updateWidgetsAccordingToState();
}

RPlugInEngineEditor::~RPlugInEngineEditor()
{
  //plugIn->removeChangeListener(this);

  automatableSliders.clear(false); 
  // don't delete objects on 'clear' as they will be deleted by the subsequent deleteAllChildren()
  // call

  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// appearance:

void RPlugInEngineEditor::loadColorScheme()
{
  String dir = FileManager::getApplicationDirectory();
  File colorFile = File( dir + String(T("/ColorScheme.xml")) );
  if( colorFile.existsAsFile() )
  {
    XmlDocument myDocument(colorFile);
    XmlElement* xmlColors = myDocument.getDocumentElement();
    String colorString;
    Colour color;
    if( xmlColors != NULL )
    {
      colorString = xmlColors->getStringAttribute(T("WebLinks"), T("ff7777ff") );
      color = Colour::fromString(colorString);
      webLink->setColour(HyperlinkButton::textColourId, color);
      delete xmlColors;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// state-management:

XmlElement* RPlugInEngineEditor::getStateAsXml(const String &stateName) const
{
  return NULL;
  /*
  if( plugIn == NULL )
    return NULL;

  return plugIn->getStateAsXml();
  */
}

bool RPlugInEngineEditor::setStateFromXml(const XmlElement &xmlState)
{
  return false;
  /*
  if( plugIn == NULL )
    return false;

  plugIn->setStateFromXml(xmlState);
  updateWidgetsAccordingToState();
  return true;
  */
}

void RPlugInEngineEditor::updateWidgetsAccordingToState()
{
  for(int i=0; i<automatableSliders.size(); i++)
    automatableSliders[i]->updateWidgetFromAssignedParameter(false);
}

void RPlugInEngineEditor::setPlugInEngineToEdit(rosic::PlugInEngine* newPlugInEngineToEdit)
{
  plugInEngine = newPlugInEngineToEdit;
}

void RPlugInEngineEditor::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();

  headline->setBounds(0, 0, getWidth(), 20);
  headline->setJustificationType(Justification::centred);
  webLink->setBounds(getWidth()-112, 0, 112-4, 20);

  x = 0;
  y = headline->getBottom();
  w = getWidth()/4;

  presetLabel->setBounds(x+4, y+4, 52, 20);

  presetPlusButton->setBounds(w-4-20, y+4, 20, 20);
  presetMinusButton->setBounds(presetPlusButton->getX()-20, y+4, 20, 20);
  presetLoadButton->setBounds(presetMinusButton->getX()-40-4, y+4, 40, 20);
  presetSaveButton->setBounds(presetLoadButton->getX()-40-4, y+4, 40, 20);

  y = presetLabel->getBottom();
  presetFileNameLabel->setBounds(presetLabel->getX(), y+4, w-8, 20);

  infoLabel->setBounds(0, getHeight()-20, 40, 20);
  infoField->setBounds(infoLabel->getRight(), getHeight()-20, getWidth()-infoLabel->getRight(),20);

  renderBackgroundImage();
}
