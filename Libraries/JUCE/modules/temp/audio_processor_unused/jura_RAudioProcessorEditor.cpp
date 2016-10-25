//#include "RAudioProcessorEditor.h"
////#include "RAudioModuleEditor.h"

RAudioProcessorEditor::RAudioProcessorEditor(RAudioProcessor* const ownerPlugIn) 
 : AudioProcessorEditor(ownerPlugIn)
{
  audioEngineEditor = new RPlugInEngineEditor(ownerPlugIn->plugInEngine);
  String headline = ownerPlugIn->getName();
  if( ownerPlugIn->isDemoVersion() )
    headline += String(" (Demo Version)");
  audioEngineEditor->setHeadline( headline );
  addAndMakeVisible( audioEngineEditor );

  //setSize(50, 50);
  //audioEngineEditor->setBounds(0, 0, getWidth(), getHeight() );

  // register ourselves with the plugIn and with the actual editor - this object will mediate the 
  // communication to allow for automaition:
  ownerPlugIn->addChangeListener(this); //??
  //synthEditor->addChangeListener(this);

  updateParametersFromFilter();
}

RAudioProcessorEditor::~RAudioProcessorEditor()
{
  getFilter()->removeChangeListener(this);
  deleteAllChildren();
}


void RAudioProcessorEditor::paint (Graphics& g)
{
  /*
  if( backgroundImage != NULL )
  {
    g.drawImage(backgroundImage, 0, 
      0, 
      getWidth(), 
      getHeight(), 
      0, 
      0, 
      backgroundImage->getWidth(), 
      backgroundImage->getHeight(), 
      false);
  }
  else
    g.fillAll(Colours::red);
  //fillRectWithDefaultBackground(g, getX(), getY(), getWidth(), getHeight() );
  */
}

//-------------------------------------------------------------------------------------------------

void RAudioProcessorEditor::changeListenerCallback (void* source)
{
  // this is the filter telling us that it's changed, so we'll update our
  // display of the time, midi message, etc.
  updateParametersFromFilter();
}

//-------------------------------------------------------------------------------------------------

void RAudioProcessorEditor::updateParametersFromFilter()
{
  RAudioProcessor* const thePlugIn = getFilter();
  //double value;
}
