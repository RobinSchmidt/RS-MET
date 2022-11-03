//#include "rojue_ColourSchemeComponent.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

ColourSchemeComponent::ColourSchemeComponent(const String& newColourSchemeComponentName) 
{
  drawWithEnclosingRectangle = true;
  setBufferedToImage(true);
}

ColourSchemeComponent::~ColourSchemeComponent()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void ColourSchemeComponent::setDescriptionField(RTextField *newDescriptionField, 
  bool callAlsoForChildColourSchemeComponents)
{
  ScopedLock scopedLock(arrayLock);
  descriptionField = newDescriptionField;
  int i;
  for(i=0; i<widgets.size(); i++)
    widgets[i]->setDescriptionField(newDescriptionField);
  for(i=0; i<widgetSets.size(); i++)
    widgetSets[i]->setDescriptionField(newDescriptionField, callAlsoForChildColourSchemeComponents);
  for(i=0; i<plots.size(); i++)
    plots[i]->setDescriptionField(newDescriptionField);
  if( callAlsoForChildColourSchemeComponents == true )
  {
    for(i=0; i<childColourSchemeComponents.size(); i++)
    {
      childColourSchemeComponents[i]->setDescriptionField(newDescriptionField, 
        callAlsoForChildColourSchemeComponents);
    }
  }
}

void ColourSchemeComponent::setEditorAppearance(int newAppearance)
{
  ScopedLock scopedLock(arrayLock);
  editorColourScheme.setAppearance(newAppearance);
  updateEmbeddedObjectsAndRepaint();
  //for(int i=0; i<childColourSchemeComponents.size(); i++)
  //  childColourSchemeComponents[i]->copyColourSettingsFrom(this);
  //updateEmbeddedObjectsAndRepaint();
}

void ColourSchemeComponent::setWidgetAppearance(int newAppearance)
{
  ScopedLock scopedLock(arrayLock);
  widgetColourScheme.setAppearance(newAppearance);
  updateEmbeddedObjectsAndRepaint();
  /*
  int i;
  for(i=0; i<widgets.size(); i++)
    widgets[i]->setColourScheme(widgetColourScheme);  
  for(i=0; i<widgetSets.size(); i++)
    widgetSets[i]->setColourScheme(widgetColourScheme);   
  repaint();
  */
}

void ColourSchemeComponent::setPlotAppearance(int newAppearance)
{
  ScopedLock scopedLock(arrayLock);
  plotColourScheme.setAppearance(newAppearance);
  updateEmbeddedObjectsAndRepaint();
  //for(int i=0; i<plots.size(); i++)
  //  plots[i]->setColourScheme(plotColourScheme); 
  //repaint();
}

void ColourSchemeComponent::setCentralHue(float newHue)
{
  ScopedLock scopedLock(arrayLock);
  widgetColourScheme.setCentralHue(newHue);
  plotColourScheme.setCentralHue(newHue);
  editorColourScheme.setCentralHue(newHue);
  updateEmbeddedObjectsAndRepaint();
}

void ColourSchemeComponent::setHueOffset(int index, float newOffset)
{
  editorColourScheme.setHueOffset(index, newOffset);
}

void ColourSchemeComponent::setSaturationMultiplier(float newMultiplier)
{
  ScopedLock scopedLock(arrayLock);
  widgetColourScheme.setSaturationMultiplier(newMultiplier);
  plotColourScheme.setSaturationMultiplier(newMultiplier);
  editorColourScheme.setSaturationMultiplier(newMultiplier);
  updateEmbeddedObjectsAndRepaint();
}

void ColourSchemeComponent::setBrightnessGamma(float newGamma)
{
  ScopedLock scopedLock(arrayLock);
  widgetColourScheme.setBrightnessGamma(newGamma);
  plotColourScheme.setBrightnessGamma(newGamma);
  editorColourScheme.setBrightnessGamma(newGamma);
  updateEmbeddedObjectsAndRepaint();
}

void ColourSchemeComponent::setGraphColours(const juce::Array<Colour> &newColours, 
                                            bool callAlsoForChildColourSchemeComponents)
{

}

void ColourSchemeComponent::copyColourSettingsFrom(const ColourSchemeComponent *componentToCopyFrom)
{
  widgetColourScheme = componentToCopyFrom->widgetColourScheme;
  plotColourScheme   = componentToCopyFrom->plotColourScheme;
  editorColourScheme = componentToCopyFrom->editorColourScheme;
  updateEmbeddedObjectsAndRepaint();
}

void ColourSchemeComponent::setWidgetsVisible(bool shouldBeVisible)
{
  ScopedLock scopedLock(arrayLock);
  for(int i=0; i<widgets.size(); i++)
  {
    Component* comp = dynamic_cast<Component*> (widgets[i]);
    if( comp != NULL )
      comp->setVisible(shouldBeVisible);
  }
  for(int i=0; i<widgetSets.size(); i++)
    widgetSets[i]->setVisible(shouldBeVisible);
}

void ColourSchemeComponent::setColourScheme(const WidgetColourScheme& newColourScheme)
{
  ScopedLock scopedLock(arrayLock);
  for(int w=0; w<widgets.size(); w++)
    widgets[w]->setColourScheme(newColourScheme);
}

void ColourSchemeComponent::setColourSchemeFromXml(const XmlElement& xml)
{
  editorColourScheme.setAppearanceFromString( xml.getStringAttribute(String("EditorAppearance")) );
  widgetColourScheme.setAppearanceFromString( xml.getStringAttribute(String("WidgetAppearance")) );
  plotColourScheme.setAppearanceFromString(   xml.getStringAttribute(String("PlotAppearance"))   );
  setCentralHue(                      (float) xml.getDoubleAttribute(String("CentralHue"))       );
  setSaturationMultiplier(            (float) xml.getDoubleAttribute(String("Saturation"))       );
  setBrightnessGamma(                 (float) xml.getDoubleAttribute(String("Gamma"), 1.0)       );

  // retrieve the hue-offsets for various parts of the editor:
  editorColourScheme.clearHueOffsets();
  bool hueOffsetFound = true;
  int  hueOffsetIndex = 1;
  while( hueOffsetFound == true )
  {
    String attributeName = String("HueOffset") + String(hueOffsetIndex);
    if( xml.hasAttribute(attributeName) )
    {
      float hueOffset = (float) xml.getDoubleAttribute(attributeName, 0.0);
      editorColourScheme.appendHueOffset(hueOffset);
      hueOffsetIndex++;
    }
    else
      hueOffsetFound = false;
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

XmlElement* ColourSchemeComponent::getColourSchemeAsXml()
{
  XmlElement *xml = new XmlElement(String("ColorScheme"));
  xml->setAttribute(String("EditorAppearance"), editorColourScheme.getAppearanceString());
  xml->setAttribute(String("WidgetAppearance"), widgetColourScheme.getAppearanceString());
  xml->setAttribute(String("PlotAppearance"),   plotColourScheme.getAppearanceString());
  xml->setAttribute(String("CentralHue"),       editorColourScheme.getCentralHue());
  xml->setAttribute(String("Saturation"),       editorColourScheme.getSaturationMultiplier());
  xml->setAttribute(String("Gamma"),            editorColourScheme.getBrightnessGamma());

  for(int i=0; i<editorColourScheme.getNumHueOffsets(); i++)
    xml->setAttribute(String("HueOffset") + String(i+1), editorColourScheme.getHueOffset(i));

  return xml;
}

//-------------------------------------------------------------------------------------------------
// others:

void ColourSchemeComponent::addWidget(RWidget *widgetToAdd, bool addAsChildComponent, 
  bool makeVisible)
{
  ScopedLock scopedLock(arrayLock);
  //widgets.addIfNotAlreadyThere(widgetToAdd);
  RAPT::rsAppendIfNotAlreadyThere(widgets, widgetToAdd);
  if( addAsChildComponent == true )
  {
    Component* comp = dynamic_cast<Component*> (widgetToAdd);
    if( comp != NULL )
    {
      addChildComponent(comp);
      if( makeVisible == true )
        comp->setVisible(true);
    }
  }
  widgetToAdd->setColourScheme(widgetColourScheme);
}

void ColourSchemeComponent::addWidgetSet(WidgetSet* widgetSetToAdd, bool addAsChildComponent, 
                                         bool makeVisible)
{
  ScopedLock scopedLock(arrayLock);
  //widgetSets.addIfNotAlreadyThere(widgetSetToAdd);
  RAPT::rsAppendIfNotAlreadyThere(widgetSets, widgetSetToAdd);
  if( addAsChildComponent == true )
  {
    //addChildComponent(widgetSetToAdd);
    addChildColourSchemeComponent(widgetSetToAdd);
    if( makeVisible == true )
      widgetSetToAdd->setVisible(true);
  }
  widgetSetToAdd->setColourScheme(widgetColourScheme);
}

void ColourSchemeComponent::addPlot(rsPlot *plotToAdd, bool addAsChildComponent, 
                                    bool makeVisible)
{
  ScopedLock scopedLock(arrayLock);
  //plots.addIfNotAlreadyThere(plotToAdd);
  RAPT::rsAppendIfNotAlreadyThere(plots, plotToAdd);
  if( addAsChildComponent == true )
  {
    addChildComponent(plotToAdd);
    if( makeVisible == true )
      plotToAdd->setVisible(true);
  }
  plotToAdd->setColourScheme(plotColourScheme);
}

void ColourSchemeComponent::addChildColourSchemeComponent(ColourSchemeComponent *editorToAdd, 
                                                          bool addAsChildComponent, 
                                                          bool makeVisible)
{
  ScopedLock scopedLock(arrayLock);
  //childColourSchemeComponents.addIfNotAlreadyThere(editorToAdd);
  RAPT::rsAppendIfNotAlreadyThere(childColourSchemeComponents, editorToAdd);
  if( addAsChildComponent == true )
  {
    addChildComponent(editorToAdd);
    if( makeVisible == true )
      editorToAdd->setVisible(true);
  }
  editorToAdd->copyColourSettingsFrom(this);
}

void ColourSchemeComponent::removeChildColourSchemeComponent(ColourSchemeComponent *childToRemove, 
                                                             bool deleteObject)
{
  ScopedLock scopedLock(arrayLock);
  removeChildComponent(childToRemove);
  //childColourSchemeComponents.removeValue(childToRemove);
  //childColourSchemeComponents.removeFirstMatchingValue(childToRemove);
  RAPT::rsRemoveFirstOccurrence(childColourSchemeComponents, childToRemove);
  if( deleteObject == true )
    delete childToRemove;
}

void ColourSchemeComponent::removeWidget(RWidget* widgetToRemove, bool removeAsChildComponent, 
  bool deleteObject)
{
  ScopedLock scopedLock(arrayLock);
  //jassert( widgets.contains(widgetToRemove) ); // trying to remove a widget that was not added before?
  //widgets.removeFirstMatchingValue(widgetToRemove);
  jassert( RAPT::rsContains(widgets, widgetToRemove) ); // trying to remove a widget that was not added before?
  RAPT::rsRemoveFirstOccurrence(widgets, widgetToRemove);
  if( removeAsChildComponent == true )
    removeChildComponent(widgetToRemove);
  if( deleteObject == true )
    delete widgetToRemove;
}

void ColourSchemeComponent::updateEmbeddedObjectsAndRepaint()
{
  ScopedLock scopedLock(arrayLock);  // added 2022/11/03

  int i;
  for(i=0; i<widgets.size(); i++)
    widgets[i]->setColourScheme(widgetColourScheme);  
  for(i=0; i<widgetSets.size(); i++)
    widgetSets[i]->setColourScheme(widgetColourScheme);  
  for(i=0; i<plots.size(); i++)
  {
    plots[i]->setPopUpColourScheme(widgetColourScheme); 
    plots[i]->setColourScheme(plotColourScheme); 
  }
  for(i=0; i<childColourSchemeComponents.size(); i++)
  {
    childColourSchemeComponents[i]->copyColourSettingsFrom(this);
  }
  repaint();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void ColourSchemeComponent::paint(Graphics &g)
{
  fillRectWithBilinearGradient(g, 0, 0, getWidth(), getHeight(),
    editorColourScheme.topLeft, editorColourScheme.topRight, editorColourScheme.bottomLeft, 
    editorColourScheme.bottomRight);
}

void ColourSchemeComponent::paintOverChildren(Graphics &g)
{
  if(drawWithEnclosingRectangle == true)
  {
    g.setColour(editorColourScheme.outline);
    g.drawRect(0, 0, getWidth(), getHeight(), 2);
  }

#if JUCE_DEBUG
  bool showSize = false; // set to true, if you want to figure out ideal gui sizes
  if(showSize)
  {
    int w = 80;
    int h = 20;
    int x = 0;
    int y = 0;
    //int x = getWidth()/2  - w/2;
    //int y = getHeight()/2 - h/2;
    g.setColour(Colours::black);
    g.fillRect(x, y, w, h);
    g.setColour(Colours::white);
    juce::String str = String(getWidth()) + " x " + String(getHeight());
    g.drawText(str, x, y, w, h, Justification::centred);
  }
#endif
}

void ColourSchemeComponent::updateWidgetsAccordingToState()
{
  ScopedLock scopedLock(arrayLock);
  for(int i = 0; i < widgets.size(); i++)
    widgets[i]->updateWidgetFromAssignedParameter(false);  
}


