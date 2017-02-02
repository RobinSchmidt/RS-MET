
//-------------------------------------------------------------------------------------------------
// construction/destruction:

Editor::Editor(const String& newEditorName) 
{
  headlineX   = 4;
  headlineY   = 4;
  closeButton = nullptr;
  resizer     = nullptr;

  setHeadlineStyle(SUB_HEADLINE);

  //setMouseClickGrabsKeyboardFocus(false); 
  // otherwise, the qwerty-keyboard in Fruity Loops (and presumably other hosts as well) won't 
  // respond anymore after clicking a plugin-GUI
}

Editor::~Editor()
{
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// setup:

void Editor::setHeadlineText(const String& newHeadlineText)
{
  headlineText = newHeadlineText;
}

void Editor::setHeadlineStyle(int newHeadlineStyle)
{
  if( newHeadlineStyle >= NO_HEADLINE && newHeadlineStyle < NUM_HEADLINE_STYLES )
    headlineStyle = newHeadlineStyle;
}

void Editor::setHeadlinePosition(int newHeadlinePosition)
{
  if( newHeadlinePosition >= TOP_CENTER && newHeadlinePosition < NUM_HEADLINE_POSITIONS )
    headlinePosition = newHeadlinePosition;

  // retrieve width and height of the headline in pixels:
  int w = boldFont16px.getTextPixelWidth(headlineText, boldFont16px.getDefaultKerning()); 
  int h = boldFont16px.getFontHeight();

  switch( headlinePosition )
  {
  case TOP_LEFT:
    {
      headlineX = 4;
      headlineY = 4;
    }
    break;
  case TOP_RIGHT:
    {
      headlineX = getWidth()-w-4;
      headlineY = 4;
    }
    break;
  case BOTTOM_LEFT:
    {
      headlineX = 4;
      headlineY = getHeight()-h-4;
    }
    break;
  case BOTTOM_RIGHT:
    {
      headlineX = getWidth()-w-4;
      headlineY = getHeight()-h-4;
    }
    break;
  }
}

void Editor::setResizable(bool shouldBeResizable)
{
  if(shouldBeResizable)
  {
    if(resizer == nullptr)
    {
      resizer = new ResizableBorderComponent(this, nullptr);
      resizer->setBorderThickness(BorderSize<int>(4, 4, 4, 4)); // 4 pixels
      resizer->setBounds(0, 0, getWidth(), getHeight());
      addAndMakeVisible(resizer);
    }
  }
  else
  {
    delete resizer;   
    resizer = nullptr;
    // it's ok to delete a child component before removing it, see 
    // Component::removeChildComponent documentation
  }
}

//-------------------------------------------------------------------------------------------------
// inquiry:

const String Editor::getHeadlineString()
{
  return headlineText;
}

int Editor::getHeadlineBottom()
{
  if( headlineStyle != NO_HEADLINE ) 
    return headlineY + boldFont16px.getFontHeight();
  else
    return 0;
}

int Editor::getHeadlineRight()
{
  return headlineX + boldFont16px.getTextPixelWidth(
    headlineText, boldFont16px.getDefaultKerning());
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void Editor::paint(Graphics &g)
{
  ColourSchemeComponent::paint(g);  
  g.setColour(editorColourScheme.outline);
  for(int i=0; i<guiLayoutRectangles.size(); i++)
  {
    fillRectWithBilinearGradient(g, guiLayoutRectangles[i], editorColourScheme.topLeft, 
      editorColourScheme.topRight, editorColourScheme.bottomLeft, editorColourScheme.bottomRight);
    g.drawRect(guiLayoutRectangles[i], 1);
  }
  if( headlineStyle != NO_HEADLINE )
    drawHeadline(g);
}

void Editor::resized()
{
  ColourSchemeComponent::resized();
  setHeadlinePosition(headlinePosition);
  if(resizer != nullptr)
    resizer->setBounds(0, 0, getWidth(), getHeight());
  //repaint();
}

//-------------------------------------------------------------------------------------------------
// others:

void Editor::addChildEditor(Editor *editorToAdd, bool addAsChildComponent, bool makeVisible)
{
  ColourSchemeComponent::addChildColourSchemeComponent(editorToAdd, addAsChildComponent, 
    makeVisible);
  childEditors.getLock().enter();
  childEditors.addIfNotAlreadyThere(editorToAdd);
  childEditors.getLock().exit();
  editorToAdd->setDescriptionField(this->getDescriptionField(), true);
}

void Editor::removeChildEditor(Editor *editorToRemove, bool deleteObject)
{
  childEditors.getLock().enter();
  childEditors.removeFirstMatchingValue(editorToRemove);
  childEditors.getLock().exit();
  ColourSchemeComponent::removeChildColourSchemeComponent(editorToRemove, deleteObject);
}

void Editor::drawHeadline(Graphics& g)
{
  if( headlineStyle == NO_HEADLINE )
    return;

  switch( headlineStyle )
  {
  case SUB_HEADLINE:   
      drawBitmapFontText(g, headlineX, headlineY, headlineText, &boldFont16px, 
        editorColourScheme.headline);
    break;
  case MAIN_HEADLINE:   
    {
      int kerning       = 3*boldFont16px.getDefaultKerning();
      int textWidth     = boldFont16px.getTextPixelWidth(headlineText, kerning);
      int w2            = textWidth/2;
      int centerX       = getWidth()/2;
      int x             = centerX - w2;

      Colour c = editorColourScheme.headlineOutline;
      drawBitmapFontText(g, x-1, headlineY-1, headlineText, &boldFont16px, c, kerning);
      drawBitmapFontText(g, x+1, headlineY+1, headlineText, &boldFont16px, c, kerning);
      drawBitmapFontText(g, x+1, headlineY-1, headlineText, &boldFont16px, c, kerning);
      drawBitmapFontText(g, x-1, headlineY+1, headlineText, &boldFont16px, c, kerning);
      drawBitmapFontText(g, x,   headlineY,   headlineText, &boldFont16px, 
        editorColourScheme.headline, kerning);
    }
    break;
  }
}

//=================================================================================================

// construction/destruction:
EditorWithStateFile::EditorWithStateFile(const String& name) : Editor(name)
{
  addChildColourSchemeComponent( stateWidgetSet = new StateLoadSaveWidgetSet() );

  //if( moduleToEdit != NULL )
  //  moduleToEdit->addStateWatcher(stateWidgetSet);
  //stateWidgetSet->setDescriptionField(infoField);

  addStateWatcher(stateWidgetSet); // the widget-set watches the state of this object...mmm...kinda circular referencing
  stateWidgetSet->stateLabel->setText(juce::String("File:"));
  stateWidgetSet->addChangeListener(this);
}

EditorWithStateFile::~EditorWithStateFile()
{
  stateWidgetSet->removeChangeListener(this);
  removeStateWatcher(stateWidgetSet);
  deleteAllChildren();
}

//-------------------------------------------------------------------------------------------------
// callbacks:

void EditorWithStateFile::resized()
{
  Editor::resized();
  stateWidgetSet->setBounds(4, getHeadlineBottom()+4, getWidth()-8, 16);
}

int EditorWithStateFile::changeListenerCallback(void* objectThatHasChanged)
{
  if( objectThatHasChanged == stateWidgetSet )
  {
    stateWidgetSet->stateFileNameLabel->setText(getStateNameWithStarIfDirty());
    updateWidgetsAccordingToState();
  }
  return 0;
}
