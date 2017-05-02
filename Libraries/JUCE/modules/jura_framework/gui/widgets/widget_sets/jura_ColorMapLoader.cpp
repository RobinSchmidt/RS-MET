ColorMapPreviewer::ColorMapPreviewer(ColorMap *mapToPreview)
{
  colorMap = mapToPreview;
  colorMap->addChangeListener(this);
}

ColorMapPreviewer::~ColorMapPreviewer()
{
  colorMap->removeChangeListener(this);
}

void ColorMapPreviewer::paint(Graphics& g)
{
  float w = (float) getWidth();
  float h = (float) getHeight();
  ColourGradient grad = colorMap->getAsColourGradient();
  grad.isRadial = false;
  if(w > h) {
    grad.point1 = Point<float>(0, 0);
    grad.point2 = Point<float>(w, 0); }
  else {
    grad.point1 = Point<float>(0, h);
    grad.point2 = Point<float>(0, 0); }
  g.setGradientFill(grad);
  g.fillAll();
}

void ColorMapPreviewer::changeListenerCallback(ChangeBroadcaster* source)
{
  if(source == colorMap)
    repaint();
}

//=================================================================================================

LoadableColorMap::LoadableColorMap()
{
  setActiveDirectory(getApplicationDirectory() + "/ColorMaps");
}

void LoadableColorMap::setStateFromXml(const XmlElement& xmlState, const juce::String& stateName,
  bool markAsClean)
{
  ColorMap::setFromXml(xmlState);
}

XmlElement* LoadableColorMap::getStateAsXml(const juce::String& stateName, bool markAsClean)
{
  return ColorMap::getAsXml();
}

//=================================================================================================

ColorMapLoader::ColorMapLoader(LoadableColorMap *mapToUpdate)
  : StateLoadSaveWidgetSet(String("ColorMapLoader"))
  , previewer(mapToUpdate)
{
  loadableColorMap = mapToUpdate;
  loadableColorMap->addStateWatcher(this);

  layout = LABEL_AND_BUTTONS_ABOVE;

  stateLabel->setText("ColorMap");
  stateLabel->setDescription("Map from values to colors");

  stateSaveButton->setVisible(false);  // later, we may want to edit and save the maps, too

  addAndMakeVisible(previewer);
}

void ColorMapLoader::stateDirtyFlagChanged(StateManager *stateManager)
{
  StateLoadSaveWidgetSet::stateDirtyFlagChanged(stateManager);
  previewer.repaint();
}

void ColorMapLoader::resized()
{
  int x = 0;
  int y = 0;
  int w = getWidth();
  int h = getHeight();
  //int x1, x2;
  if( layout == LABEL_AND_BUTTONS_ABOVE ) 
  {
    int h3 = h/3;
    stateLabel->setBounds(      0, y,  w, h3);
    x  = w-h3;
    statePlusButton->setBounds (x, y, h3, h3);
    x -= (h3-2);
    stateMinusButton->setBounds(x, y, h3, h3);
    x -= (40-2);
    stateLoadButton->setBounds (x, y, 40, h3);
    x -= (40-2);
    stateSaveButton->setBounds (x, y, 40, h3);
    x  = 0;
    y += (h3-2);
    stateFileNameLabel->setBounds(0, y, w, h3);
    y += h3;
    previewer.setBounds(0, y, w, h3);
  }
  else // ONE_LINE
  {
    jassertfalse; // not yet implemented, below is the inherited baseclass code - must be adapted
    //x2  = getWidth()-4;
    //w   = 16;
    //x1  = x2-w;
    //statePlusButton->setBounds(x1, 0, w, 16);
    //x1 -= w; 
    //stateMinusButton->setBounds(x1+2, 0, w, 16);
    //w   = 40;
    //x1  = stateMinusButton->getX()-(w-2);
    //stateLoadButton->setBounds(x1, 0, w, 16);
    //x1 -= (w-2);
    //stateSaveButton->setBounds(x1, 0, w, 16);
    //x2 = stateSaveButton->getX()+2;
    //x1 = 4;
    //w  = x2-x1;
    //stateFileNameLabel->setBounds(x1, 0, w, 16);
  }
}