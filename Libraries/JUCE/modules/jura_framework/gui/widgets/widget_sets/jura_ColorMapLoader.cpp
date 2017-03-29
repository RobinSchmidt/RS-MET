ColorMapPreviewer::ColorMapPreviewer(ColorMap *mapToPreview)
{
  colorMap = mapToPreview;
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

//=================================================================================================

ColorMapLoader::ColorMapLoader(ColorMap *mapToUpdate)
  : StateLoadSaveWidgetSet(String("ColorMapLoader"))
  , previewer(mapToUpdate)
{
  layout = LABEL_AND_BUTTONS_ABOVE;

  stateLabel->setText("ColorMap");
  stateLabel->setDescription("Map from values to colors");

  stateSaveButton->setVisible(false);  // later, we may want to edit and save the maps, too

  addAndMakeVisible(previewer);
}

void ColorMapLoader::resized()
{
  StateLoadSaveWidgetSet::resized();

  previewer.setBounds(stateFileNameLabel->getBounds());
  // hmm...maybe it would be better to show the filename and put the preview below
}