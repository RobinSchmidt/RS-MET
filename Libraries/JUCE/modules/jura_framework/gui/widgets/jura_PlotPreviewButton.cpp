PlotPreviewButton::PlotPreviewButton(const String& name, const rsPlot* plotToPreview) 
: RButton(name)
{
  plotPreviewImage = NULL;
}

PlotPreviewButton::~PlotPreviewButton()
{

}

void PlotPreviewButton::paint(Graphics &g)
{
  // not yet implemeted
  if( plotPreviewImage != NULL )
  {
    g.drawImage(*plotPreviewImage, 0, 0, getWidth(), getHeight(), 
                0, 0, plotPreviewImage->getWidth(), plotPreviewImage->getHeight());
  }
  else
  {
    RButton::paint(g);
  }

  if( getToggleState() == true )
  {
    g.setColour(Colours::blue); // use some colour from the colour scheme
    g.drawRect(0, 0, getWidth(), getHeight(), 3);
  }
}