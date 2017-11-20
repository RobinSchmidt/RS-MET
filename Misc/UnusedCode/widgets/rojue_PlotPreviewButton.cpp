#include "rojue_PlotPreviewButton.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction and static member initialization:

PlotPreviewButton::PlotPreviewButton(const String& name, const CoordinateSystemOld* plotToPreview) 
: RButton(name)
{
  plotPreviewImage = NULL;
}

PlotPreviewButton::~PlotPreviewButton()
{

}


void PlotPreviewButton::paint(Graphics &g)
{
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