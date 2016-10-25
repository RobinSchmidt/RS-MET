#include "rojue_RDraggableNumber.h"
using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction and static member initialization:

RDraggableNumber::RDraggableNumber(const String& name) : RSlider(name)
{
  //valueBox->setColour(Label::backgroundColourId, Colours::white);
  //valueBox->setColour(Label::outlineColourId,    Colours::black);
  valueOnMouseDown = 0.0;
}

RDraggableNumber::~RDraggableNumber()
{

}

//-------------------------------------------------------------------------------------------------
// others:

void RDraggableNumber::paint (Graphics& g)
{
  g.setColour(getBackgroundColour());
  g.fillRect(handleRectangle);

  // draw the outline:
  g.setColour(getOutlineColour());
  g.drawRect(handleRectangle, 2);

  // draw the value:
  String valueString = stringConversionFunction(currentValue);
  int x = getWidth() - boldFont10px.getTextPixelWidth(valueString, boldFont10px.getDefaultKerning());
  int y = handleRectangle.getY() + handleRectangle.getHeight()/2 - boldFont10px.getFontAscent()/2;
  drawBitmapFontText(g, x-4, y, valueString, &boldFont10px, getTextColour());

  // draw the name:
  x             = 4;
  Colour colour = getTextColour();
  if( layout == NAME_ABOVE )
  {
    x      = 0;
    y      = handleRectangle.getY() - boldFont10px.getFontAscent() - 2;
    //colour = colourScheme.special;
  }
  drawBitmapFontText(g, x, y, sliderName, &boldFont10px, colour);

  // gray out the slider if it's disabled:
  //if( !isEnabled() )
  //  g.fillAll(Colours::lightgrey.withAlpha(0.75f));
}

void RDraggableNumber::mouseDown (const MouseEvent& e)
{
  valueOnMouseDown = getValue();
  if( isEnabled() )
  {
    if( e.mods.isCommandDown() )
      setValue(defaultValue);
    else if( e.mods.isRightButtonDown() )
      RWidget::mouseDown(e);
  }
}

void RDraggableNumber::mouseDrag (const MouseEvent& e)
{
  double y = e.getMouseDownY() + e.getDistanceFromDragStartY();
  if( isEnabled() )
  {
    if( !e.mods.isRightButtonDown() )
    {
      double proportion  = valueToProportionOfLength(valueOnMouseDown);
      double amountToAdd = -(y/400.0);
      if( e.mods.isShiftDown() )
        amountToAdd /= 100.0; // fine adjustment with shift
      double newValue    = proportionOfLengthToValue(proportion+amountToAdd);
      setValue(newValue);
    }
  }
}

