//-------------------------------------------------------------------------------------------------
// construction/destruction:

RectangleComponent::RectangleComponent(const Colour &_fillColour, const Colour &_outlineColour, 
  int _outlineThickness)
: fillColour(_fillColour)
, outlineColour(_outlineColour)
, outlineThickness(_outlineThickness)
{
  //setOpaque(true);  // nah - this is bad for lasso-stuff
}

//-------------------------------------------------------------------------------------------------
// setup:

void RectangleComponent::setFillColour(const Colour &newColour)
{
  fillColour = newColour;
  repaint();
}

void RectangleComponent::setOutlineColour(const Colour &newColour)
{
  outlineColour = newColour;
  repaint();
}

void RectangleComponent::setOutlineThickness(int newThickness)
{
  outlineThickness = newThickness;
  repaint();
}

//-------------------------------------------------------------------------------------------------
// others:

void RectangleComponent::paint(Graphics &g)
{
  g.setColour(fillColour);
  g.fillRect(0, 0, getWidth(), getHeight());
  g.setColour(outlineColour);
  g.drawRect(0, 0, getWidth(), getHeight(), outlineThickness);
}
