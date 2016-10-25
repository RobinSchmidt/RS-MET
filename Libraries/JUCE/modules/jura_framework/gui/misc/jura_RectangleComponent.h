#ifndef jura_RectangleComponent_h
#define jura_RectangleComponent_h

/** This class is component which has the sole purpose to draw itself as an rectangle filled with 
some colour and draw an outline with another colour. */

class JUCE_API RectangleComponent : public Component
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. Initializes the colours to be used for drawing the rectangle. */
  RectangleComponent(const Colour &fillColour = Colours::white, 
    const Colour &outlineColour = Colours::black, int outlineThickness = 2);

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the colour either for the filling or for the outline. @see ColourIds */
  //virtual void setColour(const int  colourId,  const Colour &colour);

  /** Sets the colour to fill the rectangle with. */
  virtual void setFillColour(const Colour &newColour);

  /** Sets the colour for the outline. */
  virtual void setOutlineColour(const Colour &newColour);

  /** Sets the thickness of the outline in pixels. */
  virtual void setOutlineThickness(int newThickness);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Paints the rectangle. */
  virtual void paint(Graphics &g);


protected:

  Colour fillColour, outlineColour;
  int outlineThickness;

  juce_UseDebuggingNewOperator;
};

#endif  