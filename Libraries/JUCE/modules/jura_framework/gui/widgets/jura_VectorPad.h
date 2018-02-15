#ifndef jura_VectorPad_h
#define jura_VectorPad_h  


/** This is a class for controlling two parameters at once with the position of a dot in a 2D 
coordinate system.

todo: let client code restrict the usable range such that the extreme values of the parameters
do not necessarily coincide with the component edges but are somehwat more inside. This is useful
hen a coordinate system is displayed in the background whose extreme values should not be somehwat
beyond the parameter range limits. */

class JUCE_API rsVectorPad : public RWidget
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Construction/Destruction:
  
  rsVectorPad();
  virtual ~rsVectorPad();

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  virtual void assignParameterX(Parameter* newParameterX);
  
  virtual void assignParameterY(Parameter* newParameterY);

  /** Sets margins in pixels (with respect to the Components bounds) for which the min/max values
  of the parameter will be attained. For example with all margins equal to 10, the minimum value
  for x will be attained at x = 10 pixels and the maximum at x = width-10 pxiels. */
  void setMargins(double left, double right, double top, double bottom)
  {
    leftMargin   = left;
    rightMargin  = right;
    topMargin    = top;
    bottomMargin = bottom;
  }

  /** Sets the size of the dot to be drawn at the current x/y position. */
  void setDotSize(float newSize) { dotSize = newSize; }

  /** Decides whether or not the background should be filled with the background color (inherited
  for RWidget). You may want to set this to false, if you want to show a plot behind the handle. 
  When false, the background will be transparent and a plot behind the widget will show 
  through. */
  void setPaintBackground(bool shouldPaint) { paintBackground = shouldPaint; }

  //-----------------------------------------------------------------------------------------------
  // \name Callbacks

  virtual void parameterChanged(Parameter* parameterThatHasChanged) override;
  virtual void paint(Graphics& g) override;
  virtual void mouseDown(const MouseEvent& e) override;
  virtual void mouseDrag(const MouseEvent& e) override;



  
protected:

  // conversion between pixel coordinates and normalized parameter value:
  virtual double pixelToNormalizedX(double x);
  virtual double pixelToNormalizedY(double y);
  virtual double normalizedToPixelX(double x);
  virtual double normalizedToPixelY(double y);



  void setParametersFromMouseEvent(const MouseEvent& e);

  void setParametersXY(double x, double y);

  Parameter *dummyParam; // so we can use it without assigned parameters

  Parameter *paramX = nullptr, *paramY = nullptr;

  double xMin = -1, xMax = +1, yMin = -1, yMax = +1;

  double leftMargin = 0, rightMargin = 0, topMargin = 0, bottomMargin = 0;

  float dotSize = 16;
  bool paintBackground = true;

  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(rsVectorPad)
};


#endif
