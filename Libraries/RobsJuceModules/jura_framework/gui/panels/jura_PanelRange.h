#ifndef jura_PanelRange_h
#define jura_PanelRange_h

/**  This class represents the visible range for a Panel (and is used inside there). */

class PanelRange
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  PanelRange(double initMinX = -2.2, double initMaxX = +2.2,
    double initMinY = -2.2, double initMaxY = +2.2);

  /** Destructor. */
  virtual ~PanelRange();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the minimum x-value. */
  virtual void setMinX(double newMinX);

  /** Sets the maximum x-value. */
  virtual void setMaxX(double newMaxX);

  /** Sets up the minimum and maximum value for x .*/
  virtual void setRangeX(double newMinX, double newMaxX);

  /** Sets the minimum y-value. */
  virtual void setMinY(double newMinY);

  /** Sets the maximum y-value. */
  virtual void setMaxY(double newMaxY);

  /** Sets up the minimum and maximum value for y .*/
  virtual void setRangeY(double newMinY, double newMaxY);

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the minimum x-value. */
  virtual double getMinX() const { return minX; }

  /** Returns the maximum x-value. */
  virtual double getMaxX() const { return maxX; }

  /** Returns the minimum y-value. */
  virtual double getMinY() const { return minY; }

  /** Returns the maximum y-value. */
  virtual double getMaxY() const { return maxY; }

  //---------------------------------------------------------------------------------------------
  // others:

  /** Lets the range represented by the member-variables be clipped to another range such that
  the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
  virtual void clipRange(PanelRange rangeToClipTo);

  //---------------------------------------------------------------------------------------------
  // operators:

  /** Compares two PanelRanges of equality. */
  bool operator==(const PanelRange& otherRange) const
  {
    if(minX != otherRange.minX)
      return false;
    if(maxX != otherRange.maxX)
      return false;
    if(minY != otherRange.minY)
      return false;
    if(maxY != otherRange.maxY)
      return false;
    return true;
  }

  /** Compares two PanelRanges of inequality. */
  bool operator!=(const PanelRange& otherRange) const
  {
    return !(*this == otherRange);
  }

protected:

  double minX, maxX, minY, maxY;

  juce_UseDebuggingNewOperator;
};

#endif  
