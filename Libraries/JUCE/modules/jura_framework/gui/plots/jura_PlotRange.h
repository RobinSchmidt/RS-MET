#ifndef jura_PlotRange_h
#define jura_PlotRange_h

/** This class represents the visible range for a plot. */

class JUCE_API rsPlotRange
{

public:

  /** Standard Constructor. */
  rsPlotRange(double minX = -2.2, double maxX = +2.2, double minY = -2.2, double maxY = +2.2);
    // maybe use -1 and +1 for defaults, or 0 and 1

  /** Returns the minimum x-value. */
  double getMinX() const { return minX; }

  /** Returns the maximum x-value. */
  double getMaxX() const { return maxX; }

  /** Returns the minimum y-value. */
  double getMinY() const { return minY; }

  /** Returns the maximum y-value. */
  double getMaxY() const { return maxY; }

  /** Sets the minimum x-value. */
  void setMinX(double newMinX);

  /** Sets the maximum x-value. */
  void setMaxX(double newMaxX);

  /** Sets up the minimum and maximum value for x .*/
  void setRangeX(double newMinX, double newMaxX);

  /** Sets the minimum y-value. */
  void setMinY(double newMinY);

  /** Sets the maximum y-value. */
  void setMaxY(double newMaxY); 

  /** Sets up the minimum and maximum value for y .*/
  void setRangeY(double newMinY, double newMaxY);

  /** Lets the range represented by the member-variables be clipped to another range such that
  the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
  void clipRange(rsPlotRange rangeToClipTo);

  /** Compares two ranges for equality. */
  bool operator==(const rsPlotRange& r2) const  
  {
    if( r2.minX == minX && r2.maxX == maxX && r2.minY == minY && r2.maxY == maxY )
      return true;
    else
      return false;
  }

  /** Compares two ranges for inequality. */
  bool operator!=(const rsPlotRange& r2) const  
  {
    return !(*this == r2);
  }

protected:

  double minX, maxX, minY, maxY;

  juce_UseDebuggingNewOperator;
};

#endif