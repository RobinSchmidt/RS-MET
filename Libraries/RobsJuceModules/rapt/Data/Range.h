#pragma once

/** This class can be used to represent ranges that begin at some min-value and end at some
max-value. The template class must define an appropriate '-' operator (subtraction) as well as 
comparison operators.

\todo: express all comparisons internally in terms of the "<=" operator, such that the
type only needs to define this single comparison operator  */

template<class T>
class rsRange
{

public:


  /** \name Construction/Destruction */

  /** Constructor. You may pass the (initial) values for min and max - if you don't pass anything, 
  it will default to the empty range/interval from 0 to 0 (i.e. min = max = 0). */
  rsRange(T minValue = T(0), T maxValue = T(0))
  {
    rsAssert(minValue <= maxValue);
    min = minValue;
    max = maxValue;
  }


  /** \name Setup */

  /** Sets the minimum value of this range. */
  void setMin(T newMin)
  {
    rsAssert(newMin <= max);
    if(newMin <= max)
      min = newMin;
  }

  /** Sets the maximum value of this range. */
  void setMax(T newMax)
  {
    rsAssert(newMax >= min);
    if(newMax >= min)
      max = newMax;
  }

  /** Sets the min- and max-value of this range. */
  void setMinAndMax(T newMin, T newMax)
  {
    rsAssert(newMax >= newMin);
    if(newMax >= newMin)
    {
      min = newMin;
      max = newMax;
    }
  }

  /** Lets the range represented by the member-variables be clipped to another range such that
  the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
  void clipRange(const rsRange<T>& rangeToClipTo)
  {
    if(min < rangeToClipTo.getMin())
      min = rangeToClipTo.getMin();
    if(max > rangeToClipTo.getMax())
      max = rangeToClipTo.getMax();
  }


  /** \name Inquiry */

  /** Returns the minimum value of this range. */
  T getMin() const { return min; }

  /** Returns the maximum value of this range. */
  T getMax() const { return max; }

  /** Returns the size of this range (defined as difference between max and min). */
  T getSize() const { return max-min; }


  /** \name Operators */

  /** Compares two ranges for equality. */
  bool operator==(const rsRange<T>& r2) const
  {
    return r2.min == min && r2.max == max

    /*
    if(r2.min == min && r2.max == max)
      return true;
    else
      return false;
    */
  }

  /** Compares two ranges for inequality. */
  bool operator!=(const rsRange<T>& r2) const { return !(*this == r2); }

  /** Returns true, if the size (i.e. the difference between max and min) of the left operand is
  greater than the size of the right operand. */
  bool operator>(const rsRange<T>& r2) const { return getSize() > r2.getSize(); }


protected:

  /** \name Data */

  T min, max;

};


//===============================================================================================

/** This class is used to represent 2-dimensional ranges, for example minimum and maximum values 
for x and y in a coordinate system. */

template<class T>
class rsRangeXY
{

public:


  /** \name Construction/Destruction */

  /** Constructor. */
  rsRangeXY(T minX, T maxX, T minY, T maxY) : rangeX(minX, maxX), rangeY(minY, maxY) {}


  /** \name Setup */

  void setMinX(T newMinX) { rangeX.setMin(newMinX); }
  void setMaxX(T newMaxX) { rangeX.setMax(newMaxX); }
  void setMinY(T newMinY) { rangeY.setMin(newMinY); }
  void setMaxY(T newMaxY) { rangeY.setMax(newMaxY); }

  void setRangeX(T newMinX, T newMaxX) { rangeX.setMinAndMax(newMinX, newMaxX); }
  void setRangeY(T newMinY, T newMaxY) { rangeY.setMinAndMax(newMinY, newMaxY); }

  /** Sets up all 4 value at once. */
  void setRange(T newMinX, T newMaxX, T newMinY, T newMaxY)
  {
    setRangeX(newMinX, newMaxX);
    setRangeY(newMinY, newMaxY);
  }

  /** Lets the range represented by the member-variables be clipped to another range such that
  the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
  void clipRange(const rsRangeXY<T>& rangeToClipTo)
  {
    rangeX.clipTo(rangeToClipTo.rangeX);
    rangeY.clipTo(rangeToClipTo.rangeY);
  }


  /** \name Inquiry */

  T getMinX() const { return rangeX.getMin(); }
  T getMaxX() const { return rangeX.getMax(); }
  T getMinY() const { return rangeY.getMin(); }
  T getMaxY() const { return rangeY.getMax(); }


  /** \name Operators */

  /** Compares two ranges for equality. */
  bool operator==(const rsRangeXY<T>& r2) const
  {
    return r2.rangeX == rangeX && r2.rangeY == rangeY;

    /*
    if(r2.rangeX == rangeX && r2.rangeY == rangeY)
      return true;
    else
      return false;
    */
  }

  /** Compares two ranges for inequality. */
  bool operator!=(const rsRangeXY<T>& r2) const
  {
    return !(*this == r2);
  }

protected:

  /** \name Data */

  rsRange<T> rangeX, rangeY;

};

