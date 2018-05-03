#ifndef RS_RANGE_H
#define RS_RANGE_H

namespace RSLib
{

  /**

  This class can be used to represent ranges that begin at some min-value and end at some 
  max-value. The template parameter 'ValueType' must define an appropriate '-' operator 
  (subtraction) as well as comparison operators.

  \todo: express all comparisons internally in terms of the "<=" operator, such that the 
  ValueType only needs to define this single comparison operator

  */

  template<class ValueType>
  class rsRange
  {

  public:


    /** \name Construction/Destruction */

    /** Constructor. */
    rsRange(ValueType minValue, ValueType maxValue)
    {
      rsAssert(minValue <= maxValue);
      min = minValue;
      max = maxValue;
    }


    /** \name Setup */

    /** Sets the minimum value of this range. */
    void setMin(ValueType newMin)
    {
      rsAssert( newMin <= max );
      if( newMin <= max )
        min = newMin;
    }

    /** Sets the maximum value of this range. */
    void setMax(ValueType newMax)
    {
      rsAssert( newMax >= min );
      if( newMax >= min )
        max = newMax;
    }

    /** Sets the min- and max-value of this range. */
    void setMinAndMax(ValueType newMin, ValueType newMax)
    {
      rsAssert( newMax >= newMin );
      if( newMax >= newMin )
      {
        min = newMin;
        max = newMax;
      }
    }

    /** Lets the range represented by the member-variables be clipped to another range such that 
    the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
    void clipRange(const rsRange<ValueType>& rangeToClipTo)
    {
      if( min < rangeToClipTo.getMin() )
        min = rangeToClipTo.getMin();
      if( max > rangeToClipTo.getMax() )
        max = rangeToClipTo.getMax();
    }


    /** \name Inquiry */

    /** Returns the minimum value of this range. */
    ValueType getMin() const 
    { 
      return min; 
    }

    /** Returns the maximum value of this range. */
    ValueType getMax() const 
    { 
      return max; 
    }

    /** Returns the size of this range (defined as difference between max and min). */
    ValueType getSize() const 
    { 
      return max-min; 
    }


    /** \name Operators */

    /** Compares two ranges for equality. */
    bool operator==(const rsRange<ValueType>& r2) const
    {
      if( r2.min == min && r2.max == max )
        return true;
      else
        return false;
    }

    /** Compares two ranges for inequality. */
    bool operator!=(const rsRange<ValueType>& r2) const
    {
      return !(*this == r2);
    }

  protected:
        
    /** \name Data */

    ValueType min, max;

  };


  //===============================================================================================

  /**

  This class is used to represent 2-dimenstional ranges, for example minimum and maximum values for 
  x and y in a coordinate system.

  */

  template<class ValueType>
  class rsRangeXY
  {

  public:


    /** \name Construction/Destruction */

    /** Constructor. */
    rsRangeXY(ValueType minX, ValueType maxX, ValueType minY, ValueType maxY)
      : rangeX(minX, maxX), rangeY(minY, maxY)
    {

    }


    /** \name Setup */

    void setMinX(ValueType newMinX) { rangeX.setMin(newMinX); }
    void setMaxX(ValueType newMaxX) { rangeX.setMax(newMaxX); }
    void setMinY(ValueType newMinY) { rangeY.setMin(newMinY); }
    void setMaxY(ValueType newMaxY) { rangeY.setMax(newMaxY); }

    void setRangeX(ValueType newMinX, ValueType newMaxX) { rangeX.setMinAndMax(newMinX, newMaxX); }
    void setRangeY(ValueType newMinY, ValueType newMaxY) { rangeY.setMinAndMax(newMinY, newMaxY); }

    /** Sets up all 4 value at once. */
    void setRange(ValueType newMinX, ValueType newMaxX, ValueType newMinY, ValueType newMaxY)
    {
      setRangeX(newMinX, newMaxX);
      setRangeY(newMinY, newMaxY);
    }

    /** Lets the range represented by the member-variables be clipped to another range such that 
    the minima/maxima are no smaller/larger than those of the rangeToClipTo. */
    void clipRange(const rsRangeXY<ValueType>& rangeToClipTo)
    {
      rangeX.clipTo(rangeToClipTo.rangeX);
      rangeY.clipTo(rangeToClipTo.rangeY);
    }


    /** \name Inquiry */

    ValueType getMinX() const { return rangeX.getMin(); }
    ValueType getMaxX() const { return rangeX.getMax(); }
    ValueType getMinY() const { return rangeY.getMin(); }
    ValueType getMaxY() const { return rangeY.getMax(); }


    /** \name Operators */

    /** Compares two ranges for equality. */
    bool operator==(const rsRangeXY<ValueType>& r2) const
    {
      if( r2.rangeX == rangeX && r2.rangeY == rangeY )
        return true;
      else
        return false;
    }

    /** Compares two ranges for inequality. */
    bool operator!=(const rsRangeXY<ValueType>& r2) const
    {
      return !(*this == r2);
    }

  protected:

    /** \name Data */

    rsRange<ValueType> rangeX, rangeY;

  };

}

#endif
