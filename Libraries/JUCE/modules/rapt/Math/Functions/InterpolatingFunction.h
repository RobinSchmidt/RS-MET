#ifndef RAPT_INTERPOLATINGFUNCTION_H_INCLUDED
#define RAPT_INTERPOLATINGFUNCTION_H_INCLUDED

/** A class for representing a function that is defined via pairs of (x,y) datapoints. You may set
up these datatpoints and retrieve function values at arbitrary positions which are generated via
interpolating between the known datapoints .*/

template<class T>
class rsInterpolatingFunction
{

public:


  /** Adds a new datapoint at the given coordinates. */
  void addDataPoint(T x, T y);

  /** Removes the datapoint with givne index. */
  void removeDataPoint(size_t index);

  /** Moves an existing datapoint with given index to a new position. Note that this may change the
  index in the array. We keep our datapoint arrays sorted according to ascending x-values. */
  void moveDataPoint(size_t index, T newX, T newY);

  /** Returns a reference to our array of x-values. It's a constant reference because client code
  is not allowed to edit that data directly. Instead, it must use the moveDataPoint function which
  will update the datapoint and do some additional stuff. */ 
  const std::vector<T>& getValuesX() { return xValues; }

  /** Returns a (constant) reference to our array of x-values. @see getValuesX */
  const std::vector<T>& getValuesY() { return yValues; }

  /** Returns the first index i in our x-array such that x[i] > xToFind. */
  size_t firstIndexOfGreaterX(T xToFind)
  {
    for(size_t i = 0; i < xValues.size(); i++)
      if(xValues[i] > xToFind)
        return i;
    return xValues.size()-1;
  }
  // todo: use binary search with a start-index based on the previously retrieved value


  /** For internal use only... */
  T getValueLinear(T x, size_t i)
  {
    T x1 = xValues[i];
    T y1 = yValues[i];
    T x2 = xValues[i+1];
    T y2 = yValues[i+1];
    T t  = x - x1;
    T thresh = 1.e-13; // todo: use epsilon of T
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);
    return y1 + t * (y2-y1) / (x2-x1); // check this formula
  }

  /** Returns an interpolated y-value at the given value of x. */
  T getValue(T x)
  {
    //return x; // preliminary

    if(xValues.size() == 0)
      return 0;
    if(x < xValues[0])
      return xValues[0];
    size_t i = firstIndexOfGreaterX(x);
    if(i == xValues.size()-1)
      return xValues[i];
    return getValueLinear(x, i-1);
  }

protected:

  std::vector<T> xValues, yValues;

};

#endif