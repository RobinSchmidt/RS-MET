#ifndef RAPT_INTERPOLATINGFUNCTION_H_INCLUDED
#define RAPT_INTERPOLATINGFUNCTION_H_INCLUDED

/** A class for representing a node in a piecewise defined function. Each node has x- and 
y-coodinates and a shape option that determines the shape of the function segment approaching the 
node (i.e. the curve segment to the left of the node) */

template<class T>
class rsFunctionNode
{

  /** Enumeration of possible shapes. */
  enum shapes
  {
    LEFT_NEIGHBOUR,
    RIGHT_NEIGHBOUR,
    NEAREST_NEIGHBOUR,
    LINEAR,             // this is the default
    CUBIC
  };
  // maybe have log, exp, pow shapes
  // pow: y = a + b * (x+c)^d

  /** Constructor. Initializes the x,y coordinates to the given values and uses a linear shape by 
  default. */
  rsFunctionNode(T _x, T _y) : x(_x), y(_y) {}

  /** Sets the x,y coordinates of this node. */
  inline void setCoordinates(T newX, T newY) { x = newX; y = newY; }

  /** Sets the x-coordinate of this node. */
  inline void setX(T newX) { x = newX; }

  /** Sets the y-coordinate of this node. */
  inline void setY(T newY) { y = newY; }

  /** Sets the shape as one of the values in he shapes enum. This determines the shape of the line
  segment towards this node. */
  inline void setShapeType(int newType) { shapeType = newType; }

  /** Some shapes have a numeric parameter which can be set via this function. */
  inline void setShapeParameter(int newParameter) { shapeParam = newParam; }

  /** Returns the x-coordinate of this node. */
  inline T getX() { return x; }

  /** Returns the y-coordinate of this node. */
  inline T getY() { return y; }


protected:

  T x = 0, y = 0;
  int shapeType = LINEAR;
  T shapeParam = 0; // maybe this could be the derivative? of y

};

//=================================================================================================

/** A class for representing a function that is defined via pairs of (x,y) datapoints. You may set
up these datatpoints and retrieve function values at arbitrary positions which are generated via
interpolating between the known datapoints.

-maybe rename to rsNodeBasedFunction and use rsFunctionNode class
-handle endpoints by either clamping (like now), extrapolation or assuming periodicity
*/

template<class T>
class rsInterpolatingFunction
{

public:


  /** Adds a new datapoint at the given coordinates and returns the index at which it was 
  inserted. */
  size_t addDataPoint(T x, T y);

  /** Removes the datapoint with givne index. */
  void removeDataPoint(size_t index);

  /** Moves an existing datapoint with given index to a new position. Because we always keep our 
  data arrays sorted, this may change the index of the datapoint inside the array. The return value
  informs about the new index. */
  size_t moveDataPoint(size_t index, T newX, T newY);

  /** Returns a reference to our array of x-values. It's a constant reference because client code
  is not allowed to edit that data directly. Instead, it must use the moveDataPoint function which
  will update the datapoint and do some additional stuff. */ 
  const std::vector<T>& getValuesX() { return xValues; }

  /** Returns a (constant) reference to our array of x-values. @see getValuesX */
  const std::vector<T>& getValuesY() { return yValues; }


  /** Returns true, if the datapoint at index i+1 is considered to be "less than" the datapoint at 
  index i. We use this function internally to keep our arrays of values sorted. Datapoints are 
  sorted according to ascending x-values and in case of equal x-values, according to ascending 
  y-values. The caller should ensure that i <= N-2. */
  bool isNextValueLess(size_t i)
  {
    if(xValues[i+1] < xValues[i])
      return true;
    else if(xValues[i+1] > xValues[i])
      return false;
    else
      return yValues[i+1] < yValues[i]; // if x-values are equal, compare y-values
  }


  /** Returns the first index i in our x-array such that x[i] > xToFind. */
  size_t firstIndexOfGreaterX(T xToFind)
  {
    for(size_t i = 0; i < xValues.size(); i++)
      if(xValues[i] > xToFind)
        return i;
    return xValues.size();
  }
  // todo: use binary search with a start-index based on the previously retrieved value
  // see https://www.geeksforgeeks.org/binary-search/


  /** For internal use only... */
  T getValueLinear(T x, size_t i)
  {
    T x1 = xValues[i];
    T y1 = yValues[i];
    T x2 = xValues[i+1];
    T y2 = yValues[i+1];
    T thresh = T(1.e-10); // todo: use epsilon of T
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);
    return y1 + (y2-y1) * (x-x1) / (x2-x1);
  }

  /** Returns an interpolated y-value at the given value of x. */
  T getValue(T x)
  {
    if(xValues.size() == 0)
      return 0;
    if(x < xValues[0])
      return yValues[0];
    size_t i = firstIndexOfGreaterX(x);
    if(i == xValues.size())
      return yValues[i-1];
    return getValueLinear(x, i-1);
  }

protected:

  std::vector<T> xValues, yValues;

};

#endif