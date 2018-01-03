#ifndef RAPT_INTERPOLATINGFUNCTION_H_INCLUDED
#define RAPT_INTERPOLATINGFUNCTION_H_INCLUDED

/** A class for representing a node in a piecewise defined function. Each node has x- and 
y-coodinates and a shape option that determines the shape of the function segment approaching the 
node (i.e. the curve segment to the left of the node) */

template<class T>
class rsFunctionNode
{

public:

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

  /** A node is considered to be less than another, if its x-coordinate is less and in case of 
  equal x-coordinates, if its y-coordinate is less. */
  bool operator<(const rsFunctionNode& rhs) const
  {
    if(x < rhs.x)
      return true;
    else if(x > rhs.x)
      return false;
    else
      return y < rhs.y; 
  }


//protected:

  T x = 0, y = 0;
  int shapeType = LINEAR;
  T shapeParam = 0; // maybe this could be the derivative? of y

  //friend class rsNodeBasedFunction<T>;
};

//=================================================================================================

/** A class for representing a function that is defined via pairs of (x,y) datapoints. You may set
up these datatpoints and retrieve function values at arbitrary positions which are generated via
interpolating between the known datapoints.

-handle endpoints by either clamping (like now), extrapolation or assuming periodicity
*/

template<class T>
class rsNodeBasedFunction
{

public:

  /** Adds a new node at the given coordinates and returns the index at which it was inserted. */
  virtual size_t addNode(T x, T y);

  /** Removes the node with given index. */
  virtual void removeNode(size_t index);

  /** Moves an existing datapoint with given index to a new position. Because we always keep our 
  data arrays sorted, this may change the index of the datapoint inside the array. The return value
  informs about the new index. */
  virtual size_t moveNode(size_t index, T newX, T newY);

  /** Returns a reference to our array of nodes. It's a constant reference because client code
  is not allowed to edit that data directly. Instead, it must use the moveNode function which
  will update the datapoint and do some additional stuff. */ 
  const std::vector<rsFunctionNode<T>>& getNodes() { return nodes; }

  /** Returns true, if the datapoint at index i+1 is considered to be "less than" the datapoint at 
  index i. We use this function internally to keep our arrays of values sorted. Datapoints are 
  sorted according to ascending x-values and in case of equal x-values, according to ascending 
  y-values. The caller should ensure that i <= N-2. */
  bool isNextValueLess(size_t i) { return nodes[i+1] < nodes[i]; }

  /** Returns the first index i in our x-array such that x[i] > xToFind. */
  size_t firstIndexOfGreaterX(T xToFind)
  {
    for(size_t i = 0; i < nodes.size(); i++)
      if(nodes[i].x > xToFind)
        return i;
    return nodes.size();
  }
  // todo: use binary search with a start-index based on the previously retrieved value
  // see https://www.geeksforgeeks.org/binary-search/


  /** For internal use only... */
  T getValueLinear(T x, size_t i)
  {
    T x1 = nodes[i].x;
    T y1 = nodes[i].y;
    T x2 = nodes[i+1].x;
    T y2 = nodes[i+1].y;
    T thresh = T(1.e-10); // todo: use epsilon of T
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);
    return y1 + (y2-y1) * (x-x1) / (x2-x1);
  }

  /** Returns an interpolated y-value at the given value of x. */
  T getValue(T x)
  {
    if(nodes.size() == 0)
      return 0;
    if(x < nodes[0].x)
      return nodes[0].y;
    size_t i = firstIndexOfGreaterX(x);
    if(i == nodes.size())
      return nodes[i-1].y;
    return getValueLinear(x, i-1);
  }

protected:

  std::vector<rsFunctionNode<T>> nodes;

};

#endif