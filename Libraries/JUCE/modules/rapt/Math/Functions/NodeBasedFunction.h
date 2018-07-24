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
  enum  
  {
    LEFT_NEIGHBOUR,
    RIGHT_NEIGHBOUR,
    NEAREST_NEIGHBOUR,
    LINEAR,             // this is the default
    CUBIC,
    EXPONENTIAL,
    RATIONAL
  };
  // maybe have log, exp, pow shapes
  // pow: y = a + b * (x+c)^d

  /** Constructor. Initializes the x,y coordinates to the given values and uses a linear shape by 
  default. */
  //rsFunctionNode(T _x, T _y) : x(_x), y(_y) {}
  rsFunctionNode(T _x, T _y, int _shape = LINEAR, T _shapeParam = T(0)) 
    : x(_x), y(_y), shapeType(_shape), shapeParam(_shapeParam) 
  {
  
  }


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
  inline void setShapeParameter(T newParam) { shapeParam = newParam; }

  /** Returns the x-coordinate of this node. */
  inline T getX() const { return x; }

  /** Returns the y-coordinate of this node. */
  inline T getY() const { return y; }

  /** Returns the index of the shape type for this node, @see shapes.  */
  inline int getShapeType() const { return shapeType; }

  /** Returns the numerical parameter that controls the shape. */
  inline T getShapeParameter() const { return shapeParam; }

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
  
  virtual ~rsNodeBasedFunction() {}

  //-----------------------------------------------------------------------------------------------
  // \name Node manipulation

  /** Adds a new node at the given coordinates and returns the index at which it was inserted. */
  virtual size_t addNode(T x, T y);

  /** Tries to remove the node with given index and returns true, if it was actually removed (it 
  may not be removed, if the overriden isNodeRemovable function in a subclass returns false
  for the node in question). */
  virtual bool removeNode(size_t index);

  /** Moves an existing datapoint with given index to a new position. Because we always keep our 
  data arrays sorted, this may change the index of the datapoint inside the array. The return value
  informs about the new index. */
  virtual size_t moveNode(size_t index, T newX, T newY);

  /** This function is called after the node at the given index has been moved to a new x,y 
  position. It will move the node up or down in the array of nodes (by successive swaps) to keep 
  the array sorted. The return value is the new index. */
  size_t moveNodeToSortedIndex(size_t index);

  /** Sets the shape type for the node at given index. */
  virtual void setNodeShapeType(size_t index, int newType) { nodes[index].setShapeType(newType); }

  /** Sets the shape parameter for the node at given index. */
  virtual void setNodeShapeParameter(size_t index, T newParam) 
  { 
    nodes[index].setShapeParameter(newParam); 
  }

  //-----------------------------------------------------------------------------------------------
  // \name Hooks

  /** This function will be called before an attempt to remove a node and will not remove it, if 
  that function returns false. The baseclass implementation just returns true but you can override 
  it in a subclass if your subclass - for example - requires a certain minimum number of nodes. */
  virtual bool isNodeRemovable(size_t /*index*/) { return true; }

  // todo: maybe add a canNodeBeAdded(T x, T y) function. subclasses my have a maximum number of
  // nodes and/or disallow adding nodes outside a given (x,y) range

  /** This function is called after a node has been inserted or moved to a new position. You can
  override it, if you need to apply constraints to the positions of nodes, like for example that 
  the coordinates of the first and/or last node must have certain values. Because changing the 
  position of a node may change its index, too, this function should return the new index of the 
  node after the constraints have been applied. The baseclass implementation does nothing and just 
  returns the same index that you give to it. */
  virtual size_t constrainNode(size_t index) { return index; }

  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the x coordinate of the node with given index. */
  inline T getNodeX(size_t index) { return nodes[index].getX(); }

  /** Returns the y coordinate of the node with given index. */
  inline T getNodeY(size_t index) { return nodes[index].getY(); }

  inline int getNodeShapeType(size_t index) { return nodes[index].getShapeType(); }

  inline T getNodeShapeParameter(size_t index) { return nodes[index].getShapeParameter(); }


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

  //-----------------------------------------------------------------------------------------------
  // \name Output computation

  /** Returns an interpolated y-value at the given value of x. */
  T getValue(T x)
  {
    if(nodes.size() == 0)
      return 0;

    // todo: switch endpoint-mode (clamp, extrapolate, periocic, ...):
    if(x < nodes[0].x)
      return nodes[0].y;
    size_t i = firstIndexOfGreaterX(x);
    if(i == nodes.size())
      return nodes[i-1].y;

    typedef rsFunctionNode<T> FN;
    switch(nodes[i].shapeType)
    {
    //case FN::LINEAR:      return getValueLinear(x, i-1);  // redundant with default case
    case FN::EXPONENTIAL: return getValueExponential(x, i-1);
    case FN::RATIONAL:    return getValueRational(   x, i-1);
    default:              return getValueLinear(     x, i-1);
    }

    //return getValueLinear(x, i-1);
  }
  // rename to applyFunction

  /** Tries to invert the function by finding an x-value for which this function produces the 
  given y value. The implementation is based on root-finding and will succeed only, if the given y 
  is within the range of the function. If there are several possible x-values that yield the given 
  y value (as could be the case for nonmonotonic functions), the first (leftmost) x-value will be
  returned. */
  T applyInverseFunction(T y);

protected:

  /** For internal use only... */
  T getValueLinear(T x, size_t i)
  {
    T x1 = nodes[i].x;
    T y1 = nodes[i].y;
    T x2 = nodes[i+1].x;
    T y2 = nodes[i+1].y;
    //T thresh = T(1.e-10); // todo: use epsilon of T
    T thresh = RS_EPS(T);
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);
    return y1 + (y2-y1) * (x-x1) / (x2-x1); // factor into function linCurve
  }

  T getValueExponential(T x, size_t i)
  {
    T thresh = RS_EPS(T);

    T c = 0.5 * (nodes[i+1].shapeParam + 1);
    T a = 2*log((1-c)/c);
    //T a = -nodes[i+1].shapeParam;  // alpha - maybe use shapeScaler*shapeParam where shapeScaler is a 

    if(abs(a) < thresh)            // global setting for the whole function defined in this class
      return getValueLinear(x, i); // avoid div-by-zero

    // this is the same as in the linear case (factor out, if possible):
    T x1 = nodes[i].x;
    T y1 = nodes[i].y;
    T x2 = nodes[i+1].x;
    T y2 = nodes[i+1].y;
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);

    // Formulas taken from Elements of Computer Music (Moore), page 184:
    T I = (x-x1) / (x2-x1);                          // Eq 3.30
    return y1 + (y2-y1) * (1-exp(I*a)) / (1-exp(a)); // Eq 3.29 - factor into function expCurve(I, a)
  }
  // todo: interpret the shape parameter not directly as "a" - instead compute a from the condition
  // (1-exp(a/2)) / (1-exp(a)) = c  -> wolfram: solve (1-exp(a/2))/(1-exp(a)) = c for a
  // -> a = 2*log((1-c)/c) where c is the value between 0..1 of the
  // normalized transition function (between 0..1)...so, if c = 0.75, it means the normalized curve
  // goes through the point (x = 0.5,y = 0.75), with c = 0.25, it goes through (x = 0.5, y = 0.25)
  // so c is y-value in the middle of the transition function
  // -> use the same convention (shape parameter is normalized curve value at x=0.5) also for the 
  // rational mapping


  T ratCurve(T p, T a)
  {
    T ap = a*p;
    return (ap-p) / (2*ap - a - 1);
  }
  T getValueRational(T x, size_t i)
  {
    T thresh = RS_EPS(T);
    T c = 0.5 * (nodes[i+1].shapeParam + 1);
    T a = (2*c)/(2*c+1);

    // this is the same as in the linear case (factor out, if possible):
    T x1 = nodes[i].x;
    T y1 = nodes[i].y;
    T x2 = nodes[i+1].x;
    T y2 = nodes[i+1].y;
    if(fabs(x2-x1) < thresh)
      return T(0.5) * (y1+y2);

    T p = (x-x1) / (x2-x1);
    return y1 + (y2-y1) * ratCurve(p, a);
  }






  /** Simply appends a node with given coordinates and linear shape-type to our array - without 
  checking constraints or sorting. This is mainly meant to be used by subclasses to intialize the
  node array to a default state (in which case the constraints may not work as desired). */
  //void appendNode(T x, T y) { nodes.push_back(rsFunctionNode<T>(x, y)); }
  void appendNode(T x, T y, int shape = rsFunctionNode<T>::LINEAR, T shapeParam = T(0)) 
  { 
    nodes.push_back(rsFunctionNode<T>(x, y, shape, shapeParam)); 
  }



  std::vector<rsFunctionNode<T>> nodes;

};

#endif