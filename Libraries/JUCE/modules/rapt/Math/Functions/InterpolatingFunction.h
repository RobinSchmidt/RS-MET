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

  /** Moves an existing datapoint with given index to a new position. */
  void moveDataPoint(size_t index, T newX, T newY);


  T getValue(T x)
  {
    return x; // preliminary
  }

protected:

  std::vector<T> x, y;

};

#endif