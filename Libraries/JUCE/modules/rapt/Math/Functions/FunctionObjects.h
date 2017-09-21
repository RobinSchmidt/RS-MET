#ifndef RAPT_FUNCTIONOBJECTS_H_INCLUDED
#define RAPT_FUNCTIONOBJECTS_H_INCLUDED

// This file contains abstract baseclasses for function objects (a.k.a. Functors). But we use the 
// term in a loose sense in that the function-application is not necessarrily be done by the
// () operator.


/** Baseclass for mapping 2-dimensional vectors, represented as pairs of x- and y-coordinates, to 
a new location. */

template<class T>
class rsMapper2D
{

public:
	//virtual ~Mapper2D() = default;
  virtual ~Mapper2D() {}

  /** Subclasses must override this function to map the incoming xy-pair to the corresponding 
  outgoing pair. */
  virtual void map(T *x, T *y) = 0;

};

#endif
