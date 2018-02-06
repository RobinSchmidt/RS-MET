#ifndef RAPT_MAPPERS_H_INCLUDED
#define RAPT_MAPPERS_H_INCLUDED

template<class T>
class rsMapper
{
public:

  /** Subclasses need to override this function to perform the mapping. */
  virtual T map(T x) const = 0;

  // maybe implement () operator (invkoes map)
};

//=================================================================================================

/** Baseclass for mapping 2-dimensional vectors, represented as pairs of x- and y-coordinates, to 
a new location. */

template<class T>
class rsMapper2D
{

public:
  //virtual ~rsMapper2D() = default;
  virtual ~rsMapper2D() {}

  /** Subclasses must override this function to map the incoming xy-pair to the corresponding 
  outgoing pair. */
  virtual void map(T *x, T *y) const = 0;

};

//=================================================================================================

template<class T>
class rsCoordinateMapper : public rsMapper<T>
{

public:

  void setInputRange( T newMin, T newMax);
  void setOutputRange(T newMin, T newMax);
  void setLogScaled(bool shoulBeLogScaled);

  virtual T   map(T x) const override;
  virtual T unmap(T x) const;

  inline T getInMin()  const { return inMin;  }
  inline T getInMax()  const { return inMax;  }
  inline T getOutMin() const { return outMin; }
  inline T getOutMax() const { return outMax; }

protected:

  T inMin = 0, inMax = 1, outMin = 0, outMax = 1;
  bool logScaled = false;
};

//=================================================================================================

/** A class for mapping 2D coordinates, for example to convert between pixel coordinates and the
underlying model coordinates. */

template<class T>
class rsCoordinateMapper2D : public rsMapper2D<T>
{

public:

  void setInputRange( T minX, T maxX, T minY, T maxY);
  void setOutputRange(T minX, T maxX, T minY, T maxY);

  virtual void   map(T *x, T *y) const override;
  virtual void unmap(T *x, T *y) const;

  inline T   mapX(T x) const { return mapperX.  map(x); }
  inline T   mapY(T y) const { return mapperY.  map(y); }
  inline T unmapX(T x) const { return mapperX.unmap(x); }
  inline T unmapY(T y) const { return mapperY.unmap(y); }

  inline T getInMinX()  const { return mapperX.getInMin();  }
  inline T getInMaxX()  const { return mapperX.getInMax();  }
  inline T getOutMinX() const { return mapperX.getOutMin(); }
  inline T getOutMaxX() const { return mapperX.getOutMax(); }

  inline T getInMinY()  const { return mapperY.getInMin();  }
  inline T getInMaxY()  const { return mapperY.getInMax();  }
  inline T getOutMinY() const { return mapperY.getOutMin(); }
  inline T getOutMaxY() const { return mapperY.getOutMax(); }


//protected:

  rsCoordinateMapper<T> mapperX, mapperY;

};

#endif
