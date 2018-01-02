#ifndef RAPT_MAPPERS_H_INCLUDED
#define RAPT_MAPPERS_H_INCLUDED

template<class T>
class rsMapper
{
public:

  /** Subclasses need to override this function to perform the mapping. */
  virtual T map(T x) = 0;

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
  virtual void map(T *x, T *y) = 0;

};

//=================================================================================================

template<class T>
class rsCoordinateMapper : public rsMapper<T>
{

public:

  void setInputRange( T newMin, T newMax);
  void setOutputRange(T newMin, T newMax);
  void setLogScaled(bool shoulBeLogScaled);

  virtual T   map(T x) override;
  virtual T unmap(T x);

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

  virtual void   map(T *x, T *y) override;
  virtual void unmap(T *x, T *y);

  T   mapX(T x) { return mapperX.  map(x); }
  T   mapY(T y) { return mapperY.  map(y); }
  T unmapX(T x) { return mapperX.unmap(x); }
  T unmapY(T y) { return mapperY.unmap(y); }

protected:

  rsCoordinateMapper<T> mapperX, mapperY;

};

#endif
