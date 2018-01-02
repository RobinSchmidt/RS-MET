#ifndef RAPT_MAPPERS_H_INCLUDED
#define RAPT_MAPPERS_H_INCLUDED

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

/** A class for mapping 2D coordinates, for example to convert between pixel coordinates and the
underlying model coordinates. */

template<class T>
class rsCoordinateMapper2D : public rsMapper2D<T>
{

public:

  void setInputRange(T inMinX, T inMaxX, T inMinY, T inMaxY);

  void setOutputRange(T outMinX, T outMaxX, T outMinY, T outMaxY);



  virtual void map(T *x, T *y) override;

  virtual void unmap(T *x, T *y);



protected:

  void updateCoeffs();

  T minX = 0, minY = 0, maxX = 1, maxY = 1;
  bool logX, logY;



};

#endif
