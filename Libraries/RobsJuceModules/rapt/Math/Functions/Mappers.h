#ifndef RAPT_MAPPERS_H_INCLUDED
#define RAPT_MAPPERS_H_INCLUDED

template<class T>
class rsMapper
{
public:

  /** Subclasses need to override this function to perform the mapping. */
  virtual T map(T x) const = 0;

  // Maybe implement () operator (invokes map)
};




template<class T>
class rsMapperLinToExp : public rsMapper<T>
{

public:


  rsMapperLinToExp(T inMin = 0, T inMax = 1, T outMin = 1, T outMax = 2)
  {
    setRanges(inMin, inMax, outMin, outMax);
  }


  void setRanges(T inMin, T inMax, T outMin, T outMax)
  {
    this->inMin  = inMin;
    this->inMax  = inMax;
    this->outMin = outMin;
    this->outMax = outMax;

    //inScale  = T(1) / (inMax-inMin);
    argScale = log(outMax / outMin) / (inMax-inMin);


    // ToDo:  precompute 1/(inMax-inMin) and log(outMax/outMin)
  }


  T map(T x) const override
  {
    T tmp = (x - inMin);
    return outMin * std::exp(tmp * argScale);


    //T tmp = (x - inMin) / (inMax - inMin);
    //return outMin * std::exp(tmp * (log(outMax / outMin)));
  };

public:

  T inMin    = 0;
  T inMax    = 1;
  T outMin   = 1;
  T outMax   = 2;

  //T inScale  = 1;
  T argScale = 1;


  // Maybe get rid of members that arent used in the computations. Only keep inMin, exponentScale,
  // resultScale


  // See rsLinToExp(T in, T inMin, T inMax, T outMin, T outMax)


  // Maybe use log2 and exp2 instead of log and exp. That might be a bit faster. See:
  // https://cs.stackexchange.com/questions/27832/is-2x-faster-to-compute-than-expx
  // https://stackoverflow.com/questions/30222836/should-exp2-be-faster-than-exp

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
  inline bool isLogScaled() const { return logScaled; }

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
  inline bool isLogScaledX() const { return mapperX.isLogScaled(); }

  inline T getInMinY()  const { return mapperY.getInMin();  }
  inline T getInMaxY()  const { return mapperY.getInMax();  }
  inline T getOutMinY() const { return mapperY.getOutMin(); }
  inline T getOutMaxY() const { return mapperY.getOutMax(); }
  inline bool isLogScaledY() const { return mapperY.isLogScaled(); }


//protected:

  rsCoordinateMapper<T> mapperX, mapperY;

};

#endif
