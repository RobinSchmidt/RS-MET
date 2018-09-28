#ifndef RAPT_INTERPOLATINGFUNCTION_H_INCLUDED
#define RAPT_INTERPOLATINGFUNCTION_H_INCLUDED

//=================================================================================================

/** A class for interpolating between data points using various methods such as linear, cubic. 
etc. */

template<class Tx, class Ty>
class rsInterpolatingFunction
{

public:

  //rsInterpolatingFunction();

  /** Enumeration of possible interpolation modes. */
  enum  
  {
    //LEFT_NEIGHBOUR,
    //RIGHT_NEIGHBOUR,
    //NEAREST_NEIGHBOUR,
    LINEAR,             // this is the default
    CUBIC
    //QUARTIC,
    //QUINTIC,
  };

  //-----------------------------------------------------------------------------------------------
  // \name Setup:

  void setData(Tx *x, Ty *y, int numValues)
  {
    this->x = x;
    this->y = y;
    this->numValues = numValues;
    updateCoeffs(); // maybe don't directly call it here, just set a "dirty" flag and update them
                    // only when needed
  }

  void setMode(int newMode)
  {
    mode = newMode;
    updateCoeffs();
  }







  //void setPreMap(const std::function<Ty(Ty)>& newPreMap)
  void setPreMap(std::function<Ty(Ty)> newPreMap)
  {
    preMap = newPreMap;
  }

  void setPostMap(const std::function<Ty(Ty)>& newPostMap)
  {
    postMap = newPostMap;
  }

  void setPreMap(Ty(*f) (Ty x))
  {
    std::function<Ty(Ty)> tmp = f;
    setPreMap(tmp);
  }

  void setPostMap(Ty(*f) (Ty x))
  {
    std::function<Ty(Ty)> tmp = f;
    setPostMap(tmp);
  }
  // ok - this works - but needs clean up



  // void setPrePostMap
  // should allow to apply a mapping function to the y values before interpolation and another one
  // after (the post-mapping should tpyically be the inverse function of the forward mapping - for
  // example, one could use log for the pre-mapping and exp for the post mapping



  //-----------------------------------------------------------------------------------------------
  // \name Evaluation:

  void interpolate(Tx *x, Ty *y, int N, Tx *xi, Ty *yi, int Ni);

  /*
  Ty getValue(Tx x) const
  {
    return Ty(0); // preliminary
  }

  Ty operator()(Tx x) const
  {
    return getValue(x);
  }
  */


protected:

  void updateCoeffs();

  Tx *x = nullptr;
  Ty *y = nullptr;
  Ty **coeffs = nullptr;  // 1st index: datapoint, 2nd index polynomial power
  int numValues = 0;
  int mode = LINEAR;

  //std::function<Ty(Ty)> preMap  = &rsIdentity<Ty>;
  //std::function<Ty(Ty)> postMap = &rsIdentity<Ty>;
  std::function<Ty(Ty)> preMap  = nullptr;
  std::function<Ty(Ty)> postMap = nullptr;

};




// todo: drag class rsResampler over here, make it a subclass of rsInterpolatingFunction


#endif