#ifndef RAPT_INTERPOLATINGFUNCTION_H_INCLUDED
#define RAPT_INTERPOLATINGFUNCTION_H_INCLUDED

//=================================================================================================

/** A class for interpolating between data points using various methods such as linear, cubic,
etc. You can also apply mapping function before and after the interpolation - for example, you may
want to apply the interpolation in a logarithmic domain but the final result should be linearly 
scaled again. You can do this by using log and exp as pre and post mapping functions respectively.

...the implementation is still rather incomplete */

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
    CUBIC,
    NATURAL_CUBIC
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



  /** Sets the mapping function that is applied to the y-values before interpolation. */
  void setPreMap(const std::function<Ty(Ty)>& newPreMap)   { preMap = newPreMap; }

  /** Sets the mapping function that is applied to the y-values after interpolation. */
  void setPostMap(const std::function<Ty(Ty)>& newPostMap) { postMap = newPostMap; }

  /** Allows the pre-mapping function to be set via a function pointer. Pass a nullptr to reset
  the mapping function to the identity function. */
  void setPreMap(Ty(*f) (Ty x))  { setPreMap(std::function<Ty(Ty)>(f)); }

  /** Allows the post-mapping function to be set via a function pointer. Pass a nullptr to reset
  the mapping function to the identity function. */
  void setPostMap(Ty(*f) (Ty x)) { setPostMap(std::function<Ty(Ty)>(f)); }


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

  std::function<Ty(Ty)> preMap  = nullptr;
  std::function<Ty(Ty)> postMap = nullptr;

};




// todo: drag class rsResampler over here, make it a subclass of rsInterpolatingFunction


#endif