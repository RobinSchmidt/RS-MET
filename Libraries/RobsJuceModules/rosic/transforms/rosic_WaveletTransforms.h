#ifndef rosic_WaveletTransforms_h
#define rosic_WaveletTransforms_h

//// rosic includes:
//#include "../basics/rosic_HelperFunctions.h"
//#include <algorithm>

namespace rosic
{

  /** Performs a 1D Haar-transform on the passed buffer. The number of levels should be at most
  the number of times, which the length is divisible by 2 (without remainder). */
  template <class T>
  void trafoHaar1D(T *buffer, int length, int numLevels);

  /** Performs a 1D inverse Haar-transform on the passed buffer. The number of levels should be at
  most the number of times, which the length is divisible by 2 (without remainder). */
  template <class T>
  void trafoInverseHaar1D(T *buffer, int length, int numLevels);

  //===============================================================================================
  // implementation:

  template <class T>
  void trafoHaar1D(T *buffer, int length, int numLevels)
  {
    const double c = SQRT2_INV;
    //T *tmp = (T*) alloca(length*sizeof(T));
    T* tmp = new T[length]; // todo: get ird - let client pass workspace pointer
    int w  = length;
    for(int level=1; level<=numLevels; level++)
    {
      int j = w/2;
      for(int i=0; i<j; i++)
      {
        tmp[i]   = c * (buffer[2*i] + buffer[2*i+1]);
        tmp[i+j] = c * (buffer[2*i] - buffer[2*i+1]);
      }
      memcpy(buffer, tmp, w*sizeof(T));
      w /= 2;
    }
    delete[] tmp;
  }

  template <class T>
  void trafoInverseHaar1D(T *buffer, int length, int numLevels)
  {
    const double c = SQRT2_INV;
    //T *tmp = (T*) alloca(length*sizeof(T));
    T* tmp = new T[length];  // get rid
    int w  = length/numLevels;
    for(int level=1; level<=numLevels; level++)
    {
      int j = w/2;
      for(int i=0; i<j; i++)
      {
        tmp[2*i]   = c * (buffer[i] + buffer[i+j]);
        tmp[2*i+1] = c * (buffer[i] - buffer[i+j]);
      }
      memcpy(buffer, tmp, w*sizeof(T));
      w *= 2;
    }
    delete[] tmp;
  }

} // end namespace rosic

#endif // #ifndef rosic_WaveletTransforms_h
