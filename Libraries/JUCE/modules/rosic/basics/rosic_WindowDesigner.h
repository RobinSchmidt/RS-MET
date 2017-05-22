#ifndef rosic_WindowDesigner_h
#define rosic_WindowDesigner_h

//// rosic-indcludes:
//#include "../math/rosic_SpecialFunctionsReal.h"

namespace rosic
{

  /**

  This class can be used to create window-functions or apply the windows dirctly to buffers.

  */

  class WindowDesigner  
  {

  public:

    /** Enumeration of the available window types. */
    enum windowTypes
    {
      RECTANGULAR = 0,
      // TRIANGULAR // aka Bartlett, Fejer
      BLACKMAN,
      COSINE_SQUARED,
      HAMMING,
      HANN
    };

    //-------------------------------------------------------------------------------------------------------------------------------------
    // definitions of the window functions:

    /** Returns the value of a Blackman window of length N at index n. */
    static INLINE double blackmanWindow(int n, int N)
    {
      return 0.42 - 0.5*cos( 2.0*PI*n / (double) (N-1)) + 0.08*cos(4.0*PI*n / (double) (N-1)) ;
    }

    /** Returns the value of a cosine window raised to a power of length N at index n. When the power is 2 (i.e. we have a cosine squared 
    window), a succession of such windows in time, where the start of one window is one half of the window-length shifted with respect to 
    its neighbours, results in a sum of unity when the windows are added. This is property useful in overlap/add analysis/resynthesis. 
    There are other combinations for the power and the shift that also sum to some constant, for example
    power=2, shift=length/4 -> sum=2, other combinations that work: 4,4; 6,4; 6,6; 8,6 etc. Spectrally, the power parameter controls the 
    tradeoff between mainlobe width and rolloff/attenuation - the higher the power, the wider the mainlobe (due to narrower time-domain 
    extent) and the better the rolloff. */
    static INLINE double cosinePowerWindow(int n, int N, double power)
    {
      return pow( sin(PI * (double) n / (double) N), power );
    }

    /** Returns the value of a Hamming window of length N at index n. */
    static INLINE double hammingWindow(int n, int N)
    {
      return 0.54 - 0.46 * cos(2.0*PI*n / (double) (N-1));
    }

    /** Returns the value of a Hann window of length N at index n. */
    static INLINE double hannWindow(int n, int N)
    {
      return 0.5 * ( 1.0 - cos(2.0*PI*n / (double) (N-1)) );
    }

    /** Returns the value of a Kaiser window of length N at index n with form parameter beta - Kaiser windows are used for designing 
    FIR-filters of type I (odd length with left/right symmetry). */
    static INLINE double kaiserWindow(int n, int N, double beta)
    {
      int M = N-1;
      double alpha = M/2;
      double tmp   = (n-alpha)/alpha;
      return besselI0( beta * sqrt(1-tmp*tmp) ) / besselI0(beta);
    }

    // \todo static double getRequiredKaiserBeta, getRequiredKaiserLength

    //-------------------------------------------------------------------------------------------------------------------------------------
    // creating and appyling the windows:

    /** Writes a window of the desired type into the passed array. @see windowTypes */
    static void getWindow(double *window, int length, int type);

    /** Multiplies the passed buffer with a window of the desired type. @see windowTypes */
    static void applyWindow(double *buffer, int length, int type);

    /** Writes a cosine-power window into the passed array. We need a separate function for this type of window because, unlike the others, 
    it has a parameter (namely, the power)). @see windowTypes */
    static void getCosinePowerWindow(double *window, int length, double power);

    /** Writes a Kaiser window into the passed array. */
    static void getKaiserWindow(double *window, int length, double beta);

  };

} 

#endif 
