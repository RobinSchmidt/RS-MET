#ifndef rosic_VectorMixer_h
#define rosic_VectorMixer_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This is a vector mixer for four sources. You can set up the x- and y-coordiates (between -1...+1)
  and then retrieve the weighting factors for the four sources.

  */

  class VectorMixer
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    VectorMixer(); 

    /** Destructor. */
    ~VectorMixer();  

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the x- coordinate. Value is expected to be in the range -1...1. */
    void setX(double newX);

    /** Sets the y- coordinate. Value is expected to be in the range -1...1. */
    void setY(double newY);

    /** Sets the x- and y-coordinates. Values are expected to be in the range -1...1. */
    INLINE void setXY(double newX, double newY);

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the current setting of the x-coordinate, value is between -1...+1. */
    double getX() { return x; }

    /** Returns the current setting of the y-coordinate, value is between -1...+1. */
    double getY() { return y; }

    /** Assigns the weighting factors for the four signals ccording to the x-/y setting. */
    INLINE void getWeights(double *topLeft, double *topRight, double *bottomLeft, 
      double *bottomRight);

    //=============================================================================================

  protected:

    doubleA x, y;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void VectorMixer::setXY(double newX, double newY)
  {
    x = newX;
    y = newY;
  }

  INLINE void VectorMixer::getWeights(double *topLeft, double *topRight, double *bottomLeft, 
    double *bottomRight)
  {
    // calculate left/right weights:
    double tmp = 0.25 * PI * (x+1.0);
    double wl, wr;
    sinCosApprox(tmp, &wr, &wl);  // wr=sin(tmp), wl=cos(tmp)

    // calculate top/bottom weights:
    tmp = 0.25 * PI * (y+1.0);
    double wt, wb;
    sinCosApprox(tmp, &wt, &wb);  // wt=sin(tmp), wb=cos(tmp)

    // the overall 'edge' weights result from multiplication of the respective left/right, 
    // top/bottom weights:
    *topLeft     = wt*wl;
    *topRight    = wt*wr;
    *bottomLeft  = wb*wl;
    *bottomRight = wb*wr;
  }

} // end namespace rosic

#endif // rosic_VectorMixer_h
