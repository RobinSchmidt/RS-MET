#ifndef rosic_FormantRemover_h
#define rosic_FormantRemover_h

//// rosic-indcludes:
//#include "rosic_LinearPredictor.h"

namespace rosic
{

  /**

  This class removes formants from an incoming audio signal by first applying an adaptive first 
  order pre-emphasis filter, the applying an adaptive prediction error filter (as inherited from 
  LinearPredictor) and finally undoing the pre-emphasis by its inverse filter.

  */

  class FormantRemover : public LinearPredictor
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    FormantRemover(int newMaxOrder = 128);


    /** \name Audio Processing */

    /** Returns one prediction error sample and internally updates the weight vector and 
    pre-emphasis coefficient. */
    INLINE double getSample(double in); 


    /** \name Miscellaneous */

    /** Resets the state-variables and coefficient of the pre-emphasis filter and calls the 
    baseclass' method (which resets the weight vector to all zeros). */
    void reset();

  protected:

    /** \name Data */

    double coeff;
    double pastIn, pastOut;

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double FormantRemover::getSample(double in)
  {
    double predicted, error, tmp;

    // predict the input signal and establish the prediction error (which serves as pre-emphasized 
    // signal):
    predicted = coeff*pastIn;
    error     = in-predicted;

    // apply the linear prediction error filter to the pre-emphsized signal:
    tmp = LinearPredictor::getSample(error);

    // undo the pre-empasis:
    tmp     = tmp + coeff*pastOut;
    pastOut = tmp;

    // adapt the pre-emphasis coefficient:
    coeff = forgetFactor*coeff + learnRate*error*pastIn;

    // update state:
    pastIn = in;

    return tmp;
  }

}

#endif
