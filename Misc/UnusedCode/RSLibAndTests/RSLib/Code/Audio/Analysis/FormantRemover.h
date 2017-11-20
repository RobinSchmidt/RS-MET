#ifndef RS_FORMANTREMOVER_H
#define RS_FORMANTREMOVER_H

namespace RSLib
{

  /**

  This class removes formants from an incoming audio signal by first applying an adaptive first 
  order pre-emphasis filter, the applying an adaptive prediction error filter (as inherited from 
  LinearPredictor) and finally undoing the pre-emphasis by its inverse filter.

  */

  class RSLib_API rsFormantRemover : public rsLinearPredictor
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsFormantRemover(int newMaxOrder = 128);


    /** \name Audio Processing */

    /** Returns one prediction error sample and internally updates the weight vector and 
    pre-emphasis coefficient. */
    RS_INLINE double getSample(double in); 


    /** \name Misc */

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

  RS_INLINE double rsFormantRemover::getSample(double in)
  {
    double predicted, error, tmp;

    // predict the input signal and establish the prediction error (which serves as pre-emphasized 
    // signal):
    predicted = coeff*pastIn;
    error     = in-predicted;

    // apply the linear prediction error filter to the pre-emphsized signal:
    tmp = rsLinearPredictor::getSample(error);

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
