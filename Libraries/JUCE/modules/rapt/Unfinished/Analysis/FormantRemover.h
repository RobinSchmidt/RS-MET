#ifndef RAPT_FORMANTREMOVER_H
#define RAPT_FORMANTREMOVER_H

/** This class removes formants from an incoming audio signal by first applying an adaptive first
order pre-emphasis filter, the applying an adaptive prediction error filter (as inherited from
LinearPredictor) and finally undoing the pre-emphasis by its inverse filter.

...what about rsFormantPreserver...is this class still available in rosic only? */

template<class TSig, class TPar>
class rsFormantRemover : public rsLinearPredictor<TSig, TPar>
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  rsFormantRemover(int newMaxOrder = 128);


  /** \name Audio Processing */

  /** Returns one prediction error sample and internally updates the weight vector and
  pre-emphasis coefficient. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the state-variables and coefficient of the pre-emphasis filter and calls the
  baseclass' method (which resets the weight vector to all zeros). */
  void reset();

protected:

  /** \name Data */

  TSig coeff;
  TSig pastIn, pastOut;

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsFormantRemover<TSig, TPar>::getSample(TSig in)
{
  TSig predicted, error, tmp;

  // predict the input signal and establish the prediction error (which serves as pre-emphasized
  // signal, the pre-emphasis filter is a 1st order linear prediction error filter):
  predicted = coeff*pastIn;
  error     = in-predicted;

  // apply the linear prediction error filter to the pre-emphsized signal:
  tmp = rsLinearPredictor<TSig, TPar>::getSample(error);

  // undo the pre-empasis:
  tmp     = tmp + coeff*pastOut;
  pastOut = tmp;

  // adapt the pre-emphasis coefficient:
  coeff = this->forgetFactor * coeff + this->learnRate * error * pastIn;

  // update state:
  pastIn = in;

  return tmp;
}

#endif
