#ifndef rosic_FormantPreserver_h
#define rosic_FormantPreserver_h

//// rosic-indcludes:
//#include "rosic_FormantRemover.h"

namespace rosic
{

  /**

  This class removes formants from an incoming audio signal in the same way as the class
  FormantRemover from which it is inherited. Additionally, i can re-apply the removed formants at
  some later stage.

  \todo: use memmove for the shifts of the internal states

  */

  class FormantPreserver : public FormantRemover
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    FormantPreserver(int newMaxOrder = 128);

    /** Destructor.  */
    ~FormantPreserver();

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Removes the formants from an incoming audio (stereo) signal. */
    INLINE void removeFormants(double *inL, double *inR, double *outL, double *outR);

    /** Re-applies the formants which were removed by a call to removeFormants by applying the
    inverse prediction error filter (which is the all-pole model for the formants). */
    INLINE void reApplyFormants(double *inL, double *inR, double *outL, double *outR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Overrides the inherited reset()-method to reset the past outputs as well (which are stored
    in an additional member variable of this subclass). */
    void reset();

  protected:

    double preEmphPastInL, preEmphPastInR, deEmphPastOutL, deEmphPastOutR;
      // States of the pre-emphasis and de-emphasis filter.

    double *pastInputsL, *pastInputsR, *pastOutputsL, *pastOutputsR;
      // States of the predictor and model filter.

  private:

    /** This inherited method has been made private because it should not be called from outside
    anymore - use the methods removeFormants and reApplyFormants instead. */
    INLINE double getSample(double in) { return FormantRemover::getSample(in); }

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void FormantPreserver::removeFormants(double *inL, double *inR, double *outL, double *outR)
  {
    // this call adapts the weight vector and pre-emphasis coeff, the ouput value is not used
    //double dummy = getSample(*inL + *inR);

    // apply pre-emphasis filter to the input signal and update the filter's state:
    double preEmphL = *inL - coeff*preEmphPastInL;
    double preEmphR = *inR - coeff*preEmphPastInR;
    preEmphPastInL  = *inL;
    preEmphPastInR  = *inR;

    // predict signal for left and right channel with the current predictor filter:
    int i;
    double tmpL = 0.0;
    double tmpR = 0.0;
    for(i=0; i<order; i++)
    {
      tmpL += weightVector[i] * pastInputsL[i];
      tmpR += weightVector[i] * pastInputsR[i];
    }

    // shift all the values in the vector which contains the past inputs down by one position,
    // andstore the current input sample in the first position of the pastInputs-vector for the next
    // iteration:
    for(i=(order-1); i>0; i--)
    {
      pastInputsL[i] = pastInputsL[i-1];
      pastInputsR[i] = pastInputsR[i-1];
    }
    pastInputsL[0] = preEmphL;
    pastInputsR[0] = preEmphR;

    // calculate the prediction error and store it in the output slots:
    *outL = preEmphL - tmpL;
    *outR = preEmphR - tmpR;
  }

  INLINE void FormantPreserver::reApplyFormants(double *inL, double *inR, double *outL, double *outR)
  {
    // apply the inverse prediction error filter (the allpole model):
    int i;
    double tmpL = *inL;
    double tmpR = *inR;
    for(i=0; i<order; i++)
    {
      tmpL += weightVector[i] * pastOutputsL[i];
      tmpR += weightVector[i] * pastOutputsR[i];
    }

    // update the state of the allpole model filter:
    for(i=(order-1); i>0; i--)
    {
      pastOutputsL[i] = pastOutputsL[i-1];
      pastOutputsR[i] = pastOutputsR[i-1];
    }
    pastOutputsL[0] = tmpL;
    pastOutputsR[0] = tmpR;

    // apply the inverse pre-emphasis filter (de-emphasis):
    tmpL           = tmpL + coeff*deEmphPastOutL;
    tmpR           = tmpR + coeff*deEmphPastOutR;
    deEmphPastOutL = tmpL;
    deEmphPastOutR = tmpR;

    // store the result in the output slots:
    *outL = tmpL;
    *outR = tmpR;
  }

} // end namespace rosic

#endif // #ifndef rosic_FormantPreserver_h
