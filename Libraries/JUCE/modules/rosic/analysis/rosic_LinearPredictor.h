#ifndef rosic_LinearPredictor_h
#define rosic_LinearPredictor_h

#include <string.h> // for memmove
//#include <string> // for memmove

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a direct-form ('transversal') FIR (all-zero) predictor which can be used
  to estimate all-pole models or to whiten signals which have strong formants.

  */

  class LinearPredictor
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. The newMaximumOrder argument determines the maximum prediction order. */
    LinearPredictor(int newMaximumOrder = 128);

    /** Destructor.  */
    ~LinearPredictor();


    /** \name Setup */

    /** Set the sample rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the order of the prediction filter. */
    void setOrder(int newOrder);

    /** Sets the learning rate for the adaption of the weight vector. */
    void setLearnRate(double newLearnRate);

    /** Sets the forgetting rate for the adaption of the weight vector. */
    void setForgetRate(double newForgetRate);

    /** Sets the momentum term for the adaption of the weight vector - this can be seen as a
    smoother for the trajectory in weight space. */
    void setMomentum(double newMomentum);


    /** \name Inquiry */

    /** Returns a pointer to the filter coefficient vector, which is needed for implementation
    of the reconstruction(IIR)-filter. */
    double* getWeightVector() const { return weightVector; }

    /** Return the order of the prediction filter. */
    int getOrder() const { return order; }


    /** \name Audio Processing */

    /** Returns one prediction error sample and internally updates the weight vector. */
    INLINE double getSample(double in);


    /** \name Miscellaneous */

    /** Resets the weight vector, vector of past inpu samples and updateVector to all zeros. */
    void reset();

  protected:

    /** \name Data */

    int order, maxOrder;
      // The order of the prediction filter and its maximum value.

    double learnRate, forgetRate, forgetFactor, momentum;
      // The LMS adaption parameters.

    double* weightVector;
      // The weight-vector for weigthing the past inputs.

    double* pastInputs;
      // The input-vector which contains the past inputs.

    double* updateVector;
      // the vector which is added to the coefficient vector.

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double LinearPredictor::getSample(double in)
  {
    double predictedSamp, errorSamp, normalizer;
    int    i;

    // try to predict the incoming sample by means of a scalar product of the weight vector and the
    // stored past input-samples:
    predictedSamp = 0.0;
    for(i=0; i<order; i++)
      predictedSamp += weightVector[i] * pastInputs[i];

    // calculate the prediction error as the difference between the predicted and the true sample:
    errorSamp = in - predictedSamp;

    // some stability checking for debug-mode:
#ifdef _MSC_VER
    if( _isnan(errorSamp) || !_finite(errorSamp) )
      DEBUG_BREAK;
#endif

    // calculate the normalizer:
    normalizer = 0.0;
    for(i=0; i<order; i++)
      normalizer += pastInputs[i] * pastInputs[i];
    normalizer = 1.0 / (sqrt(normalizer)+0.01);

    // calculate the update vector:
    for(i=0; i<order; i++)
      updateVector[i] = normalizer*learnRate*errorSamp*pastInputs[i] + momentum*updateVector[i];

    // update the weight-vector by adding the update-vector:
    for(i=0; i<order; i++)
      weightVector[i] = forgetFactor*weightVector[i] + updateVector[i];

    // shift all the values in the vector which contains the past inputs down by one position:
    memmove(&pastInputs[1], &pastInputs[0], (order-1)*sizeof(double));
      // this is equivalent to: for(i=(order-1); i>0; i--)  pastInputs[i] = pastInputs[i-1];

    // store the current input sample in the first position of the pastInputs-vector for the next
    // iteration:
    pastInputs[0] = in;

    return errorSamp;
  }

}

#endif


