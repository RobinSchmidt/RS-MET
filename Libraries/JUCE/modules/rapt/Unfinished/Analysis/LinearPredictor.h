#ifndef RAPT_LINEARPREDICTOR_H
#define RAPT_LINEARPREDICTOR_H

/** This class implements a direct-form ('transversal') FIR (all-zero) predictor which can be used
to estimate all-pole models or to whiten signals which have strong formants. 

\todo: provid a gradient-adaptive-lattice (GAL) implementation (has easier-to-check stability and 
faster convergence). this one here is then the direct-form implementation...maybe we could also 
have a biquad-chain predictor?...or actually, it would have to be a two-zero chain..but maybe with
ARMA modeling, an actual biquad-chain predictor could be done... */

template<class TSig, class TPar>
class rsLinearPredictor
{

public:

  /** \name Construction/Destruction */

  /** Constructor. The newMaximumOrder argument determines the maximum prediction order. */
  rsLinearPredictor(int newMaximumOrder = 128);

  /** Destructor.  */
  ~rsLinearPredictor();


  /** \name Setup */

  /** Set the sample rate. */
  //void setSampleRate(TPar newSampleRate);

  /** Sets the order of the prediction filter. */
  void setOrder(int newOrder);

  /** Sets the learning rate for the adaption of the weight vector. */
  void setLearnRate(TPar newLearnRate);

  /** Sets the forgetting rate for the adaption of the weight vector. */
  void setForgetRate(TPar newForgetRate);

  /** Sets the momentum term for the adaption of the weight vector - this can be seen as a
  smoother for the trajectory in weight space. */
  void setMomentum(TPar newMomentum);


  /** \name Inquiry */

  /** Returns a pointer to the filter coefficient vector, which is needed for implementation
  of the reconstruction(IIR)-filter. */
  TSig* getWeightVector() const { return weightVector; }

  /** Return the order of the prediction filter. */
  int getOrder() const { return order; }


  /** \name Audio Processing */

  /** Returns one prediction error sample and internally updates the weight vector. */
  RS_INLINE TSig getSample(TSig in);


  /** \name Misc */

  /** Resets the weight vector, vector of past inpu samples and updateVector to all zeros. */
  void reset();

protected:

  /** \name Data */

  int order, maxOrder;
    // The order of the prediction filter and its maximum value.

  TPar learnRate, forgetRate, forgetFactor, momentum;
    // The LMS adaption parameters.

  TSig* weightVector;
    // The weight-vector for weigthing the past inputs. 

  TSig* pastInputs;
    // The input-vector which contains the past inputs.

  TSig* updateVector;
    // the vector which is added to the coefficient vector.

};

//-----------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
RS_INLINE TSig rsLinearPredictor<TSig, TPar>::getSample(TSig in)
{
  TSig predictedSamp, errorSamp, normalizer;
  int  i;

  // try to predict the incoming sample by means of a scalar product of the weight vector and the 
  // stored past input-samples:
  predictedSamp = 0.0;
  for(i = 0; i < order; i++)
    predictedSamp += weightVector[i] * pastInputs[i];

  // calculate the prediction error as the difference between the predicted and the true sample:
  errorSamp = in - predictedSamp;

  // some stability checking for debug-mode:
#ifdef _MSC_VER
  if(_isnan(errorSamp) || !_finite(errorSamp))
    RS_DEBUG_BREAK;
#endif

    // calculate the normalizer:
  normalizer = 0.0;
  for(i = 0; i < order; i++)
    normalizer += pastInputs[i] * pastInputs[i];
  normalizer = TSig(1.0) / TSig(sqrt(normalizer)+0.01);

  // calculate the update vector:
  for(i = 0; i < order; i++)
    updateVector[i] = normalizer*learnRate*errorSamp*pastInputs[i] + momentum*updateVector[i];

  // update the weight-vector by adding the update-vector:
  for(i = 0; i < order; i++)
    weightVector[i] = forgetFactor*weightVector[i] + updateVector[i];

  // shift all the values in the vector which contains the past inputs down by one position:
  memmove(&pastInputs[1], &pastInputs[0], (order-1)*sizeof(TSig));
    // this is equivalent to: for(i=(order-1); i>0; i--)  pastInputs[i] = pastInputs[i-1];

  // store the current input sample in the first position of the pastInputs-vector for the next 
  // iteration:
  pastInputs[0] = in;

  return errorSamp;
}

#endif
