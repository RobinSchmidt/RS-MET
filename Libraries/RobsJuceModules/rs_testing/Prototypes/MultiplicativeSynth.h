#pragma once


/** A prototype for a synthesis algorithm based on products of sinusoidal oscillators which I call
"multiplicative synthesis". The API and internal representation is made for convenient 
experementation in a testbed but not necessarily suitable for a realtime engine (std::vectors are 
copied and stuff). */

template<class T>
class rsMultiplicativeSynth
{

public:


  void setSampleRate(T newRate) { sampleRate = newRate; }
  void setBaseFrequency(T newFreq) { baseFreq = newFreq; }

  void setGain(T newGain) { gain = newGain; }

  void setOperatorFreqFactors(const std::vector<T>& newFactors) { opFreqFactors = newFactors; }

  void setCombinatorWeightsA(const std::vector<T>& newWeights) { cmWeightsA = newWeights; }
  void setCombinatorWeightsB(const std::vector<T>& newWeights) { cmWeightsB = newWeights; }
  void setCombinatorWeightsP(const std::vector<T>& newWeights) { cmWeightsP = newWeights; }


  int getNumPartials();


  std::vector<T> renderOutput(int numSamples);

protected:

  T getSumOfSquaresOfWeights();

  T baseFreq   = 440;
  T sampleRate = 44100;
  T gain       = 1.0;

  std::vector<T> opFreqFactors;  /**< Operator frequency multipliers.       */
  std::vector<T> cmWeightsA;     /**< Combinator weights for input A.       */
  std::vector<T> cmWeightsB;     /**< Combinator weights for input B.       */
  std::vector<T> cmWeightsP;     /**< Combinator weights for product P=A*B. */

};
