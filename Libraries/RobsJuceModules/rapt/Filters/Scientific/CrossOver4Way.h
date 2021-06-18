#ifndef RAPT_CROSSOVER4WAY_H_INCLUDED
#define RAPT_CROSSOVER4WAY_H_INCLUDED

/** This class implements a crossover filter to split the signal into several bands (at most 4). 

ToDo: rename to rsLinkwitzRileyCrossover4Way or rsLinkwitzRileyTree4, rsLinkwitzRileySplitter4 */

template<class TSig, class TPar>
class rsCrossOver4Way
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsCrossOver4Way();

  /** Destructor. */
  ~rsCrossOver4Way() = default;

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets up the sample rate. */
  void setSampleRate(TPar newSampleRate);

  /** Activates or de-activates one of the bands. */
  void setBandActive(bool shouldBeActive, int treeLevel, int indexInLevel);

  /** Sets the crossover frequency for one of the crossovers. The stage parameter refers to the 
  level in the tree (0,1,2) and the bandIndex refers to the index within the array of the 
  crossovers in the respective stage. Don't confuse this the output channel number (it is true that 
  the index in the final stage corresponds to the same output-channel-index, but that is not true 
  for all stages). */
  void setCrossoverFrequency(TPar newCrossoverFrequency, int treeLevel, int indexInLevel);

  /** Sets the slope for one of the crossovers. */
  void setSlope(int newSlope, int treeLevel, int indexInLevel);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Informs whether the band at the given treeLevel with the given index is active or not. */
  bool isBandActive(int treeLevel, int indexInLevel) const;

  /** Returns the crossover frequency of the filter-pair at the given treeLevel with the given 
  index. */
  TPar getCrossoverFrequency(int treeLevel, int indexInLevel) const;

  /** Returns the slope of the filter-pair at the given treeLevel with the given index. */
  int getSlope(int treeLevel, int indexInLevel) const;

  /** Fills the 'magnitudes' array with the magnitude response of one of the filters at the 
  frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getMagnitudeResponse(TPar* frequencies, TPar* magnitudes, int numBins, 
    int outputChannel, bool inDecibels);

  /** Returns a const reference to the first stage. */
  const rsLinkwitzRileyCrossOver<TSig, TPar>& getStage1() { return stage1; }

  /** Returns a const reference to the second stage with index i where i = 0 or 1. */
  const rsLinkwitzRileyCrossOver<TSig, TPar>& getStage2(int i) { return stage2[i]; }

  /** Returns a const reference to the compensation allpass filter for the low branch. This is the 
  allpass filter that simulates the effect of ht splitting and re-summing in the high branch. */
  const rsBiquadCascade<TSig, TPar>& getLowBranchAllpass() const
  { return lowBranchCompensationAllpass; }

  /** Returns a const reference to the compensation allpass filter for the high branch. This is the 
  allpass filter that simulates the effect of ht splitting and re-summing in the low branch. */
  const rsBiquadCascade<TSig, TPar>& getHighBranchAllpass() const
  { return highBranchCompensationAllpass; }

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Calculates an output sample frame. The inOut pointer should be a pointer to an array of 16 
  doubles representing the 8 output-bands for left and right in the following order: 
  L1, R2, L2, R2, L3, R3 and so on. The broadband input signal is expected to be passed in 
  L1, R1. If the crossover is set up to have less than 4 bands, not all output channels may be used 
  (they will be filled with zeros then). */
  inline void processSampleFrame(TSig* inOut); // use non-pointer and return value

  /** Processes one buffer of input/output samples. This is a convenience function to make it 
  easier to wrap this class into a plugin. The first index is for the channel (and is assumed to 
  range from 0...7) and the second index is the sample-frame number (assumed to range from 
  0...length-1). */
  void processBuffer(TSig** inOutBuffer, int length); // use pointer TSig* instead TSig**

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Resets the internal buffers of the filters to all zeros. */
  void resetBuffers();

  //===============================================================================================

protected:

  /** Sets up the allpass filters that compensate for the phase response of the other branch. */
  void setupCompensationAllpasses();

  // We use the crossover-object in a binary tree like structure:
  rsLinkwitzRileyCrossOver<TSig, TPar> stage1;
  rsLinkwitzRileyCrossOver<TSig, TPar> stage2[2];

  // Allpass-filters for compensating for the allpasses that results from addition of sub-branches:
  rsBiquadCascade<TSig, TPar> lowBranchCompensationAllpass, highBranchCompensationAllpass;
  // use shorter names lowBranchAllpass, high...

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TPar>
inline void rsCrossOver4Way<TSig, TPar>::processSampleFrame(TSig* inOut)
{
  // new implementation - needs test:

  TSig in = inOut[0];
  for(int i = 0; i < 4; i++)
    inOut[i] = 0.0;

  if(!stage2[0].isActive() && !stage2[1].isActive()) {      // 2 bands
    stage1.getSamplePair(&in, &inOut[0], &inOut[1]);
  }
  else if(stage2[0].isActive() && !stage2[1].isActive()) {  // 3 bands, lower band split further
    stage1.getSamplePair(&in, &inOut[0],  &inOut[2]);
    stage2[0].getSamplePair(&inOut[0], &inOut[0], &inOut[1]);
    inOut[2] = highBranchCompensationAllpass.getSampleDirect2(inOut[2]);
  }
  else if(!stage2[0].isActive() && stage2[1].isActive()) {  // 3 bands, upper band split further
    stage1.getSamplePair(&in, &inOut[0], &inOut[1]);
    stage2[1].getSamplePair(&inOut[1], &inOut[1], &inOut[2]);
    inOut[0] = lowBranchCompensationAllpass.getSampleDirect2(inOut[0]);
  }
  else  {                                                   // 4 bands
    stage1.getSamplePair(&in, &inOut[0], &inOut[2]);
    inOut[0] = lowBranchCompensationAllpass.getSampleDirect2(inOut[0]);
    inOut[2] = highBranchCompensationAllpass.getSampleDirect2(inOut[2]);
    stage2[0].getSamplePair(&inOut[0], &inOut[0], &inOut[1]);
    stage2[1].getSamplePair(&inOut[2], &inOut[2], &inOut[3]);
  }

  /*
  // needs update - multichannel should now be handled by using a multichannel signal type TSig

  double inL = inOut[0];
  double inR = inOut[1];

  for(int i = 0; i < 8; i++)
    inOut[i] = 0.0;

  if(!stage2[0].isActive() && !stage2[1].isActive())
  {
    // 2 bands:
    stage1.getSampleFrame(&inL, &inR, &inOut[0], &inOut[1], &inOut[2], &inOut[3]);
  }
  else if(stage2[0].isActive() && !stage2[1].isActive())
  {
    // 3 bands, lower band split further:
    stage1.getSampleFrame(&inL, &inR, &inOut[0], &inOut[1], &inOut[4], &inOut[5]);
    stage2[0].getSampleFrame(&inOut[0], &inOut[1], &inOut[0], &inOut[1], &inOut[2], &inOut[3]);
    highBranchCompensationAllpass.getSampleFrameDirect2(&inOut[4], &inOut[5]);
  }
  else if(!stage2[0].isActive() && stage2[1].isActive())
  {
    // 3 bands, upper band split further
    stage1.getSampleFrame(&inL, &inR, &inOut[0], &inOut[1], &inOut[2], &inOut[3]);
    stage2[1].getSampleFrame(&inOut[2], &inOut[3], &inOut[2], &inOut[3], &inOut[4], &inOut[5]);
    lowBranchCompensationAllpass.getSampleFrameDirect2(&inOut[0], &inOut[1]);
  }
  else
  {
    // 4 bands:
    stage1.getSampleFrame(&inL, &inR, &inOut[0], &inOut[1], &inOut[4], &inOut[5]);
    lowBranchCompensationAllpass.getSampleFrameDirect2(&inOut[0], &inOut[1]);
    highBranchCompensationAllpass.getSampleFrameDirect2(&inOut[4], &inOut[5]);
    stage2[0].getSampleFrame(&inOut[0], &inOut[1], &inOut[0], &inOut[1], &inOut[2], &inOut[3]);
    stage2[1].getSampleFrame(&inOut[4], &inOut[5], &inOut[4], &inOut[5], &inOut[6], &inOut[7]);
  }
  */
}

#endif