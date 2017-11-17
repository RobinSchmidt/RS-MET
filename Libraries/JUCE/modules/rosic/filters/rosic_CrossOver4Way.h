#ifndef rosic_CrossOver4Way_h
#define rosic_CrossOver4Way_h

namespace rosic
{

/** This class implements a crossover filter to split the signal into several bands (at most 4). */

class rsCrossOver4Way
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  rsCrossOver4Way();

  /** Destructor. */
  ~rsCrossOver4Way();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets up the sample rate. */
  void setSampleRate(double newSampleRate);

  /** Activates or de-activates one of the bands. */
  void setBandActive(bool shouldBeActive, int treeLevel, int indexInLevel);

  /** Sets the crossover frequency for one of the crossovers. The stage parameter refers to the 
  level in the tree (0,1,2) and the bandIndex refers to the index within the array of the 
  crossovers in the respective stage. Don't confuse this the output channel number (it is true that 
  the index in the final stage corresponds to the same output-channel-index, but that is not true 
  for all stages). */
  void setCrossoverFrequency(double newCrossoverFrequency, int treeLevel, int indexInLevel);

  /** Sets the slope for one of the crossovers. */
  void setSlope(int newSlope, int treeLevel, int indexInLevel);

  /** Switches between stereo and mono mode. */
  void setMonoMode(bool shouldBeMono);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Informs whether the band at the given treeLevel with the given index is active or not. */
  bool isBandActive(int treeLevel, int indexInLevel) const;

  /** Returns the crossover frequency of the filter-pair at the given treeLevel with the given 
  index. */
  double getCrossoverFrequency(int treeLevel, int indexInLevel) const;

  /** Returns the slope of the filter-pair at the given treeLevel with the given index. */
  int getSlope(int treeLevel, int indexInLevel) const;

  /** Fills the 'magnitudes' array with the magnitude response of one of the filters at the 
  frequencies passed in the 'frequencies' array. Both arrays are assumed to be numBins long. */
  void getMagnitudeResponse(double* frequencies, double* magnitudes, int numBins, 
    int outputChannel, bool inDecibels);

  //-----------------------------------------------------------------------------------------------
  // audio processing:

  /** Calculates an output sample frame. The inOut pointer should be a pointer to an array of 16 
  doubles representing the 8 output-bands for left and right in the following order: 
  L1, R2, L2, R2, L3, R3 and so on. The broadband input signal is expected to be passed in 
  L1, R1. If the crossover is set up to have less than 4 bands, not all output channels may be used 
  (they will be filled with zeros then). */
  INLINE void processSampleFrame(double *inOut);

  /** Processes one buffer of input/output samples. This is a convenience function to make it 
  easier to wrap this class into a plugin. The first index is for the channel (and is assumed to 
  range from 0...7) and the second index is the sample-frame number (assumed to range from 
  0...length-1). */
  void processBuffer(float **inOutBuffer, int length);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Resets the internal buffers of the filters to all zeros. */
  void resetBuffers();

  //===============================================================================================

protected:

  /** Sets up the allpass filters that compensate for the phase response of the other branch. */
  void setupCompensationAllpasses();

  // we use the crossover-object in a binary tree like structure:
  rsLinkwitzRileyCrossOverStereo stage1;
  rsLinkwitzRileyCrossOverStereo stage2[2];

  // allpass-filters for compensating for the allpasses that results from addition of sub-branches:
  rsBiquadCascadeStereo lowBranchCompensationAllpass, highBranchCompensationAllpass;

};

//-------------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsCrossOver4Way::processSampleFrame(double *inOut)
{
  double inL = inOut[0];
  double inR = inOut[1];

  for(int i=0; i<8; i++)
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
}

}

#endif
