#ifndef rosic_VectorSamplePlayer_h
#define rosic_VectorSamplePlayer_h

// rosic-indcludes:
#include "../generators/rosic_SamplePlayer.h"
#include "../others/rosic_VectorMixer.h"
#include "../modulators/rosic_LowFrequencyOscillator.h"


namespace rosic
{

  /**

  This class implements a vector synthesis source comprising four SamplePlayer objects. It is used 
  as oscillator section of the MagicCarpet synth. 

  */

  class VectorSamplePlayer 
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    VectorSamplePlayer();

    /** Destructor. */
    virtual ~VectorSamplePlayer();

    //---------------------------------------------------------------------------------------------
    // delegations:



    /** Sets up the sample-rate for the embedded sample-players and LFOs. */
    void setSampleRate(double newSampleRate);

    /** Sets up the tempo in BPM for the embedded LFOs. */
    void setBeatsPerMinute(double newBpm);

    /** Sets the current key and velocity, so the players can set up their key and velocity
    dependent parameters. */
    void setKeyAndVel(int newKey, int newVel);


    void setPlaybackFrequencyNominal(double newFrequency);

    /** Resets the embedded modules. */
    void reset();

    // maybe introduce sync setter here....

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *outL, double *outR);

    //---------------------------------------------------------------------------------------------
    // embedded objects:

    SamplePlayer samplePlayerTopLeft, samplePlayerTopRight, samplePlayerBottomLeft, 
      samplePlayerBottomRight;

    VectorMixer            vectorMixer;
    LowFrequencyOscillator xLfo, yLfo;
      // actually, these objects are not necessary per voice, but to cater for the GUI framework,
      // we have to put them here

    double topLeftGainFactor, topRightGainFactor, bottomLeftGainFactor, bottomRightGainFactor;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void VectorSamplePlayer::getSampleFrameStereo(double *outL, double *outR)
  {
    double accuL = 0.0;
    double accuR = 0.0;
    double tmpL  = 0.0;
    double tmpR  = 0.0;

    // we suppose here that some outlying class has calculated and set up the gain factors for the 
    // four sample-players, so we don't do the calculation ourselves

    samplePlayerTopLeft.getSampleFrameStereo(&tmpL, &tmpR); 
    accuL += topLeftGainFactor * tmpL;
    accuR += topLeftGainFactor * tmpR;
    samplePlayerTopRight.getSampleFrameStereo(&tmpL, &tmpR); 
    accuL += topRightGainFactor * tmpL;
    accuR += topRightGainFactor * tmpR;
    samplePlayerBottomLeft.getSampleFrameStereo(&tmpL, &tmpR); 
    accuL += bottomLeftGainFactor * tmpL;
    accuR += bottomLeftGainFactor * tmpR;
    samplePlayerBottomRight.getSampleFrameStereo(&tmpL, &tmpR); 
    accuL += bottomRightGainFactor * tmpL;
    accuR += bottomRightGainFactor * tmpR;

    *outL = accuL;
    *outR = accuR;
  }

} // end namespace rosic

#endif //  rosic_VectorSamplePlayer_h