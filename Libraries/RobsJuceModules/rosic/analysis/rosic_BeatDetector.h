#ifndef rosic_BeatDetector_h
#define rosic_BeatDetector_h

//#include "rosic_OnsetDetector.h"
//#include <limits.h>

namespace rosic
{

  /**

  This class ...needs documentation

  */

  class BeatDetector : public OnsetDetector
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Standard constrcutor. */
    BeatDetector(); 

    /** Standard constrcutor. */
    ~BeatDetector(); 

    //---------------------------------------------------------------------------------------------
    // setup:



    //---------------------------------------------------------------------------------------------
    // inquiry:




    //---------------------------------------------------------------------------------------------
    // processing:

    /** When the onset detection is done, this function should be invoked to process the onsets so 
    as to obtain estimates for the instantaneous tempo at each onset, a probability whether each 
    onset is a beat and finally make a decision whether or not the onset is considered as a 
    beat. Thus, this function assigns the members 'bpm', 'beatProbability', 'isOnset' members in
    each onset of the inherited vector of onsets. */
    void processOnsets();





    //=============================================================================================

  protected:

    /** Estimates the (possibly time-varying) tempo for each onset in the inherited vector of 
    onsets and assigns the bpm member in the respective onset accordingly. */
    void estimateTempi();

    /** Computes for each onset in the inherited vector of onsets a probability that this onset is
    a beat and assigns the beatProbability member in the respective onset accordingly - requires 
    that tempo-estimates already have been computed. */
    void computeBeatProbabilities();

    /** For each onset in the inherited vector of onsets, this function makes the decision whether
    or not the onset is considered as a beat and assigns the isBeat member in the respective onset 
    accordingly - requires that beat-probabilities already have been computed. */
    void markBeats();

    /** Wraps a value into the octave starting at xMin (inclusive) and ending at 2*xMin (exclusive) 
    by multiplying or dividing by the appropriate power of 2. */
    float wrapIntoOctave(float x, float xMin);

  };


} // end namespace rosic

#endif // rosic_BeatDetector_h
