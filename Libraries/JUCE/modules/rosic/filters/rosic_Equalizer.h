#ifndef rosic_Equalizer_h
#define rosic_Equalizer_h

//// includes from the STL:
//#include <vector>
//using std::vector;
//
//// rosic-indcludes:
//#include "rosic_TwoPoleFilter.h"

namespace rosic
{

  /**

  This class implements a series connection of an arbitrary number of parametric equalizer filters.

  \todo: implement stereo-modes: Linked, Left/Right, Mid/Side, Mono ...but maybe we should do this in its own class
  ...but the mono-switch must be handled here nonetheless


  */

  class Equalizer
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Equalizer();

    /** Destructor. */
    ~Equalizer();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Sets up the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Adds a new filter into the vector. The integer return-value informs, at which index the new band was inserted.  */
    int addBand(int newMode, double newFrequency, double newGain, double newBandwidth = 2.0*asinh(1.0/sqrt(2.0))/log(2.0));

    /** Removes a filter from the vector. The return-value informs, if there was actually a band removed (if you try to remove a 
    non-existing band it will return false */
    bool removeBand(unsigned int index);

    /** Removes all bands and returns true if there was at least one band to remove. */
    bool removeAllBands();

    /** Modifies the data of an existing band. If you try to modify a non-existent band or try to modify the data in such a way which is 
    not allowed, it will return false. */
    bool modifyBand(unsigned int index, int newMode, double newFrequency, double newGain, double newBandwidth);

    /** Changes the mode of one of the bands. */
    bool setBandMode(int index, int newMode);

    /** Changes the frequency of one of the bands. */
    bool setBandFrequency(int index, double newFrequency);

    /** Changes the gain of one of the bands. */
    bool setBandGain(int index, double newGain);

    /** Changes the bandwidth of one of the bands. */
    bool setBandBandwidth(int index, double newBadwidth);

    /** Sets the lower bandedge frequency for one of the bands by updating the bandwidth parameter (keeping the center frequency 
    constant). */
    bool setLowerBandedgeFrequency(unsigned int index, double newLowerBandedgeFrequency);

    /** Sets the upper bandedge frequency for one of the bands by updating the bandwidth parameter (keeping the center frequency 
    constant). */
    bool setUpperBandedgeFrequency(unsigned int index, double newUpperBandedgeFrequency);

    /** Changes the global gain in decibels. */
    void setGlobalGain(double newGainInDB) { globalGainFactor = dB2amp(newGainInDB); }

    /** Switches the equalizer inot mono-mode in which case only the 1st channel will be processed and the result will be copied into the 
    2nd channel as well. */
    void setMono(bool shouldBeMono) { mono = shouldBeMono; }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the number of bands. */
    int getNumBands();

    /** Returns the mode of one of the bands (0 if index is out of range). */
    int getBandMode(unsigned int index);

    /** Returns the frequency of one of the bands (0.0 if index is out of range). */
    double getBandFrequency(unsigned int index);

    /** Returns the gain of one of the bands (0.0 if index is out of range). */
    double getBandGain(unsigned int index);

    /** Returns the bandwidth of one of the bands (0.0 if index is out of range). */
    double getBandBandwidth(unsigned int index);

    /** Returns the lower bandedge frequency (defined at the half dB-gain point) of one of the bands (0.0 if index is out of range). */
    double getLowerBandedgeFrequency(unsigned int index);

    /** Returns the upper bandedge frequency (defined at the half dB-gain point) of one of the bands (0.0 if index is out of range). */
    double getUpperBandedgeFrequency(unsigned int index);

    /** Returns true if the current mode of some band supports a bandwidth parameter, false otherwise. */
    bool doesModeSupportBandwidth(unsigned int index);

    /** Returns true if the current mode of some band supports a gain parameter, false otherwise. */
    bool doesModeSupportGain(unsigned int index);

    /** Fills the 'magnitudes' array with the magnitude response of the equalizer evaluated at the frequencies passed in the 'frequencies' 
    array. Both arrays are assumed to be numBins long. */
    void getMagnitudeResponse(double* frequencies, double* magnitudes, int numBins);

    /** Returns the global gain in decibels. */
    double getGlobalGain() const { return RAPT::rsAmpToDb(globalGainFactor); }

    //-------------------------------------------------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates an output sample. */
    INLINE double getSample(double in);

    /** Calculates a stereo output sampleframe. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    /** Calculates an entire block of samples at once. */
    INLINE void processBlock(float *inOutL, float *inOutR, int numSamples);

    //-------------------------------------------------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the filters. */
    void reset();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // embedded sub-objects:

    std::vector<TwoPoleFilter> bands;   // vector of the bands

    //=====================================================================================================================================

  protected:

    double sampleRate;   
    double globalGainFactor;  
    bool   mono;

  };

  //---------------------------------------------------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE double Equalizer::getSample(double in)
  {
    double tmp = in;
    for(unsigned int i=0; i<bands.size(); i++)
      tmp = bands[i].getSample(tmp);
    return globalGainFactor*tmp;
  }

  INLINE void Equalizer::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    for(unsigned int i=0; i<bands.size(); i++)
      bands[i].getSampleFrameStereo(inOutL, inOutR);
    *inOutL *= globalGainFactor;
    *inOutR *= globalGainFactor;
  }

  INLINE void Equalizer::processBlock(float *inOutL, float *inOutR, int numSamples)
  {
    double tmpL, tmpR;
    int n;

    if( mono )
    {
      for(n=0; n<numSamples; n++)
      {
        inOutL[n] = (float) getSample( (double) inOutL[n] );
        inOutR[n] = inOutL[n];
      }
    }
    else
    {
      for(n=0; n<numSamples; n++)
      {
        tmpL = (double) inOutL[n];
        tmpR = (double) inOutR[n];
        getSampleFrameStereo(&tmpL, &tmpR);
        inOutL[n] = (float) tmpL;
        inOutR[n] = (float) tmpR;
      }
    }
  }

} 

#endif
