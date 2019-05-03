#ifndef rosic_MultiBandCompressor_h
#define rosic_MultiBandCompressor_h

namespace rosic
{


/** Baseclass for multiband effects such as multiband-compressors, -distortion. etc. */

class rsMultiBandEffect
{

public:

  /** Constructor. */
  rsMultiBandEffect();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  void setSampleRate(double newSampleRate);

  /** Adds a new band after the band at the given index. If the index is that of the last band, 
  the new band will become the last band, otherwise, it will be inserted between index and index+1.
  In the latter case, it will probably make sense, if the new frequency is somewhere between those 
  of the bands at index and index+1 - however, this class actually doesn't care about that. */
  void insertBand(int index, double splitFrequency); 

  /** Removes the band with the given index. The frequency range that was occupied by this band 
  will then either be covered by the former left or right neighbour band, depending on the second
  boolean parameter. When the first ot last band is removed, the second or second-to-last will 
  become the new first or last band. When there is only a single band, this function will have no 
  effect (you can't remove all bands - there's always at least one). */
  void removeBand(int index, bool mergeWithRightNeighbour = false);

  /** Initializes the splitter and array for splitter output. 
  
  old: After init, the MultiBandEffect will have one single band. */
  void initBands(int numberOfBands);


  void setSplitMode(int newMode);

  void setSplitFrequency(int bandIndex, double newFrequency);


  //---------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Retruns the number of bands that are currently active. */
  int getNumberOfBands() const { return splitter.getNumActiveBands(); }
  //int getNumberOfBands() const { return numBands; }

  /** Returns the upper cutoff frequency for the band with given index. */
  double getSplitFrequency(int index) const { return splitter.getSplitFrequency(index); }
    // maybe the lower (highpass) frequency would be better?

  /** Returns the magnitude response in decibels for the band with given index at the given 
  frequency. */
  double getDecibelsAt(int index, double frequency) const;

  //---------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Splits an incoming stereo signal into "numBands" bands and fills our member arrays 
  tmpL, tmpR with the individual band outputs. Subclasses are supposed to call this before 
  processing each band. */
  INLINE void split(double *inL, double *inR)
  {
    rsFloat64x2 in(*inL, *inR);
    if(splitter.getNumActiveBands() == 0)
      tmp[0] = in;
    else
      splitter.processSampleFrame(in, &tmp[0]);
  }

  /** Recombines the bands from our member arrays tmpL, tmpR into the outputs. Subclasses are 
  supposed to call this after processing each indiviudal band to produce the final recombined
  output. */
  INLINE void recombine(double *outL, double *outR)
  {
    if(getNumberOfBands() == 0) return; // bypassed, also needed to avoid access violation
    rsFloat64x2 out(0, 0);
    for(int k = 0; k < getNumberOfBands(); k++)
      out += tmp[k];
    *outL = out[0];
    *outR = out[1];
  }

  /** Returns pointer to the output of k-th band for left channel (supposed to be used by subclass
  to access individual band signals after splitting). */
  inline double* getLeft( int bandIndex) { return &tmp[bandIndex][0]; }

  /** Like getLeft, but for right channel. */
  inline double* getRight(int bandIndex) { return &tmp[bandIndex][1]; }

  /** Resets the states of the band-splitting filters. */
  void reset();

protected:


  RAPT::rsMultiBandSplitter<rsFloat64x2, double> splitter;
  std::vector<rsFloat64x2> tmp; // buffer for bandsplitter outputs
};

//=================================================================================================

/** A multiband compressor with an arbitrary number of bands and perfect reconstruction, i.e. with 
neutral settings, it will be totally transparent (up to roundoff error) and not color the signal in
any way. */

class rsMultiBandCompressor : public rsMultiBandEffect
{

public:

  //---------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsMultiBandCompressor();

  /** Destructor */
  ~rsMultiBandCompressor();

  //---------------------------------------------------------------------------------------------
  /** \name Setup */

  void setSampleRate(double newSampleRate);

  void setThreshold(int bandIndex, double newThreshold);

  void setRatio(int bandIndex, double newRatio);

  void setAttackTime(int bandIndex, double newAttackTime);

  void setReleaseTime(int bandIndex, double newReleaseTime);

  void insertBand(int index, double splitFrequency); 

  void removeBand(int index, bool mergeWithRightNeighbour = false);

  void initBands();

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Calculates one output stereo sample-frame at a time. */
  INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

  void reset();

  //===============================================================================================

protected:

  /** Clears our "compressors" array. */
  void clearCompressors();

  std::vector<Compressor*> compressors;
};

//-----------------------------------------------------------------------------------------------
// inlined functions:

INLINE void rsMultiBandCompressor::getSampleFrameStereo(double *inOutL, double *inOutR)
{
  split(inOutL, inOutR);
  for(int k = 0; k < getNumberOfBands(); k++)  // compress individual bands
    compressors[k]->getSampleFrameStereo(getLeft(k), getRight(k));
  recombine(inOutL, inOutR);
}

}

#endif 
