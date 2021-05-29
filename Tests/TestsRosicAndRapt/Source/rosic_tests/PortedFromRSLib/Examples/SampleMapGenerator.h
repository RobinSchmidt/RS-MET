#ifndef RS_SAMPLEMAPGENERATOR_H
#define RS_SAMPLEMAPGENERATOR_H

#include <rapt/rapt.h> // get rid

#include "../RSLib/Core/RSCore.h" // why that? obsolete?
using namespace RSLib;


// we have now the AudioStream class in rosic (used by the sampler engine) that provides similar 
// functionality, so maybe get rid of this:
struct rsAudioBuffer  // move to RSLib, maybe templatize on the sample-type (double/float)
{

  rsAudioBuffer()
  {
    numChannels = 0;
    numFrames   = 0;
    data        = NULL;
    dataFlat    = NULL;
  }

  ~rsAudioBuffer()
  {
    freeMemory();
  }

  void setSize(int newNumChannels, int newNumFrames)
  {
    if( newNumChannels != numChannels || newNumFrames != numFrames )
    {
      freeMemory();
      numChannels = newNumChannels;
      numFrames   = newNumFrames;
      allocateMemory();
    }
  }

  /** Returns the total number of samples in this buffer (== the number of channels times the
  number of sample-frames). */
  int getSize()
  {
    return numChannels * numFrames;
  }

  void allocateMemory()
  {
    dataFlat = new double[numChannels*numFrames];
    data     = new double*[numChannels];
    for(int i = 0; i < numChannels; i++)
      data[i] = &dataFlat[i*numFrames];
  }

  void freeMemory()
  {
    delete[] data;
    delete[] dataFlat;
  }

  int numChannels;
  int numFrames;
  double **data;     // data as 2D array, data[i][j] is the j-th frame of the i-th channel
  double *dataFlat;  // data as flat array, used internally

};


template<class T>
struct rsModalBankParameters // maybe move to RAPT, together with rsModalParameterGenerator
{
  //int numModes;    // number of modes

  T frequency;       // reference frequency in Hz - deprecate - shall be determined by the key
  T gain;            // overall gain factor
  T attack;          // attack time in seconds
  T decay;           // decay time in seconds

  std::vector<T> f;  // relative frequencies
  std::vector<T> g;  // gains
  std::vector<T> a;  // relative attack times
  std::vector<T> d;  // relative decay times
  std::vector<T> p;  // start-phases
};

template struct rsModalBankParameters<double>;
typedef rsModalBankParameters<double> rsModalBankParametersD;



//=================================================================================================

/** Baseclass for all sample-map generators. */

class SampleMapGenerator
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  SampleMapGenerator();

  /** Destructor. */
  virtual ~SampleMapGenerator();


  /** \name Setup */

  /** Sets the directory into which the samples and mapping files will be written. The soundfont
  file will end up directly in this directory and for the samples, a subdirectory will be
  created. */
  void setOutputDirectory(const rsString& newDirectory);

  /** Sets the name for the sample-map to be generated. This also determines the name for the
  sample subdirectory. */
  void setName(const rsString& newName);

  /** Sets the range of keys for which samples should be rendered. Note that the keyrange of the
  output sample map is always the full keyboard - normally less than 127 samples are rendered in
  which case the lowest and highest rendered sample will have regions which extend to lowest or
  highest key on the keyboard respectively. */
  void setKeyRangeToRender(int newLowKey, int newHighKey);



  /** \name Inquiry */

  /** Returns the estimated total size of the sample-map to be generated in bytes.  */
  //virtual long getTotalSizeEstimate() const = 0;

  /** Returns the filename that will be assigned to the sample with the given key. */
  rsString getSampleFileName(int key);

  /** Returns the relative path with respect to the sfz file of the sample for the given key. */
  rsString getSampleRelativePath(int key);

  /** Returns the relative path with respect to the sfz file of the sample for the given key. */
  rsString getAmbienceRelativePath(int key);


  /** \name Sample-map generation: */

  /** This is the single high-level function that you need to call to generate the whole
  sample-map and write it to disk. If printProgress is true, an indication of the progress will
  be written to the console output. */
  void generateSampleMap(bool printProgress);

  /** Creates the soundfont file with all the key-mapping and stuff and writes it to disk.  */
  void generateSoundFontFile();

  /** Generates all samples (for all keys) and writes them out to disk. */
  void generateAllSamples(bool printProgress);

  /** Generates a sample for a particular key and writes it out to disk. The file written will be
  normalized to 0 dBFS and the normalization factor that has been used is stored internally to be
  used later to set up a compensation gain in the sample-map file */
  void generateSampleForKey(int key);


protected:

  /** Override this to actually generate a sample for a particular key. You are supposed to write
  your result into the member "buffer". You should not normalize the data there - this will be
  taken care of by higher level functions. */
  virtual void writeUnnormalizedSampleForKeyToBuffer(int key) = 0;

  /** Genereates a string for a region in the sfz file. The root determines the rootkey of the
  sample to be used (and hence its filename) and lo and hi determine the lower and upper limit
  of the region for which this sample should be used. */
  rsString getRegionString(int root, int lo, int hi);


  /** \name Data */

  rsString outputDirectory;
  rsString sampleSubDirectory;

  rsString name;

  double gains[128];    // normalization gains for the samples
  double ambGains[128]; // normalization gains for the ambience samples

  rsAudioBuffer buffer;


  double sampleRate;  // sample-rate for the wavefiles
  int    numBits;     // number of bits for the wavefiles

  int loKey, hiKey;   // lowest and highest key for which a sample has to be generated

  bool ambience;      // determines whether ambience should be generated

  //static const int fftSize = 262144;
  static const int fftSize = 131072;
  //rsFourierTransformerRadix2 ft; // for the ambience
  RAPT::rsFourierTransformerRadix2<double> ft; // for the ambience ..// use rosic::rsFourierTransformerRadix2_D


  // \todo later provide facilities for velocity layers


};

//=================================================================================================

/** Sample-map generator using modal synthesis. */

class SampleMapGeneratorModal : public SampleMapGenerator
{

public:

  /** \name Construction/Destruction */

  /** Constructor. */
  SampleMapGeneratorModal();

  /** Destructor. */
  virtual ~SampleMapGeneratorModal();


  /** \name Setup */

  void setModalParametersForKey(int key, const rsModalBankParametersD& newParameters);

  /** Sets the level at which the rendering will be cut off. A reasonable level is -80 dB which is
  the default value. It's been made user-adjustable mainly to have a means to create shorter
  preview files. A few additional cycles will be generated for a smooth fade-out. */
  void setTruncationLevel(double newLevel)
  {
    truncationLevel = newLevel;
  }


protected:

  /** Overriden to actually do the modal synthesis. */
  virtual void writeUnnormalizedSampleForKeyToBuffer(int key);

  /** Interpolates between two sets of modal parameters given in p1, p2 according to the given
  proportion (which should be a number between 0...1. */
  rsModalBankParametersD interpolateParameters(const rsModalBankParametersD &p1,
    const rsModalBankParametersD &p2, double proportion);

  /** Returns a truncated parameter set which contains only those modes which are below the
  given cutoff frequency. */
  rsModalBankParametersD removeModesAbove(const rsModalBankParametersD &p, double cutoff);


  /** \name Data */

  rsModalBankParametersD *keyParameters[128];
    // parameters to be used for a particular key - we use a pointer array to have the option not
    // have specified parameters for each key but only a few of them (for example, for one key in
    // each octave) - whenever the pointer is NULL, an interpolation of the parameters between the
    // two next neighbours that are defined will be used. if there's only one neighbour, it's
    // parameters will be used as is (which amounts to constant extrapolation)

  //rsModalFilterBank modalFilterBank;
  RAPT::rsModalFilterBank<double, double> modalFilterBank; // use rosic::rsModalFilterBankDD

  double truncationLevel;

};


// ToDo:
// -write SampleMapGeneratorSpectralCluster, SampleMapGeneratorSinusoidal etc.
// -provide export/import functionality using a simple textfile format - maybe provide virtual
//  loadSettings/saveSettings functions
// -maybe rename writeUnnormalizedSampleForKeyToBuffer to something shorter
// -implement multithreading (we need one AudioBuffer object per thread)


#endif
