#ifndef jura_AudioFileInfo_h
#define jura_AudioFileInfo_h

//#include "../misc/rojue_MessageBoxes.h"

/** This class represents meta information about an audio file such as format, sample-rate, etc.
The file in question must be passed to the constructor, this class will then try to retrieve the
information or - if this fails - will intialize the members to default values and set the flag 
isValidAudioFile to false.

\todo: make data members protected and provide access-functions

*/

class JUCE_API AudioFileInfo
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AudioFileInfo(juce::File fileToExtractInfoFrom = juce::File::nonexistent);

  /** Destructor. */
  virtual ~AudioFileInfo();

  //-----------------------------------------------------------------------------------------------
  // setup:

  /** Sets the File to which this AudioFileInfo object belongs. */
  void setFile(const juce::File& newFile);

  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the name of the audio-file. */
  juce::String getName() const;

  /** Returns the extension of the audio-file. */
  juce::String getExtension() const;

  /** Returns the File to which this AudioFileInfo object belongs. */
  juce::File getFile() const;

  //-----------------------------------------------------------------------------------------------
  // public data members:

  /** Flag to indicate whether or not the file is a valid audio file. */
  bool isValidAudioFile;

  /** The samplerate at whcih the file was recorded/created. */
  double sampleRate;

  /** The number of samples (or sample-frames in case of multichannel files) */
  int numSamples;

  /** The number of channels in the file. */
  int numChannels;

  /** The bit-depth of the audio file. */
  int bitDepth;

protected:

  juce::File theFile;

  juce_UseDebuggingNewOperator;
};

#endif  