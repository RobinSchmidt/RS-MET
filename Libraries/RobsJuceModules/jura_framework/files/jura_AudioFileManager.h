#ifndef jura_AudioFileManager_h
#define jura_AudioFileManager_h

//#include "rojue_FileManager.h"


/** This class manages audio files...  */

class JUCE_API AudioFileManager : virtual public FileManager
{
  // why do we need virtual inheritance from FileManager? I think, it is because AudioModule 
  // derives from StateFileManager and some subclasses of it must also handle other types of files
  // such as .wav and .sfz?

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  AudioFileManager();

  /** Destructor. */
  virtual ~AudioFileManager();

  //-----------------------------------------------------------------------------------------------
  // FileManager overrides:

  /** Loads a new audio-file and informs whether or not this was operation was successful. */
  virtual bool loadFile(const juce::File& fileToLoad);

  /** Retrieves the current audio-data via getAudioData and saves it into a file.
  !!!NOT YET FUNCTIONAL!!! */
  virtual bool saveToFile(const juce::File& fileToSaveTo);

  //-----------------------------------------------------------------------------------------------
  // others:

  /** Tries to create an AudioSampleBuffer from a File and returns a pointer to it if it was
  successful (and a NULL-pointer otherwise). In the case of success, it is the callers
  responsibility to delete the AudioSampleBuffer when it is no longer needed. If something
  goes wrong, NULL will be returned and an alert-box pops up by default - by passing false as
  second parameter, the box can be supressed. */
  static AudioSampleBuffer* createAudioSampleBufferFromFile(const juce::File& theAudioFile,
    bool showAlertBoxWhenFailed = true);

  /** Same as above except that it takes a String as parameter representing the path to the File.
  The second parameter indicates whether the path should be interpreted as relative to the
  current  executable file (or .dll) - if not, it will be treated as absolute path. */
  //static AudioSampleBuffer* createAudioSampleBufferFromFile(const juce::String& filePath,
  //  bool pathIsRelativeToCurrentExecutable = true, bool showAlertBoxWhenFailed = true);

  // \todo: return a pointer to AudioFileBuffer instead of AudioSampleBuffer

  /** Saves the passed AudioSampleBuffer into a file. */
  //static void saveAudioFile(const File& fileToSaveTo, const AudioSampleBuffer &bufferToSave, 
  //  int sampleRate, int numBits, const String &format = String(T(".flac")));

protected:

  /** Creates an AudioSampleBuffer from the passed file, and - if this succeeds - calls
  setAudioData with the buffer so created. */
  virtual bool loadAudioFile(const juce::File& fileToLoad);

  /** !!!NOT YET FUNCTIONAL!!! */
  virtual bool saveAudioFile(const juce::File& fileToSaveTo, const AudioSampleBuffer& bufferToSave,
    int sampleRate, int numBits, const juce::String &format);

  /** Must be overriden by subclasses to update their internally buffered audio data according to
  the new data which come form the 'underlyingFile'. The passed buffer will not be NULL and the
  overriden function shall not delete it. */
  virtual bool setAudioData(AudioSampleBuffer* newBuffer, const juce::File& underlyingFile,
    bool markAsClean) = 0;

  ///** Must be overriden by subclasses to pass their internally buffered audio data. */
  ////virtual AudioSampleBuffer* getAudioData() = 0;

  juce_UseDebuggingNewOperator;
};

#endif  