//#include "rojue_AudioFileManager.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AudioFileManager::AudioFileManager()
{
  wildcardPatterns = String("*.wav;*.flac");
  defaultExtension = String(".flac");
  updateFileList();
}

AudioFileManager::~AudioFileManager()
{

}

//-------------------------------------------------------------------------------------------------
// static member functions:

AudioSampleBuffer* AudioFileManager::createAudioSampleBufferFromFile(const File& theAudioFile,
                                                                     bool showAlertBoxWhenFailed)
{
  WavAudioFormat     wavAudioFormat;
  FlacAudioFormat    flacAudioFormat;
  bool               fileIsReadable = false;
  //bool               success        = false;
  AudioFormatReader* reader         = NULL;

  if( !theAudioFile.existsAsFile() )
  {
    RAPT::rsError("File not found");
    if( showAlertBoxWhenFailed == true )
    {
      AlertWindow::showMessageBox(AlertWindow::WarningIcon,
        String("Alert"),
        String("Audio-File: " + theAudioFile.getFullPathName() + " not found."),
        String("OK") );
    }
    return NULL;
  }

  // check whether the file is a flac-encoded spectrum (.spec), these have to
  // treated a bit differently:
  //bool fileContainsSpectrum = false;
  //if( theAudioFile.hasFileExtension(String("spec")) )  // maybe get rid of this..
  //  fileContainsSpectrum = true;


  FileInputStream* inputStream = new FileInputStream(theAudioFile);
  if( inputStream == NULL )
  {
    RAPT::rsError("File could not be opened");
    if( showAlertBoxWhenFailed == true )
    {
      AlertWindow::showMessageBox(AlertWindow::WarningIcon,
        String("Alert"),
        String("Audio-File: " + theAudioFile.getFullPathName() + " could not be opened."),
        String("OK") );
    }
    return NULL;
  }

  if( theAudioFile.hasFileExtension(String("wav")) )
  {
    fileIsReadable = wavAudioFormat.canHandleFile(theAudioFile);
    if( fileIsReadable )
      reader = wavAudioFormat.createReaderFor(inputStream, true);
  }
  else if( theAudioFile.hasFileExtension(String("flac")) )
  {
    fileIsReadable = flacAudioFormat.canHandleFile(theAudioFile);
    if( fileIsReadable )
      reader = flacAudioFormat.createReaderFor(inputStream, true);
  }
  else if( theAudioFile.hasFileExtension(String("spec")) )
  {
    fileIsReadable = flacAudioFormat.canHandleFile(
      theAudioFile.withFileExtension(String("flac")));
    if( fileIsReadable )
      reader = flacAudioFormat.createReaderFor(inputStream, true);
  }
  if( reader == NULL )
  {
    RAPT::rsError("File not readable");
    if( showAlertBoxWhenFailed == true )
    {
      AlertWindow::showMessageBox(AlertWindow::WarningIcon,
        String("Alert"),
        String("Audio-File: " + theAudioFile.getFullPathName() + " is not readable."),
        String("OK") );
    }
    return NULL;
  }

  if( fileIsReadable )
  {
    // retrieve some info:
    int numChannels     = (int) reader->numChannels;
    int numSampleFrames = (int) reader->lengthInSamples;

    AudioSampleBuffer* buffer = new AudioSampleBuffer(numChannels, numSampleFrames);

    //buffer->readFromAudioReader(reader, 0, numSampleFrames, 0, true, true); //old
    reader->read(buffer, 0, numSampleFrames, 0, true, true);

    delete reader; // we must take care of deleting the reader
    return buffer;
  } // end of  if( fileIsReadable )
  else
  {
    RAPT::rsError("File format not supported");
    if( showAlertBoxWhenFailed == true )
    {
      AlertWindow::showMessageBox(AlertWindow::WarningIcon,
        String("Alert"),
        String("Audio-File: " + theAudioFile.getFullPathName() + " is not supported (wrong format)."),
        String("OK") );
    }
    delete reader;
    return NULL;
  }
}

/*
// old - may be deleted at some point:
AudioSampleBuffer* AudioFileManager::createAudioSampleBufferFromFile(const String& filePath,
  bool pathIsRelativeToCurrentExecutable, bool showAlertBoxWhenFailed)
{
  String fullPath;

  // complete the relative path to an absolute path if necesarry:
  if( pathIsRelativeToCurrentExecutable )
  {
    // retrieve the directory of the current executable file (or dynamic link library):
    //File   thisExeAsFile         = File::getSpecialLocation(File::currentExecutableFile);
    File   thisExeAsFile         = File::getSpecialLocation(File::currentApplicationFile);
    File   thisDirectoryAsFile   = thisExeAsFile.getParentDirectory();
    String thisDirectoryAsString = thisDirectoryAsFile.getFullPathName();

    // add the directory-string and the filePath (seperating them by a slash) to form the full
    // path:
    fullPath = thisDirectoryAsString + String("/") + filePath;
  }
  else
    fullPath = filePath; // the filePath is already absolute - we can use it as is
  
  fullPath = fullPath.replaceCharacter('\\', '/');

  return createAudioSampleBufferFromFile(File(fullPath), showAlertBoxWhenFailed);
}
*/

//-------------------------------------------------------------------------------------------------
// FileManager overrides:

bool AudioFileManager::loadFile(const File& fileToLoad)
{
  bool success = loadAudioFile(fileToLoad);
  if( success == true )
    markFileAsClean(true);
  notifyListeners();
  return success;
}

bool AudioFileManager::saveToFile(const File& fileToSaveTo)
{
  return false; // not yet functional
}

//-------------------------------------------------------------------------------------------------
// others:

bool AudioFileManager::loadAudioFile(const File& fileToLoad)
{
  AudioSampleBuffer* buffer = createAudioSampleBufferFromFile(fileToLoad);
  bool success = (buffer != NULL);
  if( success )
  {
    success = setAudioData(buffer, fileToLoad, true);
    delete buffer;
    setActiveFile(fileToLoad);
  }
  return success;
}

bool AudioFileManager::saveAudioFile(const juce::File &fileToSaveTo,
                                     const juce::AudioSampleBuffer &bufferToSave, int sampleRate,
                                     int numBits, const juce::String &format)
{
  jassertfalse;
  return false;  // not yet functional
}


// this seems to be some obsolete old code - maybe delete at some point:
/*
AudioSampleBuffer* AudioFileManager::loadAudioFile(const File& fileToLoadFrom)
{
// check whether the file is a flac-encoded spectrum (.spec), these have to
// treated a bit differently:
bool fileContainsSpectrum = false;
if( fileToLoadFrom.hasFileExtension(String(T("spec"))) )
fileContainsSpectrum = true;

FileInputStream inputStream(fileToLoadFrom);

WavAudioFormat  wavAudioFormat;
FlacAudioFormat flacAudioFormat;
bool            fileIsReadable = false;
bool            success        = false;

AudioFormatReader* reader;

if( fileToLoadFrom.hasFileExtension(String(T("wav"))) )
{
fileIsReadable = wavAudioFormat.canHandleFile(fileToLoadFrom);
if( fileIsReadable )
reader = wavAudioFormat.createReaderFor(&inputStream, false);
}
else if( fileToLoadFrom.hasFileExtension(String(T("flac"))) )
{
fileIsReadable = flacAudioFormat.canHandleFile(fileToLoadFrom);
if( fileIsReadable )
reader = flacAudioFormat.createReaderFor(&inputStream, false);
}
else if( fileToLoadFrom.hasFileExtension(String(T("spec"))) )
{
fileIsReadable = flacAudioFormat.canHandleFile(
fileToLoadFrom.withFileExtension(String(T("flac"))));
if( fileIsReadable )
reader = flacAudioFormat.createReaderFor(&inputStream, false);
}

if( fileIsReadable )
{
// retrieve some info:
int numChannels     = (int) reader->numChannels;
int numSampleFrames = (int) reader->lengthInSamples;

AudioSampleBuffer* buffer = new AudioSampleBuffer(numChannels, numSampleFrames);

buffer->readFromAudioReader(reader, 0, numSampleFrames, 0, true, true);

// remember where we are:
currentFileFullPath = fileToLoadFrom.getFullPathName();
currentDirectory    = fileToLoadFrom.getParentDirectory().getFullPathName();

// update the preset-field on the GUI
currentFileName         = fileToLoadFrom.getFileName();
currentFileNameWithStar = currentFileName + String(T("*"));

return buffer;
} // end of  if( fileIsReadable )
else
{
AudioSampleBuffer* dummyBuffer = new AudioSampleBuffer(1, 1);
dummyBuffer->clear();
return dummyBuffer;
}

}
*/




