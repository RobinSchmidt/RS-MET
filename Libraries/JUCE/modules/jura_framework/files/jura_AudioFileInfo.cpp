//#include "rojue_AudioFileInfo.h"
//using namespace rojue;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

AudioFileInfo::AudioFileInfo(File fileToExtractInfoFrom)
{
  // init the members:
  isValidAudioFile = false;
  sampleRate       = 44100.0;
  numSamples       = 0;
  numChannels      = 0;
  bitDepth         = 0;
  theFile          = fileToExtractInfoFrom;

  AudioFormatManager formatManager;
  formatManager.registerBasicFormats();
  AudioFormatReader* reader = formatManager.createReaderFor(fileToExtractInfoFrom);
  {
    if( reader != NULL )
    {
      if( reader->lengthInSamples > INT_MAX )
      {
        //showAudioFileTooLargeErrorBox();
        delete reader;
        return;
      }
      isValidAudioFile = true;
      sampleRate       = reader->sampleRate;
      numSamples       = (int) reader->lengthInSamples;
      numChannels      = reader->numChannels;
      bitDepth         = reader->bitsPerSample;
      delete reader;
    }
  }
}

AudioFileInfo::~AudioFileInfo()
{

}

//-------------------------------------------------------------------------------------------------
// setup:

void AudioFileInfo::setFile(const File& newFile)
{
  theFile = newFile;
}

//-------------------------------------------------------------------------------------------------
// inquiry:

String AudioFileInfo::getName() const
{
  return theFile.getFileName();
}

String AudioFileInfo::getExtension() const
{
  return theFile.getFileExtension();
}

File AudioFileInfo::getFile() const
{
  return theFile;
}
