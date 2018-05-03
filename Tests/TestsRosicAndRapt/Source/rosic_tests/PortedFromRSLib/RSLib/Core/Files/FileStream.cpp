using namespace RSLib;

// construction/destruction:

rsFileStream::rsFileStream(const rsString& absolutePath) : rsFile(absolutePath)
{
  fileState   = CLOSED;
  filePointer = NULL;
}

rsFileStream::~rsFileStream()
{
  closeIfOpen();
}

// setup:

void rsFileStream::setFileToUse(const rsFile &newFileToUse)
{
  closeIfOpen();
  absolutePath = newFileToUse.getAbsolutePath();
}

// opening/closing:

void rsFileStream::closeIfOpen()
{
  if( filePointer != NULL && fileState != CLOSED )
    close();
}

bool rsFileStream::close()
{
  bool success = fclose(filePointer) != 0;
  fileState    = CLOSED;
  filePointer  = NULL;
  return success;
}

// read/write:

bool rsFileStream::appendText(const rsString &textToAppend)
{
  char* tmpStr = textToAppend.getAsZeroTerminatedString();
  bool  result = appendText(tmpStr, textToAppend.getLength());
  delete tmpStr;
  return result;
}

// internal functions:

bool rsFileStream::open(const char *modeString, const int newStateIfSuccess)
{
  closeIfOpen();
  filePointer = openFile(absolutePath, modeString);
  if( filePointer == NULL )
  {
    fileState = CLOSED;
    rsError("Unable to open file");
    return false;
  }
  else
  {
    fileState = newStateIfSuccess;
    return true;
  }
}

long rsFileStream::getNumberOfBytesAfterCurrentPosition()
{
  long currentPositionInBytes = ftell(filePointer);
  bool fseekSuccess           = fseek(filePointer, 0, SEEK_END) == 0;
  long endPositionInBytes     = ftell(filePointer);
  long remainingBytes         = endPositionInBytes - currentPositionInBytes;
  fseekSuccess                = fseek(filePointer, currentPositionInBytes, SEEK_SET) == 0;
  return remainingBytes;
}
