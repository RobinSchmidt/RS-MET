using namespace RSLib;

// static data members:
//const File File::nonExistent;

// construction/destruction:

rsFile::rsFile(const rsString& absolutePath)
{
  this->absolutePath = absolutePath;
}

rsFile::~rsFile()
{

}

// inquiry:

long rsFile::getSizeInBytes() const
{
  FILE *f = openFile(absolutePath, "rb");
  if( f == NULL )
    return 0;
  else
  {
    int size;
    if( fseek(f, 0, SEEK_END) != 0 ) // fseek returns 0 in case of success
      size = 0;
    else
      size = ftell(f);               // ftell returns current position in the file
    fclose(f);
    return size;
  }
}

rsFile rsFile::withExtension(const rsString &extension) const
{
  int finalDotIndex1 = absolutePath.rsFindLastOccurrenceOf('.');
  int finalDotIndex2 = extension.rsFindLastOccurrenceOf('.');
  rsString path;
  if( finalDotIndex1 >= 0 )
    path = absolutePath.getSubString(0, finalDotIndex1-1);
  else
    path = absolutePath;
  path += rsString(".") + extension.getSubString(finalDotIndex2+1, extension.getLength()-1);
  return path;
}

// reading:

rsString rsFile::readFileAsString() const
{
  FILE *f = openFile(absolutePath,  "rb");
  if( f == NULL )
    return rsString::empty;
  else
  {
    int  length = getSizeInBytes();
    char *text  = new char[length];
    void *voidPointer  = (void*) text;
    fread(voidPointer, sizeof(char), length, f);
    fclose(f);
    rsString theString;
    theString.ensureAllocatedSize(length);
    for(int i=0; i<length; i++)
      theString.appendElement(text[i]);  // \todo: optimize - append all at once
    delete[] text;
    return theString;
  }
}

char* rsFile::readFileAsZeroTerminatedString(long &lengthExludingZero) const
{
  FILE *f = openFile(absolutePath,  "rb");
  if( f == NULL )
    return NULL;
  else
  {
    lengthExludingZero = getSizeInBytes();
    char *text         = new char[lengthExludingZero+1];
    void *voidPointer  = (void*) text;
    fread(voidPointer, sizeof(char), lengthExludingZero, f);
    fclose(f);
    text[lengthExludingZero] = '\0'; // we must append the terminating zero manually
    return text;
  }
}

// writing:

bool rsFile::writeStringToFile(const rsString &text) const
{
  FILE *f = openFile(absolutePath,  "wb");
  if( f == NULL )
    return false;
  else
  {
    int length    = text.getLength();
    char *cString = text.getAsZeroTerminatedString();
    fwrite(cString, sizeof(char), length, f);
    delete[] cString;
    fclose(f);
    return true;
  }
}

// static functions:

FILE* rsFile::openFile(const rsString& path, const char *mode)
{
  char *cString    = path.getAsZeroTerminatedString();
  FILE *fileHandle = fopen(cString, mode);
  delete[] cString;
  return fileHandle;
}
