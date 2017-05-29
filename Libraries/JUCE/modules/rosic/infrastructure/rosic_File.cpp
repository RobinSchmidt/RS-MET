//#include "rosic_File.h"
//using namespace rosic;

//-------------------------------------------------------------------------------------------------
// static data members:

const rsFile rsFile::nonExistent;

//-------------------------------------------------------------------------------------------------
// construction/destruction:

rsFile::rsFile(const rsString& absolutePath)
{
  this->absolutePath = absolutePath;
}

rsFile::~rsFile()
{

}

//-------------------------------------------------------------------------------------------------
// static member functions:



//-------------------------------------------------------------------------------------------------
// setup:



//-------------------------------------------------------------------------------------------------
// inquiry:
   
long rsFile::getSizeInBytes() const
{
  FILE *f = fopen(absolutePath.getRawString(), "r");
  if( f == NULL )
    return 0;
  else
  {
    int fseekSuccess = fseek(f, 0, SEEK_END);
    if( fseekSuccess != 0 ) // fseek returns 0 in case of success
      return 0;
    else
    {
      return ftell(f); 
      // ftell returns current position in the file which is now set to the last byte of the file 
      // due to passing 0 and SEEK_END as arguments to fseek
    }
  }
}

//-------------------------------------------------------------------------------------------------
// reading:

rsString rsFile::readFileAsString() const
{
  long numCharacters;
  char *textRaw = readFileAsZeroTerminatedString(numCharacters);
  if( textRaw == NULL )
    return rsString();
  else
  {
    rsString textAsString;
    textAsString.reserveSize(numCharacters);
    textAsString = textRaw; // this can be optimized by a method String::copyFromCharArray - avoids one String constructor call
    delete[] textRaw;
    return textAsString;
  }
}

char* rsFile::readFileAsZeroTerminatedString(long &lengthExludingZero) const
{
  FILE *f = fopen(absolutePath.getRawString(), "rt");
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

//-------------------------------------------------------------------------------------------------
// writing:

bool rsFile::appendText(const rosic::rsString &text) const
{
  if( text.containsNonPrintableCharacters() )
    return false;

  FILE *f = fopen(absolutePath.getRawString(), "wt");
  if( f == NULL )
    return false;
  else
  {
    fwrite(text.getRawString(), sizeof(char), text.getLength(), f);
    fclose(f);
    return true;
  }
}

//-------------------------------------------------------------------------------------------------
// internal functions:

