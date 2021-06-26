#ifndef RS_FILESTREAM_H
#define RS_FILESTREAM_H

namespace RSLib
{

  /**

   This subclass of File allows for opening a file and keeping it open to continuously read/write a
   stream of data from/to it. This class does not do any automatic conversions of text-files such 
   as CR/LF stuff. All files are treated as binary and all data is read and written as-is.

  */

  class RSLib_API rsFileStream : public rsFile
  {

  public:

    enum fileStates
    {
      CLOSED,
      OPEN_FOR_WRITE,
      OPEN_FOR_READ,
      OPEN_FOR_APPEND
    };


    /** \name Construction/Destruction */

    /** Constructor. Creates a file object that is associated with the given absolute path. A file 
    with such path may or may not exist. */
    rsFileStream(const rsString& absolutePath = rsString::empty);

    /** Destructor. Closes the underlying file if it is still open. */
    virtual ~rsFileStream();


    /** \name Setup */

    /** Sets up a new file to use for streaming the data. If the old one is still open when this 
    function is called, it will be closed. */
    virtual void setFileToUse(const rsFile &newFileToUse);


    /** \name Inquiry */

    /** Returns whether or not the file is open. Note that this reveals no information in which 
    mode it is open (read/write, text/binary, etc.) */
    virtual bool isOpen() const { return fileState != CLOSED; }

    /** Returns true when the file is open for writing, false otherwise. */
    virtual bool isOpenForWrite() const { return fileState == OPEN_FOR_WRITE || OPEN_FOR_APPEND; }

    /** Returns true when the file is open for reading, false otherwise. */
    virtual bool isOpenForRead() const { return fileState == OPEN_FOR_READ; }


    /** \name Opening/Closing */

    /** Tries to open the file for writing and returns if this was successful. */
    virtual bool openForWrite()  { return open("wb", OPEN_FOR_WRITE); }

    /** Tries to open the file for reading and returns if this was successful. */
    virtual bool openForRead()   { return open("rb", OPEN_FOR_READ); }

    /** Tries to open the file for appending and returns if this was successful. */
    virtual bool openForAppend() { return open("ab", OPEN_FOR_APPEND); }

    /** Closes the file in case it is open (otherwise does nothing). */
    virtual void closeIfOpen();

    /** Tries to close the file and reports whether this was successful - it will fail if the file 
    was not open. */
    virtual bool close();


    /** \name Reading */

    /** Reads text data - convenience function. @see readData */
    virtual bool readText(char *textBuffer, int textBufferLength) 
    { 
      return readData(textBuffer, textBufferLength); 
    }

    /** Reads the given number of elements from the file and returns true when the desired number 
    of elements could be read. It may be that the file ends before the desired number of 
    elements could be read - in this case, the function will read up to the end of the file and
    return false. */
    template <class DataType>
    bool readData(DataType *dataBuffer, int numElementsToRead)
    {
      if( !isOpenForRead() )
        return false;
      int elementsRead = (int)fread(dataBuffer, sizeof (DataType), numElementsToRead, filePointer);
      if (elementsRead == numElementsToRead)
        return true;
      else
        return false;
    }


    /** \name Writing */

    /** Appends text data - convenience function. @see appendData */
    virtual bool appendText(const rsString &textToAppend);

    /** Appends text data - convenience function. @see appendData */
    virtual bool appendText(const char *textBuffer, int textBufferLength) 
    { 
      return appendData(textBuffer, textBufferLength); 
    }

    /** Appends the given number of elements of type 'DataType' to the end of the file. It will 
    fail when the file is not open for writing.*/
    template <class DataType>
    bool appendData(DataType *dataBuffer, int numElementsToWrite)
    {
      if( !isOpenForWrite() )
        return false;
      fseek(filePointer, 0, SEEK_END);
      int elementsWritten = (int)fwrite(dataBuffer, sizeof (DataType), numElementsToWrite, 
        filePointer);
      if( elementsWritten == numElementsToWrite )
        return true;
      else
        return false;
    }

  protected:

    /** \name Misc */

    /** Tries to open the file in one of various modes as given in the 'modeString'. The new value 
    to which fileState member will be assigned in case of success must also be passed along. */
    virtual bool open(const char *modeString, const int newStateIfSuccess);

    /** Returns the number of bytes in the file that are after the current position in the file. */
    virtual long getNumberOfBytesAfterCurrentPosition();


    /** \name Data */

    int fileState;
    FILE *filePointer;

  };

}

#endif
