#ifndef rosic_File_h
#define rosic_File_h

//#include "../datastructures/rosic_String.h"

namespace rosic
{

/** This is a class for representing files. It is implemented as a wrapper around C file I/O
facilities.

\todo: create class FileStream that allows for opening a file and keeping it open to
continuously read/write a stream of data into it - this class here opens and closes the file for
each read/write operation ...in RSLib, there is already such a thing, try to merge the code from 
there - perhaps, we should also use (or merge) the implementation of rsString from there */

class rsFile
{

public:

  //-----------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. Creates a file object that is associated with the given absolute path. A file
  with such path may or may not exist. */
  rsFile(const rsString& absolutePath = rsString());

  /** Destructor. Closes the underlying C file in case it is still open. */
  ~rsFile();

  //-----------------------------------------------------------------------------------------------
  // setup:

  void setAbsolutePath(const rsString& newPath) { absolutePath = newPath; }

  // setRelativePathFrom


  //-----------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the absolute path of the file. */
  rsString getAbsolutePath() const { return absolutePath; }

  /** Returns the size of the file. As the return value is of type 'long' (32-bit integer), the
  maximum possible value is 2^31-1 = 2147483647 bytes which is roughly 2 gigabytes. For files
  larger than that, the return value will be undefined, so use the function with care when
  dealing with larger files. */
  long getSizeInBytes() const;

  /** Returns whether or not the file is open. Note that this reveals no information in which
  mode it is open (read/write, text/binary, etc.) */
  //bool isOpen() const { return fileIsOpen; }


  // exists, isDirectory, existsAsFile, createUniqueName, isOpen



  //-----------------------------------------------------------------------------------------------
  // reading

  /** Reads the entire file as a String and returns it. If something goes wrong, the returned
  String will be empty. */
  rsString readFileAsString() const;


  /** Reads the entire file as zero-terminated c-string and returns the pointer to the character
  array. The caller is responsible to delete it. When something goes wrong with opening the file
  for reading, it will return a NULL pointer. The reference argument is used to return the length
  of the character array excluding the terminating zero. */
  char* readFileAsZeroTerminatedString(long& lengthExludingZero) const;

  //-----------------------------------------------------------------------------------------------
  // writing:

  /** Appends the passed text at the end of the file and returns whether this was successful. It
  will fail when the passed string contains non-printable characters or when the file could not
  be openend for writing. */
  bool appendText(const rsString& text) const;


  // openForTextWriting, openForTextReading, openForBinaryWriting, openForBinaryReading, ...
  // createIfNonExistent, close


  //-----------------------------------------------------------------------------------------------
  // operators (\todo: maybe implement (in)equality, assignment, comparison):


  //---------------------------------------------------------------------------------------------
  // static functions and data members:

  /** This is an empty string. */
  //static const rsFile nonExistent; 
  // causes memory leak because of memory allocation in absolutePath member

  /** If a file with the given path/name already exists, it creates an alternative name for the
  file by appending a string with some number. If no file with the given name exists, the
  argument will be returned unchanged */
  //static String createUniqueFileName(const String& pathIncludingName);

  // separatorString, etc. ...see JUCE for inspiration

  //=============================================================================================

protected:

  rsString absolutePath;

};

} // end namespace rosic

#endif 