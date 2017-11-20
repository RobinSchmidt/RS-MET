#ifndef RS_FILE_H
#define RS_FILE_H

namespace RSLib
{

  /**

  This is a class for representing files. It is implemented as a wrapper around C file I/O 
  facilities. This class here opens and closes the file for each read/write operation - it is 
  intended mainly for reading/writing files at once. If you need to continuously stream data
  from/to a file, look at FileStream.

  */

  class RSLib_API rsFile
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. Creates a file object that is associated with the given absolute path. A file 
    with such path may or may not exist. */
    rsFile(const rsString& absolutePath = rsString::empty);

    /** Destructor. */
    virtual ~rsFile();

     
    /** \name Setup */

    virtual void setAbsolutePath(const rsString& newPath) 
    { 
      absolutePath = newPath; 
    }

    // setRelativePathFrom


    /** \name Inquiry */

    /** Returns the absolute path of the file. */
    virtual rsString getAbsolutePath() const 
    { 
      return absolutePath; 
    }

    /** Returns the size of the file. As the return value is of type 'long' (32-bit integer), the
    maximum possible value is 2^31-1 = 2147483647 bytes which is roughly 2 gigabytes. For files
    larger than that, the return value will be undefined, so use the function with care when
    dealing with larger files. */
    virtual long getSizeInBytes() const;

    // \todo exists, isDirectory, existsAsFile
    // static String getCurrentApplicationDirectory() / getCurrentSharedLibraryDirectory()

    /** Returns the absolute path of this file as String. */
    //virtual String getAbsolutePath() const { return absolutePath; }

    /** Returns a file that is the same as this but with a possibly different extension. It does 
    not matter whether or not you include the dot in the argument - the function iteself uses only 
    the portion of the passed string after the last dot.
    \todo check what happens when the original file has no extension at all
    */
    virtual rsFile withExtension(const rsString &extension) const;


    /** \name Reading */

    /** Reads the entire file as a String and returns it. If something goes wrong, the returned 
    String will be empty. */
    virtual rsString readFileAsString() const;

    /** Reads the entire file as zero-terminated c-string and returns the pointer to the character 
    array. The caller is responsible to delete it. When something goes wrong with opening the file 
    for reading, it will return a NULL pointer. The reference argument is used to return the length 
    of the character array excluding the terminating zero. */
    virtual char* readFileAsZeroTerminatedString(long &lengthExludingZero) const;


    /** \name Writing */

    /** Writes the passed string into the file and returns whether this was successful. If the file 
    exists, it will be overwritten (so be careful). It will fail when the passed string contains 
    non-printable characters or when the file could not be openend for writing. */
    virtual bool writeStringToFile(const rsString& text) const;

    // \todo createIfNonExistent
    // \todo operators (\todo: maybe implement (in)equality, assignment, comparison):


    /** \name Static Member Functions */

    /** Similar to the 'fopen' c-function but with a String as first parameter. */
    static FILE* openFile(const rsString& path, const char *mode);

    /** This represents a file that does not exist. */
    //static const File nonExistent;

    /** If a file with the given path/name already exists, it creates an alternative name for the 
    file by appending a string with some number. If no file with the given name exists, the 
    argument will be returned unchanged */
    //static String createUniqueFileName(const String& pathIncludingName);

    // separatorString, etc. ...see JUCE for inspiration

  protected:

    /** \name Data */

    rsString absolutePath;

  };

}

#endif
