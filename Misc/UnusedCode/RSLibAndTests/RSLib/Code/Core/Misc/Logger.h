#ifndef RS_LOGGER_H
#define RS_LOGGER_H

namespace RSLib
{

  /**

  This is a class for logging and printing out (into a textfile) information during runtime of an 
  application (or shared library), mainly intended for debugging purposes.

  */

  class RSLib_API rsLogger
  {

  public:

    /** Constructor. */
    rsLogger();

    /** Destructor. */
    ~rsLogger();

    /** Writes a string into the log-file, optionally appending additional information such as the
    current system time, the thread-name, etc. */
    void appendLogString(const rsString& stringToAppend, bool appendAdditionalInfo);

    /** Sets up the file into which log-strings will be written. */
    void setLogFile(const rsFile &newLogFileToUse);

  protected:

    rsFileStream logFile;

  };

}

#endif
