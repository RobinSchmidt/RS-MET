using namespace RSLib;

rsLogger::rsLogger()
{
  rsFile appFile = RSLib::getCurrentApplicationFile();
  logFile.setFileToUse(appFile.withExtension("log"));
  logFile.openForWrite();
}

rsLogger::~rsLogger()
{
  logFile.closeIfOpen();
}

void rsLogger::appendLogString(const rsString& stringToAppend, bool appendAdditionalInfo)
{
  logFile.appendText(stringToAppend);
  if( appendAdditionalInfo == true )
  {
    // append system time, thread-name, etc.
  }
}

void rsLogger::setLogFile(const rsFile &newLogFileToUse)
{
  logFile.closeIfOpen();
  logFile.setFileToUse(newLogFileToUse);
  logFile.openForWrite();
}
