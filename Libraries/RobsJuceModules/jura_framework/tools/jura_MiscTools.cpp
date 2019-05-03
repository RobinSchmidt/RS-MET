
//juce::MouseEvent createDummyMouseEvent()
//{
//  return MouseEvent(MouseInputSource(0, true), Point<int>(0, 0), ModifierKeys(0), NULL, NULL, 
//    Time::getCurrentTime(), Point<int>(0, 0), Time::getCurrentTime(), 1, false);
//}
  
void writeToLogFile(const String& logString)
{
  juce::File logFile = 
    File::getSpecialLocation(File::currentApplicationFile).withFileExtension(String("log"));
  if( !logFile.existsAsFile() )
    logFile.create();
  logFile.appendText(logString);
}

void clearLogFile()
{
  juce::File logFile = 
    File::getSpecialLocation(File::currentApplicationFile).withFileExtension(String("log"));
  if( logFile.existsAsFile() )
    logFile.deleteFile();
}

int getAvailableScreenPixelsBelow(const juce::Component* c)
{
  return c->getParentMonitorArea().getBottom() - c->getScreenBounds().getBottom();
}
