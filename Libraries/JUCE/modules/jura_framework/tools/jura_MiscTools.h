#ifndef jura_MiscTools_h
#define jura_MiscTools_h

//#define LOG_SESSION  
// Uncomment this, if you want to activate logging. If active, the macros 
// CLEAR_LOGFILE, WRITE_TO_LOGFILE(logString) will result in clearing and writing a string into a
// log file. That may be useful for debugging. If not defined, the macros will do nothing. See the
// comments of the functions writeToLogFile, clearLogFile

#ifdef LOG_SESSION
  #define WRITE_TO_LOGFILE(logString) (writeToLogFile(logString))
  #define CLEAR_LOGFILE clearLogFile()
#else
  #define WRITE_TO_LOGFILE(logString) 
  #define CLEAR_LOGFILE
#endif
  
///** Returns a dummy mouse event with some (usually meaningless) default values for the members. You
//can use this when some function requires a MouseEvent to be passed but you don't have one lying 
//around and the details of the event do not really matter. The function is needed due to lack of 
//default-constructor in class MouseEvent. */
//JUCE_API juce::MouseEvent createDummyMouseEvent();

/** Writes teh given string into a logfile which will end up in the same directory as the main
application file and will have a name like MyApplication.log */
JUCE_API void writeToLogFile(const juce::String& logString);

/** Clears the log file @see writeToLogFile. */
JUCE_API void clearLogFile();


// these are copied from RSLib - they also occur in RAPT in the same form - but we don't want to
// include RAPT already here, so we need to copy them (that violates the DRY principle, but i 
// currently don't know, how to do this any better):

#define PI 3.1415926535897932384626433832795

inline double clip(double x, double min, double max)
{
  if(x > max)
    return max;
  if(x < min)
    return min;
  return x;
}
inline double amp2dB(double amp)
{
  return 8.6858896380650365530225783783321 * log(amp);
}
inline double amp2dBWithCheck(double amp, double lowAmplitude)
{
  if( amp >= lowAmplitude )
    return amp2dB(amp);
  else
    return amp2dB(lowAmplitude);
}
inline double secondsToBeats(double timeInSeconds, double bpm)
{
  return timeInSeconds*(bpm/60.0);
}
inline double beatsToSeconds(double beat, double bpm)
{
  return (60.0/bpm)*beat;
}
inline double linToLin(double in, double inMin, double inMax, double outMin, double outMax)
{
  double tmp = (in - inMin) / (inMax - inMin);
  tmp *= (outMax - outMin);
  tmp += outMin;
  return tmp;
}
inline double linToExp(double in, double inMin, double inMax, double outMin, double outMax)
{
  double tmp = (in - inMin) / (inMax - inMin);
  return outMin * exp(tmp * (log(outMax / outMin)));
}
inline double expToLin(double in, double inMin, double inMax, double outMin, double outMax)
{
  double tmp = log(in / inMin) / log(inMax / inMin);
  return outMin + tmp * (outMax - outMin);
}


// get rid of rs prefix
// if we don't want to use inlining, the implementations must be moved to the cpp file to avoid
// linker errors ("function already defined")


#endif