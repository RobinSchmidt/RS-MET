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

/** Returns the available screen space in pixels below the given component (useful for positioning
dropdown menus. */
JUCE_API int getAvailableScreenPixelsBelow(const juce::Component* c);


// these are copied from RSLib - they also occur in RAPT in the same form - but we don't want to
// include RAPT already here, so we need to copy them (that violates the DRY principle, but i 
// currently don't know, how to do this any better):

static const double PI  = 3.1415926535897932384626433832795;
static const double INF = std::numeric_limits<double>::infinity();

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
inline int sign(double x)
{
  if(x > 0.f)
    return 1;
  else if(x < 0.f)
    return -1;
  else
    return 0;
}
inline void fillWithIndex(double *arrayToFill, int length)
{
  for(int i = 0; i < length; i++)
    arrayToFill[i] = (double)i;
}
inline void fillWithZeros(double *arrayToFill, int length)
{
  for(int i = 0; i < length; i++)
    arrayToFill[i] = 0.0;
}

inline int arrayMaxIndex(float* doubleArray, int numValues)
{
  int    maxIndex = 0;
  double maxValue = doubleArray[0];
  for(int i=0; i<numValues; i++)
  {
    if( doubleArray[i] > maxValue )
    {
      maxValue = doubleArray[i];
      maxIndex = i;
    }
  }
  return maxIndex;
}

inline int arrayMinIndex(float* doubleArray, int numValues)
{
  int    minIndex = 0;
  double minValue = doubleArray[0];
  for(int i=0; i<numValues; i++)
  {
    if( doubleArray[i] < minValue )
    {
      minValue = doubleArray[i];
      minIndex = i;
    }
  }
  return minIndex;
}

inline bool isCloseTo(double x, double targetValue, double tolerance)
{
  if(fabs(x-targetValue) <= tolerance)
    return true;
  else
    return false;
}

// some little helper/convenience functions to deal with std::vectors (move to RAPT):

using namespace std;

template<class T>
inline int size(const vector<T>& v)
{
  return (int)v.size();
}

template<class T>
inline void append(vector<T>& v, T newElement)
{
  v.push_back(newElement);
}

template<class T>
inline void remove(vector<T>& v, int index)
{
  v.erase(v.begin() + index);
}

template<class T>
inline int find(vector<T>& v, T elementToFind)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToFind)
      return i;
  return -1;
}

template<class T>
inline void removeFirstOccurrence(vector<T>& v, T elementToRemove)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToRemove){
      remove(v, i);
      return;
    }
}

template<class T>
inline bool contains(vector<T>& v, T elementToCheckFor)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToCheckFor)
      return true;
  return false;
}

template<class T>
inline void appendIfNotAlreadyThere(vector<T>& v, T newElement)
{
  if(!contains(v, newElement))
    append(v, newElement);
}




#endif