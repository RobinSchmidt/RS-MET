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

// update: RAPT is now included, so we may delete them? ...clean this up!

static const double INF = std::numeric_limits<double>::infinity(); // move to rapt

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
// ...there already are corresponding function...but they use size_t insetad of int for the indices

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
inline void insert(vector<T>& v, T newElement, int index)
{
  v.insert(v.begin() + index, newElement);
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
inline bool removeFirstOccurrence(vector<T>& v, T elementToRemove)
{
  for(int i = 0; i < size(v); i++)
    if(v[i] == elementToRemove){
      remove(v, i);
      return true; }
  return false;
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

template<class T>
inline T getAndRemoveLast(vector<T>& v)
{
  T result = v[v.size()-1];
  v.pop_back();
  return result;
}

// maybe move into file ComponentTools

/** Sets the position of the right border while keeping the left as is */
inline void setRightKeepLeft(juce::Component* c, int newRight)
{
  int x = c->getX();
  int y = c->getY();
  int w = c->getWidth();
  int h = c->getHeight();
  int oldRight = x + w;
  int newWidth = newRight - x;
  c->setBounds(x, y, newWidth, h);
}





#endif