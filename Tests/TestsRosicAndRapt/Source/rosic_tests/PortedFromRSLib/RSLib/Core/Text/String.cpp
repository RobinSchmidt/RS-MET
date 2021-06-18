//#include "String.h"
//using namespace RSLib;



namespace RSLib
{

// static data members:

const rsString rsString::empty;

// construction/destruction:

rsString::rsString()
{

}

rsString::rsString(const char *initialString)
{
  setFromZeroTerminatedString(initialString);
}

rsString::rsString(const std::string &initialString)
{
  setFromStdString(initialString);
}

rsString::rsString(const rsString &other) : rsArrayTools<char>(0)
{
  initMembers();
  copyDataFrom(other);
}

rsString::rsString(const rsArrayTools<char> &stringAsCharrsArray)
{
  initMembers();
  copyDataFrom(stringAsCharrsArray);
}

rsString::rsString(const int intValue)
{
  setFromIntValue(intValue);
}

rsString::rsString(const float floatValue)
{
  setFromDoubleValue((double) floatValue);
}

rsString::rsString(const double doubleValue)
{
  setFromDoubleValue(doubleValue);
}

rsString::~rsString()
{

}

// setup:

void rsString::setFromZeroTerminatedString(const char *cString)
{
  if( cString != NULL )
  {
    int length = (int)strlen(cString);
    ensureAllocatedSize(length);
    numUsed = length;
    for(int i = 0; i < numUsed; i++)
      elements[i] = cString[i];
  }
  else
    clear();
}

void rsString::setFromStdString(const std::string &stdString)
{
  int length = (int)stdString.length();
  if( length != 0 )
  {
    ensureAllocatedSize(length);
    numUsed = length;
    for(int i = 0; i < numUsed; i++)
      elements[i] = stdString[i];
  }
  else
    clear();
}

void rsString::setFromIntValue(const int intValue)
{
  char cString[64];
  sprintf(cString, "%d", intValue);
  setFromZeroTerminatedString(cString);
}

void rsString::setFromDoubleValue(const double doubleValue)
{
  if( doubleValue == rsInfDouble )
  {
    *this = rsString("INF");
    return;
  }
  else if( doubleValue == -rsInfDouble )
  {
    *this = rsString("-INF");
    return;
  }
  //else if( _isnan(doubleValue) )
  //else if( doubleValue == NAN )
  else if( rsIsNaN(doubleValue) )
  {
    *this = rsString("NaN");
    return;
  }

  // create temporary string long enough to hold all the characters, determine the actually 
  // required length and copy the respective part of the temporary string - the c-book says, 
  // strncpy is unreliable with respect to copying the terminating zero, so we do it manually:
  char cString[64];
  sprintf(cString, "%-.17gl", doubleValue);
  setFromZeroTerminatedString(cString);
   // http://en.wikipedia.org/wiki/Double-precision_floating-point_format says:
   // If a decimal string with at most 15 significant digits is converted to IEEE 754 double 
   // precision representation and then converted back to a string with the same number of 
   // significant digits, then the final string should match the original; and if an IEEE 754 
   // double precision is converted to a decimal string with at least 17 significant digits and 
   // then converted back to double, then the final number must match the original.

  removeAllCharactersBut("0123456789+-.eE"); 
    // strips off special characters denoting denormals (these special characters would mess
    // with back-conversion) ....still true? -> check this
}

// static member functions:

/*
int String::numberOfRequiredCharacters(int number)
{
  if( number == 0 )
    return 1;
  int numChars  = 0;
  int absValue  = abs(number);
  long long tmp = 1;             // must be able to go above INT_MAX
  while( tmp <= absValue && tmp < INT_MAX )
  {
    tmp *= 10;
    numChars++;
  }
  if(number < 0)
    numChars++; // for the minus
  return numChars;
}
*/

int rsString::compareCharacterArrays(char *left, int leftLength, char *right, int rightLength)
{
  int minLength = rsMin(leftLength, rightLength);
  for(int i = 0; i < minLength; i++)
  {
    if(      compareCharacters(left[i], right[i]) == -1 )
      return -1;
    else if( compareCharacters(left[i], right[i]) == +1 )
      return +1;
  }

  // strings match up to minLength, now the shorter one is considered smaller:
  if(      leftLength < rightLength )
    return -1;
  else if( leftLength > rightLength )
    return +1;
  else
    return 0;
}

int rsString::compareCharacters(char left, char right)
{
  if( left == right )
    return 0;

  char leftUpper  = (char) toupper(left);
  char rightUpper = (char) toupper(right);
  if( leftUpper == rightUpper )
  {
    if( left < right )       // example: 'A' < 'a'
      return -1;
    else if( left > right )  // example: 'a' > 'A'
      return +1;
    else                     // case actually already handled in if(left == right)
      return 0;
  }
  else
  {
    if( leftUpper < rightUpper )
      return -1;
    else
      return +1;
  }
}

// inquiry:

char* rsString::getAsZeroTerminatedString() const
{
  char *cString = new char[numUsed+1];
  rsCopyBuffer(elements, cString, numUsed);
  cString[numUsed] = '\0';
  return cString;
}

std::string rsString::getAsStdString() const
{
  char *cString = getAsZeroTerminatedString();
  std::string s(cString);
  delete cString;
  return s;
}

int rsString::asInt() const
{
  char* cString = getAsZeroTerminatedString();
  int   result  = atoi(cString);
  delete[] cString;
  return result;
}

double rsString::asDouble() const
{
  if( *this == rsString("INF") )
    return rsInfDouble;
  else if( *this == rsString("-INF") )
    return -rsInfDouble;
  else if( *this == rsString("NaN") )
    return rsQuietNaNDouble;
    //return NAN;

  char* cString = getAsZeroTerminatedString();

  double result = strtod(cString, NULL);  
    // atof occasionaly spits out twice the actual value for denormals with MinGW

  delete[] cString;
  return result;
}

float rsString::asFloat() const
{
  return (float) asDouble();
}

rsRange<int> rsString::findRangeEnclosedBy(const rsString& startDelimiter, 
  const rsString& endDelimiter, bool includeDelimitersInRange, int searchStart) const
{
  int startIndex = findFirstOccurrenceOf(startDelimiter, searchStart);
  int endIndex   = findFirstOccurrenceOf(endDelimiter,  startIndex+1);

  if( includeDelimitersInRange )
    endIndex += endDelimiter.getLength();
  else
    startIndex += startDelimiter.getLength();

  return rsRange<int>(startIndex, endIndex-1);
}

/*
bool String::containsNonPrintableCharacters() const
{
  for(int i = 0; i < getNumElements(); i++)
  {
    if( isprint(elements[i]) == false )
      return true;
  }
  return false;
}
*/

// manipulations:

void rsString::padToLength(int desiredLength, const char padCharacter)
{
  if( numUsed < desiredLength )
  {
    ensureAllocatedSize(desiredLength);
    for(int i = numUsed; i < desiredLength; i++)
      elements[i] = padCharacter;
    numUsed = desiredLength;
  }
}

void rsString::replaceCharacter(const char characterToReplace, const char replacementCharacter)
{
  for(int i = 0; i < numUsed; i++)
  {
    if( elements[i] == characterToReplace )
      elements[i] = replacementCharacter;
  }
}

void rsString::removeRepeatedCharacters(const char characterThatShouldNotRepeat)
{
  int read        = 0;
  int write       = 0;
  int numDropped  = 0;
  int repeatCount = 0;
  while( read < numUsed )
  {
    elements[write] = elements[read];

    if( elements[read] == characterThatShouldNotRepeat )
      repeatCount++;
    else
      repeatCount = 0;

    read++;

    if( repeatCount > 1 )
      numDropped++;
    else
      write++;
  }
  numUsed -= numDropped;
  shrinkMemoryOccupationIfAppropriate();
}


void rsString::removeLeadingCharacters(const char characterToRemoveAtBegin)
{
  int start = 0;
  while( start < numUsed && elements[start] == characterToRemoveAtBegin )
    start++;
  for(int i = start; i < numUsed; i++)
    elements[i-start] = elements[i];
  numUsed -= start;
  shrinkMemoryOccupationIfAppropriate();
}

void rsString::removeTrailingCharacters(const char characterToRemoveAtEnd)
{
  int end = numUsed-1;
  while( end >= 0 && elements[end] == characterToRemoveAtEnd )
    end--;
  numUsed = end+1;
  shrinkMemoryOccupationIfAppropriate();
}

// substring extraction:

/*
String String::getSubStringUpToCharacter(const char endCharacter, const int startIndex) const
{
  int endIndex = findFirstOccurrenceOf(endCharacter, startIndex);
  String result = getSubString(startIndex, endIndex);
  return getSubString(startIndex, endIndex);
}
*/

rsString rsString::getSubStringEnclosedBy(const rsString& startDelimiter, 
  const rsString& endDelimiter, bool includeDelimitersInrsRange, int searchStart) const
{
  rsRange<int> range = findRangeEnclosedBy(startDelimiter, endDelimiter, 
                                           includeDelimitersInrsRange, searchStart);
  return getSubString(range.getMin(), range.getMax());
}

// misc:

void rsString::printToStandardOutput() const
{
  char* cString = getAsZeroTerminatedString();
  printf("%s", cString);
  delete[] cString;
}
   
void rsString::copyToZeroTerminatedString(char *dst, int maxLengthIncludingZero) const
{
  int numToWrite = rsMin(getLength()+1, maxLengthIncludingZero);
  rsCopyBuffer(elements, dst, numToWrite-1);
  dst[numToWrite] = '\0';
}

// string-related non-member functions:

rsString rsDoubleToString2(double x)
{
  char cString[64];
  sprintf(cString, "%.2f", x);
  return rsString(cString);
}

}