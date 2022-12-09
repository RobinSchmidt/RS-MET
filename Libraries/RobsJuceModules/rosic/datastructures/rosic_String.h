#ifndef rosic_String_h
#define rosic_String_h

namespace rosic
{

/** This is a class for representing strings. It is implemented as wrapper around plain old
C strings. It contains the C-string including the terminating zero but also stores the length in
a dedicated member variable (so it doesn't have to be found using strlen each time it is needed).

\todo: get rid of the old-schoolish C-string stuff and manage the string as a plain array of
characters, maintain usedLength and allocatedLength as members, don't use the C-string
functions anymore, provide conversion functions from and to c-strings */

class rsString
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Standard constructor. Creates an empty string. */
  rsString();

  /** Constructor. Creates a string from a zero-terminated c-string. Note that when you pass the
  NULL macro, actually the constructor rsString(int) will be called with 0 as argument. */
  rsString(const char *initialString);

  /** Creates an rsString from a std::string. */
  rsString(const std::string& initialString);

  /** Copy constructor. Creates a (deep) copy of another string. */
  rsString(const rsString &other);

  /** Constructor. Creates a string from an integer number. */
  rsString(const int intValue);

  /** Constructor. Creates a string from an integer number. */
  rsString(const unsigned int uintValue);

  /** Constructor. Creates a string from a single precision floating point number. */
  rsString(const float floatValue);

  /** Constructor. Creates a string from a double precision floating point number.
  \todo: maybe include a flag 'forceExponentialNotation' with default value false. */
  rsString(const double doubleValue);

  /** Destructor. */
  ~rsString();

  //---------------------------------------------------------------------------------------------
  // setup:

  /** Reserves enough memory to store the given number of characters (excluding the terminating
  zero). If already enough memory is reserved, it will do nothing, otherwise it will re-allocate
  memory of the desired size and copy the current string into it - so this function can be used
  to make up space after the actually used memory. This can be advantageous to optimize
  concatenations. */
  void reserveSize(int numCharactersToReserve);

  /** Pads the string with the passed padding-character such that it has the desired length. In
  cases where the string has already the desired length or is even longer, the function does
  nothing. */
  void padToLength(int desiredLength, char paddingCharacter = ' ');

  /** Similar to padToLength, but appends the padding at the front of the string. */
  void prePadToLength(int desiredLength, char paddingCharacter = ' ');

  //---------------------------------------------------------------------------------------------
  // inquiry:

  /** Returns the string as raw zero terminated c-string - this is just a pointer to the
  internal character array which is maintained here, so DON'T attempt to DELETE it. */
  const char* getRawString() const { return cString; }

  /** Returns the number of characters in the string (excluding the terminating zero). */
  int getLength() const { return length; }

  /** Returns true when this string contains non-printable characters according to the isprint
  C-function.
  \todo: get rid of that legacy c-stuff, treat all characters uniformly (including '0'). */
  bool containsNonPrintableCharacters() const;

  /** Returns true when "this" string is empty, false otherwise. */
  bool isEmpty() const { return length == 0; }

  //---------------------------------------------------------------------------------------------
  // operators:

  /** Assigns one string with another one. */
  rsString& operator=(const rsString& s2)
  {
    if(reservedSize != s2.reservedSize)
    {
      delete[] cString;
      reservedSize = s2.reservedSize;
      cString      = new char[reservedSize];
    }
    length = s2.length;
    strcpy(cString, s2.cString); // use concatenateBuffers
    return *this;
  }

  /** Compares two strings of equality. Two strings are considered equal if they have the same
  length and match character by character. The reserved memory size is irrelevant for this
  comparison.  */
  bool operator==(const rsString& s2) const
  {
    if(length != s2.length)
      return false;
    else
    {
      for(int i=0; i<length; i++)
      {
        if(cString[i] != s2.cString[i])
          return false;
      }
    }
    return true;
  }

  /** Compares two strings of inequality. */
  bool operator!=(const rsString& s2) const
  {
    return !(*this == s2);
  }

  /** Checks whether the left operand is less than the right operand (using strcmp). */
  bool operator<(const rsString& s2) const
  {
    int result = strcmp(cString, s2.cString); // use compareBuffers
    return (result<0);
  }

  /** Checks whether the left operand is greater than the right operand (using strcmp). */
  bool operator>(const rsString& s2) const
  {
    int result = strcmp(cString, s2.cString); // use compareBuffers
    return (result>0);
  }

  /** Checks whether the left operand is less than or equal to the right operand (using
  strcmp). */
  bool operator<=(const rsString& s2) const
  {
    int result = strcmp(cString, s2.cString);
    return (result<=0);
  }

  /** Checks whether the left operand is greater than or equal to the right operand (using
  strcmp). */
  bool operator>=(const rsString& s2) const
  {
    int result = strcmp(cString, s2.cString);
    return (result>=0);
  }

  /** Adds another string to this string by appending it and returns the result. */
  rsString& operator+=(const rsString &s2)
  {
    reserveSize(length + s2.length);
    strcat(cString, s2.cString);
    length += s2.length;
    return *this;
  }

  /** Adds two strings by concatenating them and returns the result. */
  rsString operator+(const rsString &s2)
  {
    rsString result;
    result.reserveSize(length+s2.length);
    strcpy(result.cString, cString);
    strcat(result.cString, s2.cString);
    result.length = length+s2.length;
    return result;
  }

  //---------------------------------------------------------------------------------------------
  // conversions to numeric values:

  /** Returns the string as std::string. */
  std::string asStdString() const { return std::string(cString); }

  /** Returns the string as a integer number if conversion is possible, 0 otherwise. */
  int asInt() const { return atoi(cString); }

  /** Returns the string as double precision floating point number if conversion is possible,
  0.0 otherwise. */
  double asDouble() const;

  /** Returns the string as single precision floating point number if conversion is possible,
  0.0 otherwise. */
  float asFloat() const { return (float)atof(cString); }


  /** Returns a string where all chareacters between start and end (both inclusive) have been
  converted to lowercase. */
  rsString toLowerCase(int start, int end) const;

  /** Returns a string where all chareacters between start and end (both inclusive) have been
  converted to uppecase. */
  rsString toUpperCase(int start, int end) const;

 //---------------------------------------------------------------------------------------------
  // others:

  /** Writes the string to the standard output which is typically the console. */
  void printToStandardOutput() const;

  /** Copies the (assumed to be zero-terminated) c-string represented by the passed buffer into
  this string. */
  void readFromBuffer(char *sourceBuffer);

  /** Writes the string into the passed buffer. You should also pass the maximum number of
  characters to write including the terminating zero - if the actual string is shorter than
  that, the extra characters will be left as they are. If the string is longer, then you may
  miss a terminating zero in your buffer. */
  void writeIntoBuffer(char *targetBuffer, int maxNumCharactersToWrite) const;

  //---------------------------------------------------------------------------------------------
  // static functions and data members:

  /** Returns the number of characters that are required to represent the given integer number
  as string. */
  static int numberOfRequiredCharacters(int number);

  /** This is an empty string. - causes memory leaks? */
  //static const rsString empty;

  /** Lexicographically compares the two char-arrays and returns -1 if left < right,
  0 if left == right and +1 if left > right using the compareCharacters functions on the
  individual elements. If one array is the first section of the other, the shorter one will be
  considered to come before the longer one.

  \todo: test this function

  */
  static int compareCharacterArrays(char *left, int leftLength, char *right, int rightLength);

  /** Returns -1 if left < right, 0 if left == right and +1 if left > right according to
  alphabetical order. Uppercase letters are considered to come before lowercase letters and
  numbers are considered to come before letters. For all other characters, the order is
  determined by the ASCII standard (some may fall before numbers, some between numbers and
  characters, some between lower-/ and uppercase letters and some after lowercase letters). */
  static int compareCharacters(char left, char right);

  /** If the argument is a lowercase letter, this function returns its uppercase version.
  Otherwise, the character itself is returned. */
  static char toUpperCase(char c);


  static bool isUpperCaseLetter(const char c) { return c >= 65 && c <= 90; }

  static bool isLowerCaseLetter(const char c) { return c >= 97 && c <= 122; }

  static bool isLetter(const char c) { return isUpperCaseLetter(c) || isLowerCaseLetter(c); }

  static bool isDigit(const char c) { return c >= 48 && c <= 57; }

  static bool isLetterOrDigit(const char c) { return isLetter(c) || isDigit(c); }

  /** Creates a string from an integer number with some minimum number of characters, prepending
  whitespaces, if necessary. This is useful for alignment of strings representing numbers. It
  may also optionally include a plus sign for positive numbers to make them look more similar
  to negative numbers. \todo: test this function thouroughly */
  static rsString fromIntWithLeadingSpaces(int value, int minNumCharacters,
    bool prependPlusForPositiveNumbers = false);

/** Creates a string from a double precision float number. */
  static rsString fromDouble(double value);

  /** Creates a string with the given number of whitespace characters. */
  static rsString createWhiteSpace(int length);

  /** Creates span of repititions of the given character. */
  static rsString createSpanOfCharacters(char character, int length);

  //=============================================================================================

protected:

  /** Allocates memory area of size numCharactersToAllocateExcludingZero+1 and associates it with
  our cString member, unless the allocated memory is already exactly of the right size in which
  case nothing is done. In any event, upon return, our memory will have the right size and our
  cString pointer is associated with it. */
  void allocateMemory(int numCharactersToAllocateExcludingZero);



  char *cString;
  int  length;        // number of characters excluding the terminating zero
  int  reservedSize;  // size of the allocated memory (always greater than 'length' but not
                      // necessarily exactly length+1)

  //=============================================================================================

private:

  void initFromDoubleValue(double doubleValue);
    // used only internally by the constructors rsString(double) and rsString(float)

  /** Given the c-string "s" that represents a double precision floating point number, this
  function removes garbage characters as in .... and returns the new length of the c-string
  (excluding terminating 0) */
  int removeGarbageFromDoubleString(char *s, int length);

};

//=================================================================================================

// functions for std::string:
void removeChar(std::string& str, const char chr);
std::vector<std::string> tokenize(const std::string& str, const char splitChar);

// start: 
//  input:    index from which we start searching for the next token
//  output:   index where the next token actually starts
// length:    length of the token (output only)
void rsFindToken(const std::string& str, const std::string& seperators, int* start, int* length);
// maybe move as static function into a class rsStdStringHelpers

std::string rsGetToken(const std::string& str, size_t startIndex, const std::string& sep);


/** Replaces all occurences of oldText with newText within the subject. */
void rsReplace(std::string& subject, const std::string& oldText, const std::string& newText);
// needs unit tests



} // end namespace rosic

#endif
