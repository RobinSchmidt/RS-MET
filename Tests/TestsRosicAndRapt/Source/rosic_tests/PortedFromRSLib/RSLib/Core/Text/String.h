#ifndef RS_STRING_H
#define RS_STRING_H

namespace RSLib
{

  /**

  This is a class for representing strings. It is implemented as specialization of rs::rsArrayTools for 
  type char.

  */

  class RSLib_API rsString : public rsArrayTools<char>
  {

  public:

    /** \name Construction/Destruction */

    /** Standard constructor. Creates an empty string. */
    rsString();

    /** Constructor. Creates a string from a zero-terminated c-string. Note that when you pass the 
    NULL macro, actually the constructor String(int) will be called with 0 as argument. */
    rsString(const char *initialString);

    /** Constructor. Creates a string from a std::string object. */
    rsString(const std::string &initialString);

    /** Constructor. Creates a string from an array of characters. Although String is subclass of 
    rsArrayTools<char>, we still need to implement this. */
    rsString(const rsArrayTools<char> &stringAsCharsArray);

    /** Copy constructor. Creates a (deep) copy of another string. */
    rsString(const rsString &other);

    /** Constructor. Creates a string from an integer number. */
    rsString(const int intValue);

    /** Constructor. Creates a string from a single precision floating point number. */
    rsString(const float floatValue);

    /** Constructor. Creates a string from a double precision floating point number.
    \todo: maybe include a flag 'forceExponentialNotation' with default value false. */
    rsString(const double doubleValue);

    /** Destructor. */
    ~rsString();


    /** \name Setup */

    /** Sets this string up from a zero terminated c-string. */
    void setFromZeroTerminatedString(const char* cString);

    /** Sets this string up from a std::string object. */
    void setFromStdString(const std::string &stdString);

    /** Sets this string from an integer number. */
    void setFromIntValue(const int intValue);

    /** Sets this string from a double precision floating point number. */
    void setFromDoubleValue(const double doubleValue);


    /** \name Inquiry */

    /** Returns the a copy of string as raw zero terminated c-string. The caller is responsible 
    for eventually deleting it. */
    char* getAsZeroTerminatedString() const;

    /** Converts this rsString object into a std::string object. */
    std::string getAsStdString() const;

    /** Returns the string as a integer number if conversion is possible, 0 otherwise. */
    int asInt() const; // { return atoi(cString); }

    /** Returns the string as double precision floating point number if conversion is possible, 
    0.0 otherwise. */
    double asDouble() const;

    /** Returns the string as single precision floating point number if conversion is possible, 
    0.0 otherwise. */
    float asFloat() const;

    /** Returns the number of characters in the string (excluding the terminating zero). */
    int getLength() const { return getNumElements(); }

    /** Finds the range (start-index to end-index) of a substring that is enclosed by the given 
    start- and end-strings. The returned range optionally includes the space occupied by the 
    delimiters. If the enclosing delimiters occur multiple times, it will return the range of the 
    first occurence. If one of the delimiting strings does not occur, it will return a range with 
    -1 for the respective value of the range. */
    rsRange<int> findRangeEnclosedBy(const rsString& startDelimiter, const rsString& endDelimiter, 
      bool includeDelimitersInRange, int searchStart = 0) const;

    /** Finds all occurences of 'stringToFind' inside this string within 'searchStart' and 
    'searchEnd'. If 'searchEnd' is -1 (default) the string will be scanned until its end. */
    rsArrayTools<int> findAllOccurrencesOf(const rsString& stringToFind, const int searchStart = 0,                          
                                      int searchEnd = -1) const
    {
      if( searchEnd < 0 )
        searchEnd = numUsed-1;
      int searchLength = searchEnd-searchStart+1;
      std::vector<int> matches = RSLib::rsFindAllOccurencesOf(&elements[searchStart],              
        searchLength, stringToFind.elements, stringToFind.getLength());
      return rsArrayTools<int>(matches);
    }

    /** Returns true, if this string contains non-printable characters according to the isprint 
    C-function.  
    \todo: get rid of that legacy c-stuff, treat all characters uniformly (including '0'). */
    //bool containsNonPrintableCharacters() const;


    /** \name Manipulations */

    /** Pads the string to the desired length by appending an appropriate number of padding 
    characters (typically whitespaces). If the length of the string is already >= desiredLength,
    it will do nothing. 

    \todo: write a function: fixLength (or something) that pads or truncates the string to the 
    desired length */
    void padToLength(int desiredLength, const char padCharacter = ' ');

    /** Replaces all occurrences of 'characterToReplace' with 'replacementCharacter'. */
    void replaceCharacter(const char characterToReplace, const char replacementCharacter);

    /** Reduces all repeated occurrences of some character to just one single occurrence. */
    void removeRepeatedCharacters(const char characterThatShouldNotRepeat);

    /** Removes runs of the given character at the beginning of the string. */
    void removeLeadingCharacters(const char characterToRemoveAtBegin);

    /** Removes runs of the given character at the end of the string. */
    void removeTrailingCharacters(const char characterToRemoveAtEnd);

    /** Removes all occurences of characters which are present in the passed string. */
    void removeCharacters(const rsString& charactersToRemove) 
    { 
      removeMatchingElements(charactersToRemove); 
    }

    /** Removes all occurences of characters which are not present in the passed string. */
    void removeAllCharactersBut(const rsString& charactersToKeep) 
    { 
      keepOnlyMatchingElements(charactersToKeep); 
    }

    /** Repeats this string the given number of times. You may also pass zero in which case the 
    result will be an empty string. */
    void repeat(const int numberOfTimes) { rsArrayTools<char>::repeat(numberOfTimes); }


    /** \name Substring Extraction */

    /** Returns the substring from startIndex to endIndex (inclusive). */
    rsString getSubString(const int startIndex, const int endIndex) const
    { return rsArrayTools<char>::getSubArray(startIndex, endIndex); }

    /** Returns the substring starting at 'startIndex' and ending just before the first occurrence 
    of 'endString'. */
    //String getSubStringUpTo(const String& endString, const int startIndex) const;

    /** Returns the substring starting with 'startDelimiter' and ending with 'endDelimiter' and 
    optionally includes the delimiters themselves in the returned string. The search may start at 
    an arbitrary position inside the string as givne by searchStart. */
    rsString getSubStringEnclosedBy(const rsString& startDelimiter, const rsString& endDelimiter, 
      bool includeDelimitersInrsRange, int searchStart = 0) const;


    /** \name Operators */

    /** Checks whether the left string is lexicographically less than the right string. */
    bool operator<(const rsString& s2) const
    { return compareCharacterArrays(elements, numUsed, s2.elements, s2.numUsed) == -1; }

    /** Checks whether the left string is lexicographically greater than the right string. */
    bool operator>(const rsString& s2) const
    { return compareCharacterArrays(elements, numUsed, s2.elements, s2.numUsed) == +1; }

    /** Checks whether the left string is lexicographically less than or equal to the right 
    string. */
    bool operator<=(const rsString& s2) const
    {
      int r = compareCharacterArrays(elements, numUsed, s2.elements, s2.numUsed);
      return ( r == -1 || r == 0 );
    }

    /** Checks whether the left string is lexicographically greater than or equal to the right 
    string. */
    bool operator>=(const rsString& s2) const
    {
      int r = compareCharacterArrays(elements, numUsed, s2.elements, s2.numUsed);
      return ( r == +1 || r == 0 );
    }

    /** Adds another string to this string by appending it and returns the result. */
    rsString& operator+=(const rsString &other)
    {
      appendArray(other);
      return *this;
    }

    /** Adds two strings by concatenating them and returns the result. */
    rsString operator+(const rsString &other)
    {
      rsString result = *this;
      result.appendArray(other);
      return result;
    }


    /** \name Misc */

    /** Writes the string to the standard output which is typically the console. */
    void printToStandardOutput() const;

    /** Copies the content of this string into the character-array "dst" and appends a terminating
    zero. The second parameter determines, how many characters (including the terminating zero) may
    be written into the "dst" array. */
    void copyToZeroTerminatedString(char *dst, int maxLengthIncludingZero) const;


    /** \name Static Member Functions */

    /** Returns the number of characters that are required to represent the given integer number 
    as string. */
    //static int numberOfRequiredCharacters(int number);

    /** This is an empty string. */
    static const rsString empty;

    /** Lexicographically compares the two char-arrays and returns -1 if left < right, 0 if 
    left == right and +1 if left > right using the compareCharacters functions on the individual 
    elements. If one array is the first section of the other, the shorter one will be considered to
    come before the longer one. */
    static int compareCharacterArrays(char *left, int leftLength, char *right, int rightLength);

    /** Returns -1 if left < right, 0 if left == right and +1 if left > right according to 
    alphabetical order. Uppercase letters are considered to come before lowercase letters and 
    numbers are considered to come before letters. For all other characters, the order is 
    determined by the ASCII standard (some may fall before numbers, some between numbers and 
    characters, some between lower-/ and uppercase letters and some after lowercase letters). */
    static int compareCharacters(char left, char right);

  };

  // + operator for the case when the left hand operand is a plain character array 
  RS_INLINE rsString operator+(const char *s1, const rsString &s2)
  {
    return rsString(s1) + s2;
  }


  // string-related non-member functions:

  /** A typedef for a pointer to a function that takes a double as input and returns an 
  rsString. */
  typedef rsString(*DoubleToStringFunctionPointer) (double);

  /** Converts a double to a string with 2 decimal digits after the point. */
  rsString RSLib_API rsDoubleToString2(double x);

}

#endif
