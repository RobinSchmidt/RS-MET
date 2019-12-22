#ifndef RS_KEYVALUESTRINGARRAY_H
#define RS_KEYVALUESTRINGARRAY_H

namespace RSLib
{

  /**

  A class for representing an array of key/value pairs.
  @see rsKeyValueStringPair

  */

  class RSLib_API rsKeyValueStringArray
  {

  public:

    /** \name Setup */

    /** Adds a new string value with the given key/value to the array or updates an existing value 
    with the given key. */
    void setStringValue(const rsString& key, const rsString& value);

    /** Adds a new real number value with the given key/value to the array or updates an existing 
    value with the given key. */
    void setRealNumberValue(const rsString& key, const double value);

    /** Adds a new integer number value with the given key/value to the array or updates an 
    existing value with the given key. */
    void setIntegerValue(const rsString& key, const int value);

    /** Adds a new boolean (true/false) value with the given key/value to the array or updates an 
    existing value with the given key. */
    void setBooleanValue(const rsString& key, const bool value);

    /** Deletes all key/value pairs. */
    void clear() { keyValuePairs.clear(); }


    /** \name Inquiry */

    /** Returns the string value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    rsString getStringValue(const rsString& key, 
      const rsString& defaultReturnValue = rsString::empty);

    /** Returns the real number value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    double getRealNumberValue(const rsString& key, const double defaultReturnValue = 0.0);

    /** Returns the integer number value at the given key (if key present, otherwise returns the 
    'defaultReturnValue'). */
    int getIntegerValue(const rsString& key, const int defaultReturnValue = 0);

    /** Returns the boolean (true/false) value at the given key (if key present, otherwise returns 
    the 'defaultReturnValue'). */
    bool getBooleanValue(const rsString& key, const bool defaultReturnValue = false);

    /** Returns the array-index of the given key, -1 if the key does not exist. */
    int getIndexOfKey(const rsString& key);

    /** Returns true when the array contains a pair with the given key, false otherwise. */
    bool hasPairWithKey(const rsString& key) { return getIndexOfKey(key) != -1; }

    /** Returns the number of stored key/value pairs. */
    int getNumKeyValueStringPairs() const { return keyValuePairs.getNumElements(); }

    /** Returns the key/value pairs that is stored a the given index. The caller must make sure 
    that the index is valid - no check is made here. */
    rsKeyValueStringPair getKeyValueStringPair(int index) { return keyValuePairs[index]; }


    /** \name Operators */

    /** Compares two KeyValueStringArrays for equality. */
    bool operator==(const rsKeyValueStringArray& other) const
    {
      return keyValuePairs == other.keyValuePairs;
    }

    /** Compares two KeyValueStringArrays for inequality. */
    bool operator!=(const rsKeyValueStringArray& other) const
    {
      return !(*this == other);
    }

  protected:

    /** \name Data */

    rsArrayTools<rsKeyValueStringPair> keyValuePairs;

  };

}

#endif
