#ifndef RS_KEYVALUESTRINGPAIR_H
#define RS_KEYVALUESTRINGPAIR_H

namespace RSLib
{

  /**

  A class for representing key/value pairs where each key and value is represented by a string.

  */

  class RSLib_API rsKeyValueStringPair
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. You can pass the strings for the key and the value here or leave them empty 
    and set them later via the setters. */
    rsKeyValueStringPair(const rsString& key = rsString::empty, 
      const rsString& value = rsString::empty)
    {
      this->key   = key;
      this->value = value;
    }

    /** Destructor. */
    //~rsKeyValueStringPair();

    /** \name Setup */

    /** Sets the 'key'-string. */
    void setKey(const rsString& newKey) { key = newKey; }

    /** Sets the 'value'-string. */
    void setStringValue(const rsString& newValue) { value = newValue; }

    /** Sets the 'value'-string from a real number (of type double). */
    void setRealNumberValue(const double newValue) { value = rsString(newValue); }

    /** Sets the 'value'-string from an integer number. */
    void setIntegerValue(const int newValue) { value = rsString(newValue); }

    /** Sets the 'value'-string from a boolean (true/false) value. */
    void setBooleanValue(const bool newValue) { value = rsString((int)newValue); }


    /** \name Inquiry */

    /** Returns the 'key'-string. */
    rsString getKey() const { return key; }

    /** Returns the stored value as String. */
    rsString getStringValue() const { return value; }

    /** Returns the stored value as real number (of type double). */
    double getRealNumberValue() const { return value.asDouble(); }

    /** Returns the stored value as integer number. */
    int getIntegerValue() const { return value.asInt(); }

    /** Returns the stored value as boolean (true/false) value. */
    bool getBooleanValue() const { return (value.asInt() != 0); }

    /** Returns true when the passed key to check for is equal to this object's key, false 
    otherwise. */
    bool hasKey(const rsString& keyToCheckFor) const { return key == keyToCheckFor; }


    /** \name Operators */

    /** Compares two KeyValueStringPairs for equality. */
    bool operator==(const rsKeyValueStringPair& other) const
    {
      if( key == other.key && value == other.value )
        return true;
      else
        return false;
    }

    /** Compares two KeyValueStringPairs for inequality. */
    bool operator!=(const rsKeyValueStringPair& other) const
    {
      return !(*this == other);
    }

  protected:


    /** \name Data */

    rsString key, value;

  };

}

#endif
