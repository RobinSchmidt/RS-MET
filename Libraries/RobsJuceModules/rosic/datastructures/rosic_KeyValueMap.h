#ifndef rosic_KeyValueMap_h
#define rosic_KeyValueMap_h

//#include <stddef.h>  // for NULL macro
//#include <vector>
//#include <algorithm>
//#include <functional>

namespace rosic
{

/* Structure to hold a key/value pair where key and value can be of arbitrary type. */
template<class KeyType, class ValueType>
struct KeyValuePair
{
  KeyType   key;
  ValueType value;
};


/* Less-than comparison between two pointers to key/value pairs which returns true when the key in the first pointer is less than the
key in the second pointer. */
template<class KeyType, class ValueType>
bool keyValuePairPointerLessByKey(
  KeyValuePair<KeyType, ValueType>* pointer1, KeyValuePair<KeyType, ValueType>* pointer2)
{
  if(pointer1->key < pointer2->key)
    return true;
  else
    return false;
};

/* Less-than comparison between two pointers to key/value pairs which returns true when the value in the first pointer is less than the
value in the second pointer. */
template<class KeyType, class ValueType>
bool keyValuePairPointerLessByValue(
  KeyValuePair<KeyType, ValueType>* pointer1, KeyValuePair<KeyType, ValueType>* pointer2)
{
  if(pointer1->value < pointer2->value)
    return true;
  else
    return false;
};


/* Wraps the function keyValuePairPointerLessByKey() into a functor class, as required by the STL. */
/*
template<class KeyType, class ValueType>
struct KeyValuePairPointerLessByKey
  : std::binary_function<KeyValuePair<KeyType, ValueType>*,
                         KeyValuePair<KeyType, ValueType>*, bool>
{
  bool operator() (KeyValuePair<KeyType, ValueType>* pointer1, KeyValuePair<KeyType, ValueType>* pointer2) const
  {
    return keyValuePairPointerLessByKey(pointer1, pointer2);
  }
};
*/
// doesn't compile in c++17 - see here:
// https://stackoverflow.com/questions/33114656/replacement-for-binary-function

/* Wraps the function keyValuePairPointerLessByValue() into a functor class, as required by the STL. */
/*
template<class KeyType, class ValueType>
struct KeyValuePairPointerLessByValue
  : std::binary_function<KeyValuePair<KeyType, ValueType>*,
                         KeyValuePair<KeyType, ValueType>*, bool>
{
  bool operator() (KeyValuePair<KeyType, ValueType>* pointer1, KeyValuePair<KeyType, ValueType>* pointer2) const
  {
    return keyValuePairPointerLessByValue(pointer1, pointer2);
  }
};
*/




/**

This class implements a mapping between keys and values where both, the key and the value can be of arbitrary type. The class allows fast
translation between a key and the corresponding value and vice versa by internally maintaining two sorted arrays where one is sorted by
the key and the other is sorted by the value. This allows key/value conversion in either direction with complexitiy O(log(N)) where N is
the number of key/value pairs. This could be used, for example, for a telephone number register where one could search for the name that
corresponds to a given telephone number or vice versa with O(log(N)) complexity.

\todo: maybe rename to TwoWayAssociator

*/

template<class KeyType, class ValueType>
class KeyValueMap
{

public:

  //--------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Destructor. */
  ~KeyValueMap()
  {
    for(unsigned int i=0; i<entriesSortedByKey.size(); i++)
      delete entriesSortedByKey[i];
  }

  //--------------------------------------------------------------------------------------------------------------------------------------
  // adding elements:

  /** Inserts a key/value pair into the map. */
  void insertKeyValuePair(const KeyType key, const ValueType value)
  {
    KeyValuePair<KeyType, ValueType>* newEntry = new KeyValuePair<KeyType, ValueType>;
    newEntry->key   = key;
    newEntry->value = value;
    typename std::vector<KeyValuePair<KeyType, ValueType>*>::iterator iter;

    /* old, pre C++11:
    // insert into the sorted-by-key array at the right position:
    KeyValuePairPointerLessByKey<KeyType, ValueType> compareByKey;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), newEntry, compareByKey);
    entriesSortedByKey.insert(iter, newEntry);

    // insert into the sorted-by-value array at the right position:
    KeyValuePairPointerLessByValue<KeyType, ValueType> compareByValue;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), newEntry, compareByValue);
    entriesSortedByValue.insert(iter, newEntry);
    */

    // new, C++11:
        // insert into the sorted-by-key array at the right position:
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(),
                       newEntry, keyValuePairPointerLessByKey);
    entriesSortedByKey.insert(iter, newEntry);

    // insert into the sorted-by-value array at the right position:
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(),
                       newEntry, keyValuePairPointerLessByValue);
    entriesSortedByValue.insert(iter, newEntry);


    //keyValuePairPointerLessByKey
  }

  // todo: write mergeKeyValueMaps, ...

  //--------------------------------------------------------------------------------------------------------------------------------------
  // removing elements:

  void removePairsWithKey(const KeyType key)
  {

    //...

  }

  //--------------------------------------------------------------------------------------------------------------------------------------
  // inquiry:

  /** If a pair with the given key is present in the map, the function returns the corresponding value and assigns the reference
  parameter "wasFound" to true. If the key is not present, it returns a default value determined by ValueType's default constructor and
  assigns the reference parameter "wasFound" to false. */
  ValueType getValueForKey(KeyType key, bool& wasFound)
  {
    KeyValuePair<KeyType, ValueType>* dummyEntry = new KeyValuePair<KeyType, ValueType>;
    dummyEntry->key = key;

    /*
    // pre C++11:
    KeyValuePairPointerLessByKey<KeyType, ValueType> compareByKey;
    typename std::vector<KeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), dummyEntry, compareByKey);
    */

    // C++11:
    typename std::vector<KeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(),
                       dummyEntry, keyValuePairPointerLessByKey);



    if(iter != entriesSortedByKey.end() && (*iter)->key == key)
    {
      wasFound = true;
      delete dummyEntry;
      return (*iter)->value;
    }
    else
    {
      wasFound = false;
      delete dummyEntry;
      return ValueType();
    }
  }

  /** If a pair with the given value is present in the map, the function returns the corresponding key and assigns the reference
  parameter "wasFound" to true. If the value is not present, it returns a default key determined by KeyType's default constructor and
  assigns the reference parameter "wasFound" to false. */
  KeyType getKeyForValue(ValueType value, bool& wasFound)
  {
    KeyValuePair<KeyType, ValueType>* dummyEntry = new KeyValuePair<KeyType, ValueType>;
    dummyEntry->value = value;

    /*
    // pre C++11:
    KeyValuePairPointerLessByValue<KeyType, ValueType> compareByValue;
    typename std::vector<KeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), dummyEntry, compareByValue);
    */

    // C++11:
    typename std::vector<KeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(),
                       dummyEntry, keyValuePairPointerLessByValue);

    if(iter != entriesSortedByValue.end() && (*iter)->value == value)
    {
      wasFound = true;
      delete dummyEntry;
      return (*iter)->key;
    }
    else
    {
      wasFound = false;
      delete dummyEntry;
      KeyType result = KeyType();
      return KeyType();
    }
  }

  /** Returns the number of entries (key/value pairs) in the map. */
  int getNumEntries() const { return entriesSortedByKey.size(); }

  /** Returns true when the map is empty, false otherwise. */
  bool isEmpty() const { return getNumEntries() == 0; }

  //--------------------------------------------------------------------------------------------------------------------------------------
  // operators:

  // \todo see rsArrayTools and implement some operators...

  //=====================================================================================================================================

protected:

  std::vector<KeyValuePair<KeyType, ValueType>*> entriesSortedByKey;
  std::vector<KeyValuePair<KeyType, ValueType>*> entriesSortedByValue;
    // \todo maybe switch to rsArrayTools<KeyValuePair*> here later

  // friend decalarations to facilitate testing:
  friend bool testKeyValueMapInsert(std::string& reportString);

};

}

#endif
