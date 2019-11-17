#ifndef RS_KEYVALUEMAP_H
#define RS_KEYVALUEMAP_H

namespace RSLib
{

  /* Structure to hold a key/value pair where key and value can be of arbitrary type. */
template<class KeyType, class ValueType>
struct rsKeyValuePair
{
  KeyType   key;
  ValueType value;
};

/* Less-than comparison between two pointers to key/value pairs which returns true when the key
in the first pointer is less than the key in the second pointer. */
template<class KeyType, class ValueType>
bool keyValuePairPointerLessByKey(rsKeyValuePair<KeyType, ValueType>* pointer1,
  rsKeyValuePair<KeyType, ValueType>* pointer2)
{
  if(pointer1->key < pointer2->key)
    return true;
  else
    return false;
};

/* Wraps the function keyValuePairPointerLessByKey() into a functor class, as required by the STL. */
/*
template<class KeyType, class ValueType>
struct rsKeyValuePairPointerLessByKey
  : std::binary_function<rsKeyValuePair<KeyType, ValueType>*,
  rsKeyValuePair<KeyType, ValueType>*, bool>
{
  bool operator() (rsKeyValuePair<KeyType, ValueType>* pointer1,
    rsKeyValuePair<KeyType, ValueType>* pointer2) const
  {
    return keyValuePairPointerLessByKey(pointer1, pointer2);
  }
};
*/
// doesn't compile in c++17 - see here:
// https://stackoverflow.com/questions/33114656/replacement-for-binary-function

  /* Less-than comparison between two pointers to key/value pairs which returns true when the
  value in the first pointer is less than the value in the second pointer. */
template<class KeyType, class ValueType>
bool keyValuePairPointerLessByValue(rsKeyValuePair<KeyType, ValueType>* pointer1,
  rsKeyValuePair<KeyType, ValueType>* pointer2)
{
  if(pointer1->value < pointer2->value)
    return true;
  else
    return false;
};

/* Wraps the function keyValuePairPointerLessByValue() into a functor class, as required by
the STL. */
/*
template<class KeyType, class ValueType>
struct rsKeyValuePairPointerLessByValue
  : std::binary_function<rsKeyValuePair<KeyType, ValueType>*,
                         rsKeyValuePair<KeyType, ValueType>*, bool>
{
  bool operator() (rsKeyValuePair<KeyType, ValueType>* pointer1,
    rsKeyValuePair<KeyType, ValueType>* pointer2) const
  {
    return keyValuePairPointerLessByValue(pointer1, pointer2);
  }
};
*/

//===============================================================================================

/**

This class implements a mapping between keys and values where both, the key and the value, can be
of arbitrary type. The class allows fast translation between a key and the corresponding value
and vice versa by internally maintaining two sorted arrays where one is sorted by the key and the
other is sorted by the value. This allows key/value conversion in either direction with
complexitiy O(log(N)) where N is the number of key/value pairs. This could be used, for example,
for a telephone number register where one could search for the name that corresponds to a given
telephone number or vice versa with O(log(N)) complexity.

\todo rename to TwoWayAssociator or BidirectionalMap, implement removeByKey/removeByValue, ensure
that each key and value occurs only once

*/

template<class KeyType, class ValueType>
class rsKeyValueMap
{

public:

  /** \name Construction/Destruction */

  /** Destructor. */
  ~rsKeyValueMap()
  {
    for(unsigned int i=0; i<entriesSortedByKey.size(); i++)
      delete entriesSortedByKey[i];
  }


  /** \name Element Insertion */

  /** Inserts a key/value pair into the map. */
  void insertKeyValuePair(const KeyType key, const ValueType value)
  {
    rsKeyValuePair<KeyType, ValueType>* newEntry = new rsKeyValuePair<KeyType, ValueType>;
    newEntry->key   = key;
    newEntry->value = value;

    /*
    // for debug - to see what's wrong with the iter below:
    std::vector<int>::iterator  i1;                            // OK
    std::vector<int*>::iterator i2;                            // OK
    std::vector< rsKeyValuePair<int, int> *>::iterator i3;     // OK
    typename std::vector< rsKeyValuePair<KeyType, ValueType> *>::iterator i4; // error
    */

    typename std::vector< rsKeyValuePair<KeyType, ValueType>*>::iterator iter;
      // see http://gcc.gnu.org/ml/gcc-help/2008-01/msg00137.html why this "typename" is needed

    /*
    // pre C++11:
    // insert into the sorted-by-key array at the right position:
    rsKeyValuePairPointerLessByKey<KeyType, ValueType> compareByKey;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), newEntry,
      compareByKey);
    entriesSortedByKey.insert(iter, newEntry);

    // insert into the sorted-by-value array at the right position:
    rsKeyValuePairPointerLessByValue<KeyType, ValueType> compareByValue;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), newEntry,
      compareByValue);
    entriesSortedByValue.insert(iter, newEntry);
    */


    // insert into the sorted-by-key array at the right position:
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), newEntry,
      keyValuePairPointerLessByKey);
    entriesSortedByKey.insert(iter, newEntry);

    // insert into the sorted-by-value array at the right position:
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), newEntry,
      keyValuePairPointerLessByValue);
    entriesSortedByValue.insert(iter, newEntry);
  }

  // todo: write mergeKeyValueMaps, ...


  /** \name Element Removal */

  void removePairsWithKey(const KeyType key)
  {

    //...

  }


  /** \name Inquiry */

  /** If a pair with the given key is present in the map, the function returns the corresponding
  value and assigns the reference parameter "wasFound" to true. If the key is not present, it
  returns a default value determined by ValueType's default constructor and assigns the reference
  parameter "wasFound" to false. */
  ValueType getValueForKey(KeyType key, bool& wasFound)
  {
    rsKeyValuePair<KeyType, ValueType>* dummyEntry = new rsKeyValuePair<KeyType, ValueType>;
    dummyEntry->key = key;

    /*
    // pre C++11:
    rsKeyValuePairPointerLessByKey<KeyType, ValueType> compareByKey;
    typename std::vector<rsKeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), dummyEntry,
      compareByKey);
    */

    typename std::vector<rsKeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByKey.begin(), entriesSortedByKey.end(), dummyEntry,
      keyValuePairPointerLessByKey);

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

  /** If a pair with the given value is present in the map, the function returns the
  corresponding key and assigns the reference parameter "wasFound" to true. If the value is not
  present, it returns a default key determined by KeyType's default constructor and assigns the
  reference parameter "wasFound" to false. */
  KeyType getKeyForValue(ValueType value, bool& wasFound)
  {
    rsKeyValuePair<KeyType, ValueType>* dummyEntry = new rsKeyValuePair<KeyType, ValueType>;
    dummyEntry->value = value;

    /*
    // pre C++11:
    rsKeyValuePairPointerLessByValue<KeyType, ValueType> compareByValue;
    typename std::vector<rsKeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), dummyEntry,
      compareByValue);
    */

    typename std::vector<rsKeyValuePair<KeyType, ValueType>*>::iterator iter;
    iter = lower_bound(entriesSortedByValue.begin(), entriesSortedByValue.end(), dummyEntry,
      keyValuePairPointerLessByValue);

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
      return KeyType();
    }
  }

  /** Returns the number of entries (key/value pairs) in the map. */
  int getNumEntries() const { return entriesSortedByKey.size(); }

  /** Returns true when the map is empty, false otherwise. */
  bool isEmpty() const { return getNumEntries() == 0; }


  // \todo see rsArray and implement some operators...

protected:

  /** \name Data */

  std::vector<rsKeyValuePair<KeyType, ValueType>*> entriesSortedByKey;
  std::vector<rsKeyValuePair<KeyType, ValueType>*> entriesSortedByValue;
    // \todo maybe switch to rsArray<KeyValuePair*> here later

  // friend decalarations to facilitate testing (try to get rid of that - library code should not
  // be cluttered with stuff like that):
  friend bool testKeyValueMapInsert(std::string& reportString);

};

}

#endif
