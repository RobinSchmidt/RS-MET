#include "KeyValueMapTests.h"

bool RSLib::testKeyValueMap(std::string &reportString)
{
  std::string testName = "rsKeyValueMap";
  bool testResult = true;

  testResult &= testKeyValueMapInsert(reportString);
  testResult &= testKeyValueMapFind(reportString);

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

void fillExampleMap(RSLib::rsKeyValueMap<int, rsString> &mapToFill)
{
  mapToFill.insertKeyValuePair(4,  rsString("Four"));
  mapToFill.insertKeyValuePair(1,  rsString("One"));
  mapToFill.insertKeyValuePair(5,  rsString("Five"));
  mapToFill.insertKeyValuePair(7,  rsString("Seven"));
  mapToFill.insertKeyValuePair(2,  rsString("Two"));
  mapToFill.insertKeyValuePair(3,  rsString("Three"));
  mapToFill.insertKeyValuePair(9,  rsString("Nine"));
  mapToFill.insertKeyValuePair(12, rsString("Twelve"));
  mapToFill.insertKeyValuePair(20, rsString("Twenty"));
  mapToFill.insertKeyValuePair(11, rsString("Eleven"));
  mapToFill.insertKeyValuePair(50, rsString("Fifty"));
  mapToFill.insertKeyValuePair(80, rsString("Eighty"));
}

bool RSLib::testKeyValueMapInsert(std::string &reportString)
{
  std::string testName = "rsKeyValueMapAppend";
  bool testResult = true;

  rsKeyValueMap<int, rsString> map;
  fillExampleMap(map);

  // the two sorted arrays inside the map should now look like:
  // entriesSortedByKey:   1-One, 2-Two, 3-Three, 4-Four, 5-Five, 7-Seven, 9-Nine, 11-Eleven, 
  //                       12-Twelve, 20-Twenty, 50-Fifty, 80-Eighty
  // entriesSortedByValue: 80-Eighty, 11-Eleven, 50-Fifty, 5-Five, 4-Four, 9-Nine, 1-One, 7-Seven, 
  //                       3-Three, 12-Twelve, 20-Twenty, 2-Two
  // check, if this is the case:

  testResult &= ( map.getNumEntries() == 12 );

  testResult &= ( map.entriesSortedByKey[0]->key  ==  1 );
  testResult &= ( map.entriesSortedByKey[1]->key  ==  2 );
  testResult &= ( map.entriesSortedByKey[2]->key  ==  3 );
  testResult &= ( map.entriesSortedByKey[3]->key  ==  4 );
  testResult &= ( map.entriesSortedByKey[4]->key  ==  5 );
  testResult &= ( map.entriesSortedByKey[5]->key  ==  7 );
  testResult &= ( map.entriesSortedByKey[6]->key  ==  9 );
  testResult &= ( map.entriesSortedByKey[7]->key  == 11 );
  testResult &= ( map.entriesSortedByKey[8]->key  == 12 );
  testResult &= ( map.entriesSortedByKey[9]->key  == 20 );
  testResult &= ( map.entriesSortedByKey[10]->key == 50 );
  testResult &= ( map.entriesSortedByKey[11]->key == 80 );

  testResult &= ( map.entriesSortedByValue[0]->key  == 80 );
  testResult &= ( map.entriesSortedByValue[1]->key  == 11 );
  testResult &= ( map.entriesSortedByValue[2]->key  == 50 );
  testResult &= ( map.entriesSortedByValue[3]->key  ==  5 );
  testResult &= ( map.entriesSortedByValue[4]->key  ==  4 );
  testResult &= ( map.entriesSortedByValue[5]->key  ==  9 );
  testResult &= ( map.entriesSortedByValue[6]->key  ==  1 );
  testResult &= ( map.entriesSortedByValue[7]->key  ==  7 );
  testResult &= ( map.entriesSortedByValue[8]->key  ==  3 );
  testResult &= ( map.entriesSortedByValue[9]->key  == 12 );
  testResult &= ( map.entriesSortedByValue[10]->key == 20 );
  testResult &= ( map.entriesSortedByValue[11]->key ==  2 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}

bool RSLib::testKeyValueMapFind(std::string &reportString)
{
  std::string testName = "rsKeyValueMapFind";
  bool testResult = true;

  rsKeyValueMap<int, rsString> map;
  fillExampleMap(map);

  bool wasFound;

  testResult &= ( map.getValueForKey( 1, wasFound) == rsString("One")    );
  testResult &= ( map.getValueForKey( 2, wasFound) == rsString("Two")    );
  testResult &= ( map.getValueForKey( 3, wasFound) == rsString("Three")  );
  testResult &= ( map.getValueForKey( 4, wasFound) == rsString("Four")   );
  testResult &= ( map.getValueForKey( 5, wasFound) == rsString("Five")   );
  testResult &= ( map.getValueForKey( 6, wasFound) == rsString("")       ); // 6 not in map
  testResult &= ( map.getValueForKey( 7, wasFound) == rsString("Seven")  );
  testResult &= ( map.getValueForKey( 8, wasFound) == rsString("")       ); // 8 not in map
  testResult &= ( map.getValueForKey( 9, wasFound) == rsString("Nine")   );
  testResult &= ( map.getValueForKey(11, wasFound) == rsString("Eleven") );
  testResult &= ( map.getValueForKey(12, wasFound) == rsString("Twelve") );
  testResult &= ( map.getValueForKey(20, wasFound) == rsString("Twenty") );
  testResult &= ( map.getValueForKey(50, wasFound) == rsString("Fifty")  );
  testResult &= ( map.getValueForKey(80, wasFound) == rsString("Eighty") );

  testResult &= ( map.getKeyForValue(rsString("One"),    wasFound) == 1  );
  testResult &= ( map.getKeyForValue(rsString("Two"),    wasFound) == 2  );
  testResult &= ( map.getKeyForValue(rsString("Three"),  wasFound) == 3  );
  testResult &= ( map.getKeyForValue(rsString("Four"),   wasFound) == 4  );
  testResult &= ( map.getKeyForValue(rsString("Five"),   wasFound) == 5  );
  testResult &= ( map.getKeyForValue(rsString("Six"),    wasFound) == 0  ); // 6 not in map
  testResult &= ( map.getKeyForValue(rsString("Seven"),  wasFound) == 7  );
  testResult &= ( map.getKeyForValue(rsString("Eight"),  wasFound) == 0  ); // 8 not in map
  testResult &= ( map.getKeyForValue(rsString("Nine"),   wasFound) == 9  );
  testResult &= ( map.getKeyForValue(rsString("Eleven"), wasFound) == 11 );
  testResult &= ( map.getKeyForValue(rsString("Twelve"), wasFound) == 12 );
  testResult &= ( map.getKeyForValue(rsString("Twenty"), wasFound) == 20 );
  testResult &= ( map.getKeyForValue(rsString("Fifty"),  wasFound) == 50 );
  testResult &= ( map.getKeyForValue(rsString("Eighty"), wasFound) == 80 );

  appendTestResultToReport(reportString, testName, testResult);
  return testResult;
}


