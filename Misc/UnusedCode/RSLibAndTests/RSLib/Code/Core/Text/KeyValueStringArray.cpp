using namespace RSLib;

// setup:

void rsKeyValueStringArray::setStringValue(const rsString &key, const rsString &value)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    keyValuePairs[index].setStringValue(value);
  else
    keyValuePairs.appendElement(rsKeyValueStringPair(key, value));
}

void rsKeyValueStringArray::setRealNumberValue(const rsString &key, const double value)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    keyValuePairs[index].setRealNumberValue(value);
  else
    keyValuePairs.appendElement(rsKeyValueStringPair(key, rsString(value)));
}

void rsKeyValueStringArray::setIntegerValue(const rsString &key, const int value)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    keyValuePairs[index].setIntegerValue(value);
  else
    keyValuePairs.appendElement(rsKeyValueStringPair(key, rsString(value)));
}

void rsKeyValueStringArray::setBooleanValue(const rsString &key, const bool value)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    keyValuePairs[index].setBooleanValue(value);
  else
    keyValuePairs.appendElement(rsKeyValueStringPair(key, rsString((int)value)));
}

// inquiry:

rsString rsKeyValueStringArray::getStringValue(const rsString& key, 
  const rsString& defaultReturnValue)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    return keyValuePairs[index].getStringValue();
  else
    return defaultReturnValue;
}

double rsKeyValueStringArray::getRealNumberValue(const rsString& key, 
  const double defaultReturnValue)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    return keyValuePairs[index].getRealNumberValue();
  else
    return defaultReturnValue;
}

int rsKeyValueStringArray::getIntegerValue(const rsString& key, const int defaultReturnValue)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    return keyValuePairs[index].getIntegerValue();
  else
    return defaultReturnValue;
}

bool rsKeyValueStringArray::getBooleanValue(const rsString& key, const bool defaultReturnValue)
{
  int index = getIndexOfKey(key);
  if( index != -1 )
    return keyValuePairs[index].getBooleanValue();
  else
    return defaultReturnValue;
}

int rsKeyValueStringArray::getIndexOfKey(const rsString& key)
{
  for(int i  =0; i < keyValuePairs.getNumElements(); i++)
  {
    if( keyValuePairs[i].hasKey(key) )
      return i;
  }
  return -1;
}
