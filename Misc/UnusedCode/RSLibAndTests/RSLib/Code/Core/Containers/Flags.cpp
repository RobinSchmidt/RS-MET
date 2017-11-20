using namespace RSLib;

void rsFlags8::setFlag(int index, bool newValue)
{
  if( newValue == true )
    setFlagTrue(index);
  else
    setFlagFalse(index);
  // todo maybe this can be optimized
}

void rsFlags8::toggleFlag(int index)
{
  if( isFlagTrue(index) )
    setFlagFalse(index);
  else
    setFlagTrue(index);
  // todo maybe this can be optimized
}





rsFlagArray::rsFlagArray(rsUint64 numDesiredFlags)
{
  rsAssert( numDesiredFlags / 64 < UINT_MAX-1 );

  numFlags      = numDesiredFlags;
  numFullGroups = (rsUint32) rsMin(numFlags / 64, (rsUint64)UINT_MAX-1);
  tailLength    = (rsUint32) (numFlags - numFullGroups * 64);  
  flagGroups64  = new rsUint64[numFullGroups+1]; // +1 for the tail

  setAllFalse();
}

rsFlagArray::~rsFlagArray()
{
  delete[] flagGroups64;
}

bool rsFlagArray::areAllFlagsTrue()
{
  for(rsUint32 i = 0; i < numFullGroups; i++)
  {
    if( flagGroups64[i] != 0xFFFFFFFFFFFFFFFFULL )
      return false;
  }
  return getNumTrueTailFlags() == tailLength;
}

bool rsFlagArray::areAllFlagsFalse()
{
  for(rsUint32 i = 0; i < numFullGroups; i++)
  {
    if( flagGroups64[i] != 0 )
      return false;
  }
  return getNumTrueTailFlags() == 0;
}

rsUint64 rsFlagArray::getNumTrueFlags() const
{
  rsUint64 n = 0;
  for(rsUint32 i = 0; i < numFullGroups; i++)
    n += rsGetNumTrueBits64(flagGroups64[i]);
  n += getNumTrueTailFlags();
  return n;
}

rsUint64 rsFlagArray::getNumTrueTailFlags() const
{
  if( tailLength > 0 )
  {
    rsUint64 shift = 64 - (numFlags - numFullGroups*64);
    rsUint64 mask  = 0xFFFFFFFFFFFFFFFFULL >> shift;
    return rsGetNumTrueBits64(flagGroups64[numFullGroups] & mask);
  }
  return 0;
}
