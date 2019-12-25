#ifndef RAPT_FLAGS_H
#define RAPT_FLAGS_H

/**

This is a class to maintain a set of 8 boolean flags that can be individually accessed (to set
them and to retrieve them). The class is intended to save memory when a whole bunch of flags is
needed. Functionality wise, it's similar to the STL bitset.

\todo: write classes rslFlags16, rslFlags32, rslFlags64 in an entirely analogous way - maybe
this can be done with a template class? ...with a mask that is T::maxValue/2
or make a class rsFlagsN that can hold an arbitrary number of flags using an array of
rsFlags8 - or do both

\todo write a class rsFlagArray that contains an array of rsFlags8 and provides space
8*numArrayElements flags, access would work like that:
rsFlagArray::isFlagTrue(int index)
{
  return flagArray[index/8].isFlagTrue[index%8];
}
where the flagArray member is of type rsArrayTools<rsFlags8>. to optimize, one could perhaps use
flagArray[index>>3].isFlagTrue[index - index>>3]; -> verify this
maybe call this class rsBitArray */

class rsFlags8
{

public:


  /** \name Construction/Destruction */

  /** Standard constructor. Initializes all flags to zero. */
  rsFlags8() : flags(0)
  {

  }


  /** \name Operators */

  inline bool operator==(const rsFlags8& other) const
  {
    return flags == other.flags;
  }

  inline bool operator!=(const rsFlags8& other) const
  {
    return flags != other.flags;
  }


  /** \name Setup */

  /** Sets the flag with the given index to true. The index should be in the rnage 0...7, if
  it's greater than that, the function will do nothing. */
  inline void setFlagTrue(int index)
  {
    flags |= ((unsigned char)128 >> index);
  }

  /** Sets the flag with the given index to false. */
  inline void setFlagFalse(int index)
  {
    flags &= ~((unsigned char)128 >> index);
  }

  /** Sets all flags to true. */
  inline void setAllTrue()
  {
    flags = 255;
  }

  /** Sets all flags to false. */
  inline void setAllFalse()
  {
    flags = 0;
  }

  /** Sets the flag with the given index to the passed new boolean value.  */
  inline void setFlag(int index, bool newValue)
  {
    if( newValue == true )
      setFlagTrue(index);
    else
      setFlagFalse(index);
    // todo maybe this can be optimized
  }

  /** Toggles the flag into to the other state - a true flag becomes flase and vice versa. */
  inline void toggleFlag(int index)
  {
    if( isFlagTrue(index) )
      setFlagFalse(index);
    else
      setFlagTrue(index);
    // todo maybe this can be optimized
  }


  /** \name Inquiry */

  /** Returns whether the flag with the given index is true. The index should be in the range
  0...7, if it's greater than that, the function will return false. */
  inline bool isFlagTrue(int index) const
  {
    return ((flags << index) & (unsigned char)128) > 0;
  }

  /** Returns whether the flag with the given index is false. The index should be in the range
  0...7, if it's greater than that, the function will return false. */
  inline bool isFlagFalse(int index) const
  {
    return !isFlagTrue(index);
  }

protected:


  /** \name Data */

  rsUint8 flags;

};

//=================================================================================================

// move to a bit-twiddling file:

/** Sets the i-th bit of x to false. */
template <class T>
RS_INLINE void rsSetBitFalse(T& x, T i)
{
  x &= ~(T(1) << i);
}

/** Sets the i-th bit of x to true. */
template <class T>
RS_INLINE void  rsSetBitTrue(T& x, T i) 
{
  x |= T(1) << i;
}

/** Switches the i-th bit of x. */  
template <class T>
RS_INLINE void rsSwitchBit(T& x, T i)
{
  x ^= T(1) << i;
}

/** Fills an array of rsUint32 values with the bit pattern of a 64-bit datatype. */
template <class T>
RS_INLINE void rsGetBits64(T value, rsUint32 bits[64])
{
  // \todo convert endianness (maybe, test on different machines)
  rsUint64 *x   = reinterpret_cast<rsUint64*> (&value);
  rsUint64 mask = 0x8000000000000000ULL;
  for(int i = 0; i < 64; i++)
    bits[i] = (rsUint8) (((*x<<i) & mask) != 0);
}

/** Returns zero, if the i-th bit of x is zero, otherwise returns a value where only the i-th bit 
is set. */
template <class T>
RS_INLINE T rsGetMaskedBit(T x, T i)
{
  return x & T(1) << i;
}

/** Returns the number of true (1) bits in x. */
RS_INLINE rsUint64 rsGetNumTrueBits64(rsUint64 x)
{
  x -=  (x>>1)      & 0x5555555555555555UL;
  x  = ((x>>2)      & 0x3333333333333333UL) + (x & 0x3333333333333333UL);
  x  = ((x>>4) + x) & 0x0f0f0f0f0f0f0f0fUL;
  x *= 0x0101010101010101UL;
  return x>>56;
}

/** Returns true, iff the i-th bit of x is set ot 1 (true). */
template <class T>
RS_INLINE bool rsIsBitTrue(T x, T i)
{
  return rsGetMaskedBit(x, i) != 0;
}

/** Sets bit pattern of a 64-bit datatype from an ofarray rsUint32 values (assumed to be either 0 
or 1). */
template <class T>
RS_INLINE void rsSetBits64(T *value, rsUint32 bits[64])
{
  // \todo convert endianness (maybe, test on different machines)
  rsUint64 *x = reinterpret_cast<rsUint64*> (value);
  rsUint64 mask = 0x8000000000000000ULL;
  *x = 0;
  for(int i = 63; i >= 0; i--)
  {
    *x = *x >> 1;
    if( bits[i] != 0 )
      *x |= mask;
  }
}

//=================================================================================================
// rename to rsBitSet, use setTrue/setFalse/switch, isTrue

class rsFlagArray
{

public:


  /** \name Construction/Destruction */

  /** Standard constructor. Initializes all flags to zero. */
  inline rsFlagArray(rsUint64 numDesiredFlags);

  /** Destructor */
  inline ~rsFlagArray();





  /** \name Setup */

  /** Sets the flag with the given index to true. The index should be in the rnage 0...7, if
  it's greater than that, the function will do nothing. */
  inline void setFlagTrue(rsUint64 index)
  {
    rsSetBitTrue(flagGroups64[index/64], index%64);
  }

  /** Sets the flag with the given index to false. */
  inline void setFlagFalse(rsUint64 index)
  {
    rsSetBitFalse(flagGroups64[index/64], index%64);
  }

  /** Sets all flags to true. */
  inline void setAllTrue()
  {
    memset(flagGroups64, -1, (numFullGroups+1)*sizeof(rsUint64));
  }

  /** Sets all flags to false. */
  inline void setAllFalse()
  {
    memset(flagGroups64, 0, (numFullGroups+1)*sizeof(rsUint64));
  }

  /** Sets the flag with the given index to the passed new boolean value.  */
  //void setFlag(rsUint64 index, bool newValue);

  /** Toggles the flag into to the other state - a true flag becomes flase and vice versa. */
  //void toggleFlag(rsUint64 index);


  /** \name Inquiry */




  /** Returns whether the flag with the given index is true. */
  inline bool isFlagTrue(rsUint64 index) const
  {
    return rsIsBitTrue(flagGroups64[index/64], index%64);
  }

  /** Returns whether the flag with the given index is false. */
  inline bool isFlagFalse(rsUint64 index) const
  {
    return !isFlagTrue(index);
  }


  /** Returns true, iff all the flags are true. */
  inline bool areAllFlagsTrue();

  /** Returns true, iff all the flags are false. */
  inline bool areAllFlagsFalse();


  inline  rsUint64 getNumFlags() const
  {
    return numFlags;
  }

  /** Returns the number of flags which are true. */
  inline rsUint64 getNumTrueFlags() const;

  /** Returns the next index at which a true flag occurs, including the passed index in the
  search. If no true flags are found beyond and including the passed index, the returned index
  will be equal the number of flags (which is one unit out of range). */
  inline rsUint64 getNextTrueFlag(rsUint64 index)  const
  {
    while((index < numFlags) && (!isFlagTrue(index)))   // get rid of (int)
      index++;
    return index;
  }

protected:


  inline rsUint64 getNumTrueTailFlags() const;

  /** \name Data */

  rsUint64 numFlags;
  rsUint32 numFullGroups;
  rsUint32 tailLength;
  rsUint64 *flagGroups64;

};


inline rsFlagArray::rsFlagArray(rsUint64 numDesiredFlags)
{
  rsAssert( numDesiredFlags / 64 < UINT_MAX-1 );

  numFlags      = numDesiredFlags;
  numFullGroups = (rsUint32) rsMin(numFlags / 64, (rsUint64)UINT_MAX-1);
  tailLength    = (rsUint32) (numFlags - numFullGroups * 64);  
  flagGroups64  = new rsUint64[numFullGroups+1]; // +1 for the tail

  setAllFalse();
}

inline rsFlagArray::~rsFlagArray()
{
  delete[] flagGroups64;
}

inline bool rsFlagArray::areAllFlagsTrue()
{
  for(rsUint32 i = 0; i < numFullGroups; i++)
  {
    if( flagGroups64[i] != 0xFFFFFFFFFFFFFFFFULL )
      return false;
  }
  return getNumTrueTailFlags() == tailLength;
}

inline bool rsFlagArray::areAllFlagsFalse()
{
  for(rsUint32 i = 0; i < numFullGroups; i++)
  {
    if( flagGroups64[i] != 0 )
      return false;
  }
  return getNumTrueTailFlags() == 0;
}

inline rsUint64 rsFlagArray::getNumTrueFlags() const
{
  rsUint64 n = 0;
  for(rsUint32 i = 0; i < numFullGroups; i++)
    n += rsGetNumTrueBits64(flagGroups64[i]);
  n += getNumTrueTailFlags();
  return n;
}

inline rsUint64 rsFlagArray::getNumTrueTailFlags() const
{
  if( tailLength > 0 )
  {
    rsUint64 shift = 64 - (numFlags - numFullGroups*64);
    rsUint64 mask  = 0xFFFFFFFFFFFFFFFFULL >> shift;
    return rsGetNumTrueBits64(flagGroups64[numFullGroups] & mask);
  }
  return 0;
}

#endif